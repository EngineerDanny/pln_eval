reg.dir <- "/scratch/da2343/hmpv13_06_20_reg"
reg.dir <- "/scratch/da2343/hmpv35_06_20_reg"
reg.dir <- "/scratch/da2343/amgut1_06_20_reg"
reg.dir <- "/scratch/da2343/amgut2_06_20_reg"
reg.dir <- "/scratch/da2343/baxter_crc_06_22_reg"
reg.dir <- "/scratch/da2343/crohns_06_22_reg"
reg.dir <- "/scratch/da2343/glne007_06_22_reg"
reg.dir <- "/scratch/da2343/hmp2prot_06_22_reg"
reg.dir <- "/scratch/da2343/hmp216S_06_22_reg"
reg.dir <- "/scratch/da2343/ioral_06_22_reg"
reg.dir <- "/scratch/da2343/mixmpln_06_22_reg"
#reg.dir <- "/scratch/da2343/soilrep_06_22_reg"
#reg.dir <- "/scratch/da2343/MovingPictures_06_22_reg"
#reg.dir <- "/scratch/da2343/qa10394_06_22_reg"
#reg.dir <- "/scratch/da2343/TwinsUK_06_22_reg"

reg <- batchtools::loadRegistry(reg.dir, writeable = TRUE)
jobs.after <- batchtools::getJobTable(reg=reg)
table(jobs.after$error)
jobs.after[!is.na(error), .(error, task_id=sapply(prob.pars, "[[", "task_id"))][25:26]
ids <- jobs.after[is.na(error), job.id]

bmr = mlr3batchmark::reduceResultsBatchmark(ids, reg = reg, store_backends = FALSE)
score.dt.poisson <- bmr$score(poisson_measure)
score.dt.rmse <- bmr$score(msr("regr.rmse"))
plot_data <- merge(
  score.dt.poisson[, .(task_id, learner_id, iteration, regr.poisson_deviance)],
  score.dt.rmse[, .(task_id, learner_id, iteration, regr.rmse)],
  by = c("task_id", "learner_id", "iteration")
)

plot_data <- extract_scores_directly(ids, reg, jobs.after)

data.table::fwrite(plot_data, "/projects/genomic-ml/da2343/PLN/pln_eval/data/poisson_vs_gaussian/mixmpln.csv")