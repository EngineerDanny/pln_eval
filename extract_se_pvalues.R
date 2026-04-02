library(data.table)
base_dir <- "/projects/genomic-ml/da2343/PLN/pln_eval"
source(file.path(base_dir, "load_source.R"))

datasets <- c(
  "20170526.Castro-NallarE_2015.metaphlan_bugs_list.oralcavity_genus",
  "20170526.ChngKR_2016.metaphlan_bugs_list.skin_genus",
  "20170526.HMP_2012.metaphlan_bugs_list.nasalcavity_species",
  "crc_zeller",
  "ob_ross_genus",
  "hiv_lozupone_family",
  "nash_ob_baker_family",
  "crc_zhao_genus",
  "20190424.CosteaPI_2017.metaphlan_bugs_list.stool_genus",
  "amgut2",
  "nash_chan_family",
  "amgut1",
  "t1d_alkanani_family",
  "cdi_schubert_family",
  "2021-03-31.MehtaRS_2018.relative_abundance_genus",
  "2021-03-31.ShaoY_2019.relative_abundance_genus",
  "ob_zeevi_family",
  "otu_table_psn_v13_family",
  "diabimmune_karelia_16s_family",
  "mbqc_integrated_otus_family"
)

results <- lapply(datasets, function(ds) {
  message("Processing: ", ds)
  data_env <- new.env(parent = emptyenv())
  load(file.path(base_dir, "bmr", paste0(ds, ".RData")), envir = data_env)
  score_dt <- as.data.table(data_env$bmr$score(poisson_measure))
  score_dt[, algorithm := fcase(
    learner_id == "regr.featureless", "Featureless",
    learner_id == "regr.lasso", "GLMNet",
    learner_id == "regr.pln", "PLN"
  )]
  score_dt <- score_dt[!is.na(algorithm)]

  summary_dt <- score_dt[, .(
    mean_deviance = mean(regr.poisson_deviance),
    se_deviance   = sd(regr.poisson_deviance) / sqrt(.N)
  ), by = algorithm]

  wide <- dcast(score_dt[algorithm %in% c("GLMNet", "PLN")],
                task_id + iteration ~ algorithm,
                value.var = "regr.poisson_deviance")
  test <- t.test(wide$PLN, wide$GLMNet, paired = TRUE)

  glmnet_row <- summary_dt[algorithm == "GLMNet"]
  pln_row    <- summary_dt[algorithm == "PLN"]

  data.table(
    source_id      = ds,
    glmnet_se      = glmnet_row$se_deviance,
    pln_se         = pln_row$se_deviance,
    p_value        = test$p.value
  )
})

out <- rbindlist(results)
fwrite(out, file.path(base_dir, "out", "se_pvalues_20datasets.csv"))
cat("Done. Written to out/se_pvalues_20datasets.csv\n")
print(out)
