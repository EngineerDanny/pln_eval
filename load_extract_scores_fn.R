


# Extract scores directly - much faster approach
extract_scores_directly <- function(job_ids, reg, job_table) {
  score_list <- list()
  
  cat("Extracting scores from", length(job_ids), "jobs...\n")
  
  for(i in seq_along(job_ids)) {
    if(i %% 100 == 0) cat("Processed", i, "of", length(job_ids), "jobs\n")
    
    tryCatch({
      job_id <- job_ids[i]
      
      # Load the result
      result <- batchtools::loadResult(job_id, reg = reg)
      
      pred_obj <- mlr3::as_prediction(result$prediction$test)
      
      # Score with both measures
      poisson_score <- pred_obj$score(measures = poisson_measure)
      rmse_score <- pred_obj$score(measures = msr("regr.rmse"))
      
      # Get task_id and learner_id from job table
      job_row <- job_table[job.id == job_id]
      task_id <- job_row$prob.pars[[1]]$task_id  # Extract from prob.pars
      learner_id <- job_row$algo.pars[[1]]$learner_id  # Extract from algo.pars
      
      # You might need to extract iteration from the job structure
      # This might be in the job parameters or you might need to derive it
      iteration <- i  # Placeholder - might need adjustment
      
      # Create the row with both scores
      score_row <- data.table(
        task_id = task_id,
        learner_id = learner_id,
        iteration = iteration,
        regr.poisson_deviance = poisson_score,
        regr.rmse = rmse_score
      )
      
      score_list[[i]] <- score_row
      
    }, error = function(e) {
      cat("Error in job", job_ids[i], ":", e$message, "\n")
    })
  }
  
  # Combine all results
  rbindlist(score_list[!sapply(score_list, is.null)])
}

