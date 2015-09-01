setnames <- paste("Set", 1:n.models)


collapse_result <- function(result, value){
  return(result[[value]])
}
collapse_prediction_sets <- function(result){
  gsub(pattern = ",", replacement = " ", x = toString(result$predictors))
}
collapse_results <- function(result_set){
  control_group_Z_scores <- Reduce(cbind, lapply(result_set, collapse_result, value = "control_group_Z_scores"))
  prediction_statistics <- rbind("Z score sample" = as.numeric(sapply(result_set, collapse_result, value = "sample_Z_score")),
                                 "P value shapiro" = as.numeric(sapply(result_set, collapse_result, value = "shapiro_P_value")),
                                 "Predictor chromosomes" = sapply(result_set, collapse_prediction_sets))
  colnames(control_group_Z_scores) <- setnames
  colnames(prediction_statistics) <- setnames
  nipt_result <- list("PredictionStatistics" = prediction_statistics, "ControlZScores" = control_group_Z_scores)
  return(nipt_result)
}