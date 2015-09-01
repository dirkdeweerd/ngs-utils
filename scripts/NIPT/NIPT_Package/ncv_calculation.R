calculate_ncv_score <- function(nipt_sample, ncv_template){
  reads <- (rowSums(sumfandrautosomal(nipt_sample)))
  normalized_chromosome <- reads[ncv_template$focus_chromosome] / reads[ncv_template$denominators]
  ncv_score <- (normalized_chromosome - ncv_template$control_group_statistics[1]) / ncv_template$control_group_statistics[2]
  
  ncv_template$sample_score <- ncv_score
  ncv_template$sample_name <- nipt_sample$name  
  
  class(ncv_template)[1] <- NCV_result_class
  
  return(ncv_template)
}