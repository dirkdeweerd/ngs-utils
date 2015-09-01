ncv_template <- function(denominators, chromo_focus, nipt_sample_names, correction_status, scores, 
                         potential_denominators, statistics, type){
  new_ncv_template <- list(denominators = denominators, focus_chromosome = as.character(chromo_focus), 
                           nipt_sample_names = nipt_sample_names, correction_status = correction_status,
                           control_group_Zscores = scores, 
                           potential_denominators = potential_denominators, control_group_statistics = statistics)
  
  class(new_ncv_template) <- c(NCV_template_class, type)
  
  return(new_ncv_template)
}

