as_control_group <- function(nipt_samples, control_group_type = generic_control_group){
  if (length(unique(sapply(nipt_samples, getsamplenames))) != length(nipt_samples)){
    cat("Warning, there appear to be duplicate sample names in control group \n")
  }
  if ((length(unique(sapply(nipt_samples, getcorrectionstatus))) > 1)){
    cat("Warning, more than one correction_status in control group")
  }
  
  if ((length(unique(sapply(nipt_samples, getstrandtype))) != 1)){
     stop("More than one strand type in control group", call. = F)
  }
  control_group <- list("Samples" = nipt_samples, "correction_status" = unique(sapply(nipt_samples, getcorrectionstatus)),
                        "Description" = control_group_type)
  
  class(control_group) <- c("Control Group", unique(sapply(nipt_samples, getstrandtype)))
  return(control_group)
}