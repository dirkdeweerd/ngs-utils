mean_correction = function(m, bins_overall_mean, reads){
  cor <- bins_overall_mean / sum(m)
  reads * cor
}

sumchiscores = function(bins_list, bins_correct){
  sample_iterations <- 1:length(bins_list)
  bins_overall_mean = sum(as.numeric(unlist(bins_list))) / length(bins_list)
  bins_list_scaled  <- lapply(sample_iterations, (function(x) mean_correction(bins_list[[x]], 
                                                                              bins_overall_mean, bins_correct[[x]])))
  bins_sum_corrected = Reduce("+", bins_list_scaled)
  bins_scaled_expected = bins_sum_corrected / length(bins_list)
  bins_chi_score = bins_list_scaled
  for (i in 1:length(bins_list)) {
    bins_chi_score[[i]] = (bins_scaled_expected - bins_chi_score[[i]])^2 / bins_scaled_expected
  }
  return(Reduce("+", bins_chi_score))
}
correctsamples <- function(nipt_sample, chisumbins, degrees_of_freedom, chi_cutoff, xy_bins){
  construct_sample(autosomal_reads = lapply(X = nipt_sample$autosomal_chromosome_reads, FUN = correctbins, 
                                            chisumbins = chisumbins, degrees_of_freedom = degrees_of_freedom, 
                                            chi_cutoff = chi_cutoff), 
                   sex_reads = nipt_sample$sex_chromosome_reads, name = nipt_sample$name, 
                   correction_status_autosomal = c(nipt_sample$correction_status_autosomal, "Chi Corrected"),
                   correction_status_sex = nipt_sample$correction_status_sex)
  
}
correctsamples_include_xy <- function(nipt_sample, chisumbins, degrees_of_freedom, chi_cutoff, chisumbins_xy){
  construct_sample(autosomal_reads = lapply(X = nipt_sample$autosomal_chromosome_reads, FUN = correctbins, 
                                            chisumbins = chisumbins, degrees_of_freedom = degrees_of_freedom, 
                                            chi_cutoff = chi_cutoff), 
                   sex_reads = lapply(X = nipt_sample$sex_chromosome_reads, FUN = correctbins, chisumbins = chisumbins_xy,
                                      degrees_of_freedom = degrees_of_freedom, chi_cutoff = chi_cutoff),
                   name = nipt_sample$name, 
                   correction_status_autosomal = c(nipt_sample$correction_status_autosomal, "Chi Corrected"),
                   correction_status_sex = c(nipt_sample$correction_status_sex, "Chi Corrected"))
}
correctbins = function(bins, chisumbins, degrees_of_freedom, chi_cutoff){
  papinho <<- chisumbins
  chi_sum_bins_normalized <<- (chisumbins - degrees_of_freedom) / (sqrt( 2 * degrees_of_freedom))

  chi_sum_bins_correction_factor = as.matrix(chisumbins / degrees_of_freedom)
  index = which(chi_square_cut_off < chi_sum_bins_normalized) 
  bins[index] <- bins[index] / chi_sum_bins_correction_factor[index] 
  return(bins) 
}
chicorrect <- function(nipt_sample, control_samples, chi_cutoff = 3.5, include_XY = F){
  degrees_of_freedom  = length(control_samples$Samples) - 1
  controlbins <- lapply(X = control_samples$Samples, FUN = function(x) Reduce("+", x$autosomal_chromosome_reads))
  xy_bins <- lapply(X = control_samples$Samples, FUN = function(x) x$sex_chromosome_reads)
  chisumbins <<- sumchiscores(bins_list = controlbins, bins_correct = controlbins)
  if (include_XY == T){
    xy_controlbins <- lapply(xy_bins, function(x) Reduce("+", x))
    chisumbins_xy = sumchiscores(bins_list = controlbins, bins_correct = xy_controlbins)
    papinho <<- chisumbins_xy
    correctedcontrolsamples <- lapply(X = control_samples$Samples, FUN = correctsamples_include_xy, chisumbins = chisumbins, 
                                      degrees_of_freedom = degrees_of_freedom, chi_cutoff = chi_cutoff, 
                                      chisumbins_xy = chisumbins_xy)
    correctedsample <- correctsamples_include_xy(nipt_sample = nipt_sample, chisumbins = chisumbins, 
                                                 degrees_of_freedom = degrees_of_freedom, chi_cutoff = chi_cutoff,
                                                chisumbins_xy = chisumbins_xy)
  }
  else{
    correctedcontrolsamples <- lapply(X = control_samples$Samples, FUN = correctsamples, chisumbins = chisumbins, 
                                      degrees_of_freedom = degrees_of_freedom,
                                      chi_cutoff = chi_cutoff)
    correctedsample <- correctsamples(nipt_sample = nipt_sample, chisumbins = chisumbins, degrees_of_freedom = degrees_of_freedom,
                                      chi_cutoff = chi_cutoff)
  }
  return (list(sample = correctedsample, control_group = as_control_group(correctedcontrolsamples, control_samples$Description)))
}