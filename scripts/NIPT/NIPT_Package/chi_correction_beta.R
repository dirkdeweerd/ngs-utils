sumchiscores = function(bins_list) {
  bins_overall_mean = sum(as.numeric(unlist(bins_list))) / length(bins_list)
  mean_correction = function(m) m * bins_overall_mean / sum(m)
  bins_list_scaled = lapply(bins_list, mean_correction)
  bins_sum_corrected = Reduce("+", bins_list_scaled)
  bins_scaled_expected = bins_sum_corrected / length(bins_list)
  bins_chi_score = bins_list_scaled
  for (i in 1:length(bins_list)) {
    bins_chi_score[[i]] = (bins_scaled_expected - bins_chi_score[[i]])^2 / bins_scaled_expected
  }
  return(Reduce("+", bins_chi_score))
}
correctsamples <- function(nipt_sample, chisumbins, degrees_of_freedom){
  construct_sample(reads = lapply(nipt_sample$reads, correctbins, chisumbins = chisumbins, 
                                  degrees_of_freedom = degrees_of_freedom), name = nipt_sample$name)
}
correctbins = function(bins, chisumbins, degrees_of_freedom) {
  chi_sum_bins_normalized = (chisumbins - degrees_of_freedom) / (sqrt( 2 * degrees_of_freedom))
  chi_sum_bins_correction_factor = as.matrix(chisumbins / degrees_of_freedom)
  index = which(chi_square_cut_off < chi_sum_bins_normalized) # Variation between these bins is considered too large
  bins[index] <- bins[index] / chi_sum_bins_correction_factor[index] # Therefore, we correct them here
  return(bins) 
}
chicorrect <- function(nipt_sample, control_samples){
  degrees_of_freedom  = length(control_samples) - 1
  controlbins <- lapply(X = control_samples, FUN = function(x) Reduce("+", x$reads))
  chisumbins = sumchiscores(bins_list = controlbins)
  correctedcontrolsamples <- lapply(X = control_samples, FUN = correctsamples, chisumbins = chisumbins, degrees_of_freedom = degrees_of_freedom)
  correctedsample <- correctsamples(nipt_sample = nipt_sample, chisumbins = chisumbins, degrees_of_freedom = degrees_of_freedom)
  return (list(correctedsample, correctedcontrolsamples))
}