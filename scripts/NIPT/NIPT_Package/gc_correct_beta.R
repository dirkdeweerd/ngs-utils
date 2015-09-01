GetCorrectionFactor <- function(sample, span, gc.percentages, median_autosomals){
  fit.loess <-loess(sample ~ gc.percentages, span = span)
  fitted_values <<- cbind(fit.loess$x, median_autosomals /  fit.loess$fitted)
  return(list(median_autosomals / fit.loess$fitted, fitted_values))
}
correct.reads <- function(sample, correction.factor, indices){
  sample[indices] <- sample[indices] * correction.factor
  return(sample)
}
correct_reads_sex_chromosome <- function(bins, gc.percentages,fitted_values){
  for (chromosome in 1:nrow(bins))
  {
    gc_bins <- unname(which(gc.percentages[22+chromosome,] > 0))
    print(length(gc_bins))
    for (bin in gc_bins){
      bins[chromosome, bin] <- bins[chromosome, bin] * (fitted_values[order(abs(gc.percentages[(22+chromosome),bin] - fitted_values[,1]))[1],][2])
    }
  }
  return(bins)
}
gc.correct <- function(nipt_sample, span = 0.75, include_XY = F){
  sample <- Reduce("+", nipt_sample$autosomal_chromosome_reads)
  sex_chromosome_reads <- nipt_sample$sex_chromosome_reads
  rownames(sample) <- autosomal_chromosomes
  gc_percentages_complete <- gc.percentages
  gc_percentages_autosomal <- gc.percentages[autosomal_chromosomes, ]
  
  indices_autosomals <- (which(gc_percentages_autosomal > 0 & sample))
  median_autosomals <- median(sample[indices_autosomals])
  correction_factor_values <- GetCorrectionFactor(sample = as.vector(sample[indices_autosomals]), span = span, 
                                                  gc.percentages = as.vector(gc_percentages_autosomal[indices_autosomals]), 
                                                  median_autosomals = median_autosomals)
  correction.factor <- correction_factor_values[[1]]
  fitted_values <- correction_factor_values[[2]]
  corrected_autosomal <- lapply(nipt_sample$autosomal_chromosome_reads,  correct.reads, correction.factor = correction.factor, 
                                indices = indices_autosomals)
  if (include_XY == T){
    corrected_sex <- lapply(X = nipt_sample$sex_chromosome_reads, FUN = correct_reads_sex_chromosome, 
                            gc.percentages = gc_percentages_complete, fitted_values = fitted_values)
    corrected.sample <- construct_sample(autosomal_reads = corrected_autosomal, sex_reads = corrected_sex, 
                                         name = nipt_sample$name, correction_status_autosomal = c(nipt_sample$correction_status_autosomal, GCcorrected),
                                         correction_status_sex = c(nipt_sample$correction_status_sex, GCcorrected))
  }
  else{
  corrected.sample <- construct_sample(autosomal_reads = corrected_autosomal, sex_reads = nipt_sample$sex_chromosome_reads, 
                                       name = nipt_sample$name, correction_status_autosomal = c(nipt_sample$correction_status_autosomal, GCcorrected),
                                       correction_status_sex = c(nipt_sample$correction_status_sex))
  }
  return(corrected.sample)
}
