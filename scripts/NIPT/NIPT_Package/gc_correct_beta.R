GetCorrectionFactor <- function(sample, span, gc.percentages)
{
  fit.loess <-loess(sample ~ gc.percentages, span = span)
  return(median(sample) / fit.loess$fitted)
}
correct.reads <- function(sample, correction.factor, indices)
{
  sample[indices] <- sample[indices] * correction.factor
  return(sample)
}

gc.correct <- function(nipt_sample, span)
{
  if(is.null(attr(nipt_sample, "class"))){
    print("Object is not of type Sample")
  }
  else  UseMethod("gc.correct", nipt_sample)
}
gc.correct.SeparatedStrands <- function(nipt_sample, span = 0.75)
{
  sample <- Reduce("+", nipt_sample$reads)
  indices <- which(gc.percentages > 0 & sample > 0)
  correction.factor <- GetCorrectionFactor(sample = as.vector(sample[indices]), span = span, gc.percentages = as.vector(gc.percentages[indices]))
  corrected.sample <- construct_sample(reads = lapply(X = nipt_sample$reads, FUN= correct.reads, correction.factor, indices), name = nipt_sample$name)
  return(corrected.sample)
}

gc.correct.CombinedStrands <- function(nipt_sample, span = 0.75)
{
  sample <- Reduce("+", nipt_sample$reads)
  indices <- which(gc.percentages > 0 & sample > 0)
  correction.factor <- GetCorrectionFactor(sample = as.vector(sample[indices]), span = span, gc.percentages = as.vector(gc.percentages[indices]))
  corrected.sample <- construct_sample(reads = lapply(X = nipt_sample$reads, FUN= correct.reads, correction.factor, indices), name = nipt_sample$name)
  return(corrected.sample)
}