control.chromosomes <- c(1:12,14:17,19:20,22)

BestControlSet <-function(control.group, sample, n.of.samples)
{
  control.group.fractions <- rbind(NULL, sapply(X = control.group, FUN = ChromosomalFractionPerSample))
  sample.fractions <- t(rbind(NULL, ChromosomalFractionPerSample(bins = sample)))
  fractions.squared <- apply(X = control.group.fractions, MARGIN = 2, FUN = function(x, y) {(x - y)^2}, y = sample.fractions)
  if (nrow(fractions.squared) == 44)
  {
    control.chromosomes <- c(control.chromosomes, control.chromosomes+22)
  }
  
  sum.of.squares <- apply(fractions.squared[control.chromosomes,], 2, sum)
  sorted.scores <- order(sum.of.squares)
  
  return (control.group[sorted.scores[1:n.of.samples]])
}


