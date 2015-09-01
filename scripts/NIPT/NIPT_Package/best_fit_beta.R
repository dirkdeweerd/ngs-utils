matchcontrolgroup <-function(nipt_sample, control_group, n.of.samples)
{
  control_group_samples <- control_group$Samples
  control.group.fractions <- cbind(NULL, sapply(X = control_group_samples , FUN = chrfractions))
  sample.fractions <- sapply(list(nipt_sample), chrfractions)
  fractions.squared <- apply(X = control.group.fractions, MARGIN = 2, FUN = function(x, y) {(x - y)^2}, y = sample.fractions)
  control.chromosomes <- getcontrolchromosomes(nipt_sample)
  rownames(fractions.squared) <- c(rownames_separated_forward_autosomal, rownames_separated_reverse_autosomal)
  sum.of.squares <- apply(fractions.squared[control.chromosomes,], 2, sum)
  sorted.scores <- order(sum.of.squares)
  return(as_control_group(nipt_samples = control_group$Samples[sorted.scores[1:n.of.samples]], 
                   control_group_type = paste("Fitted to", nipt_sample$name)))
}

