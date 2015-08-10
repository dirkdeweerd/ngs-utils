control.chromosomes <- c(1:12,14:17,19:20,22)

matchcontrolgroup <- function(nipt_sample, control_group, n.of.samples)
{
  if(is.null(attr(nipt_sample, "class"))){
    print("Object is not of type Sample")
  }
  else  UseMethod("matchcontrolgroup", nipt_sample)
}

matchcontrolgroup.SeparatedStrands <-function(nipt_sample, control_group, n.of.samples)
{
  control.group.fractions <- cbind(NULL, sapply(X = control_group, FUN = chrfractions))
  sample.fractions <- chrfractions(nipt_sample = nipt_sample)
  fractions.squared <- apply(X = control.group.fractions, MARGIN = 2, FUN = function(x, y) {(x - y)^2}, y = sample.fractions)
  control.chromosomes <- c(control.chromosomes, control.chromosomes+22)
  sum.of.squares <- apply(fractions.squared[control.chromosomes,], 2, sum)
  sorted.scores <- order(sum.of.squares)
  
  return (control_group[sorted.scores[1:n.of.samples]])
}

matchcontrolgroup.CombinedStrands <-function(nipt_sample, control_group, n.of.samples)
{
  control.group.fractions <- cbind(NULL, sapply(X = control_group, FUN = chrfractions))
  sample.fractions <- chrfractions(nipt_sample = nipt_sample)
  fractions.squared <- apply(X = control.group.fractions, MARGIN = 2, FUN = function(x, y) {(x - y)^2}, y = sample.fractions)
  sum.of.squares <- apply(fractions.squared[control.chromosomes,], 2, sum)
  sorted.scores <- order(sum.of.squares)
  
  return (control_group[sorted.scores[1:n.of.samples]])
}


