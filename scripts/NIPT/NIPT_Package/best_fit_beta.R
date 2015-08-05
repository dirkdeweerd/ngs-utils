control.chromosomes <- c(1:12,14:17,19:20,22)

matchcontrolgroup <- function(sample, samplelist, n.of.samples)
{
  if(is.null(attr(sample, "class"))){
    print("Object is not of type Sample")
  }
  else  UseMethod("matchcontrolgroup", sample)
}

matchcontrolgroup.SeparatedStrands <-function(sample, samplelist, n.of.samples)
{
  control.group.fractions <- cbind(NULL, sapply(X = samplelist, FUN = chrfractions))
  sample.fractions <- chrfractions(binned.sample = sample)
  fractions.squared <- apply(X = control.group.fractions, MARGIN = 2, FUN = function(x, y) {(x - y)^2}, y = sample.fractions)
  control.chromosomes <- c(control.chromosomes, control.chromosomes+22)
  sum.of.squares <- apply(fractions.squared[control.chromosomes,], 2, sum)
  sorted.scores <- order(sum.of.squares)
  
  return (samplelist[sorted.scores[1:n.of.samples]])
}

matchcontrolgroup.CombinedStrands <-function(sample, samplelist, n.of.samples)
{
  control.group.fractions <- cbind(NULL, sapply(X = samplelist, FUN = chrfractions))
  sample.fractions <- chrfractions(binned.sample = sample)
  fractions.squared <- apply(X = control.group.fractions, MARGIN = 2, FUN = function(x, y) {(x - y)^2}, y = sample.fractions)
  sum.of.squares <- apply(fractions.squared[control.chromosomes,], 2, sum)
  sorted.scores <- order(sum.of.squares)
  
  return (samplelist[sorted.scores[1:n.of.samples]])
}


