chrfractions <- function(binned.sample)
{
  if(is.null(attr(binned.sample, "class"))){
    print("Object is not of type Sample")
  }
  else  UseMethod("chrfractions", binned.sample)
}

chrfractions.SeparatedStrands <- function(binned.sample)
{
 unlist(lapply(X = binned.sample$reads, FUN = function(x) (rowSums(x) / sum(x)) / 2))
}

chrfractions.CombinedStrands <- function(binned.sample)
{
  unlist(lapply(X = binned.sample$reads, FUN = function(x) (rowSums(x) / sum(x)) / 2))
}
