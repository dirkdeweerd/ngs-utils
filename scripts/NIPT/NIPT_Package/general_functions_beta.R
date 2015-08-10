chrfractions <- function(nipt_sample)
{
  if(is.null(attr(nipt_sample, "class"))){
    print("Object is not of type Sample")
  }
  else  UseMethod("chrfractions", nipt_sample)
}

chrfractions.SeparatedStrands <- function(nipt_sample)
{
 unlist(lapply(X = nipt_sample$reads, FUN = function(x) (rowSums(x) / sum(x)) / 2))
}

chrfractions.CombinedStrands <- function(nipt_sample)
{
  unlist(lapply(X = nipt_sample$reads, FUN = function(x) (rowSums(x) / sum(x))))
}

samplenames <- function(nipt_sample){
  nipt_sample$nam
}

collapse_results <- function(result, name){
  result[name]  
}