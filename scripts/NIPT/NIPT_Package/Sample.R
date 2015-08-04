construct.sample <- function(reads, name)
{
  new.sample <- list(reads = reads, name = name)
  if(length(reads) == 1)
  {
    type = "CombinedStrands"
  }
  
  if(length(reads) == 2)
  {
    type = "SeparatedStrands"
  }
  class(new.sample) <- c("Sample", type)
  return(new.sample)
}