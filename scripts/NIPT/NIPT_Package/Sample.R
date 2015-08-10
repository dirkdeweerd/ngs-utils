construct_sample <- function(reads, name)
{
  
  if(length(reads) == 1)
  {
    rownames(reads[[1]]) <- as.character(1:22)
    type = "CombinedStrands"
  }
  if(length(reads) == 2)
  {
    rownames(reads[[1]]) <- paste0(1:22, "F")
    rownames(reads[[2]]) <- paste0(1:22, "R")
    type = "SeparatedStrands"
  }
  
  new.sample <- list(reads = reads, name = name)
  
  class(new.sample) <- c("Sample", type)
  return(new.sample)
}