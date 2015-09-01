library(data.table)
library(Rsamtools)

.unlist <- function (x){
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}
bin.bed.sample <- function (bed.filepath)
{
  col.names = c("chrom", "startpos", "stoppos", "readname", "score", "strand")
  
  reads.data.frame <- fread(input = bed.filepath, header = F, verbose = F,
                            sep = "\t")
  setnames(x = reads.data.frame, col.names)
  bin <- bin.reads(reads.data.frame = reads.data.frame, bin.size = bin.size)
  return(bin)
}

bin.bam.sample <- function(bam.filepath, do.sort=FALSE, separate.strands=FALSE, name = basename(sub("^([^.]*).*", "\\1", bam.filepath))){
  if (do.sort == TRUE){
    temp <- sortBam(file = bam.filepath, destination = tempfile())
    bam <- scanBam(file = temp)
  }
  else{
    bam <- scanBam(bam.filepath)
  }
  bam.names <- names(bam[[1]])
  
  lst <- lapply(bam.names, function(y) .unlist(lapply(bam, "[[", y)))
  
  complete.df <- do.call("DataFrame", lst)
  
  names(complete.df) <- bam.names
  
  reads.data.frame <- data.frame(complete.df$rname, complete.df$strand, complete.df$pos)
  colnames(reads.data.frame) <- c("chrom", "strand", "startpos")
  if( separate.strands == TRUE)  {
    reads.list <- split(x = reads.data.frame, f = droplevels(reads.data.frame$strand[reads.data.frame$strand != "*"]))
  }
  else{
  reads.list <- list(reads.data.frame) 
  }
  binned_reads <- lapply(X = reads.list, FUN = bin.reads)
  testframe <<- binned_reads
  autosomal_reads <- lapply(X = binned_reads, FUN = splitchromosomes, chromosomes = autosomal_chromosomes)
  sex_reads <- lapply(X = binned_reads, FUN = splitchromosomes, chromosomes = sex_chromosomes)
  new_sample <- construct_sample(autosomal_reads = autosomal_reads, sex_reads = sex_reads, name = name, 
                                 correction_status_autosomal = Uncorrected, correction_status_sex = Uncorrected)
  
  return(new_sample)
}
bin.reads <- function(reads.data.frame){
  bin <- matrix(data = 0, nrow = 0, ncol = n.bins , dimnames = list(NULL, 1:n.bins))
  if ((length(unique(reads.data.frame$chrom[1:100]))) != 1){
    stop("BAM file appears to be unsorted", call. = F)
  }
  chromos <- rle(x = as.vector(reads.data.frame$chrom))
  min.read <- 0
  max.read <- 1
  for (chromo in 1:n_total_chromosomes){
    max.read <- (max.read -1)+ chromos[[1]][chromo] 
    reads <- sapply(X = unique(reads.data.frame[min.read:max.read,], MARGIN = 2)$startpos,  
                    FUN =  getbins, bin.size=bin.size)
    min.read <- max.read
    bins <- tabulate(reads, nbins = n.bins)
    bin <- rbind(bin, bins)
  }
  return (bin)
}

getbins <- function(pos, bin.size){
  pos %/% bin.size + 1
}