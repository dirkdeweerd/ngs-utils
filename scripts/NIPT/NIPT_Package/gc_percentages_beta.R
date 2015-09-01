dna_letters <- c("A", "T", "C", "G", "N")
chromosome_gc_percentage_list <- list()

count_bases <- function(chromosome, min_index, end_of_bin){
  genomic_interval <- chromo_current[min_index:intervals[bin]]
}


for (chromosome in 23:24){
  chromo_current <- BSgenome.Hsapiens.UCSC.hg19[[chromosome]]
  intervals <- c(seq(from = (bin.size), to = length(chromo_current), by = (bin.size), ), length(chromo_current))
  gc_percentages_matrix <- NULL
  min_index <- 1
  
  print(chromo_current[intervals[length(intervals) -1]:intervals[length(intervals)]])
  for (bin in 1:length(intervals)){
    genomic_interval <- chromo_current[min_index:intervals[bin]]
    counted_bases <- letterFrequency(x = genomic_interval, letters = dna_letters)
    gc_percentages_matrix <- rbind(gc_percentages_matrix, c(counted_bases[1], counted_bases[2], counted_bases[3], 
                                                            counted_bases[4], counted_bases[5]))
    min_index <- intervals[bin] + 1
  }
  cat("Chromosome ", bin, "Ready \n" )
  chromosome_gc_percentage_list[[chromosome]] <- gc_percentages_matrix
}
