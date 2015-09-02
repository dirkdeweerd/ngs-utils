chrfractions <- function(nipt_sample)
{
  if(is.null(attr(nipt_sample, "class"))){
    print("Object is not of type Sample")
  }
  else  UseMethod("chrfractions", nipt_sample)
}

getcontrolchromosomes <- function(nipt_sample)
{
  if(is.null(attr(nipt_sample, "class"))){
    print("Object is not of type Sample")
  }
  else  UseMethod("getcontrolchromosomes", nipt_sample)
}

getfractionscontrolgroup <- function(nipt_control_group)
{
  if(is.null(attr(nipt_control_group, "class"))){
    print("Object is not of type NIPT Control Group")
  }
  else  UseMethod("getfractionscontrolgroup", nipt_control_group)
}
getreadscontrolgroup <- function(nipt_control_group)
{
  if(is.null(attr(nipt_control_group, "class"))){
    print("Object is not of type NIPT Control Group")
  }
  else  UseMethod("getreadscontrolgroup", nipt_control_group)
}

chrfractions.SeparatedStrands <- function(nipt_sample){
 sapply(X = nipt_sample$autosomal_chromosome_reads, FUN = function(x) (rowSums(x) / sum(x)) / 2)
}

chrfractions.CombinedStrands <- function(nipt_sample){
  sapply(X = nipt_sample$autosomal_chromosome_reads, FUN = function(x) (rowSums(x) / sum(x)))
}

chrreads <- function(nipt_sample){
  sapply(X = nipt_sample$autosomal_chromosome_reads, FUN = function(x) rowSums(x)) 
}

retrieve_fractions_of_interest <- function(nipt_sample, chromo_focus, chromosomal_fracs){
  if(is.null(attr(nipt_sample, "class"))){
    print(error_message)
  }
  else UseMethod("retrieve_fractions_of_interest", nipt_sample)
}

retrieve_fractions_of_interest.CombinedStrands <- function(nipt_sample, chromo_focus, chromosomal_fracs){
  chromosomal_fracs[as.character(chromo_focus),] 
}


retrieve_fractions_of_interest.SeparatedStrands <- function(nipt_sample, chromo_focus, chromosomal_fracs){
  chromosomal_fracs[paste0(chromo_focus, "F"),] + chromosomal_fracs[paste0(chromo_focus, "R"), ]
}

setrownamesmatrix <- function(nipt_matrix){
  if (nrow(nipt_matrix) == 22){
    print(class(nipt_matrix))
    rownames(nipt_matrix) <- rownames_combined_autosomal
  }
  if (nrow(nipt_matrix) == 44){
    rownames(nipt_matrix) <- c(rownames_separated_forward_autosomal, rownames_separated_reverse_autosomal)
  }
  return (nipt_matrix)
}
getsamplenames <- function(nipt_sample){
  nipt_sample$name
}

getcorrectionstatus <- function(nipt_sample){
  nipt_sample$correction_status_autosomal_chromosomes
}

getstrandtype <- function(nipt_sample){
  class(nipt_sample)[2]
}
getcontrolchromosomes.SeparatedStrands <- function(nipt_sample){
 c(paste0(control_chromosomes, "F"), paste0(control_chromosomes, "R"))
}
getcontrolchromosomes.CombinedStrands <- function(nipt_sample){
  control_chromosomes
}
splitchromosomes <- function(read_counts, chromosomes){
  read_counts[chromosomes,]
}
appendchromosomes <- function(autosomal_reads, sex_reads){
  rbind(autosomal_reads, sex_reads)
}
getfractionscontrolgroup.SeparatedStrands <- function(nipt_control_group){
  fraction_table <- sapply(nipt_control_group$Samples, chrfractions)
  colnames(fraction_table) <- sapply(nipt_control_group$Samples, getsamplenames)
  rownames(fraction_table) <- c(rownames_separated_forward_autosomal, rownames_separated_reverse_autosomal)
  return(fraction_table)
}

getfractionscontrolgroup.CombinedStrands <- function(nipt_control_group){
  fraction_table <- sapply(nipt_control_group$Samples, chrfractions)
  colnames(fraction_table) <- sapply(nipt_control_group$Samples, getsamplenames)
  rownames(fraction_table) <- as.character(autosomal_chromosomes)
  return(fraction_table)
}

getreadscontrolgroup.SeparatedStrands <- function(nipt_control_group){
  reads_table <- sapply(nipt_control_group$Samples, chrreads)
  colnames(reads_table) <- sapply(nipt_control_group$Samples, getsamplenames)
  rownames(reads_table) <- c(rownames_separated_forward_autosomal, rownames_separated_reverse_autosomal)
  return(reads_table)
}

getreadscontrolgroup.CombinedStrands <- function(nipt_control_group){
  reads_table <- sapply(nipt_control_group$Samples, chrreads)
  colnames(reads_table) <- sapply(nipt_control_group$Samples, getsamplenames)
  rownames(reads_table) <- as.character(autosomal_chromosomes)
  return(reads_table)
}
sumfandrautosomal <- function(nipt_sample){
  summed_reads <- Reduce("+", nipt_sample$autosomal_chromosome_reads)
  rownames(summed_reads) <- autosomal_chromosomes
  return (summed_reads)
}
sumfandrsex <- function(nipt_sample){
  summed_reads <- Reduce("+", nipt_sample$sex_chromosome_reads)
  rownames(summed_reads) <- XY
  return (summed_reads)
}