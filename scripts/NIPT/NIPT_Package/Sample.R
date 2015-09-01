construct_sample <- function(autosomal_reads, sex_reads, name, correction_status_autosomal, correction_status_sex)
{
  if(length(autosomal_reads) == 1){
    rownames(autosomal_reads[[1]]) <- as.character(autosomal_chromosomes)
    rownames(sex_reads[[1]]) <- XY
    type = CombinedStrands
  }
  if(length(autosomal_reads) == 2){
    rownames(autosomal_reads[[1]]) <- rownames_separated_forward_autosomal
    rownames(autosomal_reads[[2]]) <- rownames_separated_reverse_autosomal
    
    rownames(sex_reads[[1]]) <- rownames_separated_forward_sex
    rownames(sex_reads[[2]]) <- rownames_separated_reverse_sex
    type = SeparatedStrands
  }
  if ((length(correction_status_autosomal) > 1) & Uncorrected %in% correction_status_autosomal){
    correction_status_autosomal <- correction_status_autosomal[-grep(pattern = Uncorrected, 
                                                                     x = correction_status_autosomal)]
  }
  
  if ((length(correction_status_sex) > 1) & Uncorrected %in% correction_status_sex){
    correction_status_sex <- correction_status_sex[-grep(pattern = Uncorrected, x = correction_status_sex)]
  }
  new.sample <- list(autosomal_chromosome_reads = autosomal_reads, 
                     correction_status_autosomal_chromosomes = correction_status_autosomal,
                     sex_chromosome_reads = sex_reads, correction_status_sex_chromosome = correction_status_sex,
                     name = name)
  class(new.sample) <- c(NIPT_sample_class, type)
  return(new.sample)
}