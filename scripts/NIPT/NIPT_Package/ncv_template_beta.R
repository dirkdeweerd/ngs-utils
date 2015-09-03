RetrieveSubsets <- function(candidates, number_of_elements){
  return(set_combn(candidates, m=number_of_elements))
}
CalculateVariation <- function(denominators, chromosomal_frac_control_reads, chr_focus){
  possible_denominators <- unlist(denominators)
  if (length(possible_denominators) == 1){
    mean_subset <- mean(chromosomal_frac_control_reads[chr_focus,] / chromosomal_frac_control_reads[possible_denominators,])
    sd_subset <- sd((chromosomal_frac_control_reads[chr_focus,] / (chromosomal_frac_control_reads[possible_denominators,])))
  }
  else{
    mean_subset <- mean(chromosomal_frac_control_reads[chr_focus,] / rowSums(chromosomal_frac_control_reads[possible_denominators,]))
    sd_subset <- sd((chromosomal_frac_control_reads[chr_focus,] / (rowSums(chromosomal_frac_control_reads[possible_denominators,]))))
  }
  return(list(sd_subset , mean_subset))
}
GetDenominators <- function(nipt_control_group, chr_focus, max_elements, exclude_chromosomes = NULL, include_chromosomes = NULL,
                            use_test_train_set =T, size_of_train_set = 0.6){
  exclude_chromosomes <- as.character(exclude_chromosomes)
  include_chromosomes <- as.character(include_chromosomes)
  min_variation_subset <- list()
  min_variation_vc <- NULL
  chromosomal_frac_control_reads <- sapply(lapply(nipt_control_group$Samples, sumfandrautosomal), rowSums)
  chromosomal_frac_control_reads_train <- chromosomal_frac_control_reads
  if (use_test_train_set == T){
    indices <- sample(x = 1:length(nipt_control_group$Samples), round(size_of_train_set * length(nipt_control_group$Samples)))
    nipt_control_group_train <- as_control_group(nipt_samples = nipt_control_group$Samples[indices])
    nipt_control_group <- as_control_group(nipt_samples = nipt_control_group$Samples[-indices])
    chromosomal_frac_control_reads_train <- sapply(lapply(nipt_control_group_train$Samples, sumfandrautosomal), rowSums)
    rownames(chromosomal_frac_control_reads) <- as.character(autosomal_chromosomes)
    chromosomal_frac_control_reads <- sapply(lapply(nipt_control_group$Samples, sumfandrautosomal), rowSums)
  }
  rownames(chromosomal_frac_control_reads_train) <- as.character(autosomal_chromosomes)
  rownames(chromosomal_frac_control_reads) <- as.character(autosomal_chromosomes)
  control_chromosomes <- control_chromosomes[!control_chromosomes %in% c(exclude_chromosomes, chr_focus)]
  
  if (!is.null(include_chromosomes)){
    control_chromosomes <- c(control_chromosomes, include_chromosomes)
  }
  for (n_of_elements in 1:max_elements)  {
    denominators <- as.list(RetrieveSubsets(candidates = control_chromosomes, number_of_elements = n_of_elements))
    variation_list <- lapply(FUN = CalculateVariation, X=denominators,
                             chromosomal_frac_control_reads= chromosomal_frac_control_reads_train, chr_focus=chr_focus)
    variation <- sapply(variation_list, function(x) Reduce("/", x))
    min_variation_subset[[n_of_elements]] <- as.numeric(denominators[[which.min(variation)]])
    min_variation_vc[n_of_elements] <- variation[which.min(variation)]
  }
  denominators <- min_variation_subset[[which.min(min_variation_vc)]]
  ncv_reads <- get_ncv_reads(chromosomal_frac_control_reads = chromosomal_frac_control_reads, chr_focus = chr_focus, 
                             denominators = denominators )
  scores <- ZScoresControl(denominators = denominators, ncv_reads = ncv_reads,
                           chr_focus = chr_focus, samplenames = sapply(nipt_control_group$Samples, getsamplenames))
  if (use_test_train_set == F){
    new_ncv_template <- ncv_template(denominators = denominators, chromo_focus = chr_focus, 
                                     nipt_sample_names = sapply(nipt_control_group$Samples, getsamplenames),
                                     correction_status = unique(nipt_control_group$correction_status), 
                                     scores = scores, potential_denominators = control_chromosomes,
                                     statistics = c(mean = mean(ncv_reads), SD = sd(ncv_reads), 
                                                    Shapiro_Wilk_P_value = shapiro.test(scores$NCVscore)$p.value),
                                     type = class(nipt_control_group)[2])
  }
  else{
    ncv_reads_train <- get_ncv_reads(chromosomal_frac_control_reads = chromosomal_frac_control_reads_train, chr_focus = chr_focus, 
                                     denominators = denominators )
    scores_train <- ZScoresControl(denominators = denominators, ncv_reads = ncv_reads_train,
                                   chr_focus = chr_focus, samplenames = sapply(nipt_control_group_train$Samples, getsamplenames))
    new_ncv_template <- ncv_template(denominators = denominators, chromo_focus = chr_focus, 
                                     nipt_sample_names = sapply(nipt_control_group$Samples, getsamplenames),
                                     correction_status = unique(nipt_control_group$correction_status), 
                                     scores = scores, potential_denominators = control_chromosomes,
                                     statistics = c(mean = mean(ncv_reads), SD = sd(ncv_reads), 
                                                    Shapiro_Wilk_P_value = shapiro.test(scores$NCVscore)$p.value),
                                     type = class(nipt_control_group)[2], 
                                     sample_names_train_set = sapply(nipt_control_group_train$Samples, getsamplenames),
                                     train_set_statistics = c(mean = mean(ncv_reads_train), SD = sd(ncv_reads_train), 
                                                              Shapiro_Wilk_P_value = shapiro.test(scores_train$NCVscore)$p.value),
                                     train_set_Zscores = scores_train)  
  }
  
  
  return(new_ncv_template)
}
get_ncv_reads <- function(chromosomal_frac_control_reads, chr_focus, denominators){
  ncv_reads <- chromosomal_frac_control_reads[chr_focus, ] / chromosomal_frac_control_reads[denominators,]
}
ZScoresControl <- function(denominators, ncv_reads, chr_focus, samplenames){
  scores <- as.data.frame((ncv_reads - mean(ncv_reads)) / sd(ncv_reads), row.names = samplenames)  
  colnames(scores) <- "NCVscore"
  return(scores)
}
get_ncv_statistics <- function(ncv_reads){
  return(c(mean(ncv_reads), sd(ncv_reads)))
}
