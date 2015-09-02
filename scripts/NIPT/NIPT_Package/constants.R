library(data.table)
library(Rsamtools)
library(sets)

Uncorrected <- "Uncorrected"
CombinedStrands <- "CombinedStrands"
SeparatedStrands <- 'SeparatedStrands'
GCcorrected <- "GC corrected"
ChiCorrected <- "Chi square corrected"


generic_control_group <- "General control group"
XY <- c("X", "Y")

n_total_chromosomes <- 24
n_autosomal_chromosomes <- 22
autosomal_chromosomes <- 1:22
sex_chromosomes <- 23:24 
control_chromosomes <- as.character(c(1:12,14:17,19:20,22))
chromosomes_trisomy <- c(13,18,21)

chi_square_cut_off <- 3.5
n.models <- 4
n.predictors <- 4

bin.size <- 50000
n.bins <- 4985

rownames_combined_autosomal <- as.character(autosomal_chromosomes)
rownames_combined_sex <- XY
rownames_separated_forward_autosomal <- paste0(autosomal_chromosomes, "F")
rownames_separated_reverse_autosomal <- paste0(autosomal_chromosomes, "R")
rownames_separated_forward_sex <- paste0(XY, "F")
rownames_separated_reverse_sex <- paste0(XY, "R")

NIPT_sample_class <- "NIPT Sample"
NCV_template_class <- "NCV Template"
NCV_result_class <- "NCV Result"

