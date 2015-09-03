SelectModelsFixedRegressionApproach <- function(nipt_sample, chromosomal_frac_control, chromosomes_trisomy, 
                                                frac_reads_chr_trisomy_observed, predictor.max){
  if(is.null(attr(nipt_sample, "class"))){
    print(error_message)
  }
  else UseMethod("SelectModelsFixedRegressionApproach", nipt_sample)
}

GetNextPredictor <- function(samples, frac_reads_chr_trisomy_observed, predictors, chromosomal_frac_control){
  adj.r.squares <- NULL
  chr.candidates <- rownames(samples)
  for (i in 1:length(rownames(samples)))  {
    current.pred.set <- c(predictors, chr.candidates[i])
    chr.frac.df <- t(as.matrix(rbind(chromosomal_frac_control[current.pred.set, ])))
    model <- lm(frac_reads_chr_trisomy_observed ~ chr.frac.df)
    adj.r.squares[i] <- summary(model)$adj.r.squared
  }
  adj.r.squared.ordered <- order(adj.r.squares, decreasing = TRUE)
  chromosomes <- rownames(samples)
  best.pred.set <- c(predictors, chromosomes[adj.r.squared.ordered[1]])
  best.chr.frac.df <- t(as.matrix(rbind(chromosomal_frac_control[best.pred.set,])))
  best.model <- lm(frac_reads_chr_trisomy_observed ~ best.chr.frac.df)
  predictor.adj.r.squared <- list(chromosomes[adj.r.squared.ordered[1]], summary(best.model)$adj.r.squared)
  return(predictor.adj.r.squared)
}
perform_regression <- function(nipt_sample, control_group, chromo_focus, exclude_chromosomes = NULL, include_chromosomes = NULL,
                               use_test_train_set =T, size_of_train_set = 0.6){
  if (use_test_train_set == T){
    indices <- sample(x = 1:length(control_group$Samples), round(size_of_train_set * length(control_group$Samples)))
    control_group_train <- as_control_group(nipt_samples = control_group$Samples[indices])
    control_group <- as_control_group(nipt_samples = control_group$Samples[-indices])
  }
  chromosomal_frac_control_train <<- sapply(X = control_group_train$Samples, FUN = chrfractions)
  chromosomal_frac_control_train <- setrownamesmatrix(chromosomal_frac_control_train)
  chromosomal_frac_control_test <<- sapply(X = control_group$Samples, FUN = chrfractions)
  chromosomal_frac_control_test <- setrownamesmatrix(chromosomal_frac_control_test)
  if (!is.null(exclude_chromosomes)){
    chromosomes_trisomy <- c(chromosomes_trisomy, exclude_chromosomes)
  }
  if (!is.null(include_chromosomes)){
    chromosomes_trisomy <- chromosomes_trisomy[!chromosomes_trisomy %in% include_chromosomes]
  }
  frac_reads_chr_trisomy_observed <<- retrieve_fractions_of_interest(nipt_sample = nipt_sample, chromo_focus = chromo_focus, 
                                                                    chromosomal_fracs = chromosomal_frac_control_train)
  frac_reads_chr_trisomy_observed_test <<- retrieve_fractions_of_interest(nipt_sample = nipt_sample, chromo_focus = chromo_focus, 
                                                                     chromosomal_fracs = chromosomal_frac_control_test)
  predictor.list <- SelectModelsFixedRegressionApproach(nipt_sample = nipt_sample, chromosomal_frac_control= chromosomal_frac_control_train,
                                                        chromosomes_trisomy = chromosomes_trisomy, 
                                                        frac_reads_chr_trisomy_observed = frac_reads_chr_trisomy_observed, 
                                                        predictor.max = fixed.predictors )
  prediction <- lapply(X = predictor.list, FUN = PredictTrisomy, nipt_sample = nipt_sample, chromo_focus = chromo_focus, control_group = control_group, 
                       frac_reads_chr_trisomy_observed = frac_reads_chr_trisomy_observed_test)
  result_nipt <- collapse_results(result_set = prediction)
  return (result_nipt)
}

BuildFullModel <- function(control.samples, chr.of.interest.fractions){
  papis <<- control.samples
  tonijn <<- chr.of.interest.fractions
  return (lm(chr.of.interest.fractions ~ ., data=control.samples))
}
SelectModelsFixedRegressionApproach.SeparatedStrands <- function(nipt_sample, chromosomal_frac_control, chromosomes_trisomy, 
                                                                 frac_reads_chr_trisomy_observed, predictor.max){
  predictor.list <- list()
  chr.potential.trisomic <- c(paste0(chromosomes_trisomy, "F"),paste0(chromosomes_trisomy, "R"))
  for (model in 1:n.models)  {
    predictors <- NULL
    predictor.adj.r.squares <-NULL
    predictors.complementary <- list()
    adj.r.squares <- NULL
    for (predictor in 1:n.predictors){
      potential.predictors <- chromosomal_frac_control[rownames(chromosomal_frac_control)[!rownames(chromosomal_frac_control) %in% c(chr.potential.trisomic,  unlist(predictors.complementary), unlist(predictor.list))],]
      
      predictor.adj.r.squares <- GetNextPredictor(samples = potential.predictors, frac_reads_chr_trisomy_observed=frac_reads_chr_trisomy_observed, predictors = predictors, 
                                                  chromosomal_frac_control=chromosomal_frac_control)
      predictors[predictor] <- predictor.adj.r.squares[[1]]
      predictors.complementary[[predictor]] <- unique(c(predictor.adj.r.squares[[1]], sub("F", "R", predictor.adj.r.squares[[1]]), sub("R", "F", predictor.adj.r.squares[[1]])))
      adj.r.squares[predictor] <- predictor.adj.r.squares[[2]][1]
    }
    class(predictors) <- c("Prediction Set", "SeparatedStrands")
    predictor.list[[model]] <- predictors
  }
  return(predictor.list)
}

SelectModelsFixedRegressionApproach.CombinedStrands <- function(nipt_sample, chromosomal_frac_control, chromosomes_trisomy,
                                                                frac_reads_chr_trisomy_observed, predictor.max){
  paaphoofd <<- chromosomal_frac_control
  predictor.list <- list()
  chr.potential.trisomic <- chromosomes_trisomy
  for (model in 1:n.models){
    predictors <- NULL
    predictor.adj.r.squares <-NULL
    adj.r.squares <- NULL
    for (predictor in 1:n.predictors){
      potential.predictors <- chromosomal_frac_control[rownames(chromosomal_frac_control)[!rownames(chromosomal_frac_control) %in% c(chr.potential.trisomic, unlist(predictor.list), predictors)],]
      predictor.adj.r.squares <- GetNextPredictor(samples = potential.predictors, frac_reads_chr_trisomy_observed=frac_reads_chr_trisomy_observed, predictors = predictors, 
                                                  chromosomal_frac_control=chromosomal_frac_control)
      predictors[predictor] <- predictor.adj.r.squares[[1]]
      
      adj.r.squares[predictor] <- predictor.adj.r.squares[[2]][1]
    }
    class(predictors) <- c("Prediction Set", "CombinedStrands")
    predictor.list[[model]] <- predictors
  }
  return(predictor.list)
}

PredictTrisomy <- function(predictors, nipt_sample,  chromo_focus,  control_group, frac_reads_chr_trisomy_observed){
  chromosomal_frac_control <<- sapply(X = control_group$Samples, FUN = chrfractions)
  chromosomal_frac_control <<- setrownamesmatrix(chromosomal_frac_control)
  chromosomal_frac_sample <<- chrfractions(nipt_sample = nipt_sample)
  samplereads <<- sumfandrautosomal(nipt_sample)
  vc.prac <<- 1.15 * ( 1/ sqrt(sum(samplereads[chromo_focus,])))
  print(vc.prac)
  mod <- BuildFullModel(as.data.frame(t(chromosomal_frac_control))[predictors], chr.of.interest.fractions =  frac_reads_chr_trisomy_observed)
  ratio <-  frac_reads_chr_trisomy_observed /  mod$fitted.values
  print(ratio)
  vc.theo <- sd(ratio) 
  vc <- max(c(vc.theo, vc.prac))
  sample.predicted <- predict(mod, as.data.frame(t(chromosomal_frac_sample))[predictors])
 
  sample.ratio <- retrieve_fractions_of_interest(nipt_sample = nipt_sample, chromo_focus = chromo_focus, 
                                                 chromosomal_fracs = as.matrix(chromosomal_frac_sample)) / sample.predicted  
  print(sample.ratio)

  z.score.sample <- (sample.ratio - 1) / vc
  z.score.controls <- matrix((data = as.numeric(ratio) - 1) / vc, ncol = 1, dimnames = list(sapply(X = control_group$Samples, FUN = getsamplenames), NULL))
  shapiro.regression <- shapiro.test(z.score.controls)
  
  return(list(sample_Z_score = z.score.sample, control_group_Z_scores = z.score.controls, shapiro_P_value = shapiro.regression$p.value,
              predictors = predictors))
}