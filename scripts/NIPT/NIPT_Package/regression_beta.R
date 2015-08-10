perform_regression <- function(nipt_sample, control_group, chromo_focus)
{
  if(is.null(attr(nipt_sample, "class"))){
    print("Object is not of type Sample")
  }
  else  UseMethod("perform_regression", nipt_sample)
}

GetNextPredictor <- function(samples, frac.reads.chr.trisomy.observed, predictors, chromosomal.frac.control){
  adj.r.squares <- NULL
  chr.candidates <- rownames(samples)
  for (i in 1:length(rownames(samples)))  {
    current.pred.set <- c(predictors, chr.candidates[i])
    chr.frac.df <- t(as.matrix(rbind(chromosomal.frac.control[current.pred.set, ])))
    model <- lm(frac.reads.chr.trisomy.observed ~ chr.frac.df)
    adj.r.squares[i] <- summary(model)$adj.r.squared
  }
  adj.r.squared.ordered <- order(adj.r.squares, decreasing = TRUE)
  chromosomes <- rownames(samples)
  best.pred.set <- c(predictors, chromosomes[adj.r.squared.ordered[1]])
  best.chr.frac.df <- t(as.matrix(rbind(chromosomal.frac.control[best.pred.set,])))
  best.model <- lm(frac.reads.chr.trisomy.observed ~ best.chr.frac.df)
  predictor.adj.r.squared <- list(chromosomes[adj.r.squared.ordered[1]], summary(best.model)$adj.r.squared)
  return(predictor.adj.r.squared)
}
perform_regression.SeparatedStrands <- function(nipt_sample, control_group, chromo_focus){
  chromosomal.frac.control <- sapply(X = control_group, FUN = chrfractions)
  frac.reads.chr.trisomy.observed = as.matrix(chromosomal.frac.control[chromo_focus,]) + as.matrix(chromosomal.frac.control[chromo_focus,])
  chr.potential.trisomic <- c(paste0(chromosomes.trisomy, "F"),paste0(chromosomes.trisomy, "R"))
  predictor.list <- SelectModelsFixedRegressionApproach.SeparatedStrands(chromosomal.frac.control = chromosomal.frac.control,
                                                                         chr.potential.trisomic=chr.potential.trisomic,
                                                                         frac.reads.chr.trisomy.observed=frac.reads.chr.trisomy.observed, 
                                                                         predictor.max = fixed.predictors )
  prediction <- lapply(X = predictor.list, FUN = PredictTrisomy, nipt_sample = nipt_sample, chromo_focus = 13, control_group = control_group)
  return (prediction)
}

perform_regression.CombinedStrands <- function(nipt_sample, control_group, chromo_focus){
  chromosomal_frac_control <- sapply(X = control_group, FUN = chrfractions)
  frac_reads.chr_trisomy_observed = chromosomal.frac.control[chromo_focus,]
  chr_potential_trisomic <- chromosomes.trisomy
  predictor.list.regression.fixed <- SelectModelsFixedRegressionApproach.CombinedStrands(chromosomal.frac.control = chromosomal_frac_control,
                                                                         chr.potential.trisomic=chr_potential_trisomic,
                                                                         frac.reads.chr.trisomy.observed=frac.reads.chr.trisomy.observed, 
                                                                         predictor.max = fixed.predictors )
  return (predictor.list.regression.fixed)
}

BuildFullModel <- function(control.samples, chr.of.interest.fractions){
  return (lm(chr.of.interest.fractions ~ ., data=control.samples))
}
SelectModelsFixedRegressionApproach.SeparatedStrands <- function(chromosomal.frac.control, chr.potential.trisomic, frac.reads.chr.trisomy.observed, predictor.max){
  predictor.list <- list()
  for (model in 1:n.models)  {
    predictors <- NULL
    predictor.adj.r.squares <-NULL
    predictors.complementary <- list()
    adj.r.squares <- NULL
    for (predictor in 1:n.models){
      potential.predictors <- chromosomal.frac.control[rownames(chromosomal.frac.control)[!rownames(chromosomal.frac.control) %in% c(chr.potential.trisomic,  unlist(predictors.complementary), unlist(predictor.list))],]
      
      predictor.adj.r.squares <- GetNextPredictor(samples = potential.predictors, frac.reads.chr.trisomy.observed=frac.reads.chr.trisomy.observed, predictors = predictors, 
                                                  chromosomal.frac.control=chromosomal.frac.control)
      predictors[predictor] <- predictor.adj.r.squares[[1]]
      predictors.complementary[[predictor]] <- unique(c(predictor.adj.r.squares[[1]], sub("F", "R", predictor.adj.r.squares[[1]]), sub("R", "F", predictor.adj.r.squares[[1]])))
      adj.r.squares[predictor] <- predictor.adj.r.squares[[2]][1]
    }
    class(predictors) <- c("Prediction Set", "SeparatedStrands")
    print(class(predictors))
    predictor.list[[model]] <- predictors
  }
  return(predictor.list)
}

SelectModelsFixedRegressionApproach.CombinedStrands <- function(chromosomal.frac.control, chr.potential.trisomic, frac.reads.chr.trisomy.observed, predictor.max){
  predictor.list <- list()
  for (model in 1:n.models)  {
    predictors <- NULL
    predictor.adj.r.squares <-NULL
    adj.r.squares <- NULL
    for (predictor in 1:n.predictors){
      potential.predictors <- chromosomal.frac.control[rownames(chromosomal.frac.control)[!rownames(chromosomal.frac.control) %in% c(chr.potential.trisomic, unlist(predictor.list))],]
      
      predictor.adj.r.squares <- GetNextPredictor(samples = potential.predictors, frac.reads.chr.trisomy.observed=frac.reads.chr.trisomy.observed, predictors = predictors, 
                                                  chromosomal.frac.control=chromosomal.frac.control)
      predictors[predictor] <- predictor.adj.r.squares[[1]]
      adj.r.squares[predictor] <- predictor.adj.r.squares[[2]][1]
    }
    class(predictors) <- c("Prediction Set", "CombinedStrands")
    predictor.list[[model]] <- predictors
  }
  return(predictor.list)
}