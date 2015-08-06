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
GetPredictors <- function(control_group){
  chromosomal.frac.control <- sapply(X = control_group, FUN = chrfractions)
  frac.reads.chr.trisomy.observed = as.matrix(chromosomal.frac.control[13,]) + as.matrix(chromosomal.frac.control[35,])
  chr.potential.trisomic <- c(paste0(chromosomes.trisomy, "F"),paste0(chromosomes.trisomy, "R"))
  predictor.list.regression.fixed <- SelectModelsFixedRegressionApproach(chromosomal.frac.control = chromosomal.frac.control,
                                                                         chr.potential.trisomic=chr.potential.trisomic,
                                                                         frac.reads.chr.trisomy.observed=frac.reads.chr.trisomy.observed, 
                                                                         predictor.max = fixed.predictors )
  return (predictor.list.regression.fixed)
}
PredictTrisomy <- function(predictors, chr.int,  control.samples, sample.bins.f,sample.bins.r, sample.of.interest){
  chromosomal.frac.control <- sapply(X = control_group, FUN = chrfractions)
  frac.reads.chr.trisomy.observed = as.matrix(chromosomal.frac.control[,chr.int]) + as.matrix(chromosomal.frac.control[,(chr.int + 22)])
  vc.prac <- 1.15 * (1 / sqrt(sum(sample.bins.f[chr.int,]) + sum(sample.bins.r[chr.int,])))
  mod <- BuildFullModel(control.samples[predictors], chr.of.interest.fractions =  frac.reads.chr.trisomy.observed)
  ratio <-  frac.reads.chr.trisomy.observed /  mod$fitted.values
  vc.theo <- sd(ratio) 
  vc <- max(c(vc.theo, vc.prac))
  sample.predicted <- predict(mod, sample.of.interest[predictors])
  sample.ratio <- sum(sample.of.interest[c(chr.int,(chr.int+22))]) / sample.predicted  
  z.score.sample <- (sample.ratio - 1) / vc
  z.score.controls <- (as.numeric(ratio) - 1) / vc
  shapiro.regression <- shapiro.test(z.score.controls)
}
BuildFullModel <- function(control.samples, chr.of.interest.fractions){
  return (lm(chr.of.interest.fractions ~ ., data=control.samples))
}
SelectModelsFixedRegressionApproach <- function(chromosomal.frac.control, chr.potential.trisomic, frac.reads.chr.trisomy.observed, predictor.max){
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
    predictor.list[[model]] <- predictors
  }
  return(predictor.list)
}