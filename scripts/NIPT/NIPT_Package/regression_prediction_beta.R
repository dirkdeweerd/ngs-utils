PredictTrisomy <- function(predictors, nipt_sample, chromo_focus,  control_group)
{
  if(is.null(attr(predictors, "class"))){
    print("Object is not of type Sample")
  }
  else  UseMethod("PredictTrisomy", predictors)
}

PredictTrisomy.SeparatedStrands <- function(predictors, nipt_sample,  chromo_focus,  control_group){
  chromosomal.frac.control <- sapply(X = control_group, FUN = chrfractions)
  chromosomal_frac_sample <- chrfractions(nipt_sample = nipt_sample)
  frac.reads.chr.trisomy.observed = chromosomal.frac.control[chromo_focus,] + chromosomal.frac.control[chromo_focus + 22, ]
  vc.prac <- 1.15 * ( 1/ sqrt(sum(unlist(nipt_sample$reads))))
  mod <- BuildFullModel(as.data.frame(t(chromosomal.frac.control))[predictors], chr.of.interest.fractions =  frac.reads.chr.trisomy.observed)
  ratio <-  frac.reads.chr.trisomy.observed /  mod$fitted.values
  vc.theo <- sd(ratio) 
  vc <- max(c(vc.theo, vc.prac))
  sample.predicted <- predict(mod, as.data.frame(t(chromosomal_frac_sample))[predictors])
  sample.ratio <- sum(chromosomal_frac_sample[c(chromo_focus,(chromo_focus+22))]) / sample.predicted  
  z.score.sample <- (sample.ratio - 1) / vc
  z.score.controls <- matrix((data = as.numeric(ratio) - 1) / vc, ncol = 1, dimnames = list(sapply(X = control_group, FUN = samplenames), NULL))
  print(z.score.sample)
  shapiro.regression <- shapiro.test(z.score.controls)
  return(list(sample_Z_score = z.score.sample, control_group_Z_scores = z.score.controls, shapiro_P_value = shapiro.regression$p.value,
              predictors = predictors))
}

PredictTrisomy.CombinedStrands <- function(predictors, nipt_sample,  chromo_focus,  control_group){
  chromosomal.frac.control <- sapply(X = control_group, FUN = chrfractions)
  chromosomal_frac_sample <- chrfractions(nipt_sample.sample = nipt_sample)
  frac.reads.chr.trisomy.observed = chromosomal.frac.control[chromo_focus,] 
  vc.prac <- 1.15 * ( 1/ sqrt(sum(unlist(nipt_sample$reads))))
  mod <- BuildFullModel(as.data.frame(t(chromosomal.frac.control))[predictors], chr.of.interest.fractions =  frac.reads.chr.trisomy.observed)
  ratio <-  frac.reads.chr.trisomy.observed /  mod$fitted.values
  vc.theo <- sd(ratio) 
  vc <- max(c(vc.theo, vc.prac))
  sample.predicted <- predict(mod, as.data.frame(t(chromosomal_frac_sample))[predictors])
  sample.ratio <- sum(chromosomal_frac_sample[chromo_focus]) / sample.predicted  
  z.score.sample <- (sample.ratio - 1) / vc
  z.score.controls <- (as.numeric(ratio) - 1) / vc
  print(z.score.sample)
  shapiro.regression <- shapiro.test(z.score.controls)
  return(z.score.controls)
}