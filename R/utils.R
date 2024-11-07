estimate_models <- function(df, i, treat, inst,
                            x = NULL,
                            weights = NULL,
                            formula.two.sls.y,
                            formula.two.sls.m,
                            formula.first.stage,
                            formula.two.sls.scaling = NULL,
                            formula.residuals,
                            formula.var,
                            parallel.boot = "no"){
  
  df <- df[i,]
  
  two.sls.y <- AER::ivreg(formula.two.sls.y, data = df,
                          weights = weights)
  
  two.sls.m <- AER::ivreg(formula.two.sls.m,
                          data = df,
                          weights = weights)
  
  first.stage <- lm(formula.first.stage,
                    data = df,
                    weights = weights)
  
  if(!is.null(formula.two.sls.scaling) & length(formula.two.sls.scaling) > 1){
    scaling.factors <- rep(NA, length(formula.two.sls.scaling))
    for(iv in 1:length(formula.two.sls.scaling)){
      scaling.factors[iv] <- 
        coef(AER::ivreg(formula.two.sls.scaling[[iv]],
                        data = df, weights = weights))[treat]
    }
  }
  
  # standardize controls:
  if(!is.null(x)){
    df[, x] <- apply(df[,x], 2, function(x) arm:::rescale(x, binary.inputs = "0/1"))
  }
  
  residuals.m <-  lm(formula.residuals,
                     data = df,
                     weights = weights)$residuals
  
  df$residuals.m <- residuals.m^2
  
  var.reg <- lm(formula.var,
                data = df,
                weights = weights)
  
  
  variance.effect <- coef(var.reg)[inst]

  # note that by construction of reg formula, first element is always intercept
  # and next items are instruments:
  variance.terms <- 2 * coef(var.reg)[1] + coef(var.reg)[2:(length(inst)+1)]
  
  # note that in the formula for bounds, it's
  # var(M|Z = 1) *plus* var(M|Z = 0), not minus
  
  variance.terms[variance.terms < 0] <- 0
  # since we do not take sqrt of negative numbers
  
  # update progress bar:
  if(parallel.boot == "no"){
    setTxtProgressBar(pb, rep_count)
    rep_count <<- rep_count + 1
    Sys.sleep(0.001)
  }
   
  if(!is.null(formula.two.sls.scaling) & length(formula.two.sls.scaling) > 1){
    return(
      c(
        coef(two.sls.y)[treat],
        coef(two.sls.m)[treat],
        variance.effect,
        variance.terms,
        coef(first.stage)[inst],
        scaling.factors
      )
    )
  }
  
  return(
    c(
      coef(two.sls.y)[treat],
      coef(two.sls.m)[treat],
      variance.effect,
      variance.terms,
      coef(first.stage)[inst]
    )
  )
}

compute_bounds <- function(sigma = sigma, gamma = gamma,
                           estimates = estimates){

  gamma.sigma <- expand.grid(gamma = gamma, sigma = sigma)

  # if just vector of estimates:
  if(is.null(dim(estimates))){
    out <- matrix(unlist(c(estimates["naive.estimate"] - estimates["estimate.m"]*gamma.sigma["gamma"] -
                             (sqrt(estimates["variance.term"])*gamma.sigma["sigma"])/estimates[5],
                           estimates["naive.estimate"] - estimates["estimate.m"]*gamma.sigma["gamma"] +
                             (sqrt(estimates["variance.term"])*gamma.sigma["sigma"])/estimates[5])),
                  byrow = F, nrow = 2*nrow(gamma.sigma)
    )
  } else {
    # creates (2*length(gamma)*length(sigma), nrep) matrix
    # (two bounds for each combination of gamma and sigma,
    # for each bootstrap replication)
    # ordered as param1.bound1, param2.bound1, ..., param1.bound2...:
    out <- matrix(
      unlist(
        apply(estimates, 1, function(x){
          c(x["naive.estimate"] - x["estimate.m"]*gamma.sigma["gamma"] -
              (sqrt(x["variance.term"])*gamma.sigma["sigma"])/x[5],
            x["naive.estimate"] - x["estimate.m"]*gamma.sigma["gamma"] +
              (sqrt(x["variance.term"])*gamma.sigma["sigma"])/x[5])
        })
      ), byrow = F, nrow = 2*nrow(gamma.sigma)
    )
  }
  return(out)
}
