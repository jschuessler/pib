estimate_models <- function(df, i, treat, inst,
                            x = NULL,
                            weights = NULL,
                            formula.two.sls.y,
                            formula.two.sls.m,
                            formula.first.stage,
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
  variance.term <- sum(
    coef(var.reg)[c("(Intercept)", inst)]%*%c(2, rep(1, length(inst)))
  )
  # note that in the formula for bounds, it's
  # var(M|Z = 1) *plus* var(M|Z = 0), not minus
  
  if(variance.term < 0){variance.term <- runif(1, min = 0, max = 0.05)}
  # since we do not take sqrt of negative numbers
  
  names(variance.term) <- "variance.term"
  
  # update progress bar:
  if(parallel.boot == "no"){
    setTxtProgressBar(pb, rep_count)
    rep_count <<- rep_count + 1
    Sys.sleep(0.001)
  }
  
  
  return(
    c(
      coef(two.sls.y)[treat],
      coef(two.sls.m)[treat],
      variance.effect,
      variance.term,
      coef(first.stage)[inst]
    )
  )
}

compute_bounds <- function(sigma = sigma, gamma = gamma,
                           estimates = estimates){

  gamma.sigma <- expand.grid(gamma = gamma, sigma = sigma)

  # if just vector of estimates:
  if(is.null(dim(estimates))){
    out <- matrix(unlist(c(estimates[1] - estimates[2]*gamma.sigma["gamma"] -
                             (sqrt(estimates[4])*gamma.sigma["sigma"])/estimates[5],
                           estimates[1] - estimates[2]*gamma.sigma["gamma"] +
                             (sqrt(estimates[4])*gamma.sigma["sigma"])/estimates[5])),
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
          c(x[1] - x[2]*gamma.sigma["gamma"] -
              (sqrt(x[4])*gamma.sigma["sigma"])/x[5],
            x[1] - x[2]*gamma.sigma["gamma"] +
              (sqrt(x[4])*gamma.sigma["sigma"])/x[5])
        })
      ), byrow = F, nrow = 2*nrow(gamma.sigma)
    )
  }
  return(out)
}
