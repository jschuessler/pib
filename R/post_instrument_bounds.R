#' Compute bounds from estimated components.
#'
#' @param gamma Integer or vector of integers giving lowest and highest
#'                    value for mean causal effect of post-instrument variable `M` on outcome `Y`
#' @param sigma Integer or vector of integers giving lowest and highest
#'                    value for variance of causal effect of M on Y. Must
#'                    be non-negative. Defaults to zero.
#'
#' @param n.sigma Integer: Number of distinct values of `sigma` (within bounds provided by
#'                  `sigma` argument) for which to compute bounds. Defaults to 5.
#'
#' @param df  Data frame
#' @param outcome Name of dependent variables of interest `Y`, a string
#' @param post.inst Name of post-instrument variable `M`, a string
#' @param estimates  An object from calling [estimate_sensitivity()] that
#'                    contains a vector of point estimates `t` as
#'                    well as bootstrap replicates `t0`
#' @param whichCI Determines the type of boostrap confidence interval:
#'                    `percentile` (the default), `basic`, or `both`
#' @param remove.outliers Logical. If `TRUE` (the default), outliers from bootstrapped
#'                    sampling distribution of the 2SLS regressions for `Y` and `M` are removed.
#'                    The resulting sampling distribution should better
#'                    approximate results from standard robust standard errors.
#' @param conf The confidence level of the confidence intervals. Defaults
#'                    to `0.95`
#' @details
#' The function assumes that the instrument `Z`, treatment `D`, and post-instrument
#' variable `M` are binary (in which case there is only the sensitivity parameter
#' `gamma` and `sigma` can be left at its default value zero), or that
#' `z` is binary and `D` and `M` continuous, in which case both `gamma` and
#' `sigma` are sensitivity parameters.
#'

#' @return A matrix with values of `gamma`, `sigma`, point estimates or bounds
#' of treatment effect, and confidence intervals for those point estimates
#' or bounds.
#' @references Schuessler, J., Glynn, A. N., & Rueda, M.R. (2023). Post-Instrument Bias.
#'
#' Glynn, A. N., Rueda, M. R., & Schuessler, J. (2023). Post-Instrument Bias in Linear Models. Sociological Methods & Research.
#' @examples
#' \dontrun{
#' gamma <- c(-0.02, 0.02)
#' sigma <- c(0, 0.002)
#' bounds <- post.instrument.bounds(gamma = gamma, sigma = sigma,
#' n.sigma = 2, estimates = out)}
#' @import boot
#' @export

post_instrument_bounds <- function(gamma, sigma = 0,
                                   n.sigma = 5,
                                   estimates, df,
                                   outcome, post.inst,
                                   whichCI = "percentile",
                                   remove.outliers = T,
                                   conf = 0.95){

  if(length(gamma) > 2){
    stop("gamma must be integer or vector of length 2!")
  }

  if(any(sigma < 0)){
    stop("sigma must be non-negative!")
  }else   if(length(sigma) > 2){
    stop("sigma must be integer or vector of length 2!")
  }

  gamma <- seq(gamma[1], gamma[2], length.out = 100)

  if(sum(sigma > 0)){sigma <- seq(sigma[1], sigma[2], length.out = n.sigma)}

  gamma.sigma <- expand.grid(gamma = gamma, sigma = sigma)
  n.combinations <- nrow(gamma.sigma)

  # Point estimates
  # this is 2*n.combinations x 1 (factor 2 because of lower/upper bound)
  point.estimates <- compute_bounds(gamma = gamma, sigma = sigma,
                                    estimates = estimates$t0)

  # remove bootstrap replications that
  # are outliers for at least one of the 2SLS regressions

  # create indicator for replication being outlier:
  estimates$t <- cbind(estimates$t, 0)

  if(remove.outliers){
    # flag outliers:
    estimates$t[,6][which(estimates$t[,1] %in% boxplot.stats(
      estimates$t[,1], coef = 1.5, do.conf = F)$out)] <- 1
    estimates$t[,6][which(estimates$t[,2] %in% boxplot.stats(
      estimates$t[,2], coef = 1.5, do.conf = F)$out)] <- 1
  }

  # if sigma is just scalar 0, can delete lower half
  # because lower = upper bound:
  if(sum(sigma) == 0){
    point.estimates <- matrix(point.estimates[1:n.combinations,])
    colnames(point.estimates) <- "Point Estimate"
  }

  # Confidence intervals
  # this is 2*n.combinations x R (minus outliers):
  out <- compute_bounds(gamma = gamma, sigma = sigma,
                        estimates = estimates$t[estimates$t[,6] == 0,])

  # If sigma is scalar 0, can delete lower half of CIs:
  if(sum(sigma) == 0) out <- out[1:n.combinations,]

  CI <- matrix(NA, nrow = nrow(out), ncol = 4)

  # loop over all parameter*bound combinations
  # gives 2*ncombinations x 4 matrix (2 kinds of CIs as columns)
  for(i in 1:nrow(out)){

    df <- out[i,]

    CI[i,1:2] <- boot:::basic.ci(t0 = point.estimates[i,],
                                 t = df, conf = conf)[, 4:5]
    CI[i,3:4] <- boot:::perc.ci(t = df, conf = conf)[, 4:5]
  }

  # take union of lower and upper bound CI:

  # if sigma > 0,
  # split in the middle, rbind lower part to the right
  # gives n.combinations x 8 matrix
  # then take union of lower bound CI and upper bound CI
  # yields n.combinations x 4 matrix
  if(sum(sigma) > 0) {
    CI <- cbind(CI[1:n.combinations,],
                CI[(n.combinations+1):(n.combinations*2),])

    # cols 1, 2, 5, 6, are baisc.ci, 3, 4, 7, 8, are perc.ci
    # this yields a n.combinations matrix x 4
    # cols 1 and 2 are basic CI, cols 3 and 4 are percentile CI:
    CI <- t(apply(CI, 1, function(x){
      c(
        min(x[c(1, 2, 5, 6)]), max(x[c(1, 2, 5, 6)]),
        min(x[c(3, 4, 7, 8)]), max(x[c(3, 4, 7, 8)])
      )
    }
    ))
  }

  colnames(CI) <- c("Lower Bound CI (Basic)", "Upper Bound CI (Basic)",
                    "Lower Bound CI (Percentile)", "Upper Bound CI (Percentile)")

  if(whichCI == "basic") CI <- CI[,1:2]
  if(whichCI == "percentile") CI <- CI[,3:4]

  # if sigma is greater zero, output will be bounds for point estimate
  # stack as columns
  # this creates a n.combinations x 2 matrix (lower/upper bound):
  if(sum(sigma) > 0){
    point.estimates <- cbind(point.estimates[1:n.combinations,],
                             point.estimates[(n.combinations+1):(n.combinations*2),])
    colnames(point.estimates) <- c("Lower Point Estimate",  "Upper Point Estimate")
  }

  final <- cbind(gamma.sigma, point.estimates, CI)

  return(final)

}


