#' Estimate components for sensitivity analysis.
#'
#' @param df         Data frame
#' @param R          Number of bootstrap replications
#' @param outcome    Name of dependent variables of interest `Y`, a string
#' @param treat      Name of treatment of interest `D`, a string
#' @param inst       Name of instrument `Z`, a string
#' @param post.inst  Name of post-instrument variable `M`, a string
#' @param x          Names of controls, a string or vector of strings
#' @param fe         Names of variables to be passed as fixed effects (or factors),
#'             a string or vector of strings
#' @param weights    Vector of weights to be used in all regressions
#' @param seed       Seed to make bootstrap analysis replicable. Must be set to a number
#' @param parallel   One of `c("no", "snow", "multicore")` to parallelize the bootstrap. Default is
#'             `"no"`.
#' @param ncpus      Integer: number of processes to be used in parallel operation
#' @param cl         An optional parallel or snow cluster for use if `parallel = "snow"`.
#'
#' @return An object of class `boot`, where `t0` contains point estimates for all elements of the bounds estimator, and `t` contains R bootstrap replicates thereof.
#' @references Schuessler, J., Glynn, A. N., & Rueda, M.R. (2023). Post-Instrument Bias.
#'
#' Glynn, A. N., Rueda, M. R., & Schuessler, J. (2023). Post-Instrument Bias in Linear Models. Sociological Methods & Research.
#' @examples
#' \dontrun{
#' data(carnegie_marinov)
#'y <- "new_empinxavg"
#'m <- "pwt_openk"
#'d <- "EV"
#'z <- "l2CPcol2"
#'fe <- c("ccode", "year")
#' out <- estimate_sensitivity(carnegie_marinov, R = 1501,
#'                      outcome = y, treat = d, inst = z, post.inst = m,
#'                      fe = fe,
#'                      seed = 823,
#'                      parallel = "snow",
#'                      ncpus = 4)
#'                       }
#' @import boot
#' @export
estimate_sensitivity <- function(df, R = 1501,
                                 outcome, treat, inst, post.inst,
                                 x = NULL, fe = NULL,
                                 weights = NULL,
                                 seed,
                                 parallel = "no",
                                 ncpus = getOption("boot.ncpus", 1L), cl = NULL){

  # take only complete cases:

  if(!is.null(x) & !is.null(fe)){
    df <- df[,c(outcome, treat, inst, post.inst, fe, x)]
  } else if(!is.null(x)){
    df <- df[,c(outcome, treat, inst, post.inst, x)]
  } else if(!is.null(fe)){
    df <- df[,c(outcome, treat, inst, post.inst, fe)]
  } else df <- df[,c(outcome, treat, inst, post.inst)]

  df <- df[complete.cases(df),]

  # if no controls / no FEs, convert into constants, which are effectively ignored
  # in regressions:
  if(is.null(x)) x.reg <- 1 else x.reg <- x
  if(is.null(fe)) fe.reg <- 1 else fe.reg <- paste("factor(", fe, ")", sep = "")

  # initiate progress bar:
  if(parallel == "no"){
    rep_count <<- 1

    pb <<- txtProgressBar(min = 0, max = R, initial = 0,
                          style = 3, file = "")
  }


  # create regression formulae:
  formula.two.sls.y <- as.formula(
    paste(paste(outcome, paste(c(treat, x.reg, fe.reg), collapse = " + "), sep = " ~ "),
          paste(c(inst, x.reg, fe.reg), collapse = " + "), sep = " | "))

  formula.two.sls.m <- as.formula(
    paste(paste(post.inst, paste(c(treat, x.reg, fe.reg), collapse = " + "), sep = " ~ "),
          paste(c(inst, x.reg, fe.reg), collapse = " + "), sep = " | "))

  formula.first.stage <- as.formula(
    paste(paste(d, paste(c(inst, x.reg, fe.reg), collapse = " + "), sep = " ~ ")))

  formula.residuals <- as.formula(
    paste(paste(post.inst, paste(c(inst, x.reg, fe.reg), collapse = " + "), sep = " ~ ")))

  formula.var <- as.formula(
    paste(paste("residuals.m", paste(c(inst, x.reg, fe.reg), collapse = " + "), sep = " ~ ")))


  set.seed(seed)
  bootstrap0 <- boot(data = df,
                     statistic = estimate_models,
                     R = R,
                     treat = treat,
                     inst = inst,
                     x = x,
                     weights = weights,
                     formula.two.sls.y = formula.two.sls.y,
                     formula.two.sls.m = formula.two.sls.m,
                     formula.first.stage = formula.first.stage,
                     formula.residuals = formula.residuals,
                     formula.var = formula.var,
                     parallel.boot = parallel,
                     parallel = parallel,
                     ncpus = ncpus,
                     cl = cl)

  if(parallel == "no") close(pb)

  return(bootstrap0)
}
