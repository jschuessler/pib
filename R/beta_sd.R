#' Compute standard deviation of Beta distribution.
#'
#' @param alpha,beta Integers. Need to be greater than zero. Control shape of distribution. `alpha = 4, beta = 4` gives a bell-shaped distribution
#' `alpha = 1, beta = 1` gives a Uniform distribution. `alpha = 0.5, beta = 0.5` gives
#' a high-variance distribution where very low or very high values or most likely.
#'
#' @param a,b Integers. Lower and upper bound of unit-level causal effects.
#' @details
#' Assuming that unit-level causal effects vary between `a` and `b`, this function
#' computes the implied mean effect and the implied standard deviation of the distribution of causal effects,
#' assuming a Beta distribution with parameters `alpha` and `beta`.
#'
#' The mean effect corresponds to the first sensitivity parameters. 
#' The computed standard deviation can be used to inform the magnitude of the
#' second sensitivity parameter `sigma` in [post_instrument_bounds()].
#'

#' @return A named vector with mean and standard deviation of the distribution.
#' @references Schuessler, J., Glynn, A. N., & Rueda, M.R. (2023). Post-Instrument Bias.

#' @examples
#' x <- seq(0, 1, length.out = 100)
#' plot(x, dbeta(x, shape1 = 1, shape2 = 1), type = "l",
#'ylim = c(0, 3.5),
#'xlab = expression(paste(gamma[i])),
#'ylab = "Density")
#'lines(x, dbeta(x, shape1 = 0.5, shape2 = 0.5), type = "l",
#'      lty = "dashed")
#'lines(x, dbeta(x, shape1 = 4, shape2 = 4), type = "l",
#'      lty = "dotted")
#'legend(x = "top", c("Low Variance", "Medium Variance", "High Variance"),
#'       lty = c("dotted", "solid", "dashed"))
#'
#'beta_sd(alpha = 1, beta = 1, a= 0, b = 1) # medium variance
#'beta_sd(alpha = 0.5, beta = 0.5, a = 0, b = 1) # high variance
#'beta_sd(alpha = 4, beta = 4, a = 0, b = 1) # low variance
#' @export

beta_sd <- function(alpha, beta, a, b){
  return(c(
    "mean" = (alpha*b + beta*a)/(alpha + beta),
    "standard deviation" = sqrt(
      (alpha*beta*(b - a)^2)/(
        (alpha + beta)^2*(alpha + beta + 1)
        )
      )
    )
    )
}
