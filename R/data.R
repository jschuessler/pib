#' Data from Carnegie/Marinov 2017
#'
#' A subset of data from Allison Carnegie and Nikolay Marinov (2017), "Foreign Aid, Human Rights, and Democracy Promotion:
#' Evidence from a Natural Experiment", American Journal of Political Science 61(3).
#'
#' @source <https://doi.org/10.15139/S3/X5WXLM>
"carnegie_marinov"

#' Estimated components for sensitivity analysis for Carnegie/Marinov 2017
#'
#' @examples
#' \dontrun{
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
"out_cm"

#' Data from Spenkuch/Tillmann 2017
#'
#' A subset of data from Spenkuch, Jörg L., and Philipp Tillmann.
#' "Elite influence? Religion and the electoral success of the Nazis." American Journal of Political Science 62.1 (2018): 19-36.
#'
#' @source <https://doi.org/10.7910/dvn/j2fa5r>
"spenkuch_tillmann"

#' Estimated components for sensitivity analysis for Spenkuch/Tillmann 2017
#'
#' @examples
#' \dontrun{
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

