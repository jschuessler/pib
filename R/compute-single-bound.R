#' Compute Bounds for Specific Sensitivity Parameters
#'
#' @param gamma Numeric value specifying the sensitivity parameter gamma for which to
#'             compute bounds
#' @param sigma Numeric value specifying the sensitivity parameter sigma for which to
#'             compute bounds. Can be NULL if bounds only vary with gamma.
#' @param bounds A data frame containing pre-computed bounds from [post_instrument_bounds()]
#'              with columns for gamma, sigma, point estimates, and confidence intervals
#'
#' @details
#' This function takes pre-computed bounds from [post_instrument_bounds()] and
#' computes the bounds corresponding to specific values of the sensitivity parameters gamma
#' and sigma. 
#'
#' @return A named vector containing:
#' \itemize{
#'   \item gamma: The input gamma value
#'   \item sigma: The input sigma value
#'   \item point estimate: The computed point estimate
#'   \item CI values: The computed confidence interval bounds (names match those in input bounds)
#' }
#'
#' @references
#' Schuessler, J., Glynn, A. N., & Rueda, M.R. (2023). Post-Instrument Bias.
#'
#' Glynn, A. N., Rueda, M. R., & Schuessler, J. (2023). Post-Instrument Bias in Linear Models.
#' Sociological Methods & Research.
#'
#' @examples
#' \dontrun{
#' # First compute bounds over a grid of gamma and sigma
#' bounds_grid <- post_instrument_bounds(gamma = c(-0.02, 0.02),
#'                                      sigma = c(0, 0.002),
#'                                      n.sigma = 5,
#'                                      estimates = sensitivity_estimates)
#'
#' # Then compute bounds for specific values
#' specific_bounds <- compute_single_bounds(gamma = 0.015,
#'                                         sigma = 0.001,
#'                                         bounds = bounds_grid)
#' }
#'
#' @seealso [post_instrument_bounds()]
#' @export

compute_single_bounds <- function(gamma, sigma = NULL, bounds){
  
  ci_names <- names(bounds)[grepl("CI", names(bounds))]
  
  
  out <- c("gamma" = gamma, "sigma" = sigma, 
           "point estimate" = NA, rep(NA, length(ci_names)))
  names(out) <- c("gamma", "sigma", "point estimate", ci_names)
  
  # determine gamma closest to user-specified gamma
  closest_gamma <- bounds$gamma[which.min(abs(bounds$gamma - gamma))]
  
  # determine sigma closest to user-specified sigma
  closest_sigma <-  bounds$sigma[which.min(abs(bounds$sigma - sigma))]
  
  # determine smallest and largest gamma values in bounds matrix
  largest_gamma <- bounds$gamma[which.max(bounds$gamma)]
  smallest_gamma <- bounds$gamma[which.min(bounds$gamma)] 
  
  # for point estimate
  slope <- (bounds[bounds$gamma == largest_gamma & 
                     bounds$sigma == 0, "Lower Point Estimate"] - 
              bounds[bounds$gamma == smallest_gamma & 
                       bounds$sigma == 0, "Lower Point Estimate"])/
    (largest_gamma - smallest_gamma)
  
  intercept <- bounds[bounds$gamma == largest_gamma & 
                        bounds$sigma == closest_sigma, "Lower Point Estimate"] - 
    slope*largest_gamma
  
  out["point estimate"] <- intercept + gamma*slope
  
  # for each confidence interval
  for(ci in ci_names){
    
    # take smallest and largest gamma to compute intercept and slope for that sigma
    slope <- (bounds[bounds$gamma == largest_gamma & 
                       bounds$sigma == closest_sigma, ci] - 
                bounds[bounds$gamma == smallest_gamma & 
                         bounds$sigma == closest_sigma, ci])/
      (largest_gamma - smallest_gamma)
    
    intercept <- bounds[bounds$gamma == largest_gamma & 
                          bounds$sigma == closest_sigma, ci] - 
      slope*largest_gamma
    
    # take the next-smallest or largest sigma. compute differences between the two sigmas
    # for a gamma close to user-specified gamma to compute intercept shift:
    
    if(!is.null(sigma) & length(unique(bounds$sigma)) < 2){
      stop("When sigma is specified, input bounds matrix needs to have at least two different values for sigma.")
    }
    
    unique_sigmas <- unique(bounds$sigma)
    pos <- which(unique_sigmas == closest_sigma)
    if(pos == 1) {  # if smallest
      other_sigma <- unique_sigmas[2]
    } else {
      other_sigma <- unique_sigmas[length(unique_sigmas) - 1]
    }
    
    sigmas <- c(closest_sigma, other_sigma)[order(c(closest_sigma, other_sigma))]
    
    intercept_shift <- bounds[bounds$gamma == closest_gamma & 
                                bounds$sigma == sigmas[2], ci] -
      bounds[bounds$gamma == closest_gamma & 
               bounds$sigma == sigmas[1], ci]
    
    intercept <- intercept + intercept_shift
    # use intercept and slope to predict effect for user-specified gamma, sigma
    
    out[ci] <- intercept + gamma*slope
  }
  return(out)
}



