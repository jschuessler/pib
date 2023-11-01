#' Plot estimated bounds for sensitivity analysis with respect to post-instrument covariates.
#' @param bounds           An object from [post.instrument.bounds()]
#' @param whichCI          Determines the type of boostrap confidence interval
#'                    to be plotted. Either `"percentile"` (the default) or `"basic"`
#' @param annotate.sigma   Should the confidence intervals have annotated
#'                    values of sigma? Default is `"minmax"` (annotating CI for
#'                    lowest and highest values of sigma). Alternatives can be
#'                    `"none"` or `"all"`
#' @param digits           When annotating `sigma`, this controls the rounding of
#'                    the annotated values. Default is 2.
#'
#' @param set.ylim         Logical. Should limit of Y-axis be automatically
#'                    determined? If `FALSE`, argument `ylim` has to be set.
#'                    Default is `TRUE`
#' @param set.xlim         Logical. Should limit of X-axis be automatically
#'                    determined? If `FALSE`, argument `xlim` has to be set.
#'                    Default is `TRUE`
#' @param xlim,ylim   If `set.ylim = F` or `set.ylim = F`, sets limits of axes.
#' @param horizontal.lty         Controls appearance of horizontal line at 0. Default is `"dashed"`.
#' @param ylab,xlab,... Label of Y- and X-axis, and further graphical parameters.
#' @export

plot_pib_bounds <- function(bounds, whichCI = "percentile",
                                        annotate.sigma = "minmax",
                                        digits = 2,
                                        set.ylim = T,
                                        set.xlim = T,
                                        ylab = "Effect Estimate",
                                        xlab = expression(paste("E[", gamma[i], "]")),
                                        ylim = NULL,
                                        xlim = NULL,
                                        horizontal.lty = "dashed",
                                        ...){
  if(whichCI == "percentile"){
    CInames <- c("Lower Bound CI (Percentile)",
                 "Upper Bound CI (Percentile)")
  } else if(whichCI == "basic"){
    CInames <- c("Lower Bound CI (Basic)",
                 "Upper Bound CI (Basic)")
  } else stop("whichCI must be either 'percentile' or 'basic'!")

  unique.sigma <- unique(bounds$sigma)

  if(set.ylim){
    # set ylim:
    min.lim <- min(bounds[,CInames[1]], bounds[,CInames[2]], 0)
    max.lim <- max(bounds[,CInames[1]], bounds[,CInames[2]], 0)

    lim.increase <- min(abs(c(min.lim*0.2, max.lim*0.2)))

    min.lim <- min.lim + lim.increase*sign(min.lim)
    max.lim <- max.lim + lim.increase*sign(max.lim)
    ylim <- c(min.lim, max.lim)
  }

  if(set.xlim){
    xlim <- c(min(bounds$gamma), max(bounds$gamma)*1.15)
  }

  plot(x = NULL,
       xlim = xlim,
       ylim = ylim,
       ylab = ylab,
       xlab = xlab,
       ...)

  lines(par("usr")[1:2], c(0, 0),
        lty = horizontal.lty)

  # use quadratic approximation to plot confidence intervals
  for(s in 1:length(unique.sigma)){
    model <- lm(bounds[,CInames[1]][bounds$sigma == unique.sigma[s]] ~
                  bounds$gamma[bounds$sigma == unique.sigma[s]] +
                  I(bounds$gamma[bounds$sigma == unique.sigma[s]]^2))

    lines(bounds$gamma[bounds$sigma == unique.sigma[s]],
          coef(model)[1] +
            coef(model)[2]*bounds$gamma[bounds$sigma == unique.sigma[s]] +
            coef(model)[3]*bounds$gamma[bounds$sigma == unique.sigma[s]]^2,
          lty = "solid")

    model <- lm(bounds[,CInames[2]][bounds$sigma == unique.sigma[s]] ~
                  bounds$gamma[bounds$sigma == unique.sigma[s]] +
                  I(bounds$gamma[bounds$sigma == unique.sigma[s]]^2))

    lines(bounds$gamma[bounds$sigma == unique.sigma[s]],
          coef(model)[1] +
            coef(model)[2]*bounds$gamma[bounds$sigma == unique.sigma[s]] +
            coef(model)[3]*bounds$gamma[bounds$sigma == unique.sigma[s]]^2,
          lty = "solid")

  }

  par(col = "black")

  if(annotate.sigma == "minmax" | annotate.sigma == "all"){
    # annotate min and max of sigma to CI which is closest to zero
    # these plot elements are also needed for when annotate.sigma == "all"

    x.annotation <- max(bounds$gamma) + abs(max(bounds$gamma))*0.1
    closest.to.zero <- names(which.min(abs(apply(bounds[,CInames], 2, sum))))

    text(x.annotation,
         bounds[bounds$sigma == 0 & bounds$gamma == max(bounds$gamma), closest.to.zero],
         labels = bquote(paste(sigma[gamma[i]], " = ",
                               .(round(min(bounds$sigma), digits)))))

    if(grepl("Lower", closest.to.zero)){
      text(x.annotation,
           #1/2*
           bounds[bounds$sigma == max(bounds$sigma) & bounds$gamma == max(bounds$gamma),
                  closest.to.zero] #+
           #   1/2*par()$usr[3]
           ,
           labels = bquote(paste(sigma[gamma[i]], " = ",
                                 .(round(max(bounds$sigma), digits))))
      )
    } else {
      text(x.annotation,
           #1/2*
           bounds[bounds$sigma == max(bounds$sigma) & bounds$gamma == max(bounds$gamma),
                  closest.to.zero] #+
           #   1/2*par()$usr[4]
           ,
           labels = bquote(paste(sigma[gamma[i]], " = ",
                                 .(round(max(bounds$sigma), digits))))
      )
    }
  }

  if(annotate.sigma == "all"){
    if(length(unique(bounds$sigma)) > 2){
      for(s in 2:(length(unique(bounds$sigma))-1)){

        text(x.annotation,
             bounds[bounds$sigma == unique(bounds$sigma)[s] &
                      bounds$gamma == max(bounds$gamma), closest.to.zero],
             labels = bquote(paste(sigma[gamma[i]], " = ",
                                   .(round(unique(bounds$sigma)[s], digits)))))
      }
    }
  }

  # plot point estimate (midpoint between CIs):
  lines(bounds$gamma[bounds$sigma == 0],
        c(bounds[bounds$sigma == 0, CInames[1]] +
            bounds[bounds$sigma == 0, CInames[2]])/2,
        lwd = 2)

}
