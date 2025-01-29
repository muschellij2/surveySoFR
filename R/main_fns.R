#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats coef
NULL


#' Function for doing brr
#' @param df data frame containing all variables used in initial fitting, plus survey variables (strata, PSU, weights)
#' @param fit Initial fit on the full data, pass the formula into [mgcv::gam()]
#' @param nbrr  Number of BRR samples
#' @param us Unique strata identifiers (corresponds to the unique levels of the variable specified by the strata argument)
#' @param uPSU Unique PSU identifiers (corresponds to the unique levels of the variable specified by the PSU argument)
#' @param strata Character vector of length 1 indicating the variable name for the strata identifier in the data frame used for model fitting
#' @param PSU Character vector of length 1 indicating the variable name for the PSU identifier in the data frame used for model fitting
#' @param weights Character vector of length 1 indicating the variable name for individuals' survey weights in the data frame used for model fitting
doBRR <- function(df,
                  fit,
                  nbrr,
                  us = 1:30,
                  uPSU = 1:2,
                  strata = NULL,
                  PSU = NULL,
                  weights = NULL) {
    if (!inherits(fit, "gam")) {
        stop("Regression model supplied to fit must be of class \"gam\" ")
    }
    if ("strata_PSU" %in% colnames(df)) {
        message(
            "strata_PSU column detected in dataframe used for model fitting. This column will be used for BRR."
        )
    } else {
        df$strata_PSU <- paste0(df[[strata]], "_", df[[PSU]])
    }
    n_us <- length(us)
    # get the hadamard matrix
    H_brr <- matrix(paste0(rep(us, nbrr), "_", sample(uPSU, n_us * nbrr, replace =
                                                          TRUE)), n_us, nbrr, byrow = FALSE)
    # set up empty container for storing the BRR coefficients
    coef_brr <- matrix(NA, nbrr, length(stats::coef(fit)))
    # loop over the BRR samples
    pb <- utils::txtProgressBar(min = 0,
                                max = nbrr,
                                style = 3)
    for (b in 1:nbrr) {
        # subset to the current PSU-strata sample
        data_brr <- df[df$strata_PSU %in% H_brr[, b], ]
        # double the weights in the retained sample, normalize
        weights_b <- 2 * data_brr[[weights]] / mean(data_brr[[weights]])
        data_brr$weights_b <- weights_b
        ## fit regression using the brr data
        fit_brr <- mgcv::gam(
            fit$formula,
            data = data_brr,
            family = fit$family$family,
            method = fit$method,
            weights = weights_b
        )
        # save the coefficients
        coef_brr[b, ] <- stats::coef(fit_brr)
        utils::setTxtProgressBar(pb, b)
    }
    return(list("H" = H_brr, "coefficients" = coef_brr))
}
