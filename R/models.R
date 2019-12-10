# Utility methods for model functions.

#' Title
#'
#' @param ytest Vector of true values.
#' @param ypred Vector of predicted values.
#' @param by Factor by which to split the data.
#' @param alternative See cor.test().
#'
#' @return A "pearson" data.table.
#' @export
compute_pearson <- function(ytest, ypred, by = NULL, alternative = "g") {
    if (is.null(by)) {
        pearson <- t(unlist(cor.test(ytest, ypred, alternative = alternative)))
    } else {
        # Analyse the correlation by groups
        test.split <- split(ytest, as.factor(by))
        pred.split <- split(ypred, as.factor(by))
        pearson <- t(sapply(1:length(test.split), function(r) {
            l <- cor.test(test.split[[r]], pred.split[[r]], alternative = alternative)
            l$data.name <- names(test.split[r])
            unlist(l)
        }))
    }
    pearson <- as.data.table(pearson)
    pearson[, parameter.df := as.numeric(parameter.df)]
    pearson[, estimate.cor := as.numeric(estimate.cor)]
    pearson[, p.value := as.numeric(p.value)]

    pearson
}

#' MSE
#'
#' Mean squared error of the
#'
#' @param ytest A vector of the true response values.
#' @param ypred A vector of the predicted values.
#' @param p The number of predictors used by the model.
#'
#' @return The mean squared error of the model on the hold-out fold.
#' @export
cvmse <- function(ytest, ypred) {
    mean((ytest - ypred)^2)
}

#' \ifelse{html}{\out{R<sup>2</sup>}}{\eqn{R^2}}
#'
#' @inheritParams cvmse
#'
#' @return The \ifelse{html}{\out{R<sup>2</sup>}}{\eqn{R^2}} of the model on the hold-out fold.
#' @export
cvrsq <- function(ytest, ypred) {
    SS.total      <- sum((ytest - mean(ytest))^2)
    SS.residual   <- sum((ytest - ypred)^2)

    1 - SS.residual / SS.total
}

#' Fraction of Variability
#'
#' @inheritParams cvmse
#'
#' @return The fraction of variability explained by the model in the hold-out fold.
#' @export
cvvarfrac <- function(ytest, ypred) {
    SS.total      <- sum((ytest - mean(ytest))^2)
    SS.regression <- sum((ypred - mean(ytest))^2)

    SS.regression / SS.total
}