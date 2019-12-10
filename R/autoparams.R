# TODO: validating function

#' Autoparams (S4 class)
#'
#' @slot best list.
#' @slot range list.
#' @slot max_iter numeric.
#' @slot improve_threshold numeric.
#' @slot loss function.
setClass("Autoparams",
         slots = list(
             best = "list",
             range = "list",
             max_iter = "numeric",
             improve_threshold = "numeric",
             loss = "function"))

#' Autoparams (constructor)
#'
#' @param start Named list where each element corresponds to a parameter and stores the starting value.
#' @param range Names list where each element is a vector of two elements, min and max.
#' @param loss Function of aggr; its return value is minimised during the cross-validation.
#' @param max_iter Maximum number of successful iterations.
#' @param improve_threshold Percentage of the previous loss that must be achieved.
#'
#' @return An instance of the Autoparams class.
#'
#' @export
autoparams <- function(start, range, loss = function(aggr) {1 - aggr$pred_perf_R2},
                       max_iter = 5, improve_threshold = 1) {
    new("Autoparams",
        best = start,
        range = range,
        loss = loss,
        max_iter = max_iter,
        improve_threshold = improve_threshold)
}