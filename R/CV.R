#' Cross-Validation
#'
#' @param d.train A dataset with both outcomes and predictors.
#' @param model A "model" function.
#' @param params.grid A list where each element is a vector containing all the values that are tried for that parameter.
#' @param params.auto An object of class "Autoparams" (see ?autoparams).
#' @param n.folds The number of folds.
#' @param folds.strata Factor according to which the folds are stratified.
#' @param model.strata Factor according to which the model is stratified.
#' @param fit.best Whether to fit the best model.
#' @param threads Number of threads.
#' @param verbose Whether to print what is being done.
#'
#' @return A list with table, aggr, and (if fit.best is TRUE) a fitted model.
#'
#' @export
CV <- function(d.train, model, params.grid = NULL, params.auto = NULL,
               n.folds = 5, folds.strata = NULL, model.strata = NULL,
               fit.best = FALSE, threads = 1, verbose = 1) {
    # Validate the arguments
    if (is.null(params.grid) && is.null(params.auto))
        stop("At least one of params.grid and params.auto must be specified.")
    if (!is.null(params.grid) && !is.null(params.auto))
        stop("At most one of params.grid and params.auto can be specified.")

    # Define the folds
    if (!is.null(folds.strata)) {
        folds <- cut(1:length(unique(folds.strata)), breaks = n.folds, labels = FALSE)
        folds <- sample(folds)
        folds <- rep(folds, sapply(unique(folds.strata), function(i) sum(folds.strata == i))) # FIXME: iid can be sparse
    } else {
        folds <- cut(1:nrow(d.train), breaks = n.folds, labels = FALSE)
        folds <- sample(folds)
    }

    if (!is.null(params.auto)) {
        n_params <- length(params.auto@best)
        params.auto@max_iter <- rep(params.auto@max_iter, length.out = n_params) # TODO: add warning if it doesn't recycle
        params.auto@improve_threshold <- rep(params.auto@improve_threshold, length.out = n_params) # TODO: add warning if it doesn't recycle
        loss <- params.auto@loss

        kfold.cv.table <- data.table()
        kfold.cv.aggr <- data.table()

        # assign the starting value to p
        p <- params.auto@best

        # train the first model
        kfold.cv <- parallel::mclapply(mc.cores = threads, mc.preschedule = FALSE, mc.set.seed = FALSE, 1:n.folds, function(f) {
        # kfold.cv <- lapply(1:n.folds, function(f) {
            model(d.train[folds != f], d.train[folds == f], p, model.strata)
        })
        kfold.cv.table <- rbind(kfold.cv.table, MergeFolds(kfold.cv))
        kfold.cv.aggr <- rbind(kfold.cv.aggr, AggregateFolds(kfold.cv.table[nrow(kfold.cv.table)]))
        if (verbose >= 2)
            cat("Params", unlist(p), "; loss", loss(kfold.cv.aggr[nrow(kfold.cv.aggr)]), "\n")

        # Minimise the loss for each parameter independently
        prev_loss <- loss(kfold.cv.aggr[nrow(kfold.cv.aggr)])
        for (i in 1:n_params) {
            if (is.null(params.auto@range[[i]]))
                next

            j <- 1
            max_iter <- params.auto@max_iter[i]
            imp_thresh <- params.auto@improve_threshold[i]

            # TODO: we can embed the two loops into a for loop ("right", "left") as in tuneRF.
            # go forward
            while (j <= max_iter) {
                old_p <- p[[i]]
                p[[i]] <- p[[i]] + (params.auto@range[[i]][2] - p[[i]]) / 2
                if (p[[i]] == old_p)
                    break
                kfold.cv <- parallel::mclapply(mc.cores = threads, mc.preschedule = FALSE, mc.set.seed = FALSE, 1:n.folds, function(f) {
                # kfold.cv <- lapply(1:n.folds, function(f) {
                    model(d.train[folds != f], d.train[folds == f], p, model.strata)
                })
                kfold.cv.table <- rbind(kfold.cv.table, MergeFolds(kfold.cv))
                kfold.cv.aggr <- rbind(kfold.cv.aggr, AggregateFolds(kfold.cv.table[nrow(kfold.cv.table)]))

                if (loss(kfold.cv.aggr[nrow(kfold.cv.aggr)]) < imp_thresh * prev_loss) {
                    if (verbose >= 2)
                        cat("Params", unlist(p), "; loss", loss(kfold.cv.aggr[nrow(kfold.cv.aggr)]), "< prev_loss", prev_loss, "\t:D\n")
                    prev_loss <- loss(kfold.cv.aggr[nrow(kfold.cv.aggr)])
                    params.auto@best[[i]] <- p[[i]]
                    params.auto@range[[i]][1] <- old_p
                    j <- j + 1
                } else {
                    if (verbose >= 2)
                        cat("Params", unlist(p), "; loss", loss(kfold.cv.aggr[nrow(kfold.cv.aggr)]), ">= prev_loss", prev_loss, "\t:c\n")
                    p[[i]] <- params.auto@best[[i]]
                    break
                }
            }

            # go backward
            while (j <= max_iter) {
                old_p <- p[[i]]
                p[[i]] <- p[[i]] - (p[[i]] - params.auto@range[[i]][1]) / 2
                if (p[[i]] == old_p)
                    break
                kfold.cv <- parallel::mclapply(mc.cores = threads, mc.preschedule = FALSE, mc.set.seed = FALSE, 1:n.folds, function(f) {
                # kfold.cv <- lapply(1:n.folds, function(f) {
                    model(d.train[folds != f], d.train[folds == f], p, model.strata)
                })
                kfold.cv.table <- rbind(kfold.cv.table, MergeFolds(kfold.cv))
                kfold.cv.aggr <- rbind(kfold.cv.aggr, AggregateFolds(kfold.cv.table[nrow(kfold.cv.table)]))

                if (loss(kfold.cv.aggr[nrow(kfold.cv.aggr)]) < imp_thresh * prev_loss) {
                    if (verbose >= 2)
                        cat("Params", unlist(p), "; loss", loss(kfold.cv.aggr[nrow(kfold.cv.aggr)]), "< prev_loss", prev_loss, "\t:D\n")
                    prev_loss <- loss(kfold.cv.aggr[nrow(kfold.cv.aggr)])
                    params.auto@best[[i]] <- p[[i]]
                    j <- j + 1
                } else {
                    if (verbose >= 2)
                        cat("Params", unlist(p), "; loss", loss(kfold.cv.aggr[nrow(kfold.cv.aggr)]), ">= prev_loss", prev_loss, "\t:c\n")
                    p[[i]] <- params.auto@best[[i]]
                    break
                }
            }
        }
        l <- list(table = kfold.cv.table, aggr = kfold.cv.aggr)
    }

    if (!is.null(params.grid)) {
        n_params <- length(params.grid)
        # For each fold, loop through all the parameters, then merge
        # everything at the end. This is more efficient than having the loop
        # over the parameters as the outer one.
        kfold.cv <- parallel::mclapply(mc.cores = threads, mc.preschedule = FALSE, mc.set.seed = FALSE, 1:n.folds, function(f) {
        # kfold.cv <- lapply(1:n.folds, function(f) {
            expanded.params <- do.call(expand.grid, params.grid)
            kfold.cv.tables <- apply(expanded.params, 1, function(p) {
                if (verbose >= 2)
                    cat("Fold", f, "; params", p, "\n")

                model(d.train[folds != f], d.train[folds == f], as.list(p), model.strata)
            })
            do.call("rbind", kfold.cv.tables)
        })
        kfold.cv.table <- MergeFolds(kfold.cv)
        l <- list(table = kfold.cv.table, aggr = AggregateFolds(kfold.cv.table))
    }

    if (fit.best) {
        # Grow the forest with the best parameters for this gene
        kfold.cv.params <- as.list(l$aggr[which.max(pred_perf_R2), 1:n_params]) # TODO: use loss() to find best
        kfold.cv.fit <- model(d.train, NULL, kfold.cv.params, model.strata, return.fit = TRUE)
        l <- append(l, list(fit = kfold.cv.fit))
    }

    l
}

#' Nested Cross-Validation
#'
#' @inheritParams CV
#' @param n.outer.folds Number of folds in the outer CV.
#' @param n.inner.folds Number of folds in the inner CV.
#'
#' @return A list with inner CV, table, and aggr.
#'
#' @export
NestedCV <- function(d.train, model, params.grid = NULL, params.auto = NULL, n.outer.folds = 5, n.inner.folds = 10,
                     folds.strata = NULL, model.strata = NULL, threads = 1, verbose = 1) {
    # Validate the arguments
    if (is.null(params.grid) && is.null(params.auto))
        stop("At least one of params.grid and params.auto must be specified.")
    if (!is.null(params.grid) && !is.null(params.auto))
        stop("At most one of params.grid and params.auto can be specified.")

    # Find the number of params
    if (!is.null(params.grid)) n_params <- length(params.grid)
    if (!is.null(params.auto)) n_params <- length(params.auto@best)

    # Define the outer folds
    if (!is.null(folds.strata)) {
        folds <- cut(1:length(unique(folds.strata)), breaks = n.outer.folds, labels = FALSE)
        folds <- sample(folds)
        folds <- rep(folds, sapply(unique(folds.strata), function(i) sum(folds.strata == i)))
    } else {
        folds <- cut(1:nrow(d.train), breaks = n.outer.folds, labels = FALSE)
        folds <- sample(folds)
    }

    # Threads: assign the highest possible number of cores to the outer CV
    # loop and the rest to the inner CV, without exceeding the total number
    # of cores. This is more efficient than prioritising the inner loop,
    # but the number of cores provided by the user should always be a
    # multiple of the number of outer folds, otherwise some of them are
    # wasted.
    if (threads > n.outer.folds) {
        threads.outer <- n.outer.folds
        threads.inner <- (threads - threads.outer) %/% threads.outer + 1
    } else {
        threads.inner <- 1
        threads.outer <- threads
    }

    outer.cv <- parallel::mclapply(mc.cores = threads.outer, mc.preschedule = FALSE, mc.set.seed = FALSE, 1:n.outer.folds, function(f) {
    # outer.cv <- lapply(1:n.outer.folds, function(f) {
        if (verbose)
            cat("Outer fold:", f, "\n")

        if ("n.folds" %in% names(formals(model))) {
            # model() should return a list of inner.cv and pearson
            model(d.train[folds != f], d.train[folds == f], params.grid, params.auto,
                  n.fold = n.inner.folds, folds.strata, model.strata,
                  return.fit = FALSE, threads = threads.inner, verbose = verbose - 1)
        } else {
            inner.cv <- CV(d.train[folds != f], model, params.grid, params.auto, n.folds = n.inner.folds, threads = threads.inner, verbose = verbose - 1)
            inner.cv.params <- as.list(inner.cv$aggr[which.max(pred_perf_R2), 1:n_params])
            pearson <- model(d.train[folds != f], d.train[folds == f], inner.cv.params, model.strata)
            list(inner.cv = inner.cv, pearson = pearson)
        }
    })

    # outer.cv is a list of five lists, each with two elements: inner.cv and pearson.
    # inner.cv is itself, for each fold, a list of two elements: table and aggr.
    if (any(sapply(outer.cv, function(e) class(e) == "try-error"))) {
        print(warnings())
        print(str(outer.cv))
        return(NULL)
    }

    nested.cv <- mapply(outer.cv, FUN = I)
    inner.cv <- nested.cv[1, ]
    outer.cv <- nested.cv[2, ]

    outer.cv <- lapply(outer.cv, function(f) {
        f[, 1:n_params] <- NA
        f
    })
    outer.cv.table <- MergeFolds(outer.cv)
    outer.cv.aggr_table <- AggregateFolds(outer.cv.table)

    list(inner.cv = inner.cv, table = outer.cv.table, aggr = outer.cv.aggr_table)
}

