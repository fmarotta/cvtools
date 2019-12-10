#' Merge the folds to obtain a table with all the folds.
#' Operates on a list of "pearson" objects returned by Call_RF().
#'
#' @param cv.list A list of "pearson" data.tables.
#' @param n.params The number of parameters.
#'
#' @return A single data.table with all the folds cbinded.
MergeFolds <- function(cv.list) {
    n.params <- which(names(cv.list[[1]]) == "data.name") - 1
    i <- 1 # it may happend that the first fold is useless, so we go on until we find a viable one.
    while (i <= length(cv.list) && length(cv.list[[i]]) < n.params + 6) i <- i + 1
    if (i > length(cv.list))
        print(cv.list) # DEBUG
    cv.table <- cv.list[[i]]
    names(cv.table)[(n.params + 2):(n.params + 6)] <- c("mse_1", "rsq_1", "parameter.df_1", "estimate.cor_1", "p.value_1")
    for (k in (i + 1):length(cv.list)) {
        cv.item <- cv.list[[k]]
        if (length(names(cv.item)) != n.params + 6)
            return(NULL)
        names(cv.item)[(n.params + 2):(n.params + 6)] <- c(paste0("mse_", k), paste0("rsq_", k), paste0("parameter.df_", k), paste0("estimate.cor_", k), paste0("p.value_", k))
        cv.table <- merge(cv.table, cv.item, all = TRUE)
    }
    cv.table
}

#' Compute the mean for all the folds.
#' Operates on the output of MergeFolds().
#'
#' @param cv.table The output of MergeFolds().
#' @param n.params The number of parameters.
#'
#' @return A single data.table with all the folds aggregated.
AggregateFolds <- function(cv.table) {
    n = length(cv.table)
    n.params <- which(names(cv.table) == "data.name") - 1
    cbind(
        cv.table[, 1:(n.params + 1)],
        cv.table[, .(mse_avg = rowMeans(.SD, na.rm = T),
                     mse_sd = byrow::rowSds(.SD)), .SDcols = seq(n.params + 2, n, 5)],
        cv.table[, .(rsq_avg = rowMeans(.SD, na.rm = T),
                     rsq_sd = byrow::rowSds(.SD)), .SDcols = seq(n.params + 3, n, 5)],
        cv.table[, .(rho_avg = rowMeans(.SD, na.rm = T),
                     rho_sd = byrow::rowSds(.SD),
                     pred_perf_R2 = rowMeans(.SD^2, na.rm = T)), .SDcols = seq(n.params + 5, n, 5)],
        cv.table[, .(rho_zscore = byrow::rowZscs(.SD)$zscore,
                     pred_perf_pval = byrow::rowZscs(.SD)$pvalue), .SDcols = seq(n.params + 6, n, 5)])
}

#' Compute the mean for all the genes.
#' Operates on the output of AggregateFolds().
#'
#' @param cv.aggr_table The output of AggregateFolds().
#' @param n.params The number of parameters.
#'
#' @return A single data.table with all the rows summarised.
SummariseGenes <- function(cv.aggr_table) {
    n.params <- which(names(cv.aggr_table) == "data.name") - 1
    cv.aggr_table[, .(
        prop_sign = mean(pred_perf_pval < 0.05, na.rm = TRUE),
        mean_rho_avg = mean(rho_avg, na.rm = TRUE),
        sd_rho_avg = stats::sd(rho_avg, na.rm = TRUE),
        mean_pred_perf_pval = mean(pred_perf_pval, na.rm = TRUE),
        sd_pred_perf_pval = stats::sd(pred_perf_pval, na.rm = TRUE),
        mean_rho_avg_of_sign = mean(rho_avg[pred_perf_pval < 0.05], na.rm = TRUE),
        sd_rho_avg_of_sign = stats::sd(rho_avg[pred_perf_pval < 0.05], na.rm = TRUE)
    ), by = eval(names(cv.aggr.table)[1:n.params])]
}
