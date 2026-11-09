p_adjust_global <- function(
  aggregation_strategy,
  t0,
  t_coeff,
  p
) {
  pval_temp <- global_pvalue_cpp(
    as.double(t0),
    as_stat_matrix(t_coeff, p),
    aggregation_strategy == "max"
  )

  list(
    adjusted_pvalues = rep(pval_temp, p)
  )
}
