p_adjust_twt <- function(pval, p, t0, t_coeff, aggregation_strategy) {
  # Sub- and super-level sets of the pointwise p-values, in C++ (see
  # src/adjust.cpp); the adjusted p-values keep the names of pval
  adjusted_pval <- twt_adjust_cpp(
    as.double(pval),
    as.double(t0),
    as_stat_matrix(t_coeff, p),
    aggregation_strategy == "max"
  )
  names(adjusted_pval) <- names(pval)

  list(
    adjusted_pvalues = adjusted_pval
  )
}
