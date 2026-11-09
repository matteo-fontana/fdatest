# Tests for the C++ adjustments (src/adjust.cpp) used by p_adjust_iwt(),
# p_adjust_twt() and p_adjust_global() when no mirai daemons are set: their
# results must equal those of the R implementations below (for the IWT, the R
# code that the mirai path still uses).
# jarl-ignore-file internal_function: test files legitimately access internal functions via :::

mirai::daemons(0L)

# The C++ code accumulates in long double as R does: results are identical
# when both have the same precision (the usual builds of R, including those
# where long double is double). Otherwise compare with a tolerance and skip
# the inputs full of ties, where the rounding can change a comparison.
r_digits <- if (capabilities("long.double")) .Machine$longdouble.digits else 53L
exact <- identical(as.integer(r_digits), fdatest:::ldouble_digits_cpp())
expect_same <- function(current, target, info) {
  if (exact) {
    expect_identical(current, target, info = info)
  } else {
    expect_equal(current, target, tolerance = 1e-8, info = info)
  }
}

# ---------------------------------------------------------------------------
# R implementations
# ---------------------------------------------------------------------------
iwt_reference <- function(p, pval, t0, t_coeff, n_perm, recycle, agg) {
  m <- fdatest:::iwt_pvalue_matrix_r(p, pval, t0, t_coeff, n_perm, recycle, agg)
  list(
    adjusted_pvalues = fdatest:::pval_correct_cpp(m)[1, ],
    pvalue_matrix = m
  )
}

twt_reference <- function(pval, p, t0, t_coeff, aggregation_strategy) {
  thresholds <- c(0, sort(unique(pval)), 1)
  adjusted_pval <- pval
  pval_tmp <- rep(0, p)
  for (test in seq_along(thresholds)) {
    points_1 <- which(pval <= thresholds[test])
    t0_comb <- if (aggregation_strategy == "integral") {
      sum(t0[points_1], na.rm = TRUE)
    } else {
      max(t0[points_1], na.rm = TRUE)
    }
    t_comb <- if (aggregation_strategy == "integral") {
      rowSums(t_coeff[, points_1, drop = FALSE], na.rm = TRUE)
    } else {
      apply(t_coeff[, points_1, drop = FALSE], 1, max, na.rm = TRUE)
    }
    pval_tmp[points_1] <- mean(t_comb >= t0_comb)
    adjusted_pval <- apply(rbind(adjusted_pval, pval_tmp), 2, max)

    points_2 <- which(pval > thresholds[test])
    t0_comb <- if (aggregation_strategy == "integral") {
      sum(t0[points_2], na.rm = TRUE)
    } else {
      max(t0[points_2], na.rm = TRUE)
    }
    t_comb <- if (aggregation_strategy == "integral") {
      rowSums(t_coeff[, points_2, drop = FALSE], na.rm = TRUE)
    } else {
      apply(t_coeff[, points_2, drop = FALSE], 1, max, na.rm = TRUE)
    }
    pval_tmp[points_2] <- mean(t_comb >= t0_comb)
    adjusted_pval <- apply(rbind(adjusted_pval, pval_tmp), 2, max)
  }
  list(adjusted_pvalues = adjusted_pval)
}

global_reference <- function(aggregation_strategy, t0, t_coeff, p) {
  pval_temp <- if (aggregation_strategy == "integral") {
    mean(rowSums(t_coeff) >= sum(t0))
  } else {
    mean(apply(t_coeff, 1, max) >= max(t0))
  }
  list(adjusted_pvalues = rep(pval_temp, p))
}

pointwise_pvalues <- function(t0, t_coeff) {
  t0_mat <- matrix(t0, nrow = nrow(t_coeff), ncol = ncol(t_coeff), byrow = TRUE)
  colSums(t_coeff >= t0_mat) / nrow(t_coeff)
}

# ---------------------------------------------------------------------------
# Statistics: real sign-flip tests (the first permutation is the identity, so
# T* ties T0), uniform values, discrete values (ties between sets summed in
# different orders), and special values
# ---------------------------------------------------------------------------
set.seed(1L)
cases <- list()
for (p in c(2L, 7L, 40L)) {
  for (n_perm in c(1L, 5L, 60L)) {
    x <- matrix(stats::rnorm(15L * p), 15L, p) +
      rep(seq(0, 1, length.out = p), each = 15L)
    t0 <- colMeans(x)^2
    t_coeff <- matrix(t0, nrow = 1L)
    for (b in seq_len(n_perm - 1L)) {
      signs <- sample(c(-1, 1), 15L, replace = TRUE)
      t_coeff <- rbind(t_coeff, colMeans(x * signs)^2)
    }
    label <- sprintf("p = %d, B = %d", p, n_perm)
    cases[[paste("sign flips,", label)]] <- list(t0, t_coeff, TRUE)
    cases[[paste("uniform,", label)]] <- list(
      stats::runif(p),
      matrix(stats::runif(n_perm * p), n_perm, p),
      TRUE
    )
    cases[[paste("discrete,", label)]] <- list(
      sample(0:3, p, replace = TRUE) / 7,
      matrix(sample(0:3, n_perm * p, replace = TRUE) / 7, n_perm, p),
      FALSE
    )
  }
}
t0_u <- stats::runif(12L)
t_u <- matrix(stats::runif(50L * 12L), 50L, 12L)
t0_nan <- t0_u
t0_nan[c(3L, 8L)] <- NaN
t_nan <- t_u
t_nan[c(5L, 70L, 300L)] <- NaN
t_inf <- t_u
t_inf[, 4L] <- Inf
t_inf[1:10, 9L] <- -Inf
cases[["NaN in t0"]] <- list(t0_nan, t_u, FALSE)
cases[["NaN in t_coeff"]] <- list(t0_u, t_nan, FALSE)
cases[["Inf and -Inf"]] <- list(t0_u, t_inf, FALSE)
cases[["near overflow"]] <- list(t0_u * 1e308, t_u * 1e308, FALSE)
cases[["all zero"]] <- list(rep(0, 12L), matrix(0, 50L, 12L), FALSE)
cases[["all tied"]] <- list(rep(0.5, 12L), matrix(0.5, 50L, 12L), FALSE)

for (label in names(cases)) {
  t0 <- cases[[label]][[1L]]
  t_coeff <- cases[[label]][[2L]]
  if (!exact && !cases[[label]][[3L]]) {
    next
  }
  p <- length(t0)
  n_perm <- nrow(t_coeff)
  pval <- pointwise_pvalues(t0, t_coeff)
  for (agg in c("integral", "max")) {
    for (recycle in c(FALSE, TRUE)) {
      info <- paste(label, agg, "recycle =", recycle)
      target <- suppressWarnings(
        iwt_reference(p, pval, t0, t_coeff, n_perm, recycle, agg)
      )
      current <- fdatest:::p_adjust_iwt(
        p,
        pval,
        t0,
        t_coeff,
        n_perm,
        recycle,
        agg
      )
      expect_same(current$pvalue_matrix, target$pvalue_matrix, info)
      expect_same(current$adjusted_pvalues, target$adjusted_pvalues, info)
    }
    info <- paste(label, agg)
    expect_same(
      fdatest:::p_adjust_twt(pval, p, t0, t_coeff, agg),
      suppressWarnings(twt_reference(pval, p, t0, t_coeff, agg)),
      paste("TWT,", info)
    )
    expect_same(
      fdatest:::p_adjust_global(agg, t0, t_coeff, p),
      global_reference(agg, t0, t_coeff, p),
      paste("Global,", info)
    )
  }
}

# ---------------------------------------------------------------------------
# One grid point: the IWT keeps the R computation
# ---------------------------------------------------------------------------
t0_1 <- 0.3
t_1 <- matrix(stats::runif(20L), 20L, 1L)
pval_1 <- pointwise_pvalues(t0_1, t_1)
expect_identical(
  fdatest:::p_adjust_iwt(1L, pval_1, t0_1, t_1, 20L, TRUE, "integral"),
  iwt_reference(1L, pval_1, t0_1, t_1, 20L, TRUE, "integral")
)

# ---------------------------------------------------------------------------
# The TWT keeps the names of the p-values
# ---------------------------------------------------------------------------
pval_named <- stats::setNames(pointwise_pvalues(t0_u, t_u), paste0("t", 1:12))
adj_named <- fdatest:::p_adjust_twt(pval_named, 12L, t0_u, t_u, "integral")
expect_identical(names(adj_named$adjusted_pvalues), names(pval_named))

# ---------------------------------------------------------------------------
# Permuted statistics given as a vector (a slice of the ANOVA or linear-model
# arrays with a single permutation) or as integers
# ---------------------------------------------------------------------------
t_row <- t_u[1L, , drop = FALSE]
pval_row <- pointwise_pvalues(t0_u, t_row)
expect_identical(
  fdatest:::p_adjust_twt(pval_row, 12L, t0_u, as.vector(t_row), "max"),
  fdatest:::p_adjust_twt(pval_row, 12L, t0_u, t_row, "max")
)
t_int <- matrix(sample(0:5, 30L * 6L, replace = TRUE), 30L, 6L)
t0_int <- c(1L, 4L, 2L, 5L, 0L, 3L)
expect_identical(
  fdatest:::p_adjust_global("integral", t0_int, t_int, 6L),
  global_reference("integral", t0_int, t_int, 6L)
)

# ---------------------------------------------------------------------------
# With mirai daemons the IWT uses the R computation, with the same results
# ---------------------------------------------------------------------------
pval_u <- pointwise_pvalues(t0_u, t_u)
res_cpp <- fdatest:::p_adjust_iwt(
  12L,
  pval_u,
  t0_u,
  t_u,
  50L,
  FALSE,
  "integral"
)
mirai::daemons(1L, seed = 42L, dispatcher = FALSE)
res_mirai <- fdatest:::p_adjust_iwt(
  12L,
  pval_u,
  t0_u,
  t_u,
  50L,
  FALSE,
  "integral"
)
mirai::daemons(0L)
expect_same(res_cpp, res_mirai, "mirai")

set.seed(NULL)
