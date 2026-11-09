# Tests for the C++ permutation tests (src/permutation.cpp) used by
# os_permtest() and ts_permtest() when no mirai daemons are set and the data
# are finite: with the same seed they must draw the same random numbers and
# return the same statistics as os_single_perm() and ts_single_perm(), the R
# code that the mirai path still uses.
# jarl-ignore-file internal_function: test files legitimately access internal functions via :::

mirai::daemons(0L)

# The C++ code accumulates in long double as R does: results are identical
# when both have the same precision (the usual builds of R, including those
# where long double is double). Otherwise compare the statistics with a
# tolerance. The random numbers are the same in any case.
r_digits <- if (capabilities("long.double")) .Machine$longdouble.digits else 53L
exact <- identical(as.integer(r_digits), fdatest:::ldouble_digits_cpp())

# ---------------------------------------------------------------------------
# R implementations
# ---------------------------------------------------------------------------
pointwise_pvalues <- function(t0, t_coeff) {
  t0_mat <- matrix(t0, nrow = nrow(t_coeff), ncol = ncol(t_coeff), byrow = TRUE)
  colSums(t_coeff >= t0_mat) / nrow(t_coeff)
}

os_reference <- function(coeff, n_perm) {
  t0 <- fdatest:::os_pointwise_stat(coeff)
  perms <- lapply(seq_len(n_perm), function(.x) {
    fdatest:::os_single_perm(coeff, nrow(coeff))
  })
  t_coeff <- do.call(rbind, perms)
  list(t0 = t0, t_coeff = t_coeff, pval = pointwise_pvalues(t0, t_coeff))
}

# ts_single_perm(), with drop = FALSE for unpaired tests so that one evaluation
# point stays a matrix (the C++ code handles that case)
ts_single_perm_reference <- function(
  coeff,
  n1,
  alternative,
  paired,
  standardize
) {
  if (paired) {
    return(fdatest:::ts_single_perm(
      coeff,
      n1,
      alternative,
      paired,
      standardize
    ))
  }
  coeff_perm <- coeff[sample(nrow(coeff)), , drop = FALSE]
  fdatest:::ts_pointwise_stat(coeff_perm, n1, alternative, standardize)
}

ts_reference <- function(coeff, n1, n_perm, alternative, paired, standardize) {
  t0 <- fdatest:::ts_pointwise_stat(coeff, n1, alternative, standardize)
  perms <- lapply(seq_len(n_perm), function(.x) {
    ts_single_perm_reference(coeff, n1, alternative, paired, standardize)
  })
  t_coeff <- do.call(rbind, perms)
  list(t0 = t0, t_coeff = t_coeff, pval = pointwise_pvalues(t0, t_coeff))
}

# Run both from the same seed; compare the results and the RNG state after
expect_same_perms <- function(f_cpp, f_ref, seed, info) {
  set.seed(seed)
  current <- f_cpp()
  state_current <- .Random.seed
  set.seed(seed)
  target <- f_ref()
  state_target <- .Random.seed
  expect_identical(state_current, state_target, info = paste(info, "RNG"))
  if (exact) {
    expect_identical(current, target, info = info)
  } else {
    expect_equal(current$t0, target$t0, tolerance = 1e-8, info = info)
    expect_equal(current$t_coeff, target$t_coeff, tolerance = 1e-8, info = info)
  }
}

# ---------------------------------------------------------------------------
# One sample
# ---------------------------------------------------------------------------
set.seed(1)
for (n in c(1L, 2L, 15L)) {
  for (p in c(1L, 7L)) {
    x <- matrix(rnorm(n * p), n, p) + rep(seq(0, 1, length.out = p), each = n)
    colnames(x) <- paste0("g", seq_len(p))
    xi <- matrix(sample(-3:3, n * p, replace = TRUE), n, p) * 1
    info <- sprintf("one sample n = %d, p = %d", n, p)
    expect_same_perms(
      function() fdatest:::os_permtest(x, 40L),
      function() os_reference(x, 40L),
      seed = n + p,
      info = info
    )
    expect_same_perms(
      function() fdatest:::os_permtest(xi, 40L),
      function() os_reference(xi, 40L),
      seed = 2 * n + p,
      info = paste(info, "ties")
    )
  }
}

# ---------------------------------------------------------------------------
# Two samples: paired with n2 = n1 (swaps of means), paired with n2 > n1 and
# unpaired (permuted rows), all alternatives, with and without standardizing
# ---------------------------------------------------------------------------
sizes <- list(c(2L, 2L), c(8L, 8L), c(8L, 11L), c(11L, 8L))
for (nn in sizes) {
  for (p in c(1L, 3L)) {
    n1 <- nn[1]
    n2 <- nn[2]
    x <- rbind(
      matrix(rnorm(n1 * p, mean = 0.3), n1, p),
      matrix(rnorm(n2 * p), n2, p)
    )
    for (alternative in c("two.sided", "greater", "less")) {
      for (standardize in c(FALSE, TRUE)) {
        for (paired in c(FALSE, TRUE)) {
          if (paired && n2 < n1) {
            next
          }
          info <- sprintf(
            "two samples n1 = %d, n2 = %d, p = %d, %s, standardize = %s, paired = %s",
            n1,
            n2,
            p,
            alternative,
            standardize,
            paired
          )
          expect_same_perms(
            function() {
              fdatest:::ts_permtest(
                x,
                n1,
                30L,
                alternative,
                paired,
                standardize
              )
            },
            function() {
              ts_reference(x, n1, 30L, alternative, paired, standardize)
            },
            seed = n1 + n2 + p,
            info = info
          )
        }
      }
    }
  }
}

# A large mean relative to the spread: the variance depends on cov() rounding
# its long double mean to double before summing the squared deviations
x <- 1e8 + rbind(matrix(rnorm(30, sd = 1e-3), 6, 5), matrix(rnorm(35), 7, 5))
expect_same_perms(
  function() fdatest:::ts_permtest(x, 6L, 20L, "two.sided", FALSE, TRUE),
  function() ts_reference(x, 6L, 20L, "two.sided", FALSE, TRUE),
  seed = 12,
  info = "large mean, standardized"
)

# Names, a second sample of one row (NA variance) and constant columns
x <- rbind(matrix(rnorm(15), 3, 5), matrix(rnorm(5), 1, 5))
colnames(x) <- letters[1:5]
expect_same_perms(
  function() fdatest:::ts_permtest(x, 3L, 20L, "two.sided", FALSE, TRUE),
  function() ts_reference(x, 3L, 20L, "two.sided", FALSE, TRUE),
  seed = 7,
  info = "second sample of one row"
)
xc <- cbind(1, matrix(rnorm(60), 20, 3), 0)
expect_same_perms(
  function() fdatest:::ts_permtest(xc, 10L, 20L, "less", TRUE, TRUE),
  function() ts_reference(xc, 10L, 20L, "less", TRUE, TRUE),
  seed = 8,
  info = "two samples, constant columns"
)
expect_same_perms(
  function() fdatest:::os_permtest(xc, 20L),
  function() os_reference(xc, 20L),
  seed = 9,
  info = "one sample, constant columns"
)

# ---------------------------------------------------------------------------
# Other RNG kinds, including sample(n) with sample.kind = "Rounding"
# ---------------------------------------------------------------------------
old_kind <- RNGkind()
for (kind in list(
  c("Mersenne-Twister", "Inversion", "Rounding"),
  c("Wichmann-Hill", "Inversion", "Rejection"),
  c("L'Ecuyer-CMRG", "Inversion", "Rejection")
)) {
  suppressWarnings(RNGkind(kind[1], kind[2], kind[3]))
  x <- rbind(matrix(rnorm(40), 8, 5), matrix(rnorm(45), 9, 5))
  expect_same_perms(
    function() fdatest:::ts_permtest(x, 8L, 25L, "two.sided", FALSE, FALSE),
    function() ts_reference(x, 8L, 25L, "two.sided", FALSE, FALSE),
    seed = 3,
    info = paste("unpaired,", kind[1], kind[3])
  )
  expect_same_perms(
    function() fdatest:::os_permtest(x, 25L),
    function() os_reference(x, 25L),
    seed = 4,
    info = paste("one sample,", kind[1])
  )
}
suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))

# ---------------------------------------------------------------------------
# Data that the C++ code leaves to the R code
# ---------------------------------------------------------------------------
# Non-finite values
x <- rbind(matrix(rnorm(24), 6, 4), matrix(rnorm(28), 7, 4))
x[2, 3] <- NA
expect_same_perms(
  function() fdatest:::ts_permtest(x, 6L, 20L, "greater", FALSE, FALSE),
  function() ts_reference(x, 6L, 20L, "greater", FALSE, FALSE),
  seed = 5,
  info = "two samples with NA"
)
x[2, 3] <- Inf
expect_same_perms(
  function() fdatest:::os_permtest(x, 20L),
  function() os_reference(x, 20L),
  seed = 6,
  info = "one sample with Inf"
)

# Without permutations the R code runs (and the tests fail, as before)
x <- rbind(matrix(rnorm(12), 4, 3), matrix(rnorm(9), 3, 3))
expect_null(suppressWarnings(fdatest:::os_permtest(x, 0L))$t_coeff)
expect_null(
  suppressWarnings(
    fdatest:::ts_permtest(x, 4L, 0L, "two.sided", FALSE, FALSE)
  )$t_coeff
)
expect_error(suppressWarnings(iwt1(x, n_perm = 0L)))

# Paired samples with n2 < n1 stay on the R code, where ts_single_perm() fails
# as soon as a swap reaches beyond the second sample (with 100 permutations
# and one unmatched pair, all but surely)
expect_error(fdatest:::ts_permtest(x, 4L, 100L, "two.sided", TRUE, FALSE))

# Two-sample data with an empty second sample still fail, as in the R code
expect_error(fdatest:::ts_permtest(x, nrow(x), 10L, "two.sided", FALSE, FALSE))

# ---------------------------------------------------------------------------
# Flags are read as the R code reads them: 0 and 1 work, NA fails
# ---------------------------------------------------------------------------
x <- rbind(matrix(rnorm(24), 6, 4), matrix(rnorm(28), 7, 4))
set.seed(11)
with_numbers <- fdatest:::ts_permtest(x, 6L, 20L, "less", 1, 1)
set.seed(11)
with_logicals <- fdatest:::ts_permtest(x, 6L, 20L, "less", TRUE, TRUE)
expect_identical(with_numbers, with_logicals)
expect_error(fdatest:::ts_permtest(x, 6L, 20L, "less", NA, FALSE))
expect_error(fdatest:::ts_permtest(x, 6L, 20L, "less", FALSE, NA))

# ---------------------------------------------------------------------------
# A malformed .Random.seed fails as in R, and the error leaves the C++ code
# cleanly
# ---------------------------------------------------------------------------
old_seed <- get(".Random.seed", envir = globalenv())
assign(".Random.seed", c(10403L, 1L), envir = globalenv())
expect_error(fdatest:::os_permtest(x, 5L), pattern = "Random.seed")
assign(".Random.seed", c(10403L, 1L), envir = globalenv())
expect_error(
  fdatest:::ts_permtest(x, 6L, 5L, "two.sided", FALSE, FALSE),
  pattern = "Random.seed"
)
assign(".Random.seed", old_seed, envir = globalenv())

# ---------------------------------------------------------------------------
# Numbers of permutations that the C++ code does not take go to the R code,
# which handles them as before (several numbers: the first one is used)
# ---------------------------------------------------------------------------
x <- matrix(rnorm(40), 10, 4)
expect_error(fdatest:::os_permtest(x, NA))
expect_error(fdatest:::ts_permtest(x, 5L, NA, "two.sided", FALSE, FALSE))
expect_equal(
  nrow(suppressWarnings(fdatest:::os_permtest(x, c(3L, 5L)))$t_coeff),
  3L
)
expect_equal(
  nrow(
    suppressWarnings(
      fdatest:::ts_permtest(x, 5L, c(3L, 5L), "two.sided", FALSE, FALSE)
    )$t_coeff
  ),
  3L
)

# ---------------------------------------------------------------------------
# A condition raised while the C++ loop runs (here a time limit, as an
# interrupt would be) reaches tryCatch() intact
# ---------------------------------------------------------------------------
x <- matrix(rnorm(1e5), 1e5, 1)
for (f in list(
  function() fdatest:::os_permtest(x, 1e5L),
  function() fdatest:::ts_permtest(x, 5e4L, 1e5L, "two.sided", FALSE, FALSE),
  function() fdatest:::ts_permtest(x, 5e4L, 1e5L, "two.sided", TRUE, FALSE)
)) {
  caught <- tryCatch(
    tryCatch(
      {
        setTimeLimit(elapsed = 0.2, transient = TRUE)
        f()
      },
      error = function(e) "handled"
    ),
    error = function(e) conditionMessage(e)
  )
  setTimeLimit()
  expect_identical(caught, "handled")
}
