// src/permutation.cpp
//
// Pointwise permutation tests of os_permtest() and ts_permtest() when no
// mirai daemons are set and the data are finite: the observed statistics and
// the n_perm permuted ones, with the same random numbers and the same results
// as os_single_perm() and ts_single_perm().
//
// - Signs and swaps: stats::rbinom(n, 1, 0.5) draws one uniform per value and
//   returns 1 when it is at least 0.5 (inversion branch of R's rbinom() for
//   size 1, prob 0.5), so unif_rand() >= 0.5 gives the same draws.
// - Permutations: sample(n) as R's do_sample(), with R_unif_index().
// - Means: colMeans() accumulates in long double and divides before rounding
//   to double; the code below does the same, in the same order.
// - Variances: diag(stats::cov(x)) as R's cov.c computes it (two-pass mean,
//   then the long double sum of squared deviations divided by n - 1).
//
// Ported from the fdatestcpp package (Matteo Fontana, GPL >= 3).

#include <cpp11.hpp>
#include <R_ext/Random.h>

#include "matrix.h"

#include <cmath>
#include <vector>

using namespace cpp11::literals;

namespace {

typedef long double ldouble;

// Runs draw() between GetRNGstate() and PutRNGstate(), saving the state also
// when draw() is interrupted (as an interrupted R loop keeps the draws made so
// far). Both go through cpp11::safe: they can raise R errors (a malformed
// .Random.seed, a locked binding), which must not jump over C++ frames.
template <typename F>
void with_rng_state(F&& draw) {
  cpp11::safe[GetRNGstate]();
  try {
    draw();
  } catch (cpp11::unwind_exception& e) {
    // An R condition (interrupt, time limit) is on its way to its handler,
    // with its value in CAR(e.token). cpp11::safe reuses that token for other
    // void() calls and clears its CAR on return, so keep the value across
    // PutRNGstate(): otherwise tryCatch() and try() receive NULL.
    cpp11::sexp value(CAR(e.token));
    cpp11::safe[PutRNGstate]();
    SETCAR(e.token, value);
    throw;
  } catch (...) {
    cpp11::safe[PutRNGstate]();
    throw;
  }
  cpp11::safe[PutRNGstate]();
}

// ts_pointwise_stat() for a mean difference d; alt 0 = two.sided,
// 1 = greater, 2 = less:
//   sign_diff <- sign(meandiff); sign_diff[sign_diff == -1] <- 0
//   greater: (meandiff * sign_diff)^2, less: (meandiff * (sign_diff - 1))^2
double alt_stat(double d, int alt) {
  if (alt == 0) return d * d;
  const double sd = ISNAN(d) ? d : (d > 0 ? 1.0 : 0.0);
  const double a = (alt == 1) ? d * sd : d * (sd - 1.0);
  return a * a;
}

// Variance of v as diag(stats::cov()) gives it (cov.c, use = "everything")
double r_cov_var(const double* v, int n) {
  if (n <= 1) return NA_REAL;
  ldouble sum = 0.0L;
  for (int i = 0; i < n; ++i) sum += v[i];
  ldouble tmp = sum / n;
  if (R_FINITE((double)tmp)) {
    sum = 0.0L;
    for (int i = 0; i < n; ++i) sum += (v[i] - tmp);
    tmp = tmp + sum / n;
  }
  const ldouble xm = (double)tmp;
  sum = 0.0L;
  for (int i = 0; i < n; ++i) sum += (v[i] - xm) * (v[i] - xm);
  return (double)(sum / (R_xlen_t)(n - 1));
}

}  // namespace

// One-sample sign flips: t0 = colMeans(coeff)^2 and n_perm rows of
// colMeans(coeff * signs)^2, signs <- rbinom(n, 1, 0.5) * 2 - 1.
[[cpp11::register]]
cpp11::writable::list os_permtest_cpp(cpp11::doubles_matrix<> coeff, int n_perm) {
  const int n = coeff.nrow(), p = coeff.ncol();
  const double* x = REAL(coeff.data());
  cpp11::writable::doubles t0(p);
  cpp11::sexp t_coeff = new_double_matrix(n_perm, p);
  double* pt0 = REAL(t0);
  double* ptc = REAL(t_coeff);

  // interleaved: W[2 (k n + i)] = -x[i, k], W[2 (k n + i) + 1] = x[i, k], so
  // that the sign is an index offset (negation is exact: these are the values
  // of coeff * signs)
  std::vector<double> W(2 * (size_t)n * p);
  for (size_t c = 0; c < (size_t)n * p; ++c) {
    W[2 * c] = -x[c];
    W[2 * c + 1] = x[c];
  }
  std::vector<int> idx(n);
  auto finish = [&](ldouble s) {
    const double m = (double)(s / n);
    return m * m;
  };
  // four columns in lockstep: independent long double chains, same arithmetic
  auto stats = [&](double* out, R_xlen_t stride) {
    int k = 0;
    for (; k + 3 < p; k += 4) {
      const double* w0 = W.data() + 2 * (size_t)k * n;
      const double *w1 = w0 + 2 * n, *w2 = w1 + 2 * n, *w3 = w2 + 2 * n;
      ldouble s0 = 0.0L, s1 = 0.0L, s2 = 0.0L, s3 = 0.0L;
      for (int i = 0; i < n; ++i) {
        const int j = idx[i];
        s0 += w0[j];
        s1 += w1[j];
        s2 += w2[j];
        s3 += w3[j];
      }
      out[k * stride] = finish(s0);
      out[(k + 1) * stride] = finish(s1);
      out[(k + 2) * stride] = finish(s2);
      out[(k + 3) * stride] = finish(s3);
    }
    for (; k < p; ++k) {
      const double* w = W.data() + 2 * (size_t)k * n;
      ldouble s = 0.0L;
      for (int i = 0; i < n; ++i) s += w[idx[i]];
      out[k * stride] = finish(s);
    }
  };

  for (int i = 0; i < n; ++i) idx[i] = 2 * i + 1;
  stats(pt0, 1);
  with_rng_state([&] {
    for (int b = 0; b < n_perm; ++b) {
      if ((b & 255) == 0) cpp11::check_user_interrupt();
      for (int i = 0; i < n; ++i) idx[i] = 2 * i + (unif_rand() >= 0.5);
      stats(ptc + b, n_perm);
    }
  });

  return cpp11::writable::list({"t0"_nm = t0, "t_coeff"_nm = t_coeff});
}

// Two-sample tests: coeff = rbind(data1 - mu, data2), first n1 rows sample 1.
// Paired: if_perm <- rbinom(n1, 1, 0.5) swaps rows i and n1 + i. Unpaired:
// coeff[sample(n), ]. Statistic: ts_pointwise_stat() of the permuted rows.
[[cpp11::register]]
cpp11::writable::list ts_permtest_cpp(cpp11::doubles_matrix<> coeff, int n1,
                                      int n_perm, int alt, bool paired,
                                      bool standardize) {
  const int n = coeff.nrow(), p = coeff.ncol(), n2 = n - n1;
  const double* x = REAL(coeff.data());
  cpp11::writable::doubles t0(p);
  cpp11::sexp t_coeff = new_double_matrix(n_perm, p);
  double* pt0 = REAL(t0);
  double* ptc = REAL(t_coeff);

  if (paired && !standardize && n2 == n1) {
    // Means of the two samples from interleaved values: W[2 (k n1 + i)] =
    // x[i, k] (sample 1), W[2 (k n1 + i) + 1] = x[n1 + i, k] (sample 2), so
    // that a swap is an index offset
    std::vector<double> W(2 * (size_t)n1 * p);
    for (int k = 0; k < p; ++k) {
      for (int i = 0; i < n1; ++i) {
        const size_t c = (size_t)k * n1 + i;
        W[2 * c] = x[(size_t)k * n + i];
        W[2 * c + 1] = x[(size_t)k * n + n1 + i];
      }
    }
    std::vector<int> f(n1, 0);
    auto finish = [&](ldouble s1, ldouble s2) {
      const double m1 = (double)(s1 / n1), m2 = (double)(s2 / n2);
      return alt_stat(m1 - m2, alt);
    };
    // two columns in lockstep: four independent long double chains
    auto stats = [&](double* out, R_xlen_t stride) {
      int k = 0;
      for (; k + 1 < p; k += 2) {
        const double* w = W.data() + 2 * (size_t)k * n1;
        const double* v = w + 2 * (size_t)n1;
        ldouble a1 = 0.0L, a2 = 0.0L, b1 = 0.0L, b2 = 0.0L;
        for (int i = 0; i < n1; ++i) {
          const int i1 = 2 * i + f[i], i2 = 2 * i + 1 - f[i];
          a1 += w[i1];
          a2 += w[i2];
          b1 += v[i1];
          b2 += v[i2];
        }
        out[k * stride] = finish(a1, a2);
        out[(k + 1) * stride] = finish(b1, b2);
      }
      for (; k < p; ++k) {
        const double* w = W.data() + 2 * (size_t)k * n1;
        ldouble s1 = 0.0L, s2 = 0.0L;
        for (int i = 0; i < n1; ++i) {
          s1 += w[2 * i + f[i]];
          s2 += w[2 * i + 1 - f[i]];
        }
        out[k * stride] = finish(s1, s2);
      }
    };

    stats(pt0, 1);
    with_rng_state([&] {
      for (int b = 0; b < n_perm; ++b) {
        if ((b & 255) == 0) cpp11::check_user_interrupt();
        for (int i = 0; i < n1; ++i) f[i] = unif_rand() >= 0.5;
        stats(ptc + b, n_perm);
      }
    });
  } else {
    // rows of the permuted samples: g1 (sample 1), g2 (sample 2)
    std::vector<int> g1(n1), g2(n2), perm(paired ? 0 : n), work;
    std::vector<double> v1(n1), v2(n2);
    for (int i = 0; i < n1; ++i) g1[i] = i;
    for (int i = 0; i < n2; ++i) g2[i] = n1 + i;

    auto stats = [&](double* out, R_xlen_t stride) {
      for (int k = 0; k < p; ++k) {
        const double* col = x + (size_t)k * n;
        ldouble s1 = 0.0L, s2 = 0.0L;
        for (int i = 0; i < n1; ++i) {
          v1[i] = col[g1[i]];
          s1 += v1[i];
        }
        for (int i = 0; i < n2; ++i) {
          v2[i] = col[g2[i]];
          s2 += v2[i];
        }
        const double m1 = (double)(s1 / n1), m2 = (double)(s2 / n2);
        double t = alt_stat(m1 - m2, alt);
        if (standardize) {
          // sp <- ((n1 - 1) * s1 + (n - n1 - 1) * s2) / (n - 2): R rounds each
          // product; volatile keeps the compiler from fusing them into an FMA
          volatile double a1 = (double)(n1 - 1) * r_cov_var(v1.data(), n1);
          volatile double a2 = (double)(n2 - 1) * r_cov_var(v2.data(), n2);
          const double sp = (a1 + a2) / (double)(n - 2);
          t = t / sp;
        }
        out[k * stride] = t;
      }
    };

    stats(pt0, 1);
    with_rng_state([&] {
      for (int b = 0; b < n_perm; ++b) {
        if ((b & 255) == 0) cpp11::check_user_interrupt();
        if (paired) {
          // n2 >= n1: rows n1 + i with i >= n1 are never swapped
          for (int i = 0; i < n1; ++i) {
            const bool swap = unif_rand() >= 0.5;
            g1[i] = swap ? n1 + i : i;
            g2[i] = swap ? i : n1 + i;
          }
        } else {
          // sample(n) as R's do_sample(): x[i] = i; y[i] = x[j] + 1 and
          // x[j] = x[--n] with j = R_unif_index(n)
          work.resize(n);
          for (int i = 0; i < n; ++i) work[i] = i;
          int dn = n;
          for (int i = 0; i < n; ++i) {
            const int j = (int)R_unif_index((double)dn);
            perm[i] = work[j];
            work[j] = work[--dn];
          }
          for (int i = 0; i < n1; ++i) g1[i] = perm[i];
          for (int i = 0; i < n2; ++i) g2[i] = perm[n1 + i];
        }
        stats(ptc + b, n_perm);
      }
    });
  }

  return cpp11::writable::list({"t0"_nm = t0, "t_coeff"_nm = t_coeff});
}
