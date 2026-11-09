// src/adjust.cpp
//
// Interval-wise (IWT), threshold-wise (TWT) and global adjustments of
// pointwise permutation tests, used by p_adjust_iwt(), p_adjust_twt() and
// p_adjust_global() when no mirai daemons are set.
//
// The results are identical to the R implementation: R's sum(), rowSums() and
// mean() accumulate in long double and round to double at the end, and the
// code below accumulates in long double in the same order and rounds at the
// same points. (This holds when the compiler's long double has R's long double
// precision, as on the usual builds of R; ldouble_digits_cpp() lets the tests
// check it.) The speed comes from the algorithms: O(p^2) per permutation for
// the IWT interval statistics instead of O(p^3), and O(p) for the TWT level
// sets instead of O(p^2).
//
// Ported from the fdatestcpp package (Matteo Fontana, GPL >= 3).

#include <cpp11.hpp>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

namespace {

typedef long double ldouble;

// Unit roundoff of long double, measured at run time (as R measures
// .Machine$longdouble.eps): under valgrind, or with a reduced x87 precision
// control word, long double arithmetic runs below its nominal precision.
ldouble ldouble_unit_roundoff(int* digits = nullptr) {
  volatile ldouble one = 1, u = 1, s;
  int d = 0;
  do {
    u = u / 2;
    ++d;
    s = one + u;
  } while (s != one);
  if (digits) *digits = d;
  return u;
}

// Rounding of a long double sum as R's sum() does it: totals beyond the double
// range become +/-Inf (rowSums() instead casts, as a plain (double) does).
double r_sum_round(ldouble s) {
  if (s > DBL_MAX) return R_PosInf;
  if (s < -DBL_MAX) return R_NegInf;
  return (double)s;
}

// mean() of a logical vector with cnt TRUE values out of n
double r_mean_count(int cnt, int n) {
  return (double)((ldouble)cnt / (ldouble)n);
}

// max() propagating NA, as R's max() without na.rm
double max_na(double a, double b) {
  if (ISNAN(a) || ISNAN(b)) return NA_REAL;
  return a > b ? a : b;
}

bool all_finite(const double* x, R_xlen_t n) {
  for (R_xlen_t i = 0; i < n; ++i)
    if (!R_FINITE(x[i])) return false;
  return true;
}

// ---------------------------------------------------------------------------
// IWT
//
// Interval (j, L) = grid points j, j + 1, ..., j + L - 1 (0-based, mod p when
// recycling), stored at row p - L + 1 and column j + 1 of the p-value matrix,
// as in iwt_compute_row(). Each interval statistic is accumulated from its
// start point, extending the interval one point at a time: a long double
// running sum rounded to double (what sum(..., na.rm = TRUE) does for t0 and
// rowSums(..., na.rm = TRUE) for t_coeff) or a running max. That is O(p^2) per
// permutation and reproduces R's rounding exactly.
//
// cnt[j * p + L - 1] = #{b : T*_b(j, L) >= T0(j, L)}. GENERAL = some input is
// not finite: NaN is dropped as na.rm = TRUE does (0 in sums, -Inf in max) and
// na flags comparisons that are NaN (Inf - Inf).

// Starts j..j+3 of one permutation in lockstep: four independent accumulation
// chains hide the latency of the long double add (or of max). Each chain does
// exactly the single-start arithmetic, so the results are unchanged.
template <bool MAX>
void iwt_four_starts(const double* x, const double* th, int* cnt, int p, int j,
                     bool recycle) {
  const double* x0 = x + j;
  const double *x1 = x0 + 1, *x2 = x0 + 2, *x3 = x0 + 3;
  const double* h0 = th + (size_t)j * p;
  const double *h1 = h0 + p, *h2 = h1 + p, *h3 = h2 + p;
  int* c0 = cnt + (size_t)j * p;
  int *c1 = c0 + p, *c2 = c1 + p, *c3 = c2 + p;
  // lengths common to the four starts
  const int Lc = recycle ? p : p - j - 3;
  if (MAX) {
    double m0 = R_NegInf, m1 = R_NegInf, m2 = R_NegInf, m3 = R_NegInf;
    for (int L = 0; L < Lc; ++L) {
      if (x0[L] > m0) m0 = x0[L];
      if (x1[L] > m1) m1 = x1[L];
      if (x2[L] > m2) m2 = x2[L];
      if (x3[L] > m3) m3 = x3[L];
      c0[L] += (m0 >= h0[L]);
      c1[L] += (m1 >= h1[L]);
      c2[L] += (m2 >= h2[L]);
      c3[L] += (m3 >= h3[L]);
    }
    if (!recycle) {  // starts j + 2, j + 1, j are 1, 2, 3 points longer
      for (int L = Lc; L < Lc + 1; ++L) {
        if (x2[L] > m2) m2 = x2[L];
        c2[L] += (m2 >= h2[L]);
      }
      for (int L = Lc; L < Lc + 2; ++L) {
        if (x1[L] > m1) m1 = x1[L];
        c1[L] += (m1 >= h1[L]);
      }
      for (int L = Lc; L < Lc + 3; ++L) {
        if (x0[L] > m0) m0 = x0[L];
        c0[L] += (m0 >= h0[L]);
      }
    }
  } else {
    ldouble s0 = 0.0L, s1 = 0.0L, s2 = 0.0L, s3 = 0.0L;
    for (int L = 0; L < Lc; ++L) {
      s0 += x0[L];
      s1 += x1[L];
      s2 += x2[L];
      s3 += x3[L];
      c0[L] += ((double)s0 >= h0[L]);
      c1[L] += ((double)s1 >= h1[L]);
      c2[L] += ((double)s2 >= h2[L]);
      c3[L] += ((double)s3 >= h3[L]);
    }
    if (!recycle) {
      for (int L = Lc; L < Lc + 1; ++L) {
        s2 += x2[L];
        c2[L] += ((double)s2 >= h2[L]);
      }
      for (int L = Lc; L < Lc + 2; ++L) {
        s1 += x1[L];
        c1[L] += ((double)s1 >= h1[L]);
      }
      for (int L = Lc; L < Lc + 3; ++L) {
        s0 += x0[L];
        c0[L] += ((double)s0 >= h0[L]);
      }
    }
  }
}

template <bool MAX, bool GENERAL>
void iwt_counts(const double* t0, const double* tc, int B, int p, bool recycle,
                std::vector<int>& cnt, std::vector<char>& na) {
  const double fill = MAX ? R_NegInf : 0.0;
  const size_t pp = (size_t)p * p;
  std::vector<double> x(2 * (size_t)p), th(pp, 0.0);
  cnt.assign(pp, 0);
  na.assign(GENERAL ? pp : 0, 0);

  // observed interval statistics
  for (int k = 0; k < p; ++k) {
    double v = t0[k];
    if (GENERAL && ISNAN(v)) v = fill;
    x[k] = v;
    x[k + p] = v;
  }
  for (int j = 0; j < p; ++j) {
    const int Lmax = recycle ? p : p - j;
    const double* xj = x.data() + j;
    double* thj = th.data() + (size_t)j * p;
    if (MAX) {
      double m = R_NegInf;
      for (int L = 0; L < Lmax; ++L) {
        if (xj[L] > m) m = xj[L];
        thj[L] = m;
      }
    } else {
      ldouble s = 0.0L;
      for (int L = 0; L < Lmax; ++L) {
        s += xj[L];
        thj[L] = r_sum_round(s);
      }
    }
  }

  // permuted interval statistics
  for (int b = 0; b < B; ++b) {
    if ((b & 255) == 0) cpp11::check_user_interrupt();
    for (int k = 0; k < p; ++k) {
      double v = tc[b + (R_xlen_t)k * B];
      if (GENERAL && ISNAN(v)) v = fill;
      x[k] = v;
      x[k + p] = v;
    }
    int j0 = 0;
    if (!GENERAL) {
      for (; j0 + 3 < p; j0 += 4) {
        iwt_four_starts<MAX>(x.data(), th.data(), cnt.data(), p, j0, recycle);
      }
    }
    for (int j = j0; j < p; ++j) {
      const int Lmax = recycle ? p : p - j;
      const double* xj = x.data() + j;
      const double* thj = th.data() + (size_t)j * p;
      int* cj = cnt.data() + (size_t)j * p;
      if (MAX) {
        double m = R_NegInf;
        for (int L = 0; L < Lmax; ++L) {
          if (xj[L] > m) m = xj[L];
          cj[L] += (m >= thj[L]);
        }
      } else {
        ldouble s = 0.0L;
        for (int L = 0; L < Lmax; ++L) {
          s += xj[L];
          const double v = (double)s;
          cj[L] += (v >= thj[L]);
          if (GENERAL && (ISNAN(v) || ISNAN(thj[L]))) na[(size_t)j * p + L] = 1;
        }
      }
    }
  }
}

}  // namespace

// Binary digits of long double as measured here (64 on x86_64, 53 where long
// double is double). The adjustments reproduce R exactly when this equals
// .Machine$longdouble.digits.
[[cpp11::register]]
int ldouble_digits_cpp() {
  int digits = 0;
  ldouble_unit_roundoff(&digits);
  return digits;
}

// IWT p-value matrix for intervals of length 2 to p (rows 1 to p - 1); row p
// (the pointwise p-values) is left NA for the caller to fill.
[[cpp11::register]]
cpp11::writable::doubles_matrix<> iwt_pvalue_matrix_cpp(
    cpp11::doubles t0, cpp11::doubles_matrix<> t_coeff, double n_perm,
    bool recycle, bool use_max) {
  const int p = t0.size(), B = t_coeff.nrow();
  const double* pt0 = REAL(t0);
  const double* ptc = REAL(t_coeff.data());
  const bool general = !(all_finite(pt0, p) && all_finite(ptc, (R_xlen_t)B * p));

  std::vector<int> cnt;
  std::vector<char> na;
  if (use_max) {
    if (general) {
      iwt_counts<true, true>(pt0, ptc, B, p, recycle, cnt, na);
    } else {
      iwt_counts<true, false>(pt0, ptc, B, p, recycle, cnt, na);
    }
  } else {
    if (general) {
      iwt_counts<false, true>(pt0, ptc, B, p, recycle, cnt, na);
    } else {
      iwt_counts<false, false>(pt0, ptc, B, p, recycle, cnt, na);
    }
  }

  cpp11::writable::doubles_matrix<> out(p, p);
  double* pm = REAL(out.data());
  std::fill(pm, pm + (size_t)p * p, NA_REAL);
  for (int j = 0; j < p; ++j) {
    const int Lmax = recycle ? p : p - j;
    for (int L = 2; L <= Lmax; ++L) {
      const size_t c = (size_t)j * p + (L - 1);
      pm[(size_t)(p - L) + (size_t)j * p] =
          (general && na[c]) ? NA_REAL : (double)cnt[c] / n_perm;
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// TWT
//
// For each threshold in c(0, sort(unique(pval)), 1), p_adjust_twt() tests the
// sub-level set {pval <= threshold} and the super-level set {pval >
// threshold}. With the grid points sorted by p-value (ties by index), these are
// prefixes and suffixes of the sorted order: prefixes [0, cb[h]) for h = 1..U
// and suffixes [cb[h], m) for h = 1..U - 1, where cb[h] ends the h-th group
// of equal p-values (the thresholds 0 and 1 only add empty sets or repeat
// these). The adjusted p-value of a point is the max of its own p-value and
// those of all sets containing it, NA propagating as in R's max().
//
// Set statistics. max: running max along the sorted order (exact). integral:
// R sums each set in grid order in long double; the sorted-order running sum
// differs from that only by rounding, bounded by
//   |F - s| <= (2^-53 + 2 (m - 1) u) sum|t|,
// with u the unit roundoff of long double, so the comparison with T0 is decided
// by s unless |s - F0| is within twice that bound; then the set is summed again
// in grid order (exact). The bound assumes no overflow, so sets whose sum of
// |t| exceeds DBL_MAX / 4, or whose observed sum is not finite, also take the
// exact path, as do all sets when some input is not finite. Ties between T*
// and T0 are the only ordinary cases that reach it, so the cost stays O(B p).

[[cpp11::register]]
cpp11::writable::doubles twt_adjust_cpp(cpp11::doubles pval, cpp11::doubles t0,
                                        cpp11::doubles_matrix<> t_coeff,
                                        bool use_max) {
  const int p = t0.size(), B = t_coeff.nrow();
  const double* ppv = REAL(pval);
  const double* pt0 = REAL(t0);
  const double* ptc = REAL(t_coeff.data());

  // grid points with a p-value, sorted by p-value (ties by index)
  std::vector<int> ord;
  for (int k = 0; k < p; ++k) {
    if (!ISNAN(ppv[k])) ord.push_back(k);
  }
  std::stable_sort(ord.begin(), ord.end(),
                   [&](int a, int b) { return ppv[a] < ppv[b]; });
  const int m = (int)ord.size();
  std::vector<int> rank(p, -1);
  for (int pos = 0; pos < m; ++pos) rank[ord[pos]] = pos;
  std::vector<int> cb(1, 0);
  for (int pos = 1; pos < m; ++pos) {
    if (ppv[ord[pos]] != ppv[ord[pos - 1]]) cb.push_back(pos);
  }
  if (m > 0) cb.push_back(m);
  const int U = (int)cb.size() - 1;

  // sum over the sorted positions [lo, hi) in grid order, rounded as sum()
  // (observed) or rowSums() (permuted)
  auto exact_sum = [&](const double* v, int lo, int hi, bool observed) {
    ldouble s = 0.0L;
    for (int k = 0; k < p; ++k) {
      const int r = rank[k];
      if (r >= lo && r < hi) s += v[k];
    }
    return observed ? r_sum_round(s) : (double)s;
  };

  // observed set statistics
  std::vector<double> f0P(U + 1, 0.0), f0S(U + 1, 0.0);
  if (use_max) {
    double mx = R_NegInf;
    int pos = 0;
    for (int h = 1; h <= U; ++h) {
      for (; pos < cb[h]; ++pos) {
        if (pt0[ord[pos]] > mx) mx = pt0[ord[pos]];
      }
      f0P[h] = mx;
    }
    mx = R_NegInf;
    pos = m;
    for (int h = U - 1; h >= 1; --h) {
      while (pos > cb[h]) {
        --pos;
        if (pt0[ord[pos]] > mx) mx = pt0[ord[pos]];
      }
      f0S[h] = mx;
    }
  } else {
    for (int h = 1; h <= U; ++h) f0P[h] = exact_sum(pt0, 0, cb[h], true);
    for (int h = 1; h < U; ++h) f0S[h] = exact_sum(pt0, cb[h], m, true);
  }

  const bool finite = all_finite(pt0, p) && all_finite(ptc, (R_xlen_t)B * p);
  const bool fast_sum = !use_max && finite;
  // twice the rounding bound: 2^-52 + 4 p u
  const ldouble tol_k =
      (ldouble)std::ldexp(1.0, -52) + 4 * (ldouble)p * ldouble_unit_roundoff();

  std::vector<int> cP(U + 1, 0), cS(U + 1, 0);
  std::vector<char> naP(U + 1, 0), naS(U + 1, 0);
  std::vector<double> r(p);

  // comparison of a set's sorted-order sum s (sum of |t| = a) with F0
  const ldouble big = (ldouble)DBL_MAX / 4;
  auto guarded = [&](ldouble s, ldouble a, double F0, int lo, int hi) -> int {
    if (R_FINITE(F0) && a <= big) {
      const ldouble d = s - (ldouble)F0;
      const ldouble tol = tol_k * a;
      if (d > tol) return 1;
      if (-d > tol) return 0;
      if (a == 0.0L) return 0.0 >= F0;  // every term is zero
    }
    return exact_sum(r.data(), lo, hi, false) >= F0;
  };

  for (int b = 0; b < B; ++b) {
    if ((b & 255) == 0) cpp11::check_user_interrupt();
    for (int k = 0; k < p; ++k) r[k] = ptc[b + (R_xlen_t)k * B];

    if (use_max) {
      double mx = R_NegInf;
      int pos = 0;
      for (int h = 1; h <= U; ++h) {
        for (; pos < cb[h]; ++pos) {
          if (r[ord[pos]] > mx) mx = r[ord[pos]];
        }
        cP[h] += (mx >= f0P[h]);
      }
      mx = R_NegInf;
      pos = m;
      for (int h = U - 1; h >= 1; --h) {
        while (pos > cb[h]) {
          --pos;
          if (r[ord[pos]] > mx) mx = r[ord[pos]];
        }
        cS[h] += (mx >= f0S[h]);
      }
    } else if (fast_sum) {
      ldouble s = 0.0L, a = 0.0L;
      int pos = 0;
      for (int h = 1; h <= U; ++h) {
        for (; pos < cb[h]; ++pos) {
          const double v = r[ord[pos]];
          s += v;
          a += std::fabs(v);
        }
        cP[h] += guarded(s, a, f0P[h], 0, cb[h]);
      }
      s = 0.0L;
      a = 0.0L;
      pos = m;
      for (int h = U - 1; h >= 1; --h) {
        while (pos > cb[h]) {
          --pos;
          const double v = r[ord[pos]];
          s += v;
          a += std::fabs(v);
        }
        cS[h] += guarded(s, a, f0S[h], cb[h], m);
      }
    } else {
      for (int h = 1; h <= U; ++h) {
        const double F = exact_sum(r.data(), 0, cb[h], false);
        if (ISNAN(F) || ISNAN(f0P[h])) {
          naP[h] = 1;
        } else {
          cP[h] += (F >= f0P[h]);
        }
      }
      for (int h = 1; h < U; ++h) {
        const double F = exact_sum(r.data(), cb[h], m, false);
        if (ISNAN(F) || ISNAN(f0S[h])) {
          naS[h] = 1;
        } else {
          cS[h] += (F >= f0S[h]);
        }
      }
    }
  }

  std::vector<double> pvP(U + 1, NA_REAL), pvS(U + 1, NA_REAL);
  for (int h = 1; h <= U; ++h) pvP[h] = naP[h] ? NA_REAL : r_mean_count(cP[h], B);
  for (int h = 1; h < U; ++h) pvS[h] = naS[h] ? NA_REAL : r_mean_count(cS[h], B);

  // group g (1..U) lies in the prefixes h >= g and in the suffixes h <= g - 1
  std::vector<double> preMax(U + 2, R_NegInf), sufMax(U + 1, R_NegInf);
  for (int h = U; h >= 1; --h) preMax[h] = max_na(preMax[h + 1], pvP[h]);
  for (int g = 2; g <= U; ++g) sufMax[g] = max_na(sufMax[g - 1], pvS[g - 1]);
  cpp11::writable::doubles adjusted(p);
  for (int k = 0; k < p; ++k) adjusted[k] = ppv[k];
  for (int g = 1; g <= U; ++g) {
    const double pg = ppv[ord[cb[g - 1]]];
    const double ag = max_na(max_na(pg, preMax[g]), sufMax[g]);
    for (int pos = cb[g - 1]; pos < cb[g]; ++pos) adjusted[ord[pos]] = ag;
  }
  return adjusted;
}

// ---------------------------------------------------------------------------
// Global: sum(t0) against rowSums(t_coeff), or max() against the row maxima,
// without na.rm, so any NaN makes the p-value NA.

[[cpp11::register]]
double global_pvalue_cpp(cpp11::doubles t0, cpp11::doubles_matrix<> t_coeff,
                         bool use_max) {
  const int p = t0.size(), B = t_coeff.nrow();
  const double* pt0 = REAL(t0);
  const double* ptc = REAL(t_coeff.data());

  bool bad = false;
  int c = 0;
  if (use_max) {
    double F0 = R_NegInf;
    for (int k = 0; k < p; ++k) {
      if (ISNAN(pt0[k])) {
        bad = true;
      } else if (pt0[k] > F0) {
        F0 = pt0[k];
      }
    }
    std::vector<double> mx(B, R_NegInf);
    for (int k = 0; k < p; ++k) {
      const double* col = ptc + (R_xlen_t)k * B;
      for (int b = 0; b < B; ++b) {
        if (ISNAN(col[b])) {
          bad = true;
        } else if (col[b] > mx[b]) {
          mx[b] = col[b];
        }
      }
    }
    if (!bad) {
      for (int b = 0; b < B; ++b) c += (mx[b] >= F0);
    }
  } else {
    ldouble s0 = 0.0L;
    for (int k = 0; k < p; ++k) s0 += pt0[k];
    const double F0 = r_sum_round(s0);
    std::vector<ldouble> acc(B, 0.0L);
    for (int k = 0; k < p; ++k) {
      const double* col = ptc + (R_xlen_t)k * B;
      for (int b = 0; b < B; ++b) acc[b] += col[b];
    }
    bad = ISNAN(F0);
    for (int b = 0; b < B; ++b) {
      const double F = (double)acc[b];
      if (ISNAN(F)) bad = true;
      c += (F >= F0);
    }
  }
  return bad ? NA_REAL : r_mean_count(c, B);
}
