"""
Repeated measures correlation implementation.
Based on: Bakdash JZ, Marusich LR (2017). "Repeated Measures Correlation."
"""

from __future__ import annotations

# -- Imports ---
import numpy as np
from numba import njit
from scipy import stats
import pandas as pd


# -- Numba helper functions ---
@njit
def _per_subject_sums_counts(subj_idx, x, y, K):
    """Compute per-subject sums and counts."""
    sums_x = np.zeros(K)
    sums_y = np.zeros(K)
    counts = np.zeros(K, dtype=np.int64)
    N = x.shape[0]
    for i in range(N):
        k = subj_idx[i]
        sums_x[k] += x[i]
        sums_y[k] += y[i]
        counts[k] += 1
    return sums_x, sums_y, counts


@njit
def _within_sums(subj_idx, x, y, mean_x, mean_y):
    """Compute pooled within-subject sums of squares and cross-products."""
    Sxx = 0.0
    Syy = 0.0
    Sxy = 0.0
    N = x.shape[0]
    for i in range(N):
        k = subj_idx[i]
        dx = x[i] - mean_x[k]
        dy = y[i] - mean_y[k]
        Sxx += dx * dx
        Syy += dy * dy
        Sxy += dx * dy
    return Sxx, Syy, Sxy


@njit
def _drop_nans(x, y, subject):
    """Drop entries with NaNs from x, y, and subject arrays."""
    valid_idx = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(subject)
    return x[valid_idx], y[valid_idx], subject[valid_idx]


# -- Main function ---
def rmcorr_stats(data: pd.DataFrame, x: str, y: str, subject: str, alpha: float = 0.05):
    """
    Compute repeated measures correlation statistics.

    Parameters
    ----------
    data : pandas.DataFrame
        Input dataframe containing the variables and subject identifier.
    x, y, subject : str
        Column names in ``data`` for the repeated measurements and subject IDs.
    alpha : float, default=0.05
        Significance level used for the confidence interval.

    Returns
    -------
    dict
        Dictionary with r, dof, p-value, CI95%, and post-hoc power.
    """
    x = np.asarray(data[x].values, dtype=np.float64)
    y = np.asarray(data[y].values, dtype=np.float64)
    subject = np.asarray(data[subject].values)

    x, y, subject = _drop_nans(x, y, subject)

    if x.ndim != 1 or y.ndim != 1 or subject.ndim != 1:
        raise ValueError("x, y, and subject must be 1D arrays.")
    if len(x) != len(y) or len(x) != len(subject):
        raise ValueError("x, y, and subject must have the same length.")

    unique_subjects, subj_idx = np.unique(subject, return_inverse=True)
    K = unique_subjects.shape[0]
    N = x.shape[0]

    if K < 2:
        raise ValueError("At least two subjects are required.")
    if N - K - 1 <= 0:
        raise ValueError("Insufficient data: need N - K - 1 > 0 for valid degrees of freedom.")

    # Per-subject means
    sums_x, sums_y, counts = _per_subject_sums_counts(subj_idx.astype(np.int64), x, y, K)
    mean_x = sums_x / counts
    mean_y = sums_y / counts

    # Pooled within-subject sums
    Sxx, Syy, Sxy = _within_sums(subj_idx.astype(np.int64), x, y, mean_x, mean_y)

    if Sxx == 0.0 or Syy == 0.0:
        raise ValueError("Zero within-person variance in x or y; rmcorr undefined.")

    r = Sxy / np.sqrt(Sxx * Syy)

    # Degrees of freedom
    dof = N - K - 1

    # t-statistic and p-value
    r_clamped = np.clip(r, -0.999999999, 0.999999999)
    t = r_clamped * np.sqrt(dof / (1.0 - r_clamped * r_clamped))
    pval = 2.0 * stats.t.sf(np.abs(t), dof)

    # Fisher z CI95%
    zr = np.arctanh(r_clamped)
    se = 1.0 / np.sqrt(dof)
    zcrit = stats.norm.ppf(1.0 - alpha / 2.0)
    lo = zr - zcrit * se
    hi = zr + zcrit * se
    ci95 = [np.tanh(lo), np.tanh(hi)]

    # Post-hoc power
    f2 = r_clamped ** 2 / (1.0 - r_clamped ** 2)
    ncp = np.sqrt(f2 * dof)
    crit = stats.t.ppf(1.0 - alpha / 2.0, dof)
    power = stats.nct.sf(crit, dof, ncp) + stats.nct.cdf(-crit, dof, ncp)

    return {
        "r": float(r),
        "dof": int(dof),
        "pval": float(pval),
        "CI95%": [float(ci95[0]), float(ci95[1])],
        "power": float(power),
    }


if __name__ == "__main__":  # pragma: no cover - usage example
    from pingouin import rm_corr

    rng = np.random.default_rng()
    subj = np.repeat(np.arange(5), 6)
    x = rng.normal(size=subj.size)
    y = 0.7 * (x - np.repeat([0, 1, 2, 3, 4], 6)) + rng.normal(scale=0.5, size=subj.size)
    data = pd.DataFrame({"x": x, "y": y, "subject": subj})

    res_custom = rmcorr_stats(data, "x", "y", "subject")
    print("Custom rmcorr results:")
    for k, v in res_custom.items():
        print(f"{k} : {v}")

    res_pg = rm_corr(data, x="x", y="y", subject="subject")
    print("\nPingouin results:")
    for k, v in res_pg.items():
        print(f"{k} : {v}")
