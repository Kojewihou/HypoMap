"""
This script implements the Van Elteren test, a stratified version of the Mann-Whitney U test, which is used to 
compare two independent samples while accounting for stratification in the data. The script includes two main 
functions: `VanElterenCore` for computing the U statistic and tie-corrected variance for each stratum, and 
`VanElterenTest` for aggregating the results across multiple strata and calculating the test statistic.

Overview:
- `VanElterenCore`: Computes the centered U statistic and variance for a single stratum using the Mann-Whitney U test.
- `VanElterenTest`: Applies the Van Elteren test across multiple strata and computes the effect size and p-values.

Dependencies:
- numpy: For numerical operations.
- scipy.sparse: For sparse matrix operations (not used directly but imported for potential input types).
- scanpy: For handling single-cell data structures (not used directly but assumed to be useful for input data).
- scipy.stats: For calculating chi-squared distribution and ranking data.

Functions:
- VanElterenCore(x: np.ndarray, y: np.ndarray, axis: int = 0) -> Tuple[np.ndarray, np.ndarray, int, int]: 
  Computes the U statistic and tie-corrected variance for the given input arrays x and y in a single stratum.
- VanElterenTest(X, Y, X_strata, Y_strata) -> Tuple[np.ndarray, np.ndarray]: 
  Aggregates results across strata to compute the Van Elteren test statistic, effect sizes, and p-values.

Example usage:
    fs, pvalues = VanElterenTest(X, Y, X_strata, Y_strata)
"""

from typing import Tuple
import numpy as np
import scipy.sparse as sp
import scanpy as sc
from scipy.stats import chi2
from scipy.stats._mannwhitneyu import _broadcast_concatenate, _rankdata

def VanElterenCore(x: np.ndarray, y: np.ndarray, axis: int = 0) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """
    Computes the U statistic and tie-corrected variance for a single stratum using the Mann-Whitney U test.

    Parameters
    ----------
    x : np.ndarray
        Data array for group X.
    
    y : np.ndarray
        Data array for group Y.
    
    axis : int, optional
        Axis along which the comparison is performed. Default is 0.
    
    Returns
    -------
    U : np.ndarray
        The centered and continuity-corrected Mann-Whitney U statistic.
    
    V : np.ndarray
        Tie-corrected variance of the U statistic.
    
    n1 : int
        Sample size of group X.
    
    n2 : int
        Sample size of group Y.
    """
    # Concatenate and broadcast x and y arrays for ranking
    x, y, xy = _broadcast_concatenate(x, y, axis)
    n1, n2 = x.shape[-1], y.shape[-1]

    # Compute ranks with ties
    ranks, t = _rankdata(xy, 'average', return_ties=True)
    R1 = ranks[..., :n1].sum(axis=-1)
    
    # Mann-Whitney U statistic for X
    U = R1 - n1 * (n1 + 1) / 2
    
    # Mean U
    mu = n1 * n2 / 2

    # Total sample size
    n = n1 + n2
    
    # Tie-corrected variance of U
    tie_term = (t**3 - t).sum(axis=-1)
    V = n1 * n2 / 12 * ((n + 1) - tie_term / (n * (n - 1)))
    
    # Center U about zero and apply continuity correction
    U = U - mu
    U -= np.sign(U) * 0.5

    return U, V, n1, n2

def VanElterenTest(
    X: np.ndarray,
    Y: np.ndarray,
    X_strata: np.ndarray,
    Y_strata: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Performs the Van Elteren test, a stratified version of the Mann-Whitney U test, to compare two groups while 
    accounting for stratification in the data.

    Parameters
    ----------
    X : np.ndarray
        Data array for group X (reference group).
    
    Y : np.ndarray
        Data array for group Y (comparison group).
    
    X_strata : np.ndarray
        Strata assignments for group X.
    
    Y_strata : np.ndarray
        Strata assignments for group Y.
    
    Returns
    -------
    fs : np.ndarray
        Effect sizes for each feature.
    
    pvalues : np.ndarray
        P-values for each feature based on the Van Elteren statistic.
    """
    
    unique_strata = np.unique(X_strata)
    
    Us = np.zeros(shape=(unique_strata.size, X.shape[1]))  # Centered and continuity-corrected Mann-Whitney U statistic
    Vs = np.zeros(shape=(unique_strata.size, X.shape[1]))  # Tie-corrected variance
    Ms = np.zeros(shape=(unique_strata.size, 1))  # Sample size of group X per stratum
    Ns = np.zeros(shape=(unique_strata.size, 1))  # Sample size of group Y per stratum
    
    for i, stratum in enumerate(unique_strata):
        X_stratum_mask = X_strata == stratum
        Y_stratum_mask = Y_strata == stratum

        # Skip strata with insufficient sample size
        if (sum(X_stratum_mask) < 5) or (sum(Y_stratum_mask) < 5):
            continue

        sX = X[X_stratum_mask]
        sY = Y[Y_stratum_mask]
        
        # Compute U and V for the stratum
        Us[i], Vs[i], Ms[i], Ns[i] = VanElterenCore(sX, sY, axis=0)
    
    # Calculate Van Elteren statistic and effect size
    with np.errstate(divide='ignore', invalid='ignore'):
        w = 1 / (Ms + Ns + 1)  # Locally best weighting
        VEs = np.sum(w * Us, axis=0)**2 / np.sum(w**2 * Vs, axis=0)  # Van Elteren Statistic
        fs = np.sum(w * Us, axis=0) / np.sum(w * (Ms + Ns), axis=0)  # Effect size
    
    # Compute p-values using chi-squared approximation
    pvalues = chi2.sf(VEs, df=1)
    
    return fs, pvalues
