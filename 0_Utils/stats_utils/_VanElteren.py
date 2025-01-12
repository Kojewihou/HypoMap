import numpy as np

from scipy.stats import chi2
from scipy.stats._mannwhitneyu import _broadcast_concatenate, _rankdata

def van_elteren_core(x, y, axis=0):
    x, y, xy = _broadcast_concatenate(x, y, axis)
    M, N = x.shape[-1], y.shape[-1]

    ranks, t = _rankdata(xy, 'average', return_ties=True)  
    R1 = ranks[..., :M].sum(axis=-1)
    U = R1 - M * (M + 1) / 2 # Mann-Whitney U statistic for X, symmetric to U stat for Y so works the same
    
    mu = M * N / 2  # Mean U
    
    n = M + N
    tie_term = (t**3 - t).sum(axis=-1)
    V = M * N / 12 * ((n + 1) - tie_term / (n * (n - 1))) # Tie corrected variance
    
    U = U - mu # Center U about zero
    U -= np.sign(U) * 0.5 # Continuity correction 

    return U, V, M, N

def van_elteren_test(X, Y, X_strata, Y_strata, return_effect_size=False):
    
    assert X.shape[-1] == Y.shape[-1], 'Same number of features expected.'
    
    unique_strata = np.unique(X_strata)
    n_strata = len(unique_strata)
    
    Us = np.full(shape=(n_strata, X.shape[1]), fill_value=np.nan)
    Vs = np.full(shape=(n_strata, X.shape[1]), fill_value=np.nan)
    Ms = np.full(shape=(n_strata, 1), fill_value=np.nan)
    Ns = np.full(shape=(n_strata, 1), fill_value=np.nan)

    for i, stratum in enumerate(unique_strata):
        X_stratum_mask = X_strata == stratum
        Y_stratum_mask = Y_strata == stratum
        
        if not (sum(X_stratum_mask) >= 10 and sum(Y_stratum_mask) >= 10):
            continue
        
        sX = X[X_stratum_mask]
        sY = Y[Y_stratum_mask]

        Us[i], Vs[i], Ms[i], Ns[i] = van_elteren_core(sX, sY)
    
    # Remove skipped strata from formulation
    valid_strata_mask = ~np.isnan(Ms[:,0])
    Us = Us[valid_strata_mask]
    Vs = Vs[valid_strata_mask]
    Ms = Ms[valid_strata_mask]
    Ns = Ns[valid_strata_mask]
    
    w = 1 / (Ms + Ns + 1) # locally best weighting
    # w = 1 / (Ms*Ns) # design free weighting
    with np.errstate(divide='ignore', invalid='ignore'):
        VEs = np.sum(w * Us, axis=0)**2 / np.sum(w**2 * Vs, axis=0) # Van Elteren Statistic
    
    pvalues = chi2.sf(VEs, df=1) # Van Elteren statistic can be approximated by a chi-squared distribution with 1 degree of freedom
    # valid_pvalues = ~np.isnan(pvalues) # If the variance is 0, VEs is nan therefore pvalue is nan
    # pvalues[valid_pvalues] = false_discovery_control(pvalues[valid_pvalues])
    
    if return_effect_size:
        with np.errstate(divide='ignore', invalid='ignore'):
            effect_size = np.sum(w * Us, axis=0) / np.sum(w * (Ms + Ns), axis=0) # Effect Size
        
        return effect_size, pvalues
    else:
        return pvalues