"""M0: noise covariance 3x3 per station, whitening, orientation audit (specification v2, M0).

Why. The conditional p-value of Appendix AB (Corollary C.2) and the subspace statistic assume isotropic noise per
block; real horizontal and vertical noise powers differ. Whitening v -> L^{-1/2} v makes the noise of every station
isotropic. The radiation operator must be whitened with the same matrices, otherwise the subspace no longer contains
the signal: if v = G m + z with cov(z_i) = L_i, then W v = (W G) m + W z with W = blockdiag(L_i^{-1/2}) and
cov(W z) = I. ``whiten_operators`` returns W G and its orthonormal bases.

Orientation. A misoriented horizontal pair rotates the observed P particle motion about the vertical. The audit
compares, on catalogued events, the horizontal direction of the observed P vector with the predicted ray, using the
vertical component to fix the polarity (the vertical is not affected by a horizontal misorientation).
"""
import numpy as np

from ..radiation import basis, stack_Q

CHI2_3_MEDIAN = 2.365973884375338


def tyler_shape(X, iters=100, tol=1e-9):
    """Tyler's M-estimator of scatter of centred samples X (n, 3), trace normalised to 3. Distribution-free over
    elliptical laws and with breakdown 1/3: transients (a few large samples) do not move it much."""
    X = np.asarray(X, float)
    p = X.shape[1]
    S = np.eye(p)
    for _ in range(iters):
        Si = np.linalg.inv(S)
        d = np.einsum("ni,ij,nj->n", X, Si, X)
        d = np.maximum(d, 1e-300)
        S1 = (X.T * (p / d)) @ X / len(X)
        S1 *= p / np.trace(S1)
        if np.abs(S1 - S).max() < tol:
            S = S1
            break
        S = S1
    return S


def noise_covariance(MF, mask=None, step=1):
    """(ns, 3, 3) robust noise covariance of the matched-filter output MF (ns, 3, NT) per station: Tyler shape, scale
    fixed so that the median squared Mahalanobis norm is the chi2_3 median (consistent for Gaussian noise).
    mask (ns, NT) or (NT,): samples to use (windows without detections); step: decimation (matched-filter samples
    are correlated, a step of a few samples changes nothing but speed)."""
    MF = np.asarray(MF, float)
    ns, _, nt = MF.shape
    if mask is None:
        mask = np.ones((ns, nt), bool)
    mask = np.broadcast_to(np.asarray(mask, bool), (ns, nt))
    out = np.zeros((ns, 3, 3))
    for i in range(ns):
        X = MF[i][:, mask[i]][:, ::step].T
        if len(X) < 10:
            raise ValueError(f"station {i}: {len(X)} noise samples, need >= 10")
        S = tyler_shape(X)
        d = np.einsum("ni,ij,nj->n", X, np.linalg.inv(S), X)
        out[i] = S * np.median(d) / CHI2_3_MEDIAN
    return out


def inv_sqrt(L):
    """symmetric inverse square root of each (..., 3, 3) covariance."""
    ev, V = np.linalg.eigh(L)
    return np.einsum("...ij,...j,...kj->...ik", V, 1.0 / np.sqrt(ev), V)


def whiten(MF, L):
    """whitened traces: per station, L_i^{-1/2} applied to the 3-component vector of every sample."""
    W = inv_sqrt(L)
    return np.einsum("ijk,ikt->ijt", W, MF)


def whiten_operators(G, L):
    """whitened radiation operators W G (N, 6 ns, 6), rows (phase, station, component), and their stacked bases.
    The same L_i whitens the P and the S block of station i (one noise process per station)."""
    G = np.asarray(G, float)
    N, R6, _ = G.shape
    ns = R6 // 6
    W = inv_sqrt(L)
    Gw = np.einsum("ijk,nbikm->nbijm", W, G.reshape(N, 2, ns, 3, 6)).reshape(N, R6, 6)
    return Gw, stack_Q([basis(g) for g in Gw])


def orientation_audit(v_obs, n_pred, min_h=1e-3):
    """horizontal misorientation of one station from its P vectors on catalogued events.

    v_obs (n, 3): observed P vectors (matched-filter output read at the P arrival), frame (x, y, z) of the data;
    n_pred (n, 3): predicted P polarisation (unit ray direction at the station, same frame).
    Polarity is fixed by the vertical: s = sign(v_z n_z). The angle of each event is the signed angle from the
    horizontal part of n_pred to that of s v_obs (radians, counter-clockwise in the x-y plane).
    Returns (angle, reliability r, n used): angle = circular mean, r = mean resultant length in [0, 1] (1: every event
    agrees). Events with a horizontal part below min_h x norm, or a vertical of zero, are skipped."""
    v, n = np.asarray(v_obs, float), np.asarray(n_pred, float)
    s = np.sign(v[:, 2] * n[:, 2])
    vh, nh = v[:, :2] * s[:, None], n[:, :2]
    ok = (s != 0) & (np.linalg.norm(vh, axis=1) > min_h * np.linalg.norm(v, axis=1)) & \
         (np.linalg.norm(nh, axis=1) > min_h)
    if not ok.any():
        return float("nan"), 0.0, 0
    a = np.arctan2(vh[ok, 1], vh[ok, 0]) - np.arctan2(nh[ok, 1], nh[ok, 0])
    z = np.exp(1j * a).mean()
    return float(np.angle(z)), float(np.abs(z)), int(ok.sum())


def rotate_horizontals(MF, angle):
    """undo a horizontal misorientation: rotate the (x, y) components of every sample by -angle (radians)."""
    c, s = np.cos(angle), np.sin(angle)
    out = np.array(MF, float, copy=True)
    x, y = MF[..., 0, :], MF[..., 1, :]
    out[..., 0, :] = c * x + s * y
    out[..., 1, :] = -s * x + c * y
    return out
