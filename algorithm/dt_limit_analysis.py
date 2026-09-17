# -*- coding: utf-8 -*-
"""How does dt (as a free/tunable discretization parameter) affect the S
closed-form solution, and does shrinking it recover the "ground truth"
(continuous-time, order-independent) steady state?

Computes the continuous-time steady state directly (dR/dt=0, dP/dt=0 -- no
discrete stepping, hence no move/clearance/growth *order* at all), then sweeps
dt for both markov_move_first (SIRsimulator4.m order) and
markov_clearance_first (AgentBasedModel order), tracking:
  - relative error of each discrete closed form vs the continuous truth
  - the log-log slope of that error vs dt (should be ~1, i.e. O(dt))
  - a two-point Richardson extrapolation to dt=0

Uses README.md's no-clearance test case (real 209->418 mouse connectome).

Usage:
    python dt_limit_analysis.py
"""

import numpy as np
import scipy.sparse
import scipy.sparse.linalg

import abm_s_population as sp
from compare_move_order import build_transition_matrices, markov_move_first, markov_clearance_first

V = 104.87116708452713
SPREAD_RATE = 0.02122582417807315
PKL = "../SIR_inputs/params_regionalized_voxel_weights_after_qc_0413_64pct.pkl"
GENE_DIR = "/scratch/vnathan/sir_clearance_rebuilt_conn_fixk2_03262026/derivatives/yohan_ge_filt/mr10vv0.2/"


def continuous_time_steady_state(weights_raw, distance, region_size, clearance_rate, syngene):
    """Ground truth: solve dR/dt=0, dP/dt=0 directly. No dt anywhere, no
    move/clearance/growth ordering -- everything happens simultaneously, as
    it does in the real (continuous-time) process both discretizations
    approximate."""
    N_s, N_t = weights_raw.shape
    diag_idx = np.arange(N_s)

    omega = weights_raw / weights_raw.sum(axis=1, keepdims=True)   # region->edge RATE
    omega[diag_idx, diag_idx] = 0
    omega = omega * SPREAD_RATE

    sigma = V / distance                                          # edge->region RATE
    sigma[distance == 0] = 0
    sigma[weights_raw == 0] = 0

    kappa = np.broadcast_to(clearance_rate, (N_t,))                # clearance RATE
    gamma_rs = syngene * region_size                               # growth RATE

    # P*: sigma*P = omega*R  =>  P[i,j] = omega[i,j]*R[i]/sigma[i,j]  (0 where sigma=0,
    # and since sigma=0 implies omega=0 there too -- masked consistently)
    # substitute into dR/dt=0:
    #   0 = gamma_rs - kappa*R - outflow(R) + sum_i sigma[i,j]*P[i,j]
    #     = gamma_rs - kappa*R - outflow(R) + sum_i omega[i,j]*R[i]
    outflow = np.zeros(N_t)
    outflow[:N_s] = omega.sum(axis=1)
    A_RR = scipy.sparse.diags(kappa + outflow)

    rows_i, cols_j = np.meshgrid(np.arange(N_s), np.arange(N_t), indexing="ij")
    A_from_i = scipy.sparse.coo_matrix(
        (omega.flatten(order="F"), (cols_j.flatten(order="F"), rows_i.flatten(order="F"))),
        shape=(N_t, N_s),
    )
    A_from_i_full = scipy.sparse.hstack([A_from_i, scipy.sparse.csr_matrix((N_t, N_t - N_s))])

    A = (A_RR - A_from_i_full).tocsc()
    R = scipy.sparse.linalg.spsolve(A, gamma_rs)

    # recover P from R: P[i,j] = omega[i,j]*R[i]/sigma[i,j] where sigma>0
    P = np.zeros((N_s, N_t))
    mask = sigma > 0
    P[mask] = omega[mask] * np.repeat(R[:N_s], N_t).reshape(N_s, N_t)[mask] / sigma[mask]
    return R, P


def rel(a, b):
    a = np.ravel(a); b = np.ravel(b)
    return np.max(np.abs(a - b)) / max(np.abs(b).max(), 1e-300)


def main():
    sp.Snca_dir_default = GENE_DIR
    weights_raw, distance, region_size, sources, targets, syngene = sp.load_params(PKL)
    clearance_rate = sp.load_clearance(None)

    R_cont, P_cont = continuous_time_steady_state(weights_raw, distance, region_size, clearance_rate, syngene)
    print("Continuous-time (dt-independent) ground truth computed.")
    print(f"  R range: [{R_cont.min():.4g}, {R_cont.max():.4g}]\n")

    dts = [0.1, 0.03, 0.01, 0.003, 0.001, 3e-4, 1e-4, 3e-5, 1e-5, 1e-6]
    errs_move, errs_clear = [], []
    print(f"{'dt':>10} {'move_first err':>16} {'clearance_first err':>20}")
    for dt in dts:
        weights, sconnMov = build_transition_matrices(weights_raw, distance, SPREAD_RATE, V, dt)
        sconnMov[weights == 0] = 0
        clearance = np.exp(-clearance_rate * dt)
        growth = syngene * region_size * dt

        r_move, _ = markov_move_first(weights, sconnMov, clearance, growth)
        r_clear, _ = markov_clearance_first(weights, sconnMov, clearance, growth)

        e_move, e_clear = rel(r_move, R_cont), rel(r_clear, R_cont)
        errs_move.append(e_move); errs_clear.append(e_clear)
        print(f"{dt:>10.2e} {e_move*100:>15.4f}% {e_clear*100:>19.4f}%")

    # log-log slope (should be ~1.0 for O(dt) convergence)
    log_dt = np.log(dts)
    slope_move = np.polyfit(log_dt, np.log(errs_move), 1)[0]
    slope_clear = np.polyfit(log_dt, np.log(errs_clear), 1)[0]
    print(f"\nlog-log slope (error vs dt): move_first={slope_move:.3f}, "
          f"clearance_first={slope_clear:.3f}  (1.0 = exactly first-order in dt)")

    # Richardson extrapolation to dt=0 using the two smallest dt's, assuming O(dt) error:
    # R(dt) = R_true + c*dt  =>  R_true = (dt2*R(dt1) - dt1*R(dt2)) / (dt2-dt1)
    dt1, dt2 = dts[-2], dts[-1]
    w1, s1 = build_transition_matrices(weights_raw, distance, SPREAD_RATE, V, dt1)
    s1[w1 == 0] = 0
    r1, _ = markov_clearance_first(w1, s1, np.exp(-clearance_rate*dt1), syngene*region_size*dt1)
    w2, s2 = build_transition_matrices(weights_raw, distance, SPREAD_RATE, V, dt2)
    s2[w2 == 0] = 0
    r2, _ = markov_clearance_first(w2, s2, np.exp(-clearance_rate*dt2), syngene*region_size*dt2)
    r_extrap = (dt2 * r1 - dt1 * r2) / (dt2 - dt1)
    print(f"\nRichardson extrapolation (dt={dt1:.0e},{dt2:.0e}) vs continuous truth: "
          f"rel = {rel(r_extrap, R_cont):.3e}")


if __name__ == "__main__":
    main()
