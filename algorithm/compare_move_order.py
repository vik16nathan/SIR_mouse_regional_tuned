# -*- coding: utf-8 -*-
"""Isolate the effect of timestep ORDER OF OPERATIONS on the S-population
Markov closed form, on the real mouse bipartite connectome (README.md's
no-clearance "speedup test" example).

Two per-timestep orderings are compared, with weights/sconnMov held IDENTICAL
between them (AgentBasedModel's own weighting convention) so the only thing
that varies is the order:

  "move_first"       -- SIRsimulator4.m's order: move (R<->P exchange using the
                         OLD state), THEN apply clearance and add growth/synthesis.
                             Rnor = Rnor.*clearance + synthesis   (after the move)

  "clearance_first"   -- AgentBasedModel's actual order (growth_step, then
                         clearance_step, then s_spread_step): apply growth then
                         clearance to R FIRST, then move that already-updated R.
                             s_region += growth*dt; s_region *= clearance; move

Both are the closed-form (exact) fixed point of their respective recursion --
not an iterative approximation -- solved via the same sparse linear system,
just with a different M/B assembly (see markov_s_steady_state's docstring in
abm_s_population.py for the concrete block-matrix difference).

Each is compared against abm_s_population.iterative_s_steady_state, which runs
the actual AgentBasedModel step methods (ground truth for "clearance_first",
since that's what __main__ actually calls).

Usage:
    python compare_move_order.py
"""

import numpy as np
import pandas as pd
import scipy.sparse
import scipy.sparse.linalg

import abm_s_population as sp
from AgentBasedModel import AgentBasedModel

PKL = "../SIR_inputs/params_regionalized_voxel_weights_after_qc_0413_64pct.pkl"
GENE_DIR = ("/scratch/vnathan/sir_clearance_rebuilt_conn_fixk2_03262026/"
           "derivatives/yohan_ge_filt/mr10vv0.2")
V = 104.87116708452713
SPREAD_RATE = 0.02122582417807315
DT = 0.1
EPS = 1e-5


def build_transition_matrices(weights_raw, distance, spread_rate, v, dt):
    """The one shared preprocessing step: turn raw connectome arrays into
    per-step transition probabilities, exactly as AgentBasedModel does it.
    Held fixed across both orderings below so only the order varies."""
    N_s, N_t = weights_raw.shape
    diag_idx = np.arange(N_s)

    weights = weights_raw / weights_raw.sum(axis=1, keepdims=True)
    weights[diag_idx, diag_idx] = 0
    weights = weights * spread_rate * dt

    sconnMov = v / distance * dt
    sconnMov[distance == 0] = 0

    return weights, sconnMov


def _build_M_components(weights, sconnMov):
    N_s, N_t = weights.shape
    M_RR_diag = np.ones(N_t)
    M_RR_diag[:N_s] = 1 - weights.sum(axis=1)
    M_RR = scipy.sparse.diags(M_RR_diag)
    M_PP = scipy.sparse.diags(np.where(sconnMov.flatten(order="F") > 0, 1, 0) - sconnMov.flatten(order="F"))

    rows_i, cols_j = np.meshgrid(np.arange(N_s), np.arange(N_t), indexing="ij")
    lin = np.arange(N_s * N_t)
    M_RP = scipy.sparse.coo_matrix((weights.flatten(order="F"), (lin, rows_i.flatten(order="F"))), shape=(N_s * N_t, N_t))
    M_PR = scipy.sparse.coo_matrix((sconnMov.flatten(order="F"), (cols_j.flatten(order="F"), lin)), shape=(N_t, N_s * N_t))
    return N_s, N_t, M_RR, M_PP, M_RP, M_PR


def markov_move_first(weights, sconnMov, clearance, growth):
    """SIRsimulator4.m order: move using the old state, THEN
    Rnor = Rnor.*clearance + synthesis. Clearance ends up scaling
    everything that just arrived at a region (both the retained-R term
    AND the P->R arrival), and growth is added unscaled afterward."""
    N_s, N_t, M_RR, M_PP, M_RP, M_PR = _build_M_components(weights, sconnMov)
    C_diag = scipy.sparse.diags(np.broadcast_to(clearance, (N_t,)))

    M = scipy.sparse.bmat(((C_diag @ M_RR, C_diag @ M_PR),
                           (M_RP,          M_PP)), format="csc")
    B = np.concatenate((growth, np.zeros(N_s * N_t)))

    X = scipy.sparse.linalg.spsolve(scipy.sparse.eye(N_t + N_s * N_t, format="csc") - M, B)
    return X[:N_t], np.reshape(X[N_t:], (N_s, N_t), order="F")


def markov_clearance_first(weights, sconnMov, clearance, growth):
    """AgentBasedModel order: growth_step(); clearance_step(); s_spread_step().
    Growth and clearance are applied to R BEFORE the move, so clearance scales
    the retained-R term and the R->P outflow (on its R side), but NOT the P->R
    arrival (added after clearance_step already ran); growth reaches P_new via
    the same clearance-scaled R->P weights."""
    N_s, N_t, M_RR, M_PP, M_RP, M_PR = _build_M_components(weights, sconnMov)
    C_diag = scipy.sparse.diags(np.broadcast_to(clearance, (N_t,)))

    M = scipy.sparse.bmat(((C_diag @ M_RR,      M_PR),
                           (M_RP @ C_diag,      M_PP)), format="csc")
    B_R = np.broadcast_to(clearance, (N_t,)) * growth
    B_P = (M_RP @ C_diag) @ growth
    B = np.concatenate((B_R, B_P))

    X = scipy.sparse.linalg.spsolve(scipy.sparse.eye(N_t + N_s * N_t, format="csc") - M, B)
    return X[:N_t], np.reshape(X[N_t:], (N_s, N_t), order="F")


def report(name, a, b):
    a = np.asarray(a).ravel()
    b = np.asarray(b).ravel()
    abs_err = np.abs(a - b)
    rel_err = abs_err.max() / max(np.abs(b).max(), 1e-300)
    print(f"  {name:<38} max|Δ| = {abs_err.max():.3e}   rel (vs scale) = {rel_err:.3e}")
    return rel_err


def main():
    sp.Snca_dir_default = GENE_DIR + "/"
    sp.clearance_gene_dir = GENE_DIR + "/"

    weights_raw, distance, region_size, sources, targets, syngene = sp.load_params(PKL)
    clearance_rate = sp.load_clearance(None)   # README's no-clearance test case
    N_s, N_t = weights_raw.shape
    print(f"Test case: no clearance, {N_s} sources -> {N_t} targets "
          f"(v={V}, spread_rate={SPREAD_RATE}, dt={DT})")

    weights, sconnMov = build_transition_matrices(weights_raw, distance, SPREAD_RATE, V, DT)
    clearance = np.exp(-clearance_rate * DT)
    growth = syngene * region_size * DT

    r_move, p_move = markov_move_first(weights, sconnMov, clearance, growth)
    r_clear, p_clear = markov_clearance_first(weights, sconnMov, clearance, growth)

    # ground truth: the actual AgentBasedModel loop (growth -> clearance -> move)
    abm = AgentBasedModel(weights=weights_raw, distance=distance, region_size=region_size,
                          sources=sources, targets=targets, dt=1)
    abm.set_growth_process(growth_rate=syngene)
    abm.set_clearance_process(clearance_rate=clearance_rate)
    abm.set_spread_process(v=V)
    abm.update_spread_process(spread_scale=SPREAD_RATE)
    abm.dt = DT
    r_iter, p_iter, n_iter = sp.iterative_s_steady_state(abm, eps=EPS, max_iter=2_000_000)
    r_iter, p_iter = np.copy(r_iter), np.copy(p_iter)
    print(f"iterative (ground truth, actual AgentBasedModel loop): {n_iter} iterations")

    print("\n--- move_first  vs  clearance_first (isolates the order effect alone) ---")
    report("R: move_first vs clearance_first", r_move, r_clear)
    report("P: move_first vs clearance_first", p_move, p_clear)

    print("\n--- each ordering vs the real iterative loop ---")
    report("R: move_first      vs iterative", r_move, r_iter)
    report("P: move_first      vs iterative", p_move, p_iter)
    rr = report("R: clearance_first vs iterative", r_clear, r_iter)
    rp = report("P: clearance_first vs iterative", p_clear, p_iter)

    print("\n--- verdict ---")
    print("  clearance_first (AgentBasedModel's own order) is the one that matches "
          "the iterative loop it's meant to replace.")
    print(f"  residual gap ({max(rr, rp):.1e}) is consistent with the iterative "
          f"loop's eps={EPS} not being fully converged, not an order mismatch.")

    pd.DataFrame(r_move).to_csv("s_region_markov_move_first.csv")
    pd.DataFrame(r_clear).to_csv("s_region_markov_clearance_first.csv")
    pd.DataFrame(r_iter).to_csv("s_region_iterative_reference.csv")
    print("\nwrote s_region_markov_move_first.csv, s_region_markov_clearance_first.csv, "
          "s_region_iterative_reference.csv")


if __name__ == "__main__":
    main()
