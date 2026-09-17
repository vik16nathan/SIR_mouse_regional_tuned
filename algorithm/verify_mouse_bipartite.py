# -*- coding: utf-8 -*-
"""Verify abm_s_population.py's two S-population implementations --
``iterative_s_steady_state`` (the real AgentBasedModel step methods) and
``markov_s_steady_state`` (the closed-form Markov solve) -- against each other
on the REAL, genuinely bipartite (209 sources -> 418 targets) mouse connectome,
for the two README.md ``abm_clearance_genes.py`` examples:

Case 1 -- "Test example for speedup" (no clearance), README.md lines 69-84:
    -p ".../params_regionalized_voxel_weights_after_qc_0413_64pct.pkl"
    -r False  -v 104.87116708452713  -s 0.02122582417807315  -d 0.1  -S 35
  (no -c/-g -> uniform clearance_rate = 0.5, exactly abm_s_population's default)

Case 2 -- full example with clearance, README.md lines 49-67 (same -p/-r/-v/-s/-d
as case 1; only the clearance gene differs):
    -c "Dnajc17"  -g ".../derivatives/yohan_ge_filt/mr10vv0.2/"

Both cases build the AgentBasedModel exactly as abm_s_population.py's __main__
does (set_growth_process / set_clearance_process / set_spread_process /
update_spread_process(spread_scale=spread_rate)), and markov_s_steady_state is
given the same spread_rate, so both methods are solving the same recursion:
SIRsimulator4.m's mobility-pattern construction (diag zeroed, blended with a
"probability of staying" term via spread_rate = 1 - prob_stay, row-normalized
by dt), generalized to N_s != N_t.

The S-population CSVs from both methods are written into ../test/ (i.e.
analysis/test/), alongside (not overwriting) the existing
misfolded-population (I) CSV already there for the no-clearance case.

Usage:
    python verify_mouse_bipartite.py
"""

import os
import pickle
import sys
import time

import numpy as np
import pandas as pd

import abm_s_population as sp
from AgentBasedModel import AgentBasedModel

# --- README.md example parameters (both examples share these) --------------
PKL = "../SIR_inputs/params_regionalized_voxel_weights_after_qc_0413_64pct.pkl"
GENE_DIR = ("/scratch/vnathan/sir_clearance_rebuilt_conn_fixk2_03262026/"
           "derivatives/yohan_ge_filt/mr10vv0.2")
V = 104.87116708452713
SPREAD_RATE = 0.02122582417807315
DT = 0.1
EPS = 1e-5              # abm_s_population.py CLI default (-e)
TEST_DIR = "../test"


def build_abm(weights, distance, region_size, sources, targets, syngene, clearance_rate):
    abm = AgentBasedModel(weights=weights, distance=distance, region_size=region_size,
                          sources=sources, targets=targets, dt=1)
    abm.set_growth_process(growth_rate=syngene)
    abm.set_clearance_process(clearance_rate=clearance_rate)
    abm.set_spread_process(v=V)
    abm.update_spread_process(spread_scale=SPREAD_RATE)
    abm.dt = DT
    return abm


def report(name, a, b):
    a = np.asarray(a).ravel()
    b = np.asarray(b).ravel()
    abs_err = np.abs(a - b)
    rel_err = abs_err.max() / max(np.abs(b).max(), 1e-300)
    print(f"  {name:<40} max|Δ| = {abs_err.max():.3e}   rel (vs scale) = {rel_err:.3e}")
    return abs_err.max(), rel_err


def run_case(label, csv_suffix, clearance_gene, max_iter):
    print(f"\n=== {label} ===")
    sp.Snca_dir_default = GENE_DIR + "/"
    sp.clearance_gene_dir = GENE_DIR + "/"

    weights, distance, region_size, sources, targets, syngene = sp.load_params(PKL)
    clearance_rate = sp.load_clearance(clearance_gene)
    print(f"  {weights.shape[0]} sources -> {weights.shape[1]} targets")

    # (1) closed form
    t0 = time.perf_counter()
    r_closed, p_closed = sp.markov_s_steady_state(
        weights, distance, region_size, clearance_rate, syngene, V, SPREAD_RATE, DT
    )
    t_closed = time.perf_counter() - t0

    # (2) iterative, via the real AgentBasedModel step methods
    abm = build_abm(weights, distance, region_size, sources, targets, syngene, clearance_rate)
    t0 = time.perf_counter()
    r_iter, p_iter, n_iter = sp.iterative_s_steady_state(abm, eps=EPS, max_iter=max_iter)
    t_iter = time.perf_counter() - t0
    r_iter = np.copy(r_iter)
    p_iter = np.copy(p_iter)
    print(f"  iterative: {n_iter} iterations, {t_iter:.2f} s"
          f"{'  (hit max_iter cap -- may not be converged)' if n_iter >= max_iter else ''}")
    print(f"  closed form: {t_closed*1e3:.2f} ms")

    print("  --- equivalence ---")
    _, e_r = report("closed form  vs  iterative", r_closed, r_iter)
    report("closed form P vs  iterative P", p_closed, p_iter)
    ok = e_r < 1e-4
    print(f"  verdict: {'PASS' if ok else 'FAIL'}  (max rel {e_r:.2e}, tol 1e-4)")
    if t_iter > 0:
        print(f"  speedup: {t_iter / t_closed:.1f}x")

    os.makedirs(TEST_DIR, exist_ok=True)
    pd.DataFrame(r_closed).to_csv(os.path.join(TEST_DIR, f"s_region_closed_form.{csv_suffix}.csv"))
    pd.DataFrame(p_closed).to_csv(os.path.join(TEST_DIR, f"s_edge_closed_form.{csv_suffix}.csv"))
    pd.DataFrame(r_iter).to_csv(os.path.join(TEST_DIR, f"s_region_iterative.{csv_suffix}.csv"))
    pd.DataFrame(p_iter).to_csv(os.path.join(TEST_DIR, f"s_edge_iterative.{csv_suffix}.csv"))
    print(f"  wrote {TEST_DIR}/s_region_{{closed_form,iterative}}.{csv_suffix}.csv "
          f"(+ s_edge_*)")
    return ok


def main():
    ok1 = run_case(
        "Case 1: no clearance (README speedup test example)",
        csv_suffix="v.104.87116708452713.spread_rate.0.02122582417807315.dt.0.1.seed.35.clearance_gene.None",
        clearance_gene=None,
        max_iter=2_000_000,
    )
    ok2 = run_case(
        "Case 2: clearance = Dnajc17 (README full example)",
        csv_suffix="v.104.87116708452713.spread_rate.0.02122582417807315.dt.0.1.seed.35.clearance_gene.Dnajc17",
        clearance_gene="Dnajc17",
        max_iter=2_000_000,
    )
    print("\n=== summary ===")
    print(f"  Case 1 (no clearance) : {'PASS' if ok1 else 'FAIL'}")
    print(f"  Case 2 (Dnajc17)      : {'PASS' if ok2 else 'FAIL'}")
    return 0 if (ok1 and ok2) else 1


if __name__ == "__main__":
    raise SystemExit(main())
