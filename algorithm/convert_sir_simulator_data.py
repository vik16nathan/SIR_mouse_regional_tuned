# -*- coding: utf-8 -*-
"""Convert the MATLAB reference data from

    magnesium2400/SIR_simulator @ speedup_0.2  ->  data/42regions/

into the (params, Snca) pickle format consumed by ``abm_s_population.load_params``.

``load_params`` expects:

* a ``params`` pickle -- a dict with keys
  ``weights``, ``distance``, ``region_size``, ``sources``, ``targets``
* an ``Snca.pkl`` -- a DataFrame indexed by the region labels in ``sources``
  with a ``"Snca"`` column

The MATLAB pipeline (main.m) also uses GBA for the clearance term, so a
matching ``GBA.pkl`` is written as well.

Nothing here applies the SIRsimulator4.m pre-processing (the prob_stay mix, the
row-normalisation, zscore/normcdf).  Exactly like the mouse connectome pickle,
``weights``/``distance`` hold the raw structural connectivity (strength / length)
with a zero diagonal; the pre-processing is the caller's job.

Usage:
    python convert_sir_simulator_data.py \
        --repo /scratch/vnathan/SIR_simulator_mehul \
        --parcellation 42regions --sc sc35 \
        --out-dir ./sir_simulator_ref
"""

import argparse
import os
import pickle

import numpy as np
import pandas as pd
import scipy.io as sio


def convert(repo, parcellation, sc, out_dir):
    data_dir = os.path.join(repo, "data", parcellation)
    ws = sio.loadmat(os.path.join(data_dir, "workspace.mat"))
    scmat = sio.loadmat(os.path.join(data_dir, f"{sc}.mat"))

    sconn_den = np.array(scmat["sconnDen"], dtype=float)
    sconn_len = np.array(scmat["sconnLen"], dtype=float)
    n_regions = sconn_den.shape[0]

    # match SIRsimulator4.m lines 35-36: zero the self connections
    np.fill_diagonal(sconn_den, 0.0)
    np.fill_diagonal(sconn_len, 0.0)

    roi_size = np.array(ws["ROIsize"], dtype=float).ravel()
    snca = np.array(ws["SNCA"], dtype=float).ravel()
    gba = np.array(ws["GBA"], dtype=float).ravel()

    labels = [f"region_{i}" for i in range(n_regions)]

    params = {
        "weights": sconn_den,          # structural connectivity strength (raw)
        "distance": sconn_len,         # structural connectivity length   (raw)
        "region_size": roi_size,       # voxel counts per region
        "sources": labels,
        "targets": labels,
    }

    os.makedirs(out_dir, exist_ok=True)
    params_path = os.path.join(out_dir, f"params_{parcellation}_{sc}.pkl")
    snca_path = os.path.join(out_dir, "Snca.pkl")
    gba_path = os.path.join(out_dir, "GBA.pkl")

    with open(params_path, "wb") as f:
        pickle.dump(params, f)
    pd.DataFrame({"Snca": snca}, index=labels).to_pickle(snca_path)
    pd.DataFrame({"GBA": gba}, index=labels).to_pickle(gba_path)

    print(f"wrote {params_path}")
    print(f"wrote {snca_path}")
    print(f"wrote {gba_path}")
    print(f"  n_regions            = {n_regions}")
    print(f"  sconnDen nnz off-diag = {int(np.count_nonzero(sconn_den))}")
    print(f"  sconnLen nnz off-diag = {int(np.count_nonzero(sconn_len))}")
    return params_path, snca_path, gba_path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--repo", default="/scratch/vnathan/SIR_simulator_mehul",
                   help="clone of magnesium2400/SIR_simulator (speedup_0.2)")
    p.add_argument("--parcellation", default="42regions",
                   choices=["42regions", "65regions", "119regions"])
    p.add_argument("--sc", default="sc35",
                   help="structural-connectivity file stem (42regions: "
                        "sc25/sc30/sc35/sc40; others: sc)")
    p.add_argument("--out-dir", default="./sir_simulator_ref")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    convert(args.repo, args.parcellation, args.sc, args.out_dir)
