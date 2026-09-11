# -*- coding: utf-8 -*-
"""This script runs the agent-based Susceptible-Infected-Removed (SIR) Model.
Authors:
    Vikram Nathan, Ying-Qiu Zheng, Shady Rahayel
    
For running the model, run:
    
    python abm.py --retro True --speed 10 --spreading-rate 0.01 --time 1000
        --delta-t 0.1 --seed -1 --seed-amount 1
        
--retro True specifies a retrograde spreading
--speed is the spreading speed of agents in edges
--spreading-rate is the probability of staying inside a region
--time is the spreading time of agents
--delta-t is the size of timesteps
--seed is an integer that refers to the list of regions listed
alphabetically from the Allen Mouse Brain Atlas (see params_nature_retro.pickle)
CP = 35, ACB = 3, and CA1 = 24
--seed-amount is the initial injected amount of infected agents
        
This generates arrays containing the number of normal and infected agents
at each iteration for every region of the Allen Mouse Brain Atlas.
The distribution of normal agents can be found in .s_region_history
The distribution of infected agents can be found in .i_region_history

"""

import sys
import numpy as np
import pandas as pd
import pickle
sys.path.insert(1, './SIR_mouse/model/')
from AgentBasedModel import AgentBasedModel
from scipy.stats import zscore, norm
import scipy.sparse
from tqdm import tqdm
import argparse
import time


SIR_input_dir="../../derivatives/SIR_inputs/"
Snca_dir_default="../../derivatives/yohan_ge_filt/mr10vv0.2/"
DEFAULT_PARAMS_ANTERO=SIR_input_dir+"params_regionalized_voxel_weights_after_qc_0413_64pct.pkl"
DEFAULT_PARAMS_RETRO=SIR_input_dir+"params_regionalized_voxel_weights_after_qc_0413_64pct_retro.pkl"

def parse_arguments():
    parser = argparse.ArgumentParser() 
    parser.add_argument(
        "-g", "--GEdir", default=None, dest="clearance_gene_dir",
         nargs='?', help="Clearance Gene Directory"
    )

    parser.add_argument(
        "-o", "--output-dir", dest="output_dir",
        nargs='?', help="Result Output Directory"
    )
    
    parser.add_argument(
        "-r", "--retro", default=True, dest="retro",
        type=str2bool, nargs='?', help="Retrograde spreading (True) "
    )
    parser.add_argument(
        "-v", "--speed", default=10, dest="v", nargs='?',
        help="Spreading speed", type=float
    )
    parser.add_argument(
        "-s", "--spreading-rate", default=0.01, type=float,
        nargs='?', help="Spreading rate", dest="spread_rate"
    )
    parser.add_argument(
        "-t", "--time", default=1000, type=int, nargs='?',
        dest="total_time", help="Total spreading time"
    )
    parser.add_argument(
        "-d", "--delta-t", default=0.1, type=float, nargs='?',
        dest="dt", help="Size of time increment"
    )
    parser.add_argument(
        "-S", "--seed", default=35, type=int, nargs='?',
        dest="seed", help="Simulated seeding site of misfolded alpha-syn"
        # injecting into the CP; CP = 35; CA1 = 24
    )
    parser.add_argument(
        "-i", "--seed-amount", default="1", type=float,
        dest="injection_amount", help="Seeding amount of misfolded alpha-syn"
    )
    parser.add_argument(
        "-c", "--clearance", default=None, dest="clearance_gene", type=none_or_str,
        help="Specify the gene modulating clearance (omit, or pass 'None', for uniform clearance)"
    )
    parser.add_argument(
        "-k1", "--k1", default=None, type=float, dest="k1_atrophy", nargs='?',
        help="k1_atrophy"
    )
    parser.add_argument(
        "-k2", "--k2", default=None, type=float, dest="k2_atrophy", nargs='?',
        help="k2_atrophy"
    )
    parser.add_argument(
        "-p", "--params", default=None,
        dest="params", nargs='?', help="Specify the .pkl file with the connectome "
        "(default: most recently updated regionalized-voxel-weights connectome, "
        "chosen based on --retro)"
    )
    parser.add_argument(
        "-e", "--epsilon", default=1e-5, type=float, nargs='?', #Shady's default was 1e-7 but this was waaaay too slow
        dest="eps", help="Convergence threshold"
    )
    parser.add_argument(
        "-x", "--suffix", default="",
        dest="suffix", nargs='?', help="Descriptive suffix for simulation output"
    )
    args = parser.parse_args()
    if args.params is None:
        args.params = DEFAULT_PARAMS_RETRO if args.retro else DEFAULT_PARAMS_ANTERO
    return args


def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def none_or_str(v):
    if v is None or v.lower() == 'none':
        return None
    return v


def load_params(params_file):
 
    # load homogeneous values
    # choose the rate from 0.1 to 0.9
    #homorate = 0.5
    #syngene = np.transpose(np.full((1,209),homorate))
    with open(params_file, "rb") as f:
        params = pickle.load(f)
        # if retro is False, load data from 'params_nature.pickle'; anterograde spreading
    #print(params)
    #FILTER TO GET RID OF INVALID REGIONS!
    weights = params['weights']
    distance = params['distance']
    region_size = params['region_size']
    sources = params['sources']
    targets = params['targets']

   
    #Load Yohan's SNCA
    ge = pd.read_pickle(Snca_dir_default+'Snca.pkl')
    syngene = ge.loc[sources, "Snca"] #can modify to be another gene if needed... 
    syngene = norm.cdf(zscore(syngene))
    
    return (
            weights, distance, np.append(region_size, region_size),
            sources, targets, np.append(syngene, syngene)
        )

def load_clearance(clearance_gene=None): 
    if clearance_gene is None:
        return 0.5 # default clearance rate
        # if clearance_gene is not specified (i.e., it is None), the function returns the default clearance rate of 0.5
    else:
        #CHANGE THIS: preprocessed .csv file, need to COPY to both hemispheres and get rid of regions
        ge = pd.read_csv(clearance_gene_dir+clearance_gene+'.csv',index_col=0)
        epr=ge.loc[:,clearance_gene]
        #epr = ge.loc[sources, clearance_gene] 
        # #everything is in the right order - see SIR_mouse_exploration.ipynb on cicws        
        epr = np.append(epr, epr)
        # appends the data to itself

        return norm.cdf(zscore(epr.flatten()))
        # normalizes and transforms the data, and returns the resulting array as a cumulative distribution function


if __name__ == "__main__":
    # run ABM
    # read arguments
    args = parse_arguments()

    retro = args.retro
    #injection_site = args.injection_site
    v = args.v
    spread_rate = args.spread_rate
    dt = args.dt
    seed = args.seed
    injection_amount = args.injection_amount
    total_time = args.total_time
    clearance_gene = args.clearance_gene
    k1_atrophy = args.k1_atrophy
    k2_atrophy = args.k2_atrophy
    params_file = args.params
    suffix = args.suffix
    eps = args.eps

    output_dir=args.output_dir
    clearance_gene_dir=args.clearance_gene_dir

    weights, distance, region_size, sources, targets, syngene = load_params(params_file)
    #print(weights.dtype)
    #weights = weights/np.max(weights) ##prevent overflow errors
    clearance_rate = load_clearance(clearance_gene)
    #with open('snca_norm.pickle', 'wb') as f:
    #    pickle.dump(syngene, f)
    #clearance_rate = load_clearance(clearance_gene) 
    #CHANGE THIS: add a filter step to get rid of regions that aren't represented in the clearance gene expression file
    abm = AgentBasedModel(
        weights=weights, distance=distance, region_size=region_size,
        sources=sources, targets=targets, dt=1
    )
    # reads input arguments passed to the script through the command line
        # such as the spread rate, the injection amount, and the total time. 

    abm.set_growth_process(growth_rate=syngene)
    abm.set_clearance_process(clearance_rate=clearance_rate)
    abm.set_spread_process(v=v)
    abm.update_spread_process(spread_scale=spread_rate) #CHANGE THIS: look into v scale
    #CHANGE THIS: look into (lack??) of update_growth_process, update_clearance_process, update_trans_process

    # calls several functions to load the model parameters, 
        # the clearance rate of proteins, and to set up the growth, clearance, and spreading processes of the ABM

    ###########SHOULD BE ABLE TO GET RID OF ENTIRE S_SPREAD_STEP BASED ON CLOSED-FORM SOLUTION###############
    # growth process
    #print("Begin protein growth process....")
    #start_time = time.time()
    #for t in range(30000):
    #    prev = np.copy(abm.s_region)
    #    abm.growth_step()
    #    abm.clearance_step()
    #    abm.s_spread_step()
    #    if np.where(np.abs(prev - abm.s_region) / abm.s_region > eps, 1, 0).sum() == 0:
    #        #print("entered")
    #        break
    #abm.dt = 0.1
    #for t in range(500000000):
    #    prev = np.copy(abm.s_region)
    #    abm.growth_step()
    #    abm.clearance_step()
    #    abm.s_spread_step()
    #    if np.where(np.abs(prev - abm.s_region) / abm.s_region > eps, 1, 0).sum() == 0:
    #        break
    #abm.dt = dt
    #for t in range(1000000000):
    #    prev = np.copy(abm.s_region)
    #    abm.growth_step()
    #    abm.clearance_step()
    #    abm.s_spread_step()
    #    if np.where(np.abs(prev - abm.s_region) / abm.s_region > eps, 1, 0).sum() == 0:
    #        break
    #stop_time = time.time()
    #elapsed_time = stop_time - start_time
    #print(f"Protein growth time: {elapsed_time:.4f} seconds")
    # runs the protein growth process in three stages
    # In each stage, the ABM model computes the protein growth, clearance, and spreading steps until the growth process stops, 
    # which is detected by comparing the changes in protein concentrations between time steps

    ###########################################################################################################
    #####################TODO: ###########################

    sconnMov = v / distance * dt
    sconnMov[distance == 0] = 0                  # SIRsimulator4.m l.39: sconnMov(sconnLen == 0) = 0
    N = len(sources)
    print("Calculate S protein steady state...")
    ###1. Build the Markov Transition Components
    ###To make this efficient and prevent memory explosion, we use sparse matrices.

    ###[R -> R]: Proteins that stay in regions (1 minus sum of outward weights)
    M_RR = np.diag(1 - weights.sum(axis=1))

    ###[P -> P]: Proteins that stay in paths (1 minus outward path connections)
    M_PP = np.diag(np.where(sconnMov.flatten(order="F") > 0, 1, 0) - sconnMov.flatten(order="F"))

    ###For moving between R and P, we map the 2D matrix indices to the 1D vector
    list_1_n = list(range(0, N))
    rows_i, cols_j = np.meshgrid(list_1_n, list_1_n, indexing="ij")

    ###[R -> P]: R(i) sends proteins to P(i,j) based on weights(i,j)
    ###rows_i(:) gives the source region 'i' for every linear index in P
    list_1_n2 = list(range(0, N**2))
    M_RP = scipy.sparse.coo_matrix((weights.flatten(order="F"), (list_1_n2, rows_i.flatten(order="F"))), shape=(N**2, N)).toarray()

    ###[P -> R]: P(i,j) sends proteins to R(j) based on sconnMov(i,j)
    ###cols_j(:) gives the destination region 'j' for every linear index in P
    M_PR = scipy.sparse.coo_matrix((sconnMov.flatten(order="F"), (cols_j.flatten(order="F"), list_1_n2)), shape=(N, N**2)).toarray()

    ###Combine into a single Movement matrix with clearance in the regions afterwards:
    clearance = np.exp(-clearance_rate * dt)

    C_diag = np.diag(clearance)
    M = np.vstack((np.hstack((C_diag @ M_RR, C_diag @ M_PR)),
                  np.hstack((M_RP, M_PP))))

    ###2. Synthesis
    synthesis_rate = syngene
    synthesis = (synthesis_rate * region_size) * dt
    B = np.concatenate((synthesis, np.zeros(N**2)))

    ###3. Run the Markov Loop
    print('normal alpha synuclein growth (Markov Analytic Solution)')
    X, resid, rank, s = np.linalg.lstsq((np.eye(N + N**2) - M), B, rcond=None)

    ###4. Unpack results back to original shapes
    Rnor0 = X[0:N]
    Pnor0 = np.reshape(X[N:N + N**2], (N, N), order="F")

    ####Assign S values to ABM

    abm.s_region = Rnor0
    abm.s_edge = Pnor0

    pd.DataFrame(abm.s_region).to_csv("s_region_closed_form.csv")
    pd.DataFrame(abm.s_edge).to_csv("s_edge_closed_form.csv")