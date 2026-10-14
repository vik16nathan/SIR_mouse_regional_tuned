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
import scipy.sparse.linalg
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
    parser.add_argument(
        "--s-method", default="iterative", choices=["iterative", "closed_form"],
        dest="s_method", help="How to compute the S (normal alpha-syn) steady "
        "state before injection: the original iterative convergence loop, or "
        "the closed-form Markov solution (see abm_s_population.py)"
    )
    parser.add_argument(
        "--s-dt", default=1e-6, type=float, dest="s_dt",
        help="Timestep used by the closed-form S solve (--s-method closed_form only); "
        "independent of -d/--delta-t, which controls the I spreading simulation"
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


def compute_s_iterative(abm, dt, eps=1e-5):
    """Original 3-stage iterative convergence of S to steady state (unchanged
    from the inline logic this used to be): dt=1 (abm's construction default)
    for a coarse first pass, then dt=0.1, then the caller's own dt -- the same
    dt that is later used for the I spreading simulation.
    """
    for t in range(30000):
        prev = np.copy(abm.s_region)
        abm.growth_step()
        abm.clearance_step()
        abm.s_spread_step()
        if np.where(np.abs(prev - abm.s_region) / abm.s_region > eps, 1, 0).sum() == 0:
            break
    abm.dt = 0.1
    for t in range(500000000):
        prev = np.copy(abm.s_region)
        abm.growth_step()
        abm.clearance_step()
        abm.s_spread_step()
        if np.where(np.abs(prev - abm.s_region) / abm.s_region > eps, 1, 0).sum() == 0:
            break
    abm.dt = dt
    for t in range(1000000000):
        prev = np.copy(abm.s_region)
        abm.growth_step()
        abm.clearance_step()
        abm.s_spread_step()
        if np.where(np.abs(prev - abm.s_region) / abm.s_region > eps, 1, 0).sum() == 0:
            break
    return abm


def markov_s_steady_state(weights, distance, region_size, clearance_rate, syngene, v, spread_rate, dt):
    """Closed-form steady state of the S (normal alpha-syn) growth process --
    ported verbatim from abm_s_population.py in this same project, see that
    file's docstring for the full derivation. Returns the same (s_region,
    s_edge) pair that compute_s_iterative leaves on `abm` after its
    convergence loop.
    """
    sconnMov = v / distance * dt
    sconnMov[distance == 0] = 0
    N_s = weights.shape[0]
    N_t = weights.shape[1]

    diag_idx = np.arange(N_s)
    weights = weights / weights.sum(axis=1, keepdims=True)
    weights[diag_idx, diag_idx] = 0
    weights = weights * spread_rate * dt

    print("Calculate S protein steady state (closed form)...")
    M_RR_diag = np.ones(N_t)
    M_RR_diag[:N_s] = 1 - weights.sum(axis=1)
    M_RR = scipy.sparse.diags(M_RR_diag)

    M_PP = scipy.sparse.diags(np.where(sconnMov.flatten(order="F") > 0, 1, 0) - sconnMov.flatten(order="F"))

    list_1_n_s = list(range(0, N_s))
    list_1_n_t = list(range(0, N_t))
    rows_i, cols_j = np.meshgrid(list_1_n_s, list_1_n_t, indexing="ij")

    list_1_n2 = list(range(0, N_s * N_t))
    M_RP = scipy.sparse.coo_matrix((weights.flatten(order="F"), (list_1_n2, rows_i.flatten(order="F"))), shape=(N_s * N_t, N_t))
    M_PR = scipy.sparse.coo_matrix((sconnMov.flatten(order="F"), (cols_j.flatten(order="F"), list_1_n2)), shape=(N_t, N_s * N_t))

    clearance = np.broadcast_to(np.exp(-clearance_rate * dt), (N_t,))
    C_diag = scipy.sparse.diags(clearance)
    M = scipy.sparse.bmat(((C_diag @ M_RR,      M_PR),
                           (M_RP @ C_diag,      M_PP)), format="csc")

    synthesis_rate = syngene
    growth = (synthesis_rate * region_size) * dt
    B_R = (C_diag @ M_RR) @ growth
    B_P = (M_RP @ C_diag) @ growth
    B = np.concatenate((B_R, B_P))

    X = scipy.sparse.linalg.spsolve(scipy.sparse.eye(N_t + N_s * N_t, format="csc") - M, B)

    Rnor0 = X[0:N_t]
    Pnor0 = np.reshape(X[N_t:N_t + N_s * N_t], (N_s, N_t), order="F")
    return Rnor0, Pnor0


def compute_s_closed_form(abm, weights, distance, region_size, clearance_rate, syngene, v, spread_rate, s_dt):
    """Set abm.s_region/abm.s_edge to the closed-form Markov steady state,
    bypassing compute_s_iterative's convergence loop entirely. s_dt is
    independent of the dt used for I spreading.
    """
    s_region, s_edge = markov_s_steady_state(
        weights, distance, region_size, clearance_rate, syngene, v, spread_rate, s_dt
    )
    abm.s_region[:] = s_region
    abm.s_edge[:] = s_edge
    return abm


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
    eps = args.eps
    s_method_tag = "s_iterative" if args.s_method == "iterative" else f"s_closed_form.s_dt.{args.s_dt:g}"
    suffix = (args.suffix + "." if args.suffix else "") + s_method_tag

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

    # growth process: either the original iterative convergence loop, or the
    # closed-form Markov steady state (--s-method closed_form, timestep --s-dt,
    # independent of -d/--delta-t which only controls I spreading below)
    print(f"Begin protein growth process (s_method={args.s_method})....")
    start_time = time.time()
    if args.s_method == "closed_form":
        abm = compute_s_closed_form(
            abm, weights, distance, region_size, clearance_rate, syngene,
            v, spread_rate, args.s_dt
        )
        abm.dt = dt
    else:
        abm = compute_s_iterative(abm, dt, eps)
    stop_time = time.time()
    elapsed_time = stop_time - start_time
    print(f"Protein growth time: {elapsed_time:.4f} seconds")

    # spread process
    print("Begin protein spreading process...")
    start_time = time.time()
    print("Inject infectious proteins into {}...".format(abm.targets[seed]))
    abm.injection(seed=seed, amount=injection_amount)

    abm.preallocate_history(total_time)
    for t in tqdm(range(total_time)):
        prev = np.copy(abm.i_region)
        abm.growth_step()
        abm.clearance_step()
        abm.trans_step()
        abm.s_spread_step()
        abm.i_spread_step()
        abm.record_history_region()
        abm.record_history_edge()
    abm.finalize_history()

    # runs the protein spreading process for the specified total time.
    # injects infectious proteins into the ABM's target region,
    # simulates the protein growth, clearance, and spreading steps,
    # and records the resulting protein concentrations at each time step.

    stop_time = time.time()
    elapsed_time = stop_time - start_time
    print(f"Protein spreading time: {elapsed_time:.4f} seconds")
    # saves the results of the ABM simulation in a pickle file

    # generates arrays containing the number of normal and infected agents at each iteration - LOOK AT ZHENG PAPER
    # for every region of the Allen Mouse Brain Atlas
        # The distribution of normal agents can be found in .s_region_history
        # The distribution of infected agents can be found in .i_region_history

        # i_region_history = the number of infected agents in every region for every time step),
        # s_region_history = the number of susceptible agents in every region for every time step)
        # i_edge_history = the number of infected agents in every edge for every time step),
        # s_edge_history = the number of susceptibles agents in every edge for every time step)


    start_time = time.time()
    infected_region_time=abm.i_region_history
    #infected_edge=results.i_edge_history
    normal_region_time=abm.s_region_history
    #normal_edge=results.s_edge_history
    dt=abm.dt
    # size of timesteps

    # Total number of proteins per regions per timestep
    total_proteins_region_time = normal_region_time + infected_region_time

    # Ratio of infected over total proteins per region at every timestep
    infected_proteins_ratio_region_time = infected_region_time / total_proteins_region_time
    #ratio(Rmis_all<1) = 0; % remove possible NaNs... from MATLAB code

    if k1_atrophy is None and k2_atrophy is None:
        simulated=infected_region_time.T

        abm_filename = "abm_spread_v.{0}.spread_rate.{1}.dt.{2}.seed.{3}.injection_amount.{4}.clearance_gene.{5}." \
                .format(v, spread_rate, dt, seed, injection_amount, clearance_gene) + "." + suffix
        np.savetxt(output_dir+abm_filename + '.csv', simulated, delimiter=',', fmt='%s')
        sys.exit(0)


    abm_filename = "abm_spread_v.{0}.spread_rate.{1}.dt.{2}.seed.{3}.injection_amount.{4}.clearance_gene.{5}.k1.{6}.k2.{7}" \
        .format(v, spread_rate, dt, seed, injection_amount, clearance_gene, k1_atrophy, k2_atrophy) + "." + suffix

    if retro is True:
        abm_filename = "retro_" + abm_filename

    abm_filename_pickle = abm_filename + ".pickle"

    with open(output_dir+abm_filename_pickle, "wb") as f:
        pickle.dump(abm, f)

    ## Normalized measure of connectivity strength for every target regions
    # step 1: Calculate the sum of connectivity strength for each source region
    sum_strength = np.sum(weights, axis=1)
        # weights is defined as a matrix with dimensions N_regions by N_regions; 209 x 418
        # np.sum(weights, axis=1); sum across every row to get total connectivity strength across all injection sites (target regions) for each region

    # step 2: Calculate the ratio of weights to the sum of strength for each source region
    ratio_weights = weights / np.repeat(sum_strength[:, np.newaxis], len(targets), axis=1)
        # np.repeat(sum_strength[:, np.newaxis], len(targets), axis=1)
        # takes the sum of the connectivity strength for each source regions across every target region repeated N_target_regions times

    # # Print the resulting ratio weights
    # print(ratio_weights.shape) (209,418)
    #         # from this ratio, the resulting weights matrix provides a normalized measure of connectivity strength
    #         # that takes into account the relative contribution of each source region to the overall connectivity with the target regions.

    k1 = k1_atrophy # 0.5
    #k2 = 1 - k1
    k2 = k2_atrophy # 0.5

    #############################BUG FIX########################
    # input weights of deafferentation (scaled by structural connectivity)

    # Create the new matrix with dimensions (full_target, full_target)
    new_ratio_weights = np.zeros((len(targets), len(targets)))
    # Copy ipsi source --> full target matrix for the top half of the full target --> full target conn. matrix
    new_ratio_weights[:len(sources), :] = ratio_weights

    # Copy ipsi source --> contra target (RHS) to contra source --> ipsi target (LHS) (assume symmetry)
    new_ratio_weights[len(sources):, :len(sources)] = ratio_weights[:, len(sources):]

    # Copy ipsi source --> ipsi target (RHS) to contra source --> contra target (LHS) (assume symmetry)
    new_ratio_weights[len(sources):, len(sources):] = ratio_weights[:, :len(sources)]


    # neuronal loss caused by lack of input from neighboring regions
    ratio_cum = np.dot(new_ratio_weights, (1 - np.exp(-infected_proteins_ratio_region_time.T * dt)))
            # apply a decay or scaling factor to the weights, modelling time-dependent changes or
            # attenuation in connectivity strengths as a result of atrophy/deaff for example
    # print(ratio_cum.shape) #(418, 1000)

    # simulate a backward shift of one time step
    ratio_cum = np.hstack([np.zeros((len(targets), 1)), ratio_cum[:, :-1]])
            # shifts the ratio_cum matrix one time step back by
                    # appending a column of zeros at the beginning and removing the last column of ratio_cum.
            # The result is that each column in ratio_cum is shifted one position to the right,
                    # with a column of zeros added at the leftmost position.
    # print(ratio_cum.shape) #(418, 1000)

    ratio_cum = (k2 * ratio_cum) + (k1 * (1 - np.exp(-infected_proteins_ratio_region_time.T * dt)))
    # update the ratio_cum matrix by applying a weighted combination of the shifted ratio_cum and a term involving k1, ratio, and dt.
        # k2 * ratio_cum: This term scales the shifted ratio_cum matrix by a factor k2.
        # k1 * (1 - exp(-ratio * dt)): calculates the exponential decay factor 1 - exp(-ratio * dt) and scales it by k1.
        # The resulting ratio_cum matrix is obtained by adding these two terms together.
    # print(ratio_cum.shape) #(418, 1000)

    # add all the increments across each dt
    simulated_atrophy = np.cumsum(ratio_cum, axis=1)
        # np.cumsum(ratio_cum, axis=1) computes the cumulative sum along the second axis (axis=1) of the ratio_cum array.
        # The resulting array simulated_atrophy will have the same shape as ratio_cum, and
        # each element will contain the cumulative sum of the corresponding elements in ratio_cum.
    # print(simulated_atrophy)
    # print(simulated_atrophy.shape) #(418, 1000)

    # Save variables to a file
    with open(output_dir+'saved_atrophy_' + abm_filename_pickle, 'wb') as f:
        pickle.dump((simulated_atrophy, infected_proteins_ratio_region_time, targets, dt), f)

    #save to a csv file

    #import rpy2.robjects as robjects
    np.savetxt(output_dir+abm_filename + '.csv', simulated_atrophy, delimiter=',', fmt='%s')
    stop_time = time.time()
    elapsed_time = stop_time - start_time
    print(f"Regional atrophy map post time: {elapsed_time:.4f} seconds")
