import numpy as np
import pandas as pd
import os
from pathlib import Path

os.chdir(".")

clearance_gene_dir="../../derivatives/yohan_ge_filt/mr10vv0.2/"
clearance_list_path=clearance_gene_dir+"gene_names_mr10vv0.2.csv"
####Atrophy/pSyn maps to consider####
atrophy_maps_dir="../../preprocessed/steph_janice_regional_jacobians/"
connectome_params_dir="../../derivatives/SIR_inputs/"
pSyn_dir="../../preprocessed/shady_ihc/"
sir_commands_dir="../sir_command_files/baselines_max_corr/"
results_dir="../sir_result_csvs/baselines/"

cp_hipp_atrophy_maps_full={'CP':[
    "Time_1cp_rgn_t_stats_full_hemiHuPff_hemiPBS.pkl",
    "Time_1cp_rgn_t_stats_full_hemiMsPff_hemiPBS.pkl",
    "Time_2cp_rgn_t_stats_full_hemiHuPff_hemiPBS.pkl",
    "Time_2cp_rgn_t_stats_full_hemiMsPff_hemiPBS.pkl",
    "Time_3cp_rgn_t_stats_full_hemiHuPff_hemiPBS.pkl",
    "Time_3cp_rgn_t_stats_full_hemiMsPff_hemiPBS.pkl",
    "Time_4cp_rgn_t_stats_full_hemiHuPff_hemiPBS.pkl",
    "Time_4cp_rgn_t_stats_full_hemiMsPff_hemiPBS.pkl"
    ], 

    'DG':[
    "Time_1hipp_rgn_t_stats_hemiHuPff_hemiPBS.pkl",
    "Time_2hipp_rgn_t_stats_hemiHuPff_hemiPBS.pkl",
    "Time_3hipp_rgn_t_stats_hemiHuPff_hemiPBS.pkl",
    "Time_4hipp_rgn_t_stats_hemiHuPff_hemiPBS.pkl"
    ]
}

cp_hipp_atrophy_maps_comparison = {'CP':[
    "Time_3cp_rgn_t_stats_full_hemiHuPff_hemiPBS.pkl",
    "Time_3cp_rgn_t_stats_full_hemiMsPff_hemiPBS.pkl" ##strongest atrophy correlation from Tullo et al., 2025 (bioRxiv)
    ], 

    'DG':[
    "Time_3hipp_rgn_t_stats_hemiHuPff_hemiPBS.pkl"
    ]
}

###epsilon for convergence
eps=1e-5

###dictionary to convert from epi --> "source number" in the list of sources in the .pkl file 
epi_num_dict = {'CP':35, 'DG':40, 'CA1':24}
cp_hipp_pSyn_maps_full={'CP':[
    "Total_Pathology_0.5_MPICP_injection.pkl",
    "Total_Pathology_12_MPICP_injection.pkl",
    "Total_Pathology_18_MPICP_injection.pkl",
    "Total_Pathology_1_MPICP_injection.pkl",
    "Total_Pathology_24_MPICP_injection.pkl",
    "Total_Pathology_3_MPICP_injection.pkl",
    "Total_Pathology_6_MPICP_injection.pkl"
    ], 
    
    'CA1':[
    "Total_Pathology_0.5_MPIHIP_injection.pkl",
    "Total_Pathology_12_MPIHIP_injection.pkl",
    "Total_Pathology_18_MPIHIP_injection.pkl",
    "Total_Pathology_1_MPIHIP_injection.pkl",
    "Total_Pathology_24_MPIHIP_injection.pkl",
    "Total_Pathology_3_MPIHIP_injection.pkl",
    "Total_Pathology_6_MPIHIP_injection.pkl"
    ]}

cp_hipp_pSyn_maps_comparison = {'CP':[
    "Total_Pathology_12_MPICP_injection.pkl",
    "Total_Pathology_24_MPICP_injection.pkl"
    ], 
    
    'CA1':[
    "Total_Pathology_24_MPIHIP_injection.pkl"
    ]}

connectome_params_list=["params_nature_yohan_ccfv3.pkl", 
                        "params_homogeneous_weights_after_qc_0413_64pct.pkl", 
                        "params_regionalized_voxel_weights_after_qc_0413_64pct.pkl"]

connectome_params_list_retro=["params_nature_retro_yohan_ccfv3.pkl", 
                        "params_homogeneous_weights_after_qc_0413_64pct_retro.pkl", 
                        "params_regionalized_voxel_weights_after_qc_0413_64pct_retro.pkl"  ]

connectome_params_suffixes=["original_oh","rebuilt_homogeneous","rebuilt_rvm"]

connectome_params_dict={'antero':connectome_params_list, 'retro':connectome_params_list_retro}

#######come up with a filtered list to only include the rebuilt RVM connectome as the cleanest/most robust connectome#####
connectome_params_list_rvm=["params_regionalized_voxel_weights_after_qc_0413_64pct.pkl"]

connectome_params_list_retro_rvm=["params_regionalized_voxel_weights_after_qc_0413_64pct_retro.pkl"]

connectome_params_suffixes_rvm=["rebuilt_rvm"]

connectome_params_dict_rvm={'antero':connectome_params_list_rvm, 'retro':connectome_params_list_retro_rvm}


if __name__ == "__main__":

    clearance_list_df = pd.read_csv(clearance_list_path)
    clearance_list = np.array(clearance_list_df['gene']).flatten()
    print(clearance_list)

    #######baselines - no clearance#########
    num_reps_baseline = 40

    ########################################
    ###########  ATROPHY  ##################
    ########################################

    for epicentre in ['CP', 'DG']:
        epi = str(epi_num_dict[epicentre])

        for atrophy_map in cp_hipp_atrophy_maps_comparison[epicentre]:

            atrophy_map_path = atrophy_maps_dir + atrophy_map
            atrophy_map_stem = Path(atrophy_map).stem

            for direction in ['antero', 'retro']:

                for i, conn_suff in enumerate(connectome_params_suffixes):

                    params_file = connectome_params_dict[direction][i]
                    params_file_path = connectome_params_dir + params_file

                    command_filename = (
                        sir_commands_dir
                        + f"{atrophy_map_stem}_{direction}_{conn_suff}_baseline.txt"
                    )

                    with open(command_filename, "w") as f:

                        for rep in range(num_reps_baseline):

                            cmd = (
                                f"python3 ../algorithm/abm_optuna_general.py "
                                f"-a True "
                                f"-g {clearance_gene_dir} "
                                f"-m {atrophy_map_path} "
                                f"-p {params_file_path} "
                                f"-r {'True' if direction == 'retro' else 'False'} "
                                f"-t 1000 "
                                f"-d 0.1 "
                                f"-e {eps} "
                                f"-S {epi} "
                                f"-x {atrophy_map_stem}_{rep}_{conn_suff} "
                                f">> {results_dir}{atrophy_map_stem}_{conn_suff}_{direction}.csv"
                                f"\n"
                            )

                            f.write(cmd)


    ########################################
    ############  pSyn  ####################
    ########################################

    for epicentre in ['CP', 'CA1']:
        epi = str(epi_num_dict[epicentre])

        for pSyn_map in cp_hipp_pSyn_maps_comparison[epicentre]:

            pSyn_map_path = pSyn_dir + pSyn_map
            pSyn_map_stem = Path(pSyn_map).stem

            for direction in ['antero', 'retro']:

                for i, conn_suff in enumerate(connectome_params_suffixes):

                    params_file = connectome_params_dict[direction][i]
                    params_file_path = connectome_params_dir + params_file

                    command_filename = (
                        sir_commands_dir
                        + f"{pSyn_map_stem}_{direction}_{conn_suff}_baseline.txt"
                    )

                    with open(command_filename, "w") as f:

                        for rep in range(num_reps_baseline):

                            cmd = (
                                f"python3 ../algorithm/abm_optuna_general.py "
                                f"-a False "
                                f"-g {clearance_gene_dir} "
                                f"-m {pSyn_map_path} "
                                f"-p {params_file_path} "
                                f"-r {'True' if direction == 'retro' else 'False'} "
                                f"-t 1000 "
                                f"-d 0.1 "
                                f"-e {eps} "
                                f"-S {epi} "
                                f"-x {pSyn_map_stem}_{rep}_{conn_suff} "
                                f">> {results_dir}{pSyn_map_stem}_{conn_suff}_{direction}.csv"
                                f"\n"
                            )

                            f.write(cmd)
    

    #############################################
    ########baselines - SINGLE connectome########
    #############################################


    ###########  ATROPHY  ##################
    
    for epicentre in ['CP', 'DG']:
        epi = str(epi_num_dict[epicentre])
        conn_suff="rebuilt_rvm"
        for atrophy_map in cp_hipp_atrophy_maps_full[epicentre]:

            atrophy_map_path = atrophy_maps_dir + atrophy_map
            atrophy_map_stem = Path(atrophy_map).stem

            for direction in ['antero', 'retro']:
                params_file = connectome_params_dict_rvm[direction][0] ##sole element
                params_file_path = connectome_params_dir + params_file

                command_filename = (
                    sir_commands_dir
                    + f"{atrophy_map_stem}_{direction}_{conn_suff}_baseline.txt"
                )

                with open(command_filename, "w") as f:

                    for rep in range(num_reps_baseline):

                        cmd = (
                            f"python3 ../algorithm/abm_optuna_general.py "
                            f"-a True "
                            f"-g {clearance_gene_dir} "
                            f"-m {atrophy_map_path} "
                            f"-p {params_file_path} "
                            f"-r {'True' if direction == 'retro' else 'False'} "
                            f"-t 1000 "
                            f"-d 0.1 "
                            f"-e {eps} "
                            f"-S {epi} "
                            f"-x {atrophy_map_stem}_{rep}_{conn_suff} "
                            f">> {results_dir}{atrophy_map_stem}_{conn_suff}_{direction}.csv"
                            f"\n"
                        )

                        f.write(cmd)

    

    for epicentre in ['CP', 'CA1']:
        epi = str(epi_num_dict[epicentre])
        conn_suff="rebuilt_rvm"

        for pSyn_map in cp_hipp_pSyn_maps_full[epicentre]:

            pSyn_map_path = pSyn_dir + pSyn_map
            pSyn_map_stem = Path(pSyn_map).stem

            for direction in ['antero', 'retro']:

                params_file = connectome_params_dict_rvm[direction][0]
                params_file_path = connectome_params_dir + params_file

                command_filename = (
                    sir_commands_dir
                    + f"{pSyn_map_stem}_{direction}_{conn_suff}_baseline.txt"
                )

                with open(command_filename, "w") as f:

                    for rep in range(num_reps_baseline):

                        cmd = (
                            f"python3 ../algorithm/abm_optuna_general.py "
                            f"-a False "
                            f"-g {clearance_gene_dir} "
                            f"-m {pSyn_map_path} "
                            f"-p {params_file_path} "
                            f"-r {'True' if direction == 'retro' else 'False'} "
                            f"-t 1000 "
                            f"-d 0.1 "
                            f"-e {eps} "
                            f"-S {epi} "
                            f"-x {pSyn_map_stem}_{rep}_{conn_suff} "
                            f">> {results_dir}{pSyn_map_stem}_{conn_suff}_{direction}.csv"
                            f"\n"
                        )

                        f.write(cmd)
        

    ##############CLEARANCE GENES - fix direction to be retro for hipp; antero for CP###############
    ##############allows for comparison####################
    sir_commands_dir="../sir_command_files/clearance/"
    results_dir="../sir_result_csvs/clearance/"

    for epicentre in ['CP', 'DG']:
        epi = str(epi_num_dict[epicentre])
        conn_suff="rebuilt_rvm"
        if epicentre == 'CP':
            direction = "antero"
        else:
            direction = "retro"
        for atrophy_map in cp_hipp_atrophy_maps_comparison[epicentre]:

            atrophy_map_path = atrophy_maps_dir + atrophy_map
            atrophy_map_stem = Path(atrophy_map).stem

            ###make path for all clearance results
            atrophy_results_dir = results_dir+"atrophy/"+atrophy_map_stem+"/"+direction+"/"
            Path(atrophy_results_dir).mkdir(parents=True, exist_ok=True)

            params_file = connectome_params_dict_rvm[direction][0] ##sole element
            params_file_path = connectome_params_dir + params_file

            command_filename = (
                sir_commands_dir
                + f"{atrophy_map_stem}_{direction}_{conn_suff}_clearance.txt"
            )

            ##load clearance genes

            with open(command_filename, "w") as f:

                for gene in clearance_list:

                    cmd = (
                        f"python3 ../algorithm/abm_optuna_general.py "
                        f"-a True "
                        f"-g {clearance_gene_dir} "
                        f"-c {gene} "
                        f"-m {atrophy_map_path} "
                        f"-p {params_file_path} "
                        f"-r {'True' if direction == 'retro' else 'False'} "
                        f"-t 1000 "
                        f"-d 0.1 "
                        f"-e {eps} "
                        f"-S {epi} "
                        f">> {atrophy_results_dir}{gene}.csv"
                        f"\n"
                    )

                    f.write(cmd)

    for epicentre in ['CP', 'CA1']:
        epi = str(epi_num_dict[epicentre])
        conn_suff="rebuilt_rvm"
        
        if epicentre == 'CP':
            direction = "antero"
        else:
            direction = "retro"

        for pSyn_map in cp_hipp_pSyn_maps_comparison[epicentre]:


            pSyn_map_path = pSyn_dir + pSyn_map
            pSyn_map_stem = Path(pSyn_map).stem

            ###make path for all clearance results
            pSyn_results_dir = results_dir+"pSyn/"+pSyn_map_stem+"/"+direction+"/"
            Path(pSyn_results_dir).mkdir(parents=True, exist_ok=True)

            params_file = connectome_params_dict_rvm[direction][0]
            params_file_path = connectome_params_dir + params_file

            command_filename = (
                sir_commands_dir
                + f"{pSyn_map_stem}_{direction}_{conn_suff}_clearance.txt"
            )

            with open(command_filename, "w") as f:

                for gene in clearance_list:

                    cmd = (
                        f"python3 ../algorithm/abm_optuna_general.py "
                        f"-a False "
                        f"-g {clearance_gene_dir} "
                        f"-c {gene} "
                        f"-m {pSyn_map_path} "
                        f"-p {params_file_path} "
                        f"-r {'True' if direction == 'retro' else 'False'} "
                        f"-t 1000 "
                        f"-d 0.1 "
                        f"-e {eps} "
                        f"-S {epi} "
                        f">> {pSyn_results_dir}{gene}.csv"
                        f"\n"
                    )

                    f.write(cmd)

        