
import numpy as np
import pandas as pd
import os
from pathlib import Path

os.chdir(".")
num_reps_baseline = 40
###dictionary to convert from epicentre --> number within .pkl file sources
epi_num_dict = {'CP':35, 'DG':40, 'CA1':24}

clearance_gene_dir="../../derivatives/yohan_ge_filt/mr10vv0.2/"
clearance_list_path=clearance_gene_dir+"gene_names_mr10vv0.2.csv"
sir_commands_dir="../sir_command_files/clearance/"
results_dir="../sir_result_csvs/clearance/"
atrophy_maps_dir="../../preprocessed/steph_janice_regional_jacobians/"
connectome_params_dir="../../derivatives/SIR_inputs/"
pSyn_dir="../../preprocessed/shady_ihc/"

cp_hipp_atrophy_maps_comparison = {'CP':[
    "Time_3cp_rgn_t_stats_full_hemiHuPff_hemiPBS.pkl",
    "Time_3cp_rgn_t_stats_full_hemiMsPff_hemiPBS.pkl" ##strongest atrophy correlation from Tullo et al., 2025 (bioRxiv)
    ], 

    'DG':[
    "Time_3hipp_rgn_t_stats_hemiHuPff_hemiPBS.pkl"
    ]
}

cp_hipp_pSyn_maps_comparison = {'CP':[
    "Total_Pathology_12_MPICP_injection.pkl",
    "Total_Pathology_24_MPICP_injection.pkl"
    ], 
    
    'CA1':[
    "Total_Pathology_24_MPIHIP_injection.pkl"
    ]}

###epsilon for convergence
eps=1e-5

###list of connectome parameters
connectome_params_list_rvm=["params_regionalized_voxel_weights_after_qc_0413_64pct.pkl"]

connectome_params_list_retro_rvm=["params_regionalized_voxel_weights_after_qc_0413_64pct_retro.pkl"]

connectome_params_suffixes_rvm=["rebuilt_rvm"]

connectome_params_dict_rvm={'antero':connectome_params_list_rvm, 'retro':connectome_params_list_retro_rvm}


if __name__ == "__main__":

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
                + f"{atrophy_map_stem}_{direction}_{conn_suff}_clearance_top40.txt"
            )

            ##load top40 clearance genes
            top40_clearance_genes= np.array(pd.read_csv(atrophy_results_dir+"top40_genes.csv")).flatten()


            with open(command_filename, "w") as f:

                for gene in top40_clearance_genes:

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
                        f"-x {atrophy_map_stem}_{gene}_{conn_suff} "
                        f">> {atrophy_results_dir}{gene}.csv" ## note: will rewrite results for a given gene!!! our merge script account for this.
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
                + f"{pSyn_map_stem}_{direction}_{conn_suff}_clearance_top40.txt"
            )


            top40_clearance_genes= np.array(pd.read_csv(pSyn_results_dir+"top40_genes.csv")).flatten()

            with open(command_filename, "w") as f:

                for gene in top40_clearance_genes:

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
                        f"-x {pSyn_map_stem}_{gene}_{conn_suff} "
                        f">> {pSyn_results_dir}{gene}.csv"
                        f"\n"
                    )

                    f.write(cmd)
        

    
    #####################################################################################
    #######FOR THE TOP GENE: REPEAT 40 TIMES (for parameter curvilinear plot)############
    ######################################################################################

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
                + f"{atrophy_map_stem}_{direction}_{conn_suff}_clearance_topGene_40reps.txt"
            )

            ##load top40 clearance genes
            top40_clearance_genes= np.array(pd.read_csv(atrophy_results_dir+"top40_genes.csv")).flatten()
            top_gene = top40_clearance_genes[0]
            print(top_gene)
            gene = top_gene

            with open(command_filename, "w") as f:

                for rep in range(num_reps_baseline):

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
                        f"-x {atrophy_map_stem}_{gene}_{rep}_{conn_suff} "
                        f">> {atrophy_results_dir}{gene}.csv" ## note: will rewrite results for a given gene!!! our merge script account for this.
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
                + f"{pSyn_map_stem}_{direction}_{conn_suff}_clearance_topGene_40reps.txt"
            )


            ##load top40 clearance genes
            top40_clearance_genes= np.array(pd.read_csv(pSyn_results_dir+"top40_genes.csv")).flatten()
            top_gene = top40_clearance_genes[0]
            print(top_gene)
            gene = top_gene

            with open(command_filename, "w") as f:

                for rep in range(num_reps_baseline):

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
                        f"-x {pSyn_map_stem}_{gene}_{rep}_{conn_suff} "
                        f">> {pSyn_results_dir}{gene}.csv"
                        f"\n"
                    )

                    f.write(cmd)