import numpy as np
import pandas as pd
import pickle

conn_dir = "../derivatives/regionalized_connectomes/"
params_dir = "./SIR_mouse/model/"
output_dir="../derivatives/SIR_inputs/"

##maximum value of connection strength in the connectome from Oh et al., 2014;
##prevents overflow error when training models in Optuna, especially with larger v parameter

#MAX_VAL_SCONN=
def load_filter_params(params_dir):
    old_weights_antero_path=params_dir+'params_nature.pickle'
    old_weights_retro_path=params_dir+'params_nature_retro.pickle'
    old_snca_path =params_dir+'snca_nature.pickle'
    num_reg = 213 #for MRI (had to discard SUBd and SUBv)
    with open(old_weights_antero_path, "rb") as f:
                params = pickle.load(f)
                # if retro is False, load data from 'params_nature.pickle'; anterograde spreading
                
    weights = params['sconn']

    with open(old_weights_retro_path, "rb") as f:
                params = pickle.load(f)
                # if retro is False, load data from 'params_nature.pickle'; anterograde spreading
                
    weights_retro = params['sconn']

    with open(old_snca_path, "rb") as f:
        snca = pickle.load(f)

    sources = params['source_list']
    targets = params['target_list'] 

    ###confirm that retrograde spreading is just the transpose
    weights_ipsi = weights[:, :num_reg]
    weights_contra = weights[:, num_reg:]

    weights_retro_ipsi = weights_retro[:, :num_reg]
    weights_retro_contra = weights_retro[:, num_reg:]

    ###assert that transpose operation for "retrograde" spreading is justified
    np.all(weights_ipsi.T == weights_retro_ipsi), np.all(weights_contra.T == weights_retro_contra)

    ###dorsal and ventral subiculum are no longer in CCFv3, 2020 version 
    ###MEPO and RM are ADDITIONAL REGIONS that are at the midline/too small to resolve 200 um gene expression
    source_indices_filt = ~np.isin(sources, ["SUBd", "SUBv", "MEPO", "RM"])
    target_indices_filt = ~np.isin(targets, ["SUBd", "SUBv", "MEPO", "RM"])
    
    ###get rid of regions without "child" labels; weren't processed in Yohan's GE filtering pipeline b/c they're too small

    sdist_filt = params['sdist'][source_indices_filt,:][:, target_indices_filt]
    roisize_filt = params['roisize'][source_indices_filt]
    source_list_filt = np.array(sources)[source_indices_filt]
    target_list_filt = np.array(targets)[target_indices_filt]

    weights_filt = weights[source_indices_filt,:][:, target_indices_filt]
    weights_bin = np.where(weights_filt > 0, 1, 0)
    sparsity_weights_filt = 1- np.sum(weights_bin)/(weights_bin.shape[0]*weights_bin.shape[1]) ##ALL connections, none filtered (64% sparsity)

    return sdist_filt, roisize_filt, source_list_filt, target_list_filt, sparsity_weights_filt

def process_regionalized_connectome(conversion_table, conn_list, prefix, prefix_thresh, output_dir):
    prefix = output_prefixes[i]
    prefix_thresh = output_prefixes_64pct[i]

    ##automatically skips SUBd and SUBv (not in most recent CCFv3 ontology)
    ##load ipsi/contra
    conn_ipsi = pd.read_csv(conn_list[0], index_col=0) 
    conn_ipsi.columns = np.array(conn_ipsi.columns).astype(int)
    conn_contra = pd.read_csv(conn_list[1], index_col=0)
    conn_contra.columns = np.array(conn_contra.columns).astype(int)
    
    
    source_numbers = np.array([conversion_table.loc[np.where(conversion_table.loc[:,"abbreviation"] == s)[0][0],"structure ID"] for s in source_list_filt])
    target_numbers = np.array([conversion_table.loc[np.where(conversion_table.loc[:,"abbreviation"] == t)[0][0],"structure ID"] for t in target_list_filt])

    half = int(len(target_numbers)/2)
    print(target_numbers[:int(len(target_numbers)/2)] == target_numbers[int(len(target_numbers)/2):])

    conn_ipsi_ordered=conn_ipsi.loc[source_numbers, target_numbers[:int(len(target_numbers)/2)]]
    conn_contra_ordered=conn_contra.loc[source_numbers, target_numbers[int(len(target_numbers)/2):]]

    assert np.all(conn_ipsi_ordered.index.values == source_numbers)
    assert np.all(conn_ipsi_ordered.columns.values == target_numbers[:half])

    assert np.all(conn_contra_ordered.index.values == source_numbers)
    assert np.all(conn_contra_ordered.columns.values == target_numbers[half:])

    ##create anterograde/retrograde spreading matrices (n x 2n, ipsi + contra)
    conn_antero = np.hstack((np.array(conn_ipsi_ordered), np.array(conn_contra_ordered)))
    conn_retro = np.hstack((np.array(conn_ipsi_ordered).T, np.array(conn_contra_ordered).T))

    ##percentile filtering (ensure same level of sparsity as original connectome)
    pct_thresh = np.percentile(conn_antero, sparsity_weights_filt*100)

    conn_antero_thresh = np.where(conn_antero > pct_thresh, conn_antero, 0)
    conn_retro_thresh = np.where(conn_retro > pct_thresh, conn_retro, 0)

    weights_bin = np.where(conn_retro_thresh > 0, 1, 0)
    weights_bin
    print(1 - np.sum(weights_bin)/(weights_bin.shape[0]*weights_bin.shape[1])) ##ALL connections, none filtered (64% sparsity)

    ###write out .pkl files, keeping all other params the same as the original .pkl files for SIR modelling
    with open(output_dir+prefix+".pkl", 'wb') as f: #dropped regions with low resolution in CCFv3 template
        pickle.dump({'weights':conn_antero, 'distance':sdist_filt, 
                'region_size':roisize_filt, 'sources':source_list_filt, 'targets':target_list_filt}, f)
    
    ##repeat for retrograde
    with open(output_dir+prefix+"_retro.pkl", 'wb') as f: #dropped regions with low resolution in CCFv3 template
        pickle.dump({'weights':conn_retro, 'distance':sdist_filt, 
                'region_size':roisize_filt, 'sources':source_list_filt, 'targets':target_list_filt}, f)
    
    ###########REPEAT FOR THRESHOLDED CONNECTOMES######################
    with open(output_dir+prefix_thresh+".pkl", 'wb') as f: #dropped regions with low resolution in CCFv3 template
        pickle.dump({'weights':conn_antero_thresh, 'distance':sdist_filt, 
                'region_size':roisize_filt, 'sources':source_list_filt, 'targets':target_list_filt}, f)
    
    ##repeat for retrograde
    with open(output_dir+prefix_thresh+"_retro.pkl", 'wb') as f: #dropped regions with low resolution in CCFv3 template
        pickle.dump({'weights':conn_retro_thresh, 'distance':sdist_filt, 
                'region_size':roisize_filt, 'sources':source_list_filt, 'targets':target_list_filt}, f)
    return

if __name__ == "__main__":
    sdist_filt, roisize_filt, source_list_filt, target_list_filt, sparsity_weights_filt = load_filter_params(params_dir)
    oh_conn_ipsi_old=conn_dir+"oh_conn_ipsi_old_211.csv"
    oh_conn_contra_old=conn_dir+"oh_conn_contra_old_211_filled_in.csv"
    oh_conn_ipsi_contra_old = [oh_conn_ipsi_old, oh_conn_contra_old]

    oh_conn_ipsi_new=conn_dir+"oh_conn_ipsi_new_211.csv"
    oh_conn_contra_new=conn_dir+"oh_conn_contra_new_211_filled_in.csv"
    oh_conn_ipsi_contra_new = [oh_conn_ipsi_new, oh_conn_contra_new]

    knox_conn_ipsi_old=conn_dir+'knox_conn_ipsi_old_211.csv'
    knox_conn_contra_old=conn_dir+'knox_conn_contra_old_211_filled_in.csv'
    knox_conn_ipsi_contra_old=[knox_conn_ipsi_old, knox_conn_contra_old]

    knox_conn_ipsi_new=conn_dir+"knox_conn_ipsi_new_211.csv"
    knox_conn_contra_new=conn_dir+"knox_conn_contra_new_211_filled_in.csv"
    knox_conn_ipsi_contra_new = [knox_conn_ipsi_new, knox_conn_contra_new]

    connectome_lists = [oh_conn_ipsi_contra_old, oh_conn_ipsi_contra_new, knox_conn_ipsi_contra_old, knox_conn_ipsi_contra_new]

    ###load conversion table
    allen_rgn_conversion_table_path="../preprocessed/allen_template_inputs/allen_ccfv3_tree_wang_2020_s2.xlsx"
    conversion_table = pd.read_excel(allen_rgn_conversion_table_path, header=1)

    ##suffixes to write out processed .pkl
    output_prefixes=["params_homogeneous_weights_old","params_homogeneous_weights_after_qc_0127", "params_regionalized_voxel_weights_old", "params_regionalized_voxel_weights_after_qc_0127"]
    output_prefixes_64pct=[output_prefixes[i] + "_64pct" for i in range(len(output_prefixes))]
    i=0
    for conn_list in connectome_lists:
        prefix = output_prefixes[i]
        prefix_thresh = output_prefixes_64pct[i]
        process_regionalized_connectome(conversion_table,
                conn_list, prefix, prefix_thresh, output_dir)
        i+=1
    
