## SIR_mouse_regional_tuned

The following repository contains the scripts to run a parameter-tuned version of the SIR model from Rahayel et al., 2022. The goal of these additions is to help find the parameters and clearance genes that best predict the network-spreading of alpha synuclein. In doing so, we will explore how these parameters differ when predicting IHC-derived aSyn accumulation vs. MRI-derived atrophy, and predicting the spreading of aSyn from a PD-like caudoputamen epicentre vs. a hippocampal positive control (see Tullo et al., 2024 for a description of the experiments).

Each subdirectory contains a separate README file for the corresponding scripts. All raw data accessed within ```data_processing``` can be made available upon request.

We recommend creating a set of jobs using the scripts provided in ```batch_run``` and then submitting these jobs using ```qbatch```. 

See ```pyproject.toml``` for a list of python package dependencies.

TODO: use `uvr` to manage R package dependencies.


### Order of operations

* ```data_processing```
* ```algorithm```
* ```batch_run```
* ```visualization```

### How to run the SIR model 

*A priori* definitions of a single spreading direction and clearance gene (or `-c None` for uniform clearance) are required for both of these scripts. 

* ```abm_optuna_general.py```: uses *Optuna* Bayesian parameter tuning to gridsearch the numerical parameter space. (for example runs across different clearance gene/output configurations, see sir_command_files). Additional parameters added (not shared with ```abm_clearance_genes.py```): 
  * m: output map (either IHC or MRI-derived regionalized, whole-brain data)
  * a: atrophy (True or False)
  * x: suffix (for outputted Optuna .csv file that shows how parameter tuning affects simulated vs. empirical correlation in each optimization trial)
 

* ```abm_clearance_genes.py```: runs the SIR model simulations for a single parameter configuration (usually after optimization)
  * t: number of timesteps
  * p: path for `params` file that packages SIR inputs/outputs (see ```data_processing```)
  * g: path for clearance gene data
  * o: output directory for simulated I agents/atrophy
  * r: retrograde (True/False): for labelling purposes
  * e: epsilon determining numerical convergence for initializing S fraction (TODO: REMOVE DEPENDENCY, use closed-form solution)
  * d: delta-t for discretizing diff-eq solution (equations unchanged from Rahayel et al., 2022)
  * S: numerical index of the epicentre region (35: CP, 40: DG, 24: CA1)
  * c: clearance gene name (pulls from clearance gene directory in g)
  * v: velocity (along edges) (unitless; 0-500)
  * s: spreading rate (probability of protein exiting region and entering edge)
  * i: injection amount of aSyn (unitless)
  * k1: presynaptic neuronal loss (0-1): ONLY use when predicting atrophy 
  * k2: postsynaptic neuronal loss (deafferentation) parameter (0-1): ONLY use when predicting atrophy 

Example usage (within the ```batch_run``` directory): 

````
 python3 ../algorithm/abm_clearance_genes.py \
  -t 1000 \
  -p "../../derivatives/SIR_inputs/params_regionalized_voxel_weights_after_qc_0413_64pct.pkl" \
  -g "../../derivatives/yohan_ge_filt/mr10vv0.2/" \
  -o "../simulations/" \
  -r "False" \
  -e 1e-05 \
  -d 0.1 \
  -S 35 \
  -c "Dnajc17" \
  -v  104.87116708452713 \
  -s 0.02122582417807315 \
  -i 97.7903636896988 \
  -k1 0.03066890459365682 \
  -k2 0.42327816662213974
````

