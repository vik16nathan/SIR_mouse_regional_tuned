#!/bin/bash

#GOAL: simulate SIR for the top clearance gene/parameter configurations for each pathological dataset


###ATROPHY####

##CP, MsPFF vs. PBS, anterograde###

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


##DG, HuPFF vs. PBS, retrograde###
python3 ../algorithm/abm_clearance_genes.py \
  -t 1000 \
  -p "../../derivatives/SIR_inputs/params_regionalized_voxel_weights_after_qc_0413_64pct_retro.pkl" \
  -g "../../derivatives/yohan_ge_filt/mr10vv0.2/" \
  -o "../simulations/" \
  -r "True" \
  -e 1e-05 \
  -d 0.1 \
  -S 40 \
  -c "Twist2" \
  -v  1.0862247470116264 \
  -s 0.027270129876801334 \
  -i 35.5593507054226 \
  -k1 0.7296383313819099 \
  -k2 0.059192558911420556

###pSyn###
##CP, 24 mpi, anterograde###
python3 ../algorithm/abm_clearance_genes.py \
  -t 1000 \
  -p "../../derivatives/SIR_inputs/params_regionalized_voxel_weights_after_qc_0413_64pct.pkl" \
  -g "../../derivatives/yohan_ge_filt/mr10vv0.2/" \
  -o "../simulations/" \
  -r "False" \
  -e 1e-05 \
  -d 0.1 \
  -S 35 \
  -c "Zmat4" \
  -v 2.723990344784835 \
  -s 0.05812310477637486 \
  -i 82.23688396384566

##CA1, 24 mpi, retrograde###
python3 ../algorithm/abm_clearance_genes.py \
  -t 1000 \
  -p "../../derivatives/SIR_inputs/params_regionalized_voxel_weights_after_qc_0413_64pct_retro.pkl" \
  -g "../../derivatives/yohan_ge_filt/mr10vv0.2/" \
  -o "../simulations/" \
  -r "True" \
  -e 1e-05 \
  -d 0.1 \
  -S 24 \
  -c "Ctnnbip1" \
  -v  0.29188965403566514 \
  -s 0.0006938150874085919 \
  -i 54.177718195890705

