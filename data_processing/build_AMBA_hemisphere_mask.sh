#!/bin/bash

module load minc-toolkit-v2
module load ANTs


##define original mask and annotations (same labels for RH/LH)
annotations="../preprocessed/allen_template_inputs/AMBA_25um_resampled_50um_int.mnc"
allenmask="../preprocessed/allen_template_inputs/mask_50um.mnc"

# Find size of image
mincinfo $allenmask

# Make image half size; starting x coordinate is from mincinfo start + (length/2)*step
mincresample -step 0.05 0.05 0.05 -start 0 -7.9 -2.7 -nelements 114 264 160 $allenmask /tmp/test.mnc -clobber

# Fill half-size image with 1's
minccalc -expression '1' /tmp/test.mnc /tmp/test2.mnc -clobber

# Reshape all-ones image back to original size
antsApplyTransforms -d 3 -i /tmp/test2.mnc -r $allenmask -o /tmp/test3.mnc

# Multiply labels by this to zero out half
minccalc -expression 'A[0]*A[1]' $annotations /tmp/test3.mnc /tmp/half-labels.mnc -clobber

mv /tmp/half-labels.mnc ../preprocessed/allen_template_inputs/AMBA_25um_resampled_50um_int_RH.mnc