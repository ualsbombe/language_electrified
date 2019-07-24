#! /bin/bash

## export variables
export SUBJECTS_DIR=/home/lau/analyses/fieldtrip_datasets/language_MEG/data/freesurfer
export SUBJECT=sub-V1002_space-CTF

nii_path=/home/lau/analyses/fieldtrip_datasets/language_MEG/data/$SUBJECT/anat
cd $nii_path
filename=${SUBJECT}_T1w.nii
recon-all -subjid $SUBJECT -i $filename -openmp 32

#! /bin/bash

export SUBJECTS_DIR=/home/lau/analyses/fieldtrip_datasets/language_MEG/data/freesurfer
export SUBJECT=sub-V1002_space-CTF

recon-all -subjid $SUBJECT -all -notal-check -openmp 32
