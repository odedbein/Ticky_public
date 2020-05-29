   #!/bin/bash -e -u
# Author: Alexa Tompary
# List of relative file paths starting in each subject directory
echo 'running globals'
PROJECT_DIR=/Volumes/Oded/Bein/TickyReanalysis
SUBJECTS_DIR=SubData
FSL_TEMPLATES=$PROJECT_DIR/fsl_templates
#folders of each subject
DATA_DIR=data
DESIGN_DIR=design
#PREPROC_DIR=analysis/preproc
#TASK_DIR=analysis/task
#SM_DIR=analysis/sm_sm6
OUTPUT_DIR=analysis_files
Tmaps_DIR=single_trial_Tmaps
ROI_DIR=ROIs
MPRAGE_DIR=mprage

#KERNEL_SIZE=8 # size of ROI kernel
#SMOOTH_SIZE=6 # size of smoothing kernel

#master_run=localizer_run1 # align all ROIs and statistical maps to this run's native space
#standard=$PROJECT_DIR/design/standard

#ALL_SUBJECTS=`ls -1d subjects/* | cut -c 10-12`