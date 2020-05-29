#!/bin/sh

#  Afni2NII_multiple_subs.sh
#  
#
#  Created by Oded Bein on 4/14/16.
#

#DEFENITIONS
PROJECT_DIR=/Volumes/Oded/Bein/TickyReanalysis
analysis_dir=UnivariateGLM
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

#/Volumes/Oded/Bein/TickyReanalysis/scripts/PreprocessingAndModel/step0_TickyReanalysis_globals.sh
declare -a subjects=(AB AD AK AR AT BW CR DH DM EB JA JD JG JM JR JW JR JW KZ LD SB YE)
#finished: AB LD SB DM CR BW DH JA AT AD EB AK; Joanna: JW KZ YE AR are currently running
#seahorse: JD JG JM JR are currently running
#notes for participant JG: had 119 volumes in the
#LOOP FOR ALL SUBJECTS
for subj in "${subjects[@]}"
do

echo "analyzing subject $subj"

subj_input_dir=$PROJECT_DIR/$SUBJECTS_DIR/$subj/$DATA_DIR

subj_output_dir=$PROJECT_DIR/$SUBJECTS_DIR/$subj/$analysis_dir
if [ ! -d $subj_output_dir ]; then
mkdir $subj_output_dir
fi

subj_output_dir=$PROJECT_DIR/$SUBJECTS_DIR/$subj/$analysis_dir/AllTrials
if [ ! -d $subj_output_dir ]; then
mkdir $subj_output_dir
fi


subj_design_dir=$subj_output_dir/$DESIGN_DIR
if [ ! -d $subj_design_dir ]; then
mkdir $subj_design_dir
fi


subj_reg_dir=$PROJECT_DIR/$SUBJECTS_DIR/$subj/regressorsUnivariate/All

    for r in {1..10}; do
    echo "creating a model for $subj run ${r}"

    #create the design for the cue
    cat $FSL_TEMPLATES/univarAll.fsf \
    | sed "s|subj_output_dir|$subj_output_dir/run${r}|g" \
    | sed "s|curr_sess_data|$subj_input_dir/all_mtc_ps_run${r}|g" \
    | sed "s|reg_dir|$subj_reg_dir|g" \
    | sed "s|runX|run${r}|g" \
    > $subj_design_dir/univar_run${r}.fsf

    #run feat -  model
#f=`echo Cue$trial`
    if [ ! -d "$subj_output_dir/run${r}.feat" ]; then

    echo "running FEAT model on run${r}"
    feat $subj_design_dir/univar_run${r}.fsf &
    sleep 60
    $PROJECT_DIR/scripts/PreprocessingAndModel/wait-for-feat.sh $subj_output_dir/run${r}.feat
    #ln -s $subj_input_dir/${f}.feat/reg/ $subj_output_dir/${f}.feat/reg
    fi

    #copy and clear files
    sleep 20
#    mkdir $subj_output_dir/Cue$((trial)).feat
#    cp $subj_output_dir/Cue$((trial))temp.feat/design.fsf $subj_output_dir/Cue$((trial)).feat/design.fsf
#    cp $subj_output_dir/Cue$((trial))temp.feat/report_log.html $subj_output_dir/Cue$((trial)).feat/report_log.html
#    cp $subj_output_dir/Cue$((trial))temp.feat/design.png $subj_output_dir/Cue$((trial)).feat/design.png
#    cp $subj_output_dir/Cue$((trial))temp.feat/stats/pe1.nii.gz $subj_output_dir/Cue$((trial)).feat/pe1.nii.gz
#    cp $subj_output_dir/Cue$((trial))temp.feat/stats/tstat1.nii.gz $subj_Tmaps_dir/tstat_Cue$((trial)).nii.gz
#
#    rm -r $subj_output_dir/Cue$((trial))temp.feat

    done #runs loop

done #subjects loop



