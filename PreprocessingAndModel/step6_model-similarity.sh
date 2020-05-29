#!/bin/sh

#
#  Created by Oded Bein on 4/14/16.
#

#DEFENITIONS
/Volumes/Oded/Bein/TickyReanalysis/scripts/TickyReanalysis_globals.sh
declare -a subjects=(AB AD AK AR AT BW CR DH DM EB JA JD JG JM JR JW JR JW KZ LD SB YE)
#finished: AB LD SB DM CR BW DH JA AT AD EB AK; Joanna: JW KZ YE AR are currently running
#seahorse: JD JG JM JR are currently running
#notes for participant JG: had 119 volumes in one of the scans
#LOOP FOR ALL SUBJECTS
for subj in "${subjects[@]}"
do

echo "analyzing subject $subj"
subj_design_dir=$PROJECT_DIR/$SUBJECTS_DIR/$subj/$DESIGN_DIR
if [ ! -d $subj_design_dir ]; then
mkdir $subj_design_dir
fi
subj_input_dir=$PROJECT_DIR/$SUBJECTS_DIR/$subj/$DATA_DIR

subj_output_dir=$PROJECT_DIR/$SUBJECTS_DIR/$subj/$OUTPUT_DIR
if [ ! -d $subj_output_dir ]; then
mkdir $subj_output_dir
fi
subj_Tmaps_dir=$PROJECT_DIR/$SUBJECTS_DIR/$subj/$Tmaps_DIR
if [ ! -d $subj_Tmaps_dir ]; then
mkdir $subj_Tmaps_dir
fi
subj_reg_dir=$PROJECT_DIR/$SUBJECTS_DIR/$subj/regressors

    for r in {7..8}; do
    echo "creating single-trial encoding models $subj run ${r}"
    for t in {1..27}; do # number of trials in each run
    trial=$(((r-1)*27+t))
    #create the design for the cue
    cat $FSL_TEMPLATES/single_trial_cue_model.fsf \
    | sed "s|Output_single_model|$subj_output_dir/Cue$((trial))temp|g" \
    | sed "s|curr_sess_data|$subj_input_dir/all_mtc_ps_run${r}|g" \
    | sed "s|reg_dir|$subj_reg_dir|g" \
    | sed "s|NumModel|$((trial))Model|g" \
    | sed "s|CueNum|Cue$trial|g" \
    > $subj_design_dir/sim_run${r}_Cue$trial.fsf

    #run feat - cue model
    f=`echo Cue$trial`
    if [ ! -d "$subj_output_dir/${f}temp.feat" ]; then

    echo "running FEAT pattern similarity model on ${f}"
    feat $subj_design_dir/sim_run${r}_${f}.fsf &
    sleep 60
    $PROJECT_DIR/scripts/wait-for-feat.sh $subj_output_dir/${f}temp.feat
    #ln -s $subj_input_dir/${f}.feat/reg/ $subj_output_dir/${f}.feat/reg
    fi

    #copy and clear files
    sleep 20
    mkdir $subj_output_dir/Cue$((trial)).feat
    cp $subj_output_dir/Cue$((trial))temp.feat/design.fsf $subj_output_dir/Cue$((trial)).feat/design.fsf
    cp $subj_output_dir/Cue$((trial))temp.feat/report_log.html $subj_output_dir/Cue$((trial)).feat/report_log.html
    cp $subj_output_dir/Cue$((trial))temp.feat/design.png $subj_output_dir/Cue$((trial)).feat/design.png
    cp $subj_output_dir/Cue$((trial))temp.feat/stats/pe1.nii.gz $subj_output_dir/Cue$((trial)).feat/pe1.nii.gz
    cp $subj_output_dir/Cue$((trial))temp.feat/stats/tstat1.nii.gz $subj_Tmaps_dir/tstat_Cue$((trial)).nii.gz

    rm -r $subj_output_dir/Cue$((trial))temp.feat

    #create the design for the image
    cat $FSL_TEMPLATES/single_trial_image_model.fsf \
    | sed "s|Output_single_model|$subj_output_dir/Im$((trial))temp|g" \
    | sed "s|curr_sess_data|$subj_input_dir/all_mtc_ps_run${r}|g" \
    | sed "s|reg_dir|$subj_reg_dir|g" \
    | sed "s|NumModel|$((trial))Model|g" \
    | sed "s|ImNum|Im$trial|g" \
    > $subj_design_dir/sim_run${r}_Im$trial.fsf

    #run feat - image model
    f=`echo Im$trial`
    if [ ! -d "$subj_output_dir/${f}temp.feat" ]; then

    echo "running FEAT pattern similarity model on ${f}"
    feat $subj_design_dir/sim_run${r}_${f}.fsf &
    sleep 60
    $PROJECT_DIR/scripts/wait-for-feat.sh $subj_output_dir/${f}temp.feat
    #ln -s $subj_input_dir/${f}.feat/reg/ $subj_output_dir/${f}.feat/reg
    fi

    #copy and clear files
    sleep 20
    mkdir $subj_output_dir/Im$((trial)).feat
    cp $subj_output_dir/Im$((trial))temp.feat/design.fsf $subj_output_dir/Im$((trial)).feat/design.fsf
    cp $subj_output_dir/Im$((trial))temp.feat/report_log.html $subj_output_dir/Im$((trial)).feat/report_log.html
    cp $subj_output_dir/Im$((trial))temp.feat/design.png $subj_output_dir/Im$((trial)).feat/design.png
    cp $subj_output_dir/Im$((trial))temp.feat/stats/pe1.nii.gz $subj_output_dir/Im$((trial)).feat/pe1.nii.gz
    cp $subj_output_dir/Im$((trial))temp.feat/stats/tstat1.nii.gz $subj_Tmaps_dir/tstat_Im$((trial)).nii.gz

    rm -r $subj_output_dir/Im$((trial))temp.feat

    done #trials loop
    done #runs loop

done #subjects loop



