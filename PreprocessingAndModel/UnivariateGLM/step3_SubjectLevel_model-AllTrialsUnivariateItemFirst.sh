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

#subjects that did the layout first:
#AB AD AT BW EB JA JG   JW KZ YE

#subjects that did the item first:
#AR CR DH DM JD     JM JR SB


#LD: lay/item/lay/item/lay/item/lay/item
#AK: lay/item/lay/item/lay/item/lay

#/Volumes/Oded/Bein/TickyReanalysis/scripts/PreprocessingAndModel/step0_TickyReanalysis_globals.sh
declare -a subjects=(AB AD AK AR AT BW CR DH DM EB JA JD JG JM JR JW JR JW KZ LD SB YE)
#finished: AB LD SB DM CR BW DH JA AT AD EB AK; Joanna: JW KZ YE AR are currently running
#seahorse: JD JG JM JR are currently running
#notes for participant JG: had 119 volumes in the
#LOOP FOR ALL SUBJECTS
for subj in "${subjects[@]}"
do

echo "analyzing subject $subj"

subj_output_dir=$PROJECT_DIR/$SUBJECTS_DIR/$subj/$analysis_dir/AllTrials
subj_input_dir=$subj_output_dir

subj_design_dir=$subj_output_dir/$DESIGN_DIR

## no need to run this part bc I did it for the analysis collapsed on tasks:

##subj_reg_dir=$PROJECT_DIR/$SUBJECTS_DIR/$subj/regressorsUnivariate/All
#
##fsl is dumb in the sense that it doesn't run high level analysis without registration. To go around that, we do registration, but make sure it doesn't do any transforamtion by having the identity matrix as the transformation matrx. To prevent from doing any interpulation, we make the stanrdard ref image the mean_func image. for ref on this proceedure, see Mumdford's video: https://www.youtube.com/watch?v=U3tG7JMEf7M
#
##copy the necessary files for the dummy registration:
#for r in {1..10}; do
#
#subj_reg_dir=$subj_output_dir/run${r}.feat/reg
#if [ ! -d $subj_reg_dir ]; then
#mkdir $subj_reg_dir
#fi
#
##copy the identity matrix from fsl's general folder:
#cp $FSLDIR/etc/flirtsch/ident.mat $subj_reg_dir/example_func2standard.mat
#
##copy the mean_func image from the specific run to be the ref image for the dummy registration
#cp $subj_output_dir/run${r}.feat/mean_func.nii.gz $subj_reg_dir/standard.nii.gz
#
#done #runs loop

echo "creating a model for $subj"

#create the design file:
cat $FSL_TEMPLATES/univarAllSubjLevelItemFirst.fsf \
| sed "s|subjX|$subj|g" \
> $subj_design_dir/univar_SubjLevelItemFirst.fsf

#run feat -  model
if [ ! -d "$subj_output_dir/SubjLevelTasks.feat" ]; then

echo "running FEAT model on $subj"
feat $subj_design_dir/univar_SubjLevelItemFirst.fsf &
sleep 60
$PROJECT_DIR/scripts/PreprocessingAndModel/wait-for-feat.sh $subj_output_dir/SubjLevelTasks.gfeat
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


done #subjects loop



