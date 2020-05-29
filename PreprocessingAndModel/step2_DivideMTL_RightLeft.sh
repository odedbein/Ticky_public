#!/bin/sh

#  DivideMTL_RightLeft.sh
#  
#
#  Created by Oded Bein on 10/21/16.
#
proj_dir='/Volumes/Oded/Bein/TickyReanalysis'
declare -a subjects=(AB AD AK AR AT BW CR DH DM EB JA JD JG JM JR JW JR JW KZ LD SB YE)
scripts_dir=$proj_dir/scripts
Rside=$scripts_dir/Rside.nii.gz
Lside=$scripts_dir/Lside.nii.gz
declare -a reg=(Para Peri Ent)

#LOOP FOR ALL SUBJECTS - create the nii files
for subj in "${subjects[@]}"
do
subj_dir=$proj_dir/SubData/$subj
subj_roi_dir=$subj_dir/ROIs
for roi in "${reg[@]}"
do
fslmaths $subj_roi_dir/$roi.nii.gz -mul $Rside $subj_roi_dir/r$roi.nii.gz
fslmaths $subj_roi_dir/$roi.nii.gz -mul $Lside $subj_roi_dir/l$roi.nii.gz
done
done





