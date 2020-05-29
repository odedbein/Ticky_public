#!/bin/sh

#  Entorhinal_dropout.sh
#  
#
#  Created by Oded Bein on 7/19/17.
#
#some definitions:
declare -a subjects=(AB AD AK AR AT BW CR DH DM EB JA JD JG JM JR JW KZ LD SB YE)
#declare -a subjects=(AD AK AR AT BW CR DH DM EB JA JD JG JM JR JW JR JW KZ LD SB YE)
proj_dir='/Volumes/Oded/Bein/TickyReanalysis'


#01. copy all files pre-epi masking, and convert them to nii files
hd_dir=/Volumes/3/Ticky/MTL/SubData

#LOOP FOR ALL SUBJECTS - create the nii files
for subj in "${subjects[@]}"
do
subj_Afni_dir=$proj_dir/AfniFiles/$subj
subj_dir=$proj_dir/SubData/$subj
subj_roi_dir=$subj_dir/ROIs

file_name=$hd_dir/${subj}/ROIs/parahip/Full_Cort_${subj}_EPI_al+orig.BRIK
if [ ! -f $file_name ]; then
echo 'missing EPI_al for subject' $subj
else
new_file_name=$subj_Afni_dir/Full_Cort_${subj}_EPI_al+orig.BRIK
cp $file_name $new_file_name
file_name=$hd_dir/${subj}/ROIs/parahip/Full_Cort_${subj}_EPI_al+orig.HEAD
new_file_name=$subj_Afni_dir/Full_Cort_${subj}_EPI_al+orig.HEAD
cp $file_name $new_file_name
nii_file_name=$subj_roi_dir/Full_Cort_EPI_al.nii.gz

echo 'copying EPI_al and converting to nii for subject' $subj
3dAFNItoNIFTI -prefix $nii_file_name $new_file_name
fi
done

#02. extract only the entorhinal from the mask, sperately for right and left:
for subj in "${subjects[@]}"
do
subj_dir=$proj_dir/SubData/$subj
subj_roi_dir=$subj_dir/ROIs

fslmaths $subj_roi_dir/Full_Cort_EPI_al -thr 6 -uthr 6 -bin $subj_roi_dir/lEnt_pre_epi_mask
fslmaths $subj_roi_dir/Full_Cort_EPI_al -thr 3 -uthr 3 -bin $subj_roi_dir/rEnt_pre_epi_mask
done


#02.1. for some subjects, ent was "3" for right and left, so divide:
declare -a subjects=(JM JR JW KZ LD SB)
Rside=$proj_dir/scripts/PreprocessingAndModel/Rside.nii.gz
Lside=$proj_dir/scripts/PreprocessingAndModel/Lside.nii.gz

for subj in "${subjects[@]}"
do
subj_dir=$proj_dir/SubData/$subj
subj_roi_dir=$subj_dir/ROIs
fslmaths $subj_roi_dir/Full_Cort_EPI_al -thr 3 -uthr 3 -bin $subj_roi_dir/Ent_pre_epi_mask
fslmaths $subj_roi_dir/Ent_pre_epi_mask -mul $Rside $subj_roi_dir/rEnt_pre_epi_mask
fslmaths $subj_roi_dir/Ent_pre_epi_mask -mul $Lside $subj_roi_dir/lEnt_pre_epi_mask
done

for subj in "${subjects[@]}"
do
subj_dir=$proj_dir/SubData/$subj
subj_roi_dir=$subj_dir/ROIs
fslview $subj_roi_dir/Full_Cort_EPI_al $subj_roi_dir/lEnt_pre_epi_mask -l Red $subj_roi_dir/rEnt_pre_epi_mask -l Green &
done


#03. write information to a text file
#lEnt:
#copy this line to the command line, then press ctrl+d to close the file, then run the loop:
cat > $proj_dir/Ent_dropout/lEnt_pre_epi_mask_size.txt

for subj in "${subjects[@]}"
do
subj_dir=$proj_dir/SubData/$subj
subj_roi_dir=$subj_dir/ROIs
cluster -i $subj_roi_dir/lEnt_pre_epi_mask -t 1 >> $proj_dir/dropout_analysis/lEnt_dropout/lEnt_pre_epi_mask_size.txt
done

#copy this line to the command line, then press ctrl+d to close the file, then run the loop:
cat > $proj_dir/Ent_dropout/lEnt_size.txt
for subj in "${subjects[@]}"
do
subj_dir=$proj_dir/SubData/$subj
subj_roi_dir=$subj_dir/ROIs
cluster -i $subj_roi_dir/lEnt -t 1 >> $proj_dir/dropout_analysis/lEnt_dropout/lEnt_size.txt
done

#rEnt
#copy this line to the command line, then press ctrl+d to close the file, then run the loop:
cat > $proj_dir/dropout_analysis/rEnt_dropout/rEnt_pre_epi_mask_size.txt

for subj in "${subjects[@]}"
do
subj_dir=$proj_dir/SubData/$subj
subj_roi_dir=$subj_dir/ROIs
cluster -i $subj_roi_dir/rEnt_pre_epi_mask -t 1 >> $proj_dir/dropout_analysis/rEnt_dropout/rEnt_pre_epi_mask_size.txt
done

#copy this line to the command line, then press ctrl+d to close the file, then run the loop:
cat > $proj_dir/dropout_analysis/rEnt_dropout/rEnt_size.txt
for subj in "${subjects[@]}"
do
subj_dir=$proj_dir/SubData/$subj
subj_roi_dir=$subj_dir/ROIs
cluster -i $subj_roi_dir/rEnt -t 1 >> $proj_dir/dropout_analysis/rEnt_dropout/rEnt_size.txt
done

