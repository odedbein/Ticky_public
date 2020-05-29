#!/bin/sh

#  Afni2NII_multiple_subs.sh
#  
#
#  Created by Oded Bein on 4/14/16.
#

#DEFENITIONS

proj_dir='/Volumes/Oded/Bein/TickyReanalysis';
#declare -a subjects=(AB AD AK AR AT BW CR DH DM EB JA JD JG JM JR JW JR JW KZ LD SB YE)
declare -a subjects=(AD AK AR AT BW CR DH DM EB JA JD JG JM JR JW JR JW KZ LD SB YE)

#LOOP FOR ALL SUBJECTS - create the nii files
for subj in "${subjects[@]}"
do
    subj_Afni_dir=$proj_dir/AfniFiles/$subj
    subj_dir=$proj_dir/SubData/$subj
    if [ ! -d $subj_dir ]; then
    mkdir $subj_dir;
    fi
    subj_roi_dir=$subj_dir/mprage
    if [ ! -d $subj_roi_dir ]; then
    mkdir $subj_roi_dir;
    fi
flag_file=true
file_name=$subj_Afni_dir/anat_oblique_susan_al_ns+orig.BRIK
if [ ! -f $file_name ]; then
    file_name=$subj_Afni_dir/anat_oblique_susan_alepi_ns+orig.BRIK
    if [ ! -f $file_name ]; then
    echo 'missing file for subject' $subj
    flag_file=false
    fi
fi

if [ $flag_file ]; then
    new_file_name=$subj_roi_dir/anat_oblique_susan_al_ns.nii.gz
    3dAFNItoNIFTI -prefix $new_file_name $file_name
    fslreorient2std $new_file_name $subj_roi_dir/anat_oblique_susan_al_ns_std
fi

done

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

#
#file_name=$subj_Afni_dir/mean_EPI+orig.BRIK
#new_file_name=$subj_roi_dir/mean_EPI.nii.gz
#3dAFNItoNIFTI -prefix $new_file_name $file_name
#fslreorient2std $new_file_name $subj_roi_dir/mean_EPI_std

#
#file_name=$subj_Afni_dir/Para+orig.BRIK
#new_file_name=$subj_roi_dir/Para.nii.gz
#3dAFNItoNIFTI -prefix $new_file_name $file_name
#
#file_name=$subj_Afni_dir/Ent+orig.BRIK
#new_file_name=$subj_roi_dir/Ent.nii.gz
#3dAFNItoNIFTI -prefix $new_file_name $file_name
#
#file_name=$subj_Afni_dir/Peri+orig.BRIK
#new_file_name=$subj_roi_dir/Peri.nii.gz
#3dAFNItoNIFTI -prefix $new_file_name $file_name
#
#file_name=$subj_Afni_dir/Non_MTL+orig.BRIK
#new_file_name=$subj_roi_dir/Non_MTL.nii.gz
#3dAFNItoNIFTI -prefix $new_file_name $file_name
#
#    file_name=$subj_Afni_dir/$subj'_all_mtc_ps+orig.BRIK'
#    new_file_name=$subj_dir/all_mtc_ps.nii.gz
#    3dAFNItoNIFTI -prefix $new_file_name $file_name
#
#    file_name=$subj_Afni_dir/'Hip15_oblique_al+orig.BRIK'
#    new_file_name=$subj_roi_dir/Hip15_oblique_al.nii.gz
#    3dAFNItoNIFTI -prefix $new_file_name $file_name
#
#    file_name=$subj_Afni_dir/'Ant_Post15+orig.BRIK'
#    new_file_name=$subj_roi_dir/Ant_Post15.nii.gz
#    3dAFNItoNIFTI -prefix $new_file_name $file_name
#done

for subj in "${subjects[@]}"
do

subj_dir=$proj_dir/SubData/$subj
subj_roi_dir=$subj_dir/ROIs
file_name=$subj_roi_dir/mean_EPI.nii.gz
new_file_name=$subj_dir/mprage/mean_EPI.nii.gz
mv $file_name $new_file_name

file_name=$subj_roi_dir/mean_EPI_std.nii.gz
new_file_name=$subj_dir/mprage/mean_EPI_std.nii.gz
mv $file_name $new_file_name
done

for subj in "${subjects[@]}"
do
subj_dir=$proj_dir/SubData/$subj/mprage
file_name=$subj_dir/anat_oblique_susan.nii.gz
fslreorient2std $file_name $subj_dir/anat_oblique_susan_std
done



#LOOP FOR ALL SUBJECTS - create the
for subj in "${subjects[@]}"
do

    subj_dir=$proj_dir/SubData/$subj
    subj_roi_dir=$subj_dir/ROIs
    fslmaths $subj_roi_dir/Hip15_oblique_al.nii.gz -bin $subj_roi_dir/All_Hip15.nii.gz
    fslmaths $subj_roi_dir/All_Hip15.nii.gz -mul $subj_roi_dir/Ant_Post15.nii.gz $subj_roi_dir/Hip15_AntMidPost.nii.gz

done




