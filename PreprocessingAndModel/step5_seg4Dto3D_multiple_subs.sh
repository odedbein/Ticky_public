#!/bin/sh

#  Afni2NII_multiple_subs.sh
#  
#
#  Created by Oded Bein on 4/14/16.
#

#DEFENITIONS

proj_dir='/Volumes/Oded/Bein/TickyReanalysis';
declare -a subjects=(AB AD AK AR AT BW CR DH DM EB JA JD JG JM JR JW JR JW KZ LD SB YE)


#LOOP FOR ALL SUBJECTS - create the nii files
for subj in "${subjects[@]}"
do
    subj_dir=$proj_dir/SubData/$subj
    if [ ! -d $subj_dir/data ]; then
    mkdir $subj_dir/data;
    fi
    echo "segmenting to sessions subject $subj"
    file_name=$subj_dir/all_mtc_ps.nii.gz
    fslsplit $file_name $subj_dir/data/all_mtc_ps

    ls $subj_dir/data > $subj_dir/file_names;

    for sess in {1..10}; do
        files_names3d=
        echo "session $sess"
        first_scan=$((120*(sess-1)+1))
        last_scan=$((first_scan+119))
        for line in $(seq $first_scan $last_scan); do
            files_names3d+="$subj_dir/data/";
            files_names3d+=$(sed -n "${line}p" "$subj_dir/file_names");
            files_names3d+=" ";
        done

        fslmerge -t $subj_dir/data/all_mtc_ps_run$sess $files_names3d
        rm $files_names3d
    done
done


#Participant JG had a missing vlume at the end of the first scan, so need to group the volumes differently.
#sess 1:
sess=1
files_names3d=
echo "session $sess"
first_scan=$((120*(sess-1)+1))
last_scan=$((first_scan+118))
for line in $(seq $first_scan $last_scan); do
files_names3d+="$subj_dir/data/";
files_names3d+=$(sed -n "${line}p" "$subj_dir/file_names");
files_names3d+=" ";
done

fslmerge -t $subj_dir/data/all_mtc_ps_run$sess $files_names3d
rm $files_names3d

#sess 2-10:
for sess in {2..10}; do
files_names3d=
echo "session $sess"
first_scan=$((120*(sess-1)))
last_scan=$((first_scan+119))
for line in $(seq $first_scan $last_scan); do
files_names3d+="$subj_dir/data/";
files_names3d+=$(sed -n "${line}p" "$subj_dir/file_names");
files_names3d+=" ";
done

fslmerge -t $subj_dir/data/all_mtc_ps_run$sess $files_names3d
rm $files_names3d
done




