# Ticky_public

These are the scripts I used to analyze the data for Bein et al., 2020 (accepted), NatComms: *"Mnemonic prediciton errors bias hippocampal states".
I did not clean them well enough, and also - this project is from the begining of my PhD. So the code is messy. Feel free to contact at oded.bein@nyu.edu
Relevant data is abvailable on: https://osf.io/re2wd/
* note that this is a reanalysis of data previously published in Duncan et al., 2012, Hippocampus: *Evidence for Area CA1 as a Match/Mismatch Detector: A High-Resolution fMRI Study of the Human Hippocampus".
* Some of the scripts here were used to prepare the data (e.g., convert the data to nifti, break rois to hemispheres etc.). I put them here just for documentation, but as the data that I uploaded to osf is already prepared for analysis, these scripts are unnecessary. They are marked with:  "Already done, unecessary now."
* I used participants folders with initials as folder names. In the osf data base, I re-named the folders to be s01-s020. But the scripts are with the subjects initials. So if one downloads the data and want to use these scripts to analyse these data - change the subject names.
* in the ms, subject 20 was removed from analysis due to entorhinal dropout. I uploaded their data, because no reason not to use if analyzing other ROIs.

## PreprocessingAndModel folder:

step0_TickyReanalysis_globals: definition file, to be run before starting the analysis/called by other scripts

step1_Afni2NII_multiple_subs: I converted the files to nifti from KD's Afni format. Already done, unecessary now.
step2_DivideMTL_RightLeft: divide MTL regions to right and left hemisphere. Already done, unecessary now.
step3_remake_behavioral_files: for some participants, I needed to re-create their behavior files to be in the standard format because they were in a different format in KD's folder. Already done, unecessary now.

step4_create_regressors_files: Creates the regressors for the LSS and the univariate GML models. In the osf repository I uploaded the regressors.

step5_seg4Dto3D_multiple_subs: I think I had the data all runs concatinated. This breaks the data to runs.  Already done, unecessary now.

step6_model-similarity: runs the LSS models to obtain the single-trial t-stat. Each trial had a cue part and an image part. Creates one t-stat for the cue and one for the image per each trial.

wait_for_feat: a script that waits for feat to end running. The modeling code and other codes use it.

