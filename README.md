# Ticky_public

These are the scripts I used to analyze the data for Bein et al., 2020 (accepted), NatComms: *"Mnemonic prediciton errors bias hippocampal states".
I did not clean them well enough, and also - this project is from the begining of my PhD. So the code is messy. Feel free to contact at oded.bein@nyu.edu
Relevant data is abvailable on: https://osf.io/re2wd/
* note that this is a reanalysis of data previously published in Duncan et al., 2012, Hippocampus: *Evidence for Area CA1 as a Match/Mismatch Detector: A High-Resolution fMRI Study of the Human Hippocampus".
* Some of the scripts here were used to prepare the data (e.g., convert the data to nifti, break rois to hemispheres etc.). I put them here just for documentation, but as the data that I uploaded to osf is already prepared for analysis, these scripts are unnecessary. They are marked with:  "Already done, unecessary now."
* I used participants folders with initials as folder names. In the osf data base, I re-named the folders to be s01-s020. But the scripts are with the subjects initials. So if one downloads the data and want to use these scripts to analyse these data - change the subject names.
* in the ms, subject 20 was removed from analysis due to entorhinal dropout. I uploaded their data, because no reason not to use if analyzing other ROIs.
* numbers refer to the order that things should be done, there is some flexibilty though...

## 1. PreprocessingAndModel folder:

step0_TickyReanalysis_globals: definition file, to be run before starting the analysis/called by other scripts

step1_Afni2NII_multiple_subs: I converted the files to nifti from KD's Afni format. Already done, unecessary now.
step2_DivideMTL_RightLeft: divide MTL regions to right and left hemisphere. Already done, unecessary now.
step3_remake_behavioral_files: for some participants, I needed to re-create their behavior files to be in the standard format because they were in a different format in KD's folder. Already done, unecessary now.

step4_create_regressors_files: Creates the regressors for the LSS and the univariate GML models. In the osf repository I uploaded the regressors.

step5_seg4Dto3D_multiple_subs: I think I had the data all runs concatinated. This breaks the data to runs.  Already done, unecessary now.

step6_model-similarity: runs the LSS models to obtain the single-trial t-stat. Each trial had a cue part and an image part. Creates one t-stat for the cue and one for the image per each trial.

entorhinal_dropout: the script I used to calculate entorhinal dropout - calculate how many voxels from entorhinal anatomical mask (registered to epi) are there in the epi data - which is masked for signal dropout

wait_for_feat: a script that waits for feat to end running. The modeling code and other codes use it.

### 1.1. UnivariateGLM

This analysis was done after I did the single-trial models, to do the voxel selection. I ran it because I wanted a model that is as independent as possible from the single-trial t-stats I used for RSA/connnectivity. Note that the univariate effects reported in the paper are from averaging single-trials (see below). The reasoning was that the main aim was to control for univaraite effects in the RSA/connecitvity analysis. So I wanted to use univariate activation from the same data. If I remember correctly, the pattern was very similar in the two approaches.

Some steps in the analysis that were done before and are identical and should not be done again (e.g., converting from afni to nifti files) are not in this folder. for these steps, see the mother dir of Preprocessing and model


#contrasts 10.23.17:
1-9: the single conditions
1: 0-changes
9: 4-changes
10: 1-changes
11: 2-changes
12: 3-changes

#in the subject level analysis, for each cope, cope1 is layout, cope2 is item, cope3 averages both

#AK and LD subject level analysis ran through the gui - different number of scans

scripts:
step1-create regressors_files: termed step1 here bc I ran it after I did the single trials
step2_FirstLevel_model-AllTrialsUnivariate: runs the first level (run level) univariate analysis
step3_SubjectLevel_model-AllTrialsUnivariateItemFirst: run the subject level feat analysis, 
step3_SubjectLevel_model-AllTrialsUnivariateLayoutFirst: run the subject level feat analysis
*these step3 scripts also prepare the files for dummy registration - see notes in the scripts*
*since all analyses were ROI, never ran the group level*

## 1.2 fsl templates folder:
The way the models run is that they take an fsf template file, and change it to be specific to each participant, saves it for each participant, and run feat calling the subject-specific files. These are the templates, and their names are self-explanatory.

## 2. GetSingleTrialsData:
The file in this folder grabs the single trial t-stats from the nifti files per roi and make a matlab structure with the data, per trial and per voxel. It creates a reg_data**.mat file that is placed in each participant's data folder.
The get_CueIm_data is for both the hippocampus and main MTL rois.
* it uses spm_read_vol - I'd now replace that to niftiread function - but pay attention to x/y/z coordinates becasue different functions read the nifti files differently. As long as the command for the ROIs is the same as for the uploading of the data itself - or that the function read nifti files the same way, you're good.

## 3. connectivity

these files grab the reg_data**.mat file, calculate connectivity, and creates a matlab structure with the group results of all rois. Per pairs of rois, it has the connectivity value for each participant, in each task and number of changes.

connectivity_allRegs_SeparateTaskNumChanges: that's the main script, it calculates the connectivity btw regions as above, that's what reported in the paper.

## 4. RSA

There are two RSA analyses in the paper, prediction strength, and predicion error. They have two different folders.

## 4.1. prediction strength

This analysis generally computes similarity between a cue and the match (0-changes images) of that cue, vs. the match image of other cues.
CueMatchImagePredictionRSACompOtherRooms: that's the main script that computes the prediction strength. outputs a matlab structure with the prediction strength per roi, per participant, per task. Also creates an R data file, to be used later on for stats. This analysis is reported for CA1 in the supplementary.

CueMatchImagePredictionRSACompOtherRooms_selectVoxels: same as above, but with voxel selection - in the paper we adopted voxel selection for RSA. this is reported in the main text in prediction strength.

## 4.2. prediction error

This analysis computes the similarity between the cue and the image.
CueImageVsOtherRoomsRSA_SeprateTaskNumChanges: that's the main script that computes the analysis. In the paper we reported RSA with voxel selection, the results of this script, i.e., without voxel selection, are reported in the supplementary.

CueImageVsOtherRoomsRSA_SeprateTaskNumChanges_selectVoxels: same as above, with voxel selection - reported in the main text.

## 5. voxel selection (Univariate proper GLM)
In the manusctipt we did voxel selection based on activity level in voxels. To do that, I used the average activation from a GLM modelling average task activation across trials (for the different bins), than averaged across the 5 levels of changes. I call it "proper GLM" bc it's the traditional way of doing univariate analysis based on task, while for this experiment, the reported univariate was done by averaging across single-trials t-stats (see below, this was bc the main aim was controlling for our other effects).

get_Im_data: gets the data per participant from the nifti files to a matlab stracture

ActivityNumChanges_SeparateTask_allregs: This script also outputs a structure "ActiveVoxels" that has the active voxels within each participant that will be taken for the voxel selection for RSA (high 2thirds of the voxels). 

## 6. univariate analysis (based on average single trials)

Note that the univariate effects reported in the paper are from averaging single-trials. The reasoning was that the main aim was to control for univaraite effects in the RSA/connecitvity analysis. These are the scripts that do this.

In the main text, we report the control for univariate ananlysis for the connectivity finding. In the supplementary, we report the values of the univariate analysis, as well as controling for RSA results.

AverageActivityNumChanges_SeparateTask_cue_image: takes the single trial data and computes the average univariate activation. Note that the single-trial data should have been pulled already from the nifti files (see files in the GetSingleTrialsData folder). This was used to control for univariate activation in the connectivity results (main text), and these are data reported in the Supplementary Figure.

AverageActivity_cue_image_SelVoxProperGLM: in the main text, we reported RSA using voxel selection. Then, we need to do voxel selection to control for univariate activation (these control analyses are reported in the supplementary) - this computes the univariate activation with voxel selection.

## 7. all group level analyses  

Mostly, the analyses scripts above created a data structure with all subjects and all bins/conditions etc - ready to run group level stats and plots. This folder does that for all analyses. I like having all the group level stuff in one folder, across various analyses, because it's easy for me to look at things in parallel, analyse across analyses.

*some scripts have the option of looking by accuracy - we did that in response to reviewer, we didn't find any differences based on accuracy, so the results in the paper are based on all items, regardless of accuracy in the task
*I started the group analyses in matlab, then shifted to R, that's why the scripts are a mix of the two, and there's some redundency.

*In some group analysis plots I calculate SEM based on some within-participant measures. While this is useful for looking at the data, for reporting, SEM should be done the standard way. Since most plots were taken from R eventually anyway, I didn't bother changing it.

*The prediction strength and prediction error were reported in the main text with voxel-selection. We also report without - note that the same group-level matlab scripts provide stats for both - just with uploading a different data structure.

*for convenience, I also placed in this folder txt/Excel files that are needed for the group analyses in R.

TickyDataR.xlsx: includes spreadsheets with data that is later uploaded in the GroupImaging.R scripts. Whenever needed, I explain below how the data for this file was created.

##### behavior:
Behavior_acc_RTs_SeparateTaskNumChanges: this runs the behavioral analysis and creates a group structure.

plot_Behavior_SeparateTask_NumChanges: plots some things and does stats. I shifted to R along the way, but kept this script. It is used to produce the data structure in an easy way for data analysis in R. *note that plots and stats were eventually done and reported using R - see script below.

Ticky_behavior.Rmd: basic stats and plotting of the behavior.

##### connectivity:

plot_region_connectivity_SeparateTask_NumChanges: this script does stats and plots, I shifted to R along the way, but kept this script. it's pretty much unnecessary now. 
makeRdataStructureHemiRegANOVA: makes the data file for R to run the three-way region by changes by hemisphere ANOVA
makeRdataStructureConUnivar: makes the data file for each pair of regions for connectivity, including univariate activation and RT/accuracy models. Control for univariate is reported, RT/accuracy models was done in response to reviewers. eventually we took the actual RT/accuracy rates per participant, but this works the same.

##### prediction strength:
Con_CueMatchImVsOtherRoomsPredictabilityRSA_AcrossSubCorrel: this computes the prediction strength, as well as the correlation with connectivity. 

CA1_prediction_CA1_ent_con.xlsx: control analyses for the prediction strength and correlation were done with R (GroupImaging.R, see below). The values for the R analyses are loaded from this excel file.
Here is where I took the values from:
The activation values for control for perdiction strength: the script that runs the prediction strength create these values (CueMatchImagePredictionRSACompOtherRooms_selectVoxels, see above). I took the data from AvPerTaskNo0Changes, which means that each task is averaged first, then in the excel file I average across tasks when needed. This is for the cues and the match images. if you check intact_image_act.all_items you'll see 6 columns, they match the 6 columns in the excel file that correspond to the images (columns L-Q).
The activation values for control for the connectivity I copy pasted from the data after script that computes the univariate activation based on single-trials (AverageActivityNumChanges_SeparateTask_cue_image, see above).

##### prediction error:
plot_region_CueImageOtherImRSA_SeparateTask_NumChanges: This script runs the group-level analysis and plot stuff for the PE analysis. results and stats are based on that script, and the bars in Figure 3 in the ms is from this file.
*this script also has the option of looking at similarity subtracting other images, a control that is reported in Supplementary Note 5.

The data for the control analysis for univariate (see GroupImaging.R) was done by copy-paste to the relevant excel sheet, from running the script above and the relevant univariate analysis (single-trials).

##### univariate
plot_region_activity_SeparatedTask_NumberChanges: does some plots (not that supp figure is from the Groupimaging script), and run stats - stats reported in the paper are based on this file 

##### GroupImaging.R script:
Includes a lot of the control analyses that I've done in the paper, as well as the code that plots most of the graphs. See the script itself for the specific plots/stats.

