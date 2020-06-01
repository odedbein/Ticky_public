function [ResultsCueImSim, ResultsCueImSimItemsCount, ResultsCueImSimOnlyNum, ResultsCueImSimAllTrials]=CueMatchImagePredictionRSACompOtherRooms_selectVoxels(ActiveVoxels,engram,remove_outliers)
%this script computes the "prediction" similarity: similarity between the
%cue and the 0-changes image. I.e., the extent to which participants
%reinstated the correct image during the cue. it compares that to the cue
%vs. all other scenes. Thus, cannot do it for the 0-changes condition -
%because that is comparing withn/across trials


% the reported analysis removes same block trials (RemoveSameBlock=1), obv. FishTrans=1
% 

% to select voxels, you need the activeVoxels structure:
%load('/Volumes/data/Bein/TickyReanalysis/results/Univar_properGLM_withVoxelSelection.mat')

% in the paper I didn't remove outliers, tried that for the re-submission -
% didn't change much I think.
warning('off','all')
FishTrans=1;
RemoveSameBlock=1; %i built it in the analysis that I only compute the control across block, so take the main analysis across block only as well.

if RemoveSameBlock
    display('analysis excludes cue of the same block as the 0-changes image, aside from in the 0-changes condition');
else
    display('analysis DOES NOT exclude cue of the same block as the 0-changes image, aside from in the 0-changes condition');
end

if remove_outliers
    disp('remove outliers is based on Remove same block and only same-task cue/match image trials.');
    disp('This means that a lot same-block corr will be removed as if they are outleirs.');
    disp('If want to look at analysis that does not remove same block for some reason - don''t remove outliers.');
    
    out_ttl='_outliers_removed';
else
    display('analysis DOES NOT remove outleirs');
    out_ttl='';
end

med_ttl='HighThird';
glm_ttl='ProperGLM';
mm_ttl='Av5cond';

if engram
    mydir='/data/Bein';
else
    mydir='/Volumes/data/Bein';
end
proj_dir=fullfile(mydir,'TickyReanalysis');
rmpath('/Volumes/data/Bein/fMRI_course/AnalysisScripts');
results_fname=sprintf('CueMatchImPredCompOtherRoomsSelVox_%s_%s_%s%s.mat',glm_ttl,mm_ttl, med_ttl,out_ttl);
results_fname=fullfile(proj_dir,'results','rsa_prediction',results_fname);

subjects={'AB'; 'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD';'JG'; 'JM'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};

subj_dir=fullfile(proj_dir,'SubData');
trial_types={...
    'R0I0';%1
    'R0I1';%2
    'R0I2';%3
    'R1I0';%4
    'R1I1';%5
    'R1I2';%6
    'R2I0';%7
    'R2I1';%8
    'R2I2'...%9
    };
num_comparisons=4*2;%1-4 changes * lay/item

reg_names={'lCA1'};

ResultsCueImSim={};
ResultsCueImSimOnlyNum={};
ResultsCueImSimItemsCount={};
ResultsCueImSimItemsCount.CueNoChange={};

ResultsCueImSimAllTrials={};
%for each reegion:
%1. all items
%2. all items without the no-change
%3. all items without same block ones

%prepare the header:
for reg=1:numel(reg_names)
    ResultsCueImSim.(reg_names{reg}).all_items{1,1}='subjects';
    ResultsCueImSim.(reg_names{reg}).all_items(1,(2:num_comparisons+1))={'lay: 1-change','lay: 2-changes','lay: 3-changes','lay: 4-changes','item: 1-change','item: 2-changes','item: 3-changes','item: 4-changes'};
    ResultsCueImSim.(reg_names{reg}).Acc=ResultsCueImSim.(reg_names{reg}).all_items;
    ResultsCueImSim.(reg_names{reg}).NonAcc=ResultsCueImSim.(reg_names{reg}).all_items;
end

for subj=1:numel(subjects)
    if strcmp(char(subjects(subj)),'LD')
        num_trials=216;
    elseif strcmp(char(subjects(subj)),'AK')
        num_trials=189;
    else
        num_trials=270;
    end
    fprintf('analyzing subj %s\n',char(subjects(subj)));
    %load the data:
    load(fullfile(subj_dir,char(subjects(subj)),'data','reg_data.mat'),'reg_data');
    
    
    %load the behavioral data, to know which items go to which condition:
    behav_filename=fullfile(subj_dir,char(subjects(subj)),sprintf('output_subject_%s.txt',subjects{subj}));
    if strcmp(char(subjects(subj)),'AR')
        subj_behavior=textscan(fopen(behav_filename),'%d %d %d %.1f %s %s %d %s %.1f %s');
    else
        subj_behavior=textscan(fopen(behav_filename),'%d %d %d %.1f %s %d %s %.1f %s %s %.3f %.3f %d');
    end
    
    %gett accuracy, trial type, stimuli
    acc=subj_behavior{3}(1:num_trials);
    type=subj_behavior{5}(1:num_trials);
    if strcmp(char(subjects(subj)),'AR')
        stim=subj_behavior{6}(1:num_trials);
    else
        stim=subj_behavior{10}(1:num_trials);
    end
    block=subj_behavior{1}(1:num_trials);
    
    if strcmp(char(subjects(subj)),'AR')
        task=subj_behavior{8}(1:num_trials);
    else
        task=subj_behavior{7}(1:num_trials);
    end
    
    if strcmp(char(subjects(subj)),'LD') %need to cut the first two sessions - brain data only from session 3.
        acc=subj_behavior{3}(((270/10*2)+1):end);
        type=subj_behavior{5}(((270/10*2)+1):end);
        task=subj_behavior{7}(((270/10*2)+1):end);
        stim=subj_behavior{10}(((270/10*2)+1):end);
        block=subj_behavior{1}(((270/10*2)+1):end);
    end
    if isempty(acc)
        fprintf('shoot.. no accurcay column, check behavioral file\n');
    end
    
    %edit the stim to extract just the room's name
    for s=1:num_trials
        curr_room=stim{s};
        stim{s}=curr_room(7:12);
    end
    
    find_types=zeros(num_trials,numel(trial_types));
    for cond=1:numel(trial_types)
        for i=1:num_trials
            find_types(i,cond)=strcmp(type{i},trial_types(cond));
        end
    end
    
    %binarise the task vector, easier for later: 1 - lay task, 0 - item task
    task_binary=zeros(num_trials,1);
    for i=1:num_trials
        task_binary(i)=strcmp(task{i},'Lay');
    end
    %check for problems:
    if ~any(task_binary)
        fprintf('shoot.. all task is zero - check behavioral file\n');
    end
    
    conds=zeros(num_trials,num_comparisons);
    conds_temp=[];
    %number of changes, regardless of accuracy:
    conds_temp(:,1)=find_types(:,1);%no changes
    conds_temp(:,2)=find_types(:,2)+find_types(:,4);%1 change
    conds_temp(:,3)=find_types(:,3)+find_types(:,5)+find_types(:,7);%2 changes
    conds_temp(:,4)=find_types(:,6)+find_types(:,8);%3 changes
    conds_temp(:,5)=find_types(:,9);
    %layout task
    %conds(:,1)= (conds_temp(:,1) & task_binary==1); %no changes
    conds(:,1)= (conds_temp(:,2) & task_binary==1);%1 change
    conds(:,2)= (conds_temp(:,3) & task_binary==1);%2 changes
    conds(:,3)= (conds_temp(:,4) & task_binary==1);%3 changes
    conds(:,4)= (conds_temp(:,5) & task_binary==1);%4 changes
    
    %item task
    %conds(:,6)= (conds_temp(:,1) & task_binary==0); %no changes
    conds(:,5)= (conds_temp(:,2) & task_binary==0);%1 change
    conds(:,6)= (conds_temp(:,3) & task_binary==0);%2 changes
    conds(:,7)= (conds_temp(:,4) & task_binary==0);%3 changes
    conds(:,8)= (conds_temp(:,5) & task_binary==0);%4 changes
    
    %build the info matrix:
    %col 1: trial number of the relevant no-change images (col 1),
    %col 2: task of the no-changes image (lay-1,item-0)
    %col 3: are the current item and the image of the same task? (1-yes, 0-no)
    %col 4: block # of the no-changes image
    %col 5: are the current item and the image of the same block? (1-yes, 0-no)
    %col 6: is the current item a no-change item
    noChngeMat=zeros(num_trials,6);
    for i=1:num_trials
        cue=stim{i}(1:2);
        %now search for the 0-changes condition
        noChange=[cue 'R0I0'];
        for nc=1:num_trials
            if strcmp(stim{nc},noChange)
                noChngeMat(i,1)=nc; %trial of the 0-changes
                noChngeMat(i,2)=task_binary(nc); %was the 0-changes in the layout/task
                noChngeMat(i,3)=double(task_binary(i)==task_binary(nc)); %was the 0-changes in the same task as the current one
                noChngeMat(i,4)=block(nc);
                noChngeMat(i,5)=double(block(i)==block(nc));%this trial is in the same block as the 0-changes trial
                noChngeMat(i,6)=double(i==nc);%this trial is the 0-changes trial
                break
            end
        end
    end
    all_0changes=unique(noChngeMat(:,1));
    %remove in case that there are zeros:
    all_0changes=all_0changes(all_0changes~=0);
    ResultsCueImSimItemsCount.all_items{subj,1}=subjects{subj};
    ResultsCueImSimItemsCount.all_items(subj,(2:num_comparisons+1))=num2cell(sum(conds));
    
    %prepare the matrix with only accuracte responses
    conds_accurate=conds;
    conds_accurate(~acc,:)=0;
    %save the number of items per each subject
    ResultsCueImSimItemsCount.Acc{subj,1}=subjects{subj};
    ResultsCueImSimItemsCount.Acc(subj,(2:num_comparisons+1))=num2cell(sum(conds_accurate));
    
    %prepare the matrix with only NON accuracte responses
    conds_non_accurate=conds;
    conds_non_accurate(acc==1,:)=0;
    %save the number of items per each subject
    ResultsCueImSimItemsCount.NonAcc{subj,1}=subjects{subj};
    ResultsCueImSimItemsCount.NonAcc(subj,(2:num_comparisons+1))=num2cell(sum(conds_non_accurate));
    
    %now get the data for each region and calculate the correlations:
    for reg=1:numel(reg_names)
        
        if strcmp(char(subjects(subj)),'EB') && (strcmp(reg_names{reg},'lacc') || strcmp(reg_names{reg},'acc') || strcmp(reg_names{reg},'racc'))
            %display('skipping EB lacc');
            %skip this region - subj do not have a lacc region
        else
            
            %put the subj_name in the reuslts structure:
            ResultsCueImSim.(reg_names{reg}).all_items(subj+1,1)=subjects(subj);
            ResultsCueImSim.(reg_names{reg}).Acc(subj+1,1)=subjects(subj);
            ResultsCueImSim.(reg_names{reg}).NonAcc(subj+1,1)=subjects(subj);
            
            cue_data=reg_data.(reg_names{reg}).cue;
            if size(cue_data,1) > 10 %only do if more than 10 voxels - otherwise there's no selection of voxels... and the analysis is meaningless
                
                %select voxels:
                curr_voxels=ActiveVoxels.((reg_names{reg})).(mm_ttl).(med_ttl).(subjects{subj});
                cue_data=cue_data(curr_voxels,:);
                %check the data for no nans:
                if ~isempty(find(isnan(cue_data)))
                    fprintf('cue data has nans\n')
                end
                image_data=reg_data.(reg_names{reg}).Image;
                image_data=image_data(curr_voxels,:);
                %check the data for no nans:
                if ~isempty(find(isnan(image_data)))
                    fprintf('image data has nans\n')
                end
                
                if size(cue_data,1) > 10 %only do the analysis for more than 10 voxels:
                    %% clean the data:
                    %compute the mean and SD to be able to exclude outliers
                    %since everything is computed across runs - exclude based on
                    %the entire sample of images.
                    %compute the correlations between all cues and their 0-changes
                    %images- get a measure of trial-by-trial prediction-all items
                    if remove_outliers
                        all_correl=[];
                        for i=1:num_trials
                            mi=noChngeMat(i,1); %this is the trial of the match image
                            if mi~=0 %a match trial exists
                                if (task_binary(i)==task_binary(mi)) && (block(i)~=block(mi)) %match image is in the same task,
                                    %and not in the same block: include this correlation, and the others
                                    %add the correlation to the match image:
                                    all_correl=[all_correl;corr(cue_data(:,i),image_data(:,noChngeMat(i,1)))];
                                    curr_0changes=[];
                                    for checkTask=1:length(all_0changes)
                                        ct=all_0changes(checkTask);
                                        if (task_binary(i)==task_binary(ct)) && (block(i)~=block(ct)) %this cue and the 0-changes image are of the same task
                                            curr_0changes=[curr_0changes; ct];
                                        end
                                    end
                                    curr_0changes(curr_0changes==noChngeMat(i,1))=[];
                                    all_correl=[all_correl; corr(cue_data(:,i),image_data(:,curr_0changes))']; %that is the average of similarity to all the other scenes.
                                end
                            end
                        end
                        %compute mean and SD :
                        mean_correl=nanmean(all_correl);
                        sd_correl=nanstd(all_correl);
                    end
                    
                    %% now for real compute, and remove outliers
                    %compute the correlations between all cues and their 0-changes
                    %images- get a measure of trial-by-trial prediction-all items
                    
                    %also include: univariate activation per cue
                    %match images activation will be later
                    reg_correl=nan(num_trials,2);
                    cue_act=nan(num_trials,1);
                    image_act=nan(num_trials,2);
                    for i=1:num_trials
                        if noChngeMat(i,1)==0
                            reg_correl(i,:)=nan;
                        else
                            reg_correl(i,1)=corr(cue_data(:,i),image_data(:,noChngeMat(i,1)));
                            curr_0changes=[];
                            for checkTask=1:length(all_0changes)
                                ct=all_0changes(checkTask);
                                if (task_binary(i)==task_binary(ct)) && (block(i)~=block(ct)) %this cue and the 0-changes image are of the same task
                                    curr_0changes=[curr_0changes; ct];
                                end
                            end
                            curr_0changes(curr_0changes==noChngeMat(i,1))=[];
                            reg_correl(i,2)=nanmean(corr(cue_data(:,i),image_data(:,curr_0changes))); %that is the average of similarity to all the other scenes.
                            cue_act(i)=nanmean(cue_data(:,i));
                            image_act(i,1)=nanmean(image_data(:,noChngeMat(i,1)));
                            image_act(i,2)=nanmean(nanmean(image_data(:,curr_0changes))); %mean across images and across voxels
                        end
                    end
                    
                    %same for accuracte responses only:
                    %4/14/19 - the way I did it before was effectively taking all of the
                    %other scenes, even for accurate trials. The idea was that accuracy for the 0-changes image is not important
                    %since this is the template. This is still done (AccCueOnly)
                    %what I decided to check was to also take only 0-changes scenes
                    %in which participants were accurate, both for the "same" and
                    %for "others":
                    %note that this cannot be done for inaccurate responses, since
                    %participants were pretty much accurate in all of them, so the
                    %non acc is taking all of them..
                    reg_correl_acc=nan(num_trials,2);
                    cue_act_acc=nan(num_trials,1);
                    image_act_acc=nan(num_trials,2);
                    for i=1:num_trials
                        if (acc(i)==1) && (noChngeMat(i,1)~=0) && (acc(noChngeMat(i,1))==1) %participant was correct in both the cue and the 0-changes
                            reg_correl_acc(i,1)=corr(cue_data(:,i),image_data(:,noChngeMat(i,1)));
                            curr_0changes=[];
                            for checkTask=1:length(all_0changes)
                                ct=all_0changes(checkTask);
                                if (task_binary(i)==task_binary(ct)) && (block(i)~=block(ct)) && (acc(ct)==1)
                                    %this cue and the 0-changes image are of the same task and accuracy for the other one was 1
                                    %so add the index to the list of curr_0changes:
                                    curr_0changes=[curr_0changes; ct];
                                end
                            end
                            %take out the one that is the same one:
                            curr_0changes(curr_0changes==noChngeMat(i,1))=[];
                            %average correlation
                            reg_correl_acc(i,2)=nanmean(corr(cue_data(:,i),image_data(:,curr_0changes))); %that is the average of similarity to all the other scenes.
                            cue_act_acc(i)=nanmean(cue_data(:,i));
                            image_act_acc(i,1)=nanmean(image_data(:,noChngeMat(i,1)));
                            image_act_acc(i,2)=nanmean(nanmean(image_data(:,curr_0changes))); %mean across images and across voxels
                        end
                    end
                    
                    if remove_outliers
                        reg_correl(reg_correl > (mean_correl + 3*sd_correl))=nan;
                        reg_correl(reg_correl < (mean_correl - 3*sd_correl))=nan;
                        
                        reg_correl_acc(reg_correl_acc > (mean_correl + 3*sd_correl))=nan;
                        reg_correl_acc(reg_correl_acc < (mean_correl - 3*sd_correl))=nan;
                    end
                    
                    %add the difference score:
                    reg_correl=[reg_correl reg_correl(:,1)-reg_correl(:,2)];
                    reg_correl_acc=[reg_correl_acc reg_correl_acc(:,1)-reg_correl_acc(:,2)];
                    image_act=[image_act image_act(:,1)-image_act(:,2)];
                    image_act_acc=[image_act_acc image_act_acc(:,1)-image_act_acc(:,2)];
                    
                    % SEE NOTE ABOVE ABOUT OUTLIER REMOVAL!!!!
                    %because some subjects have different number of trials, need to put
                    %in a different structrue per participant:
                    ResultsCueImSimAllTrials.All.(reg_names{reg}).all_items.(subjects{subj})=reg_correl;
                    ResultsCueImSimAllTrials.All.(reg_names{reg}).AccCueAnd0Changes.(subjects{subj})=reg_correl_acc(acc==1,:);
                    ResultsCueImSimAllTrials.All.(reg_names{reg}).AccCueOnly.(subjects{subj})=reg_correl(acc==1,:);
                    ResultsCueImSimAllTrials.All.(reg_names{reg}).NonAcc.(subjects{subj})=reg_correl(acc==0,:);
                    %now remove the no-change condition:
                    curr_reg_correl=reg_correl(noChngeMat(:,6)==0,:);
                    curr_acc=acc(noChngeMat(:,6)==0);
                    ResultsCueImSimAllTrials.NoChangeRemoved.(reg_names{reg}).all_items.(subjects{subj})=curr_reg_correl;
                    ResultsCueImSimAllTrials.NoChangeRemoved.(reg_names{reg}).AccCueOnly.(subjects{subj})=curr_reg_correl(curr_acc==1,:);
                    ResultsCueImSimAllTrials.NoChangeRemoved.(reg_names{reg}).NonAcc.(subjects{subj})=curr_reg_correl(curr_acc==0,:);
                    curr_reg_correl=reg_correl_acc(noChngeMat(:,6)==0,:);
                    ResultsCueImSimAllTrials.NoChangeRemoved.(reg_names{reg}).AccCueAnd0Changes.(subjects{subj})=curr_reg_correl(curr_acc==1,:);
                    %now without trials in which the cue and the no-change are in the same block (becasue presumably cannot similarity calculate within the same block - note that this removes also the 0-changes condition):
                    curr_reg_correl=reg_correl(noChngeMat(:,5)==0,:);
                    curr_acc=acc(noChngeMat(:,5)==0);
                    ResultsCueImSimAllTrials.SameBlockRemoved.(reg_names{reg}).all_items.(subjects{subj})=curr_reg_correl;
                    ResultsCueImSimAllTrials.SameBlockRemoved.(reg_names{reg}).AccCueOnly.(subjects{subj})=curr_reg_correl(curr_acc==1,:);
                    ResultsCueImSimAllTrials.SameBlockRemoved.(reg_names{reg}).NonAcc.(subjects{subj})=curr_reg_correl(curr_acc==0,:);
                    curr_reg_correl=reg_correl_acc(noChngeMat(:,5)==0,:);
                    ResultsCueImSimAllTrials.SameBlockRemoved.(reg_names{reg}).AccCueAnd0Changes.(subjects{subj})=curr_reg_correl(curr_acc==1,:);
                    
                    %for single trials R analysis:
                    %1. keep track on which task
                    %2. remove same block
                    %3. Take only same task
                    %4. fisher transform
                    curr_trials=(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==2); %take only same task, and remove same block
                    %if we only have 0/1 that means that the intercept would be the average diff in the
                    %0 task. So need to have 2 variables, one for each task
                    curr_tasks_lay=task_binary(curr_trials);
                    curr_tasks_item=~task_binary(curr_trials);
                    curr_block=double(block(curr_trials));
                    if FishTrans
                        curr_reg_correl=atanh(reg_correl(curr_trials,1:3)); %don't take the univar
                    else
                        curr_reg_correl=reg_correl(curr_trials,1:3);
                    end
                    %put it all together:
                    ResultsCueImSimAllTrials.TrialsForR.(reg_names{reg}).all_items.(subjects{subj})=[curr_reg_correl curr_tasks_lay curr_tasks_item curr_block];
                    
                    
                    %% now compute and average for each condition:
                    
                    
                    %this analysis is dones seprately for each task. Here, can devide
                    %based on whether the task of the cue, the no-changes image, or
                    %both. I chose to do it by a match - if both are of the lay/item
                    %task - I did the other scenes for only matched in task, so if
                    %I change how I calculate the main analysis, I need to change
                    %how I do the control.
                    
                    %%average all items together, collapse on conditions - take out
                    %%the 0-changes condition, because this is not a fair
                    %%comparison - it's within trial, and all controls are across
                    %%trials.
                    if RemoveSameBlock
                        if FishTrans
                            ResultsCueImSimOnlyNum.AvAllTrialsNo0Changes.(reg_names{reg}).all_items(subj,:)=atanh(nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==2,:)));
                            ResultsCueImSimOnlyNum.AvAllTrialsNo0Changes.(reg_names{reg}).AccCueOnly(subj,:)=atanh(nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (acc==1))==3,:)));
                            ResultsCueImSimOnlyNum.AvAllTrialsNo0Changes.(reg_names{reg}).AccCueAnd0Changes(subj,:)=atanh(nanmean(reg_correl_acc(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (acc==1))==3,:)));
                            ResultsCueImSimOnlyNum.AvAllTrialsNo0Changes.(reg_names{reg}).NonAcc(subj,:)=atanh(nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (acc==0))==3,:)));
                            
                        else
                            ResultsCueImSimOnlyNum.AvAllTrialsNo0Changes.(reg_names{reg}).all_items(subj,:)=nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==2,:));
                            ResultsCueImSimOnlyNum.AvAllTrialsNo0Changes.(reg_names{reg}).AccCueOnly(subj,:)=nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (acc==1))==3,:));
                            ResultsCueImSimOnlyNum.AvAllTrialsNo0Changes.(reg_names{reg}).AccCueAnd0Changes(subj,:)=nanmean(reg_correl_acc(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (acc==1))==3,:));
                            ResultsCueImSimOnlyNum.AvAllTrialsNo0Changes.(reg_names{reg}).NonAcc(subj,:)=nanmean(reg_correl_acc(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (acc==0))==3,:));
                            
                        end
                    else %didn't remove same block - so need to remove 0-changes condition  -currently doesn't have the stringent acc
                        if FishTrans
                            ResultsCueImSimOnlyNum.AvAllTrialsNo0Changes.(reg_names{reg}).all_items(subj,:)=atanh(nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,6)==0))==2,:)));
                            ResultsCueImSimOnlyNum.AvAllTrialsNo0Changes.(reg_names{reg}).AccCueOnly(subj,:)=atanh(nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,6)==0) + (acc==1))==3,:)));
                        else
                            ResultsCueImSimOnlyNum.AvAllTrialsNo0Changes.(reg_names{reg}).all_items(subj,:)=nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,6)==0))==2,:));
                            ResultsCueImSimOnlyNum.AvAllTrialsNo0Changes.(reg_names{reg}).AccCueOnly(subj,:)=nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,6)==0) + (acc==1))==3,:));
                        end
                    end
                    
                    %average based on task:
                    col=0;
                    for t=1:-1:0
                        if RemoveSameBlock
                            if FishTrans
                                ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).all_items(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (task_binary==t))==3,:)));
                                ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).AccCueOnly(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (acc==1) + (task_binary==t))==4,:)));
                                ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).AccCueAnd0Changes(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl_acc(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (acc==1) + (task_binary==t))==4,:)));
                                ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).nonAcc(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (acc==0) + (task_binary==t))==4,:)));
                                
                            else
                                ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).all_items(subj,(col*3)+1:(col*3)+3)=nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (task_binary==t))==3,:));
                                ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).AccCueOnly(subj,(col*3)+1:(col*3)+3)=nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (acc==1) + (task_binary==t))==4,:));
                                ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).AccCueAnd0Changes(subj,(col*3)+1:(col*3)+3)=nanmean(reg_correl_acc(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (acc==1) + (task_binary==t))==4,:));
                                ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).NonAcc(subj,(col*3)+1:(col*3)+3)=nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (acc==0) + (task_binary==t))==4,:));
                                
                            end
                            %cue and image activation:
                            ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).cue_act.all_items(subj,(col+1))=nanmean(cue_act(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (task_binary==t))==3));
                            ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).intact_image_act.all_items(subj,(col*3)+1:(col*3)+3)=nanmean(image_act(((noChngeMat(:,3)==1) + (noChngeMat(:,5)==0) + (task_binary==t))==3,:));
                        else %didn't remove same block - so need to remove 0-changes condition
                            if FishTrans
                                ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).all_items(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,6)==0) + (task_binary==t))==3,:)));
                                ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).AccCueOnly(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,6)==0) + (acc==1) + (task_binary==t))==4,:)));
                            else
                                ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).all_items(subj,(col*3)+1:(col*3)+3)=nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,6)==0) + (task_binary==t))==3,:));
                                ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_names{reg}).AccCueOnly(subj,(col*3)+1:(col*3)+3)=nanmean(reg_correl(((noChngeMat(:,3)==1) + (noChngeMat(:,6)==0) + (acc==1) + (task_binary==t))==4,:));
                            end
                        end
                        col=col+1;
                    end
                    
                    %%average based on conditions:
                    for c=1:size(conds,2) %for each condition
                        col=c-1;
                        if RemoveSameBlock
                            if reg==1 %collect the number of items
                                %note that this may be inaccurate for subjects like AK
                                %that have missing values because the 0-changes is
                                %in the blocks that he misses - these are the "0"
                                %in the first column of noChngeMat
                                % however,as long as I have noChngeMat(:,3)==1 as a
                                % condition - it's okay, becasue for every missing
                                % trial, noChngeMat(:,3)==0, so these trials will
                                % be removed
                                ResultsCueImSimItemsCount.CueNoChange.all_items(subj,c)=sum(((conds(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3);
                                ResultsCueImSimItemsCount.CueNoChange.AccCueOnly(subj,c)=sum(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3);
                                ResultsCueImSimItemsCount.CueNoChange.NonAcc(subj,c)=sum(((conds_non_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3);
                                %ResultsCueImSimItemsCount.CueNoChange.AccCueOnly(subj,c)=sum(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3);
                                
                            end
                            if FishTrans
                                %fisher transform:
                                ResultsCueImSimOnlyNum.(reg_names{reg}).all_items(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl(((conds(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                ResultsCueImSimOnlyNum.(reg_names{reg}).AccCueOnly(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                ResultsCueImSimOnlyNum.(reg_names{reg}).AccCueAnd0Changes(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl_acc(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                ResultsCueImSimOnlyNum.(reg_names{reg}).NonAcc(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl(((conds_non_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                
                                %ResultsCueImSim.(reg_names{reg}).all_items{subj+1,1+((col*3)+1:(col*3)+3)}=atanh(nanmean(reg_correl(((conds(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                %ResultsCueImSim.(reg_names{reg}).AccCueOnly{subj+1,1+((col*3)+1:(col*3)+3)}=atanh(nanmean(reg_correl(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                %ResultsCueImSim.(reg_names{reg}).AccCueAnd0Changes{subj+1,1+((col*3)+1:(col*3)+3)}=atanh(nanmean(reg_correl_acc(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                %ResultsCueImSim.(reg_names{reg}).NonAcc{subj+1,1+((col*3)+1:(col*3)+3)}=atanh(nanmean(reg_correl(((conds_non_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                            else
                                ResultsCueImSimOnlyNum.(reg_names{reg}).all_items(subj,(col*3)+1:(col*3)+3)=(nanmean(reg_correl(((conds(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                ResultsCueImSimOnlyNum.(reg_names{reg}).AccCueOnly(subj,(col*3)+1:(col*3)+3)=(nanmean(reg_correl(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                ResultsCueImSimOnlyNum.(reg_names{reg}).AccCueAnd0Changes(subj,(col*3)+1:(col*3)+3)=(nanmean(reg_correl_acc(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                ResultsCueImSimOnlyNum.(reg_names{reg}).NonAcc(subj,(col*3)+1:(col*3)+3)=(nanmean(reg_correl(((conds_non_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                
                                %ResultsCueImSim.(reg_names{reg}).all_items{subj+1,1+((col*3)+1:(col*3)+3)}=(nanmean(reg_correl(((conds(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                %ResultsCueImSim.(reg_names{reg}).AccCueOnly{subj+1,1+((col*3)+1:(col*3)+3)}=(nanmean(reg_correl(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                %ResultsCueImSim.(reg_names{reg}).NonAcc{subj+1,1+((col*3)+1:(col*3)+3)}=(nanmean(reg_correl(((conds_non_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                %ResultsCueImSim.(reg_names{reg}).AccCueAnd0Changes{subj+1,1+((col*3)+1:(col*3)+3)}=(nanmean(reg_correl_acc(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1) + (noChngeMat(:,5)==0))==3,:)));
                                
                            end
                        else %don't remove same block - that's an issue - you need to remove. accuracy is not updated
                            if reg==1 %collect the number of items
                                %note that this may be inaccurate for subjects like AK
                                %that have missing values because the 0-changes is
                                %in the blocks that he misses - these are the "0"
                                %in the first column of noChngeMat
                                % however,as long as I have noChngeMat(:,3)==1 as a
                                % condition - it's okay, becasue for every missing
                                % trial, noChngeMat(:,3)==0, so these trials will
                                % be removed
                                
                                ResultsCueImSimItemsCount.CueNoChange.all_items(subj,c)=sum(((conds(:,c)==1) + (noChngeMat(:,3)==1))==2);
                                ResultsCueImSimItemsCount.CueNoChange.Acc(subj,c)=sum(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1))==2);
                                ResultsCueImSimItemsCount.CueNoChange.NonAcc(subj,c)=sum(((conds_non_accurate(:,c)==1) + (noChngeMat(:,3)==1))==2);
                            end
                            if FishTrans
                                %fisher transform:
                                ResultsCueImSimOnlyNum.(reg_names{reg}).all_items(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl(((conds(:,c)==1) + (noChngeMat(:,3)==1))==2,:)));
                                ResultsCueImSimOnlyNum.(reg_names{reg}).Acc(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1))==2,:)));
                                ResultsCueImSimOnlyNum.(reg_names{reg}).NonAcc(subj,(col*3)+1:(col*3)+3)=atanh(nanmean(reg_correl(((conds_non_accurate(:,c)==1) + (noChngeMat(:,3)==1))==2,:)));
                                
                                %ResultsCueImSim.(reg_names{reg}).all_items{subj+1,1+((col*3)+1:(col*3)+3)}=atanh(nanmean(reg_correl(((conds(:,c)==1) + (noChngeMat(:,3)==1))==2,:)));
                                %ResultsCueImSim.(reg_names{reg}).Acc{subj+1,1+((col*3)+1:(col*3)+3)}=atanh(nanmean(reg_correl(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1))==2,:)));
                                %ResultsCueImSim.(reg_names{reg}).NonAcc{subj+1,1+((col*3)+1:(col*3)+3)}=atanh(nanmean(reg_correl(((conds_non_accurate(:,c)==1) + (noChngeMat(:,3)==1))==2,:)));
                            else
                                ResultsCueImSimOnlyNum.(reg_names{reg}).all_items(subj,(col*3)+1:(col*3)+3)=(nanmean(reg_correl(((conds(:,c)==1) + (noChngeMat(:,3)==1))==2,:)));
                                ResultsCueImSimOnlyNum.(reg_names{reg}).Acc(subj,(col*3)+1:(col*3)+3)=(nanmean(reg_correl(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1))==2,:)));
                                ResultsCueImSimOnlyNum.(reg_names{reg}).NonAcc(subj,(col*3)+1:(col*3)+3)=(nanmean(reg_correl(((conds_non_accurate(:,c)==1) + (noChngeMat(:,3)==1))==2,:)));
                                
                                %ResultsCueImSim.(reg_names{reg}).all_items{subj+1,1+((col*3)+1:(col*3)+3)}=(nanmean(reg_correl(((conds(:,c)==1) + (noChngeMat(:,3)==1))==2,:)));
                                %ResultsCueImSim.(reg_names{reg}).Acc{subj+1,1+((col*3)+1:(col*3)+3)}=(nanmean(reg_correl(((conds_accurate(:,c)==1) + (noChngeMat(:,3)==1))==2,:)));
                                %ResultsCueImSim.(reg_names{reg}).NonAcc{subj+1,1+((col*3)+1:(col*3)+3)}=(nanmean(reg_correl(((conds_non_accurate(:,c)==1) + (noChngeMat(:,3)==1))==2,:)));
                            end
                        end
                    end
                end %ends the if for > 10 voxels - to run the analysis
            end %ends the if for > 10 voxels - to select voxels
        end %ends the if for subject EB
    end %ends the regions loop
    
end %ends the subjects loop
save(results_fname,'ResultsCueImSim', 'ResultsCueImSimItemsCount', 'ResultsCueImSimOnlyNum', 'ResultsCueImSimAllTrials');
end





