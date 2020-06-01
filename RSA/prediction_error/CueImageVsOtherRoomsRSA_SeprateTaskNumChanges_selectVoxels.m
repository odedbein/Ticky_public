function [ResultsCueImSim, ResultsCueImSimItemsCount, ResultsCueImSimOnlyNum]=CueImageVsOtherRoomsRSA_SeprateTaskNumChanges_selectVoxels(ActiveVoxels,engram,remove_outliers)
% to select voxels, you need the activeVoxels structure:
%load('/Volumes/data/Bein/TickyReanalysis/results/Univar_properGLM_withVoxelSelection.mat')

%med_more: choose whether based on median or on more than.. currently have
%more t>05
% 1: median, 0: more than..

warning('off','all')

if engram
    mydir='/data/Bein';
else
    mydir='/Volumes/data/Bein';
end

med_ttl='HighThird';
glm_ttl='ProperGLM';
mm_ttl='Av5cond';


if remove_outliers
    disp('analysis removes outleirs');
    out_ttl='_outliers_removed';
else
    disp('analysis DOES NOT remove outleirs');
    out_ttl='';
end

proj_dir=fullfile(mydir,'TickyReanalysis');
rmpath('/Volumes/data/Bein/fMRI_course/AnalysisScripts');
results_fname=sprintf('CueImageVsOtherRoomsExcludeSameRoom_Ftrans_%s%s%s%s.mat',glm_ttl,mm_ttl,med_ttl,out_ttl);
results_fname=fullfile(proj_dir,'results','rsa_PE',results_fname);
%THIS HAS ALL OF THEM: subjects={'AB'; 'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD'; 'JG'; 'JM'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};
%JG has only 1199 images instead of 1200 - take care of it later
%subjects={'JA';'JG';'JW';'YE'; 'AR'}; - the ones that I editted their
%behavioral files
%LD and AK - only 216 trials - and the entire first block is no responses - maybe
%the data is from runs 2-8?
subjects={'AB'; 'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD';'JG'; 'JM'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};

subj_dir=fullfile(proj_dir,'SubData');
region_mat_dir='regions_mat'; %a directory to store the matfiles of each region
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
num_comparisons=5*2;%0-4 changes * lay/item

reg_names={'lCA1'};

ResultsCueImSim={};
ResultsCueImSimOnlyNum={};
ResultsCueImSimItemsCount={};

%prepare the header:
for reg=1:numel(reg_names)
    ResultsCueImSim.(reg_names{reg}).same.all_items{1,1}='subjects';
    ResultsCueImSim.(reg_names{reg}).same.all_items(1,(2:num_comparisons+1))={'lay: 0-changes','lay: 1-change','lay: 2-changes','lay: 3-changes','lay: 4-changes','item: 0-changes','item: 1-change','item: 2-changes','item: 3-changes','item: 4-changes'};
    ResultsCueImSim.(reg_names{reg}).same.Acc=ResultsCueImSim.(reg_names{reg}).same.all_items;
    ResultsCueImSim.(reg_names{reg}).same.NonAcc=ResultsCueImSim.(reg_names{reg}).same.all_items;
    
    ResultsCueImSim.(reg_names{reg}).same_others.all_items=ResultsCueImSim.(reg_names{reg}).same.all_items;
    ResultsCueImSim.(reg_names{reg}).same_others.Acc=ResultsCueImSim.(reg_names{reg}).same.all_items;
    ResultsCueImSim.(reg_names{reg}).same_others.NonAcc=ResultsCueImSim.(reg_names{reg}).same.all_items;
    
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
    
    %gett the accuracy, trial type and stim
    acc=subj_behavior{3}(1:num_trials);
    type=subj_behavior{5}(1:num_trials);
    if strcmp(char(subjects(subj)),'AR')
        task=subj_behavior{8}(1:num_trials);
        stim=subj_behavior{6}(1:num_trials);
    else
        task=subj_behavior{7}(1:num_trials);
        stim=subj_behavior{10}(1:num_trials);
    end
    %block number
    block=subj_behavior{1}(1:num_trials);
    
    if strcmp(char(subjects(subj)),'LD') %need to cut the first two sessions - brain data only from session 3.
        acc=subj_behavior{3}(((270/10*2)+1):end);
        type=subj_behavior{5}(((270/10*2)+1):end);
        task=subj_behavior{7}(((270/10*2)+1):end);
        stim=subj_behavior{10}(((270/10*2)+1):end);
        block=subj_behavior{1}(((270/10*2)+1):end);
    end
    
    %edit the stim to extract just the room's name
    for s=1:num_trials
        curr_room=stim{s};
        stim{s}=curr_room(7:12);
    end
    
    if isempty(acc)
        fprintf('shoot.. no accurcay column, check behavioral file\n');
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
    conds(:,1)= (conds_temp(:,1) & task_binary==1); %no changes
    conds(:,2)= (conds_temp(:,2) & task_binary==1);%1 change
    conds(:,3)= (conds_temp(:,3) & task_binary==1);%2 changes
    conds(:,4)= (conds_temp(:,4) & task_binary==1);%3 changes
    conds(:,5)= (conds_temp(:,5) & task_binary==1);%4 changes
    
    %item task
    conds(:,6)= (conds_temp(:,1) & task_binary==0); %no changes
    conds(:,7)= (conds_temp(:,2) & task_binary==0);%1 change
    conds(:,8)= (conds_temp(:,3) & task_binary==0);%2 changes
    conds(:,9)= (conds_temp(:,4) & task_binary==0);%3 changes
    conds(:,10)= (conds_temp(:,5) & task_binary==0);%4 changes
    
    
    %save the number of items per each subject
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
        %reg=17;
        % display(sprintf('subject %s, region %s \n',char(subjects(subj)),reg_names{reg}));
        if strcmp(char(subjects(subj)),'EB') && (strcmp(reg_names{reg},'lacc') || strcmp(reg_names{reg},'acc') || strcmp(reg_names{reg},'racc'))
            %display('skipping EB lacc');
            %skip this region - subj do not have a lacc region
        else
            %put the subj_name in the reuslts structure:
            ResultsCueImSim.(reg_names{reg}).same.all_items(subj+1,1)=subjects(subj);
            ResultsCueImSim.(reg_names{reg}).same.Acc(subj+1,1)=subjects(subj);
            ResultsCueImSim.(reg_names{reg}).same.NonAcc(subj+1,1)=subjects(subj);
            ResultsCueImSim.(reg_names{reg}).same_others.all_items(subj+1,1)=subjects(subj);
            ResultsCueImSim.(reg_names{reg}).same_others.Acc(subj+1,1)=subjects(subj);
            ResultsCueImSim.(reg_names{reg}).same_others.NonAcc(subj+1,1)=subjects(subj);
            
            cue_data=reg_data.(reg_names{reg}).cue;
            if size(cue_data,1) > 10 %only do if more than 10 voxels - otherwise there's no selection of voxels... and the analysis is meaningless
                %select voxels:
                curr_voxels=ActiveVoxels.((reg_names{reg})).(mm_ttl).(med_ttl).(subjects{subj});
                cue_data=cue_data(curr_voxels,:);
                %size(cue_data,1)
                %check the data for no nans:
                if ~isempty(find(isnan(cue_data)))
                    fprintf('cue data has nans\n')
                end
                image_data=reg_data.(reg_names{reg}).Image;
                image_data=image_data(curr_voxels,:);
                %check the data for no nans:
                if ~isempty(find(isnan(image_data)))
                    fprintf('cue data has nans\n')
                end
                
                %compute correlations between all cues and images:
                reg_correl=nan(num_trials,3);
                if size(cue_data,1) > 10 %only do the analysis for more than 10 voxels:
                    %% clean the data:
                    %it doesn't make sense to remove outliers for
                    %same and others together - because same is within trial, and
                    %others is across... so I do it separately. I'm sure there's a
                    %more elegant way to do it, but  to group all of the "other"
                    %correlations - I just run the loop twice...
                    if remove_outliers
                        all_others_final=[];
                        %this is just to gather all others:
                        for i=1:num_trials
                            curr_type=find_types(i,:)==1;
                            curr_task=task_binary(i);
                            all_others=find((find_types(:,curr_type) + (task_binary==curr_task))==2); %items that are of the same task and trial type
                            all_others(all_others==i)=[];%don't take current one.
                            %remove same block and same item:
                            cue=stim{i}(1:2); %this is the item of the current cue-probe
                            for oi=1:length(all_others)
                                %same block
                                if block(all_others(oi))==block(i)
                                    all_others(oi)=nan;
                                else %if removed for same block, no need to worry about item - and it wouldn't work bc Nan
                                    %same item:
                                    other=stim{all_others(oi)}(1:2);
                                    if strcmp(cue,other)
                                        all_others(oi)=nan;
                                    end
                                end
                            end
                            
                            all_others=all_others(~isnan(all_others));
                            other_corr=[];
                            for oi=1:length(all_others)
                                other_corr=[other_corr; corr(cue_data(:,i),image_data(:,all_others(oi)))];
                            end
                            all_others_final=[all_others_final;other_corr];
                        end
                        
                        %compute mean and SD for others:
                        mean_others=nanmean(all_others_final);
                        sd_others=nanstd(all_others_final);
                    end
                    
                    %% now actually compute the same and other correlations:
                    for i=1:num_trials
                        %compute the fisher transformed correlation value
                        %corr=corrcoef(cue_data(:,i),image_data(:,i));
                        reg_correl(i,1)=corr(cue_data(:,i),image_data(:,i));
                        curr_type=find_types(i,:)==1;
                        curr_task=task_binary(i);
                        %curr_block=floor((i-1)/27)+1;
                        all_others=find((find_types(:,curr_type) + (task_binary==curr_task))==2); %items that are of the same task and trial type
                        all_others(all_others==i)=[];%don't take current one.
                        %remove same block and same item:
                        cue=stim{i}(1:2); %this is the item of the current cue-probe
                        for oi=1:length(all_others)
                            %same block
                            if block(all_others(oi))==block(i)
                                all_others(oi)=nan;
                            else %if removed for same block, no need to worry about item - and it wouldn't work bc Nan
                                %same item:
                                other=stim{all_others(oi)}(1:2);
                                if strcmp(cue,other)
                                    all_others(oi)=nan;
                                end
                            end
                        end
                        
                        all_others=all_others(~isnan(all_others));
                        other_corr=[];
                        for oi=1:length(all_others)
                            other_corr=[other_corr corr(cue_data(:,i),image_data(:,all_others(oi)))];
                        end
                        
                        if remove_outliers
                            other_corr(other_corr >  (mean_others + 3*sd_others))=[];
                            other_corr(other_corr <  (mean_others - 3*sd_others))=[];
                        end
                        reg_correl(i,2)=nanmean(other_corr);
                    end
                    
                    %clean the data: it doesn't make sense to remove outliers for
                    %same and others together - because same is within trial, and
                    %others is across... so I do it separately - so here is the
                    %cleaning of the same trials:
                    if remove_outliers
                        Mcorrel=nanmean(reg_correl(:,1));
                        SDcorrel=nanstd(reg_correl(:,1));
                        reg_correl(reg_correl(:,1)>(Mcorrel+(3*SDcorrel)),:)=nan;
                        reg_correl(reg_correl(:,1)<(Mcorrel-(3*SDcorrel)),:)=nan;
                    end
                    
                    %compute the difference:
                    reg_correl(:,3)=reg_correl(:,1)-reg_correl(:,2);
                    
                    %% make the R data structure: include condition, match_mis contrast, task, block
                    num_changes=nan(size(reg_correl,1),1);
                    for c=1:size(conds_temp,2) %for each condition
                        num_changes(conds_temp(:,c)==1)=(c-1);
                    end
                    match_mis=num_changes>0;
                    
                    if any(isnan(num_changes))
                        fprintf('script is bad, nans in num_changes');
                    end
                    
                    %prep the mat:
                    curr_mat=[atanh(reg_correl) num_changes match_mis task_binary double(block)];
                    curr_mat(isnan(curr_mat(:,1)),:)=[];
                    ResultsCueImSimOnlyNum.(reg_names{reg}).AllTrialsR.(subjects{subj})=curr_mat;
                    
                    
                    %% average based on conditions:
                    for c=1:size(conds,2) %for each condition
                        %only same:
                        ResultsCueImSimOnlyNum.(reg_names{reg}).same.all_items(subj,c)=atanh(nanmean(reg_correl(conds(:,c)==1,1)));
                        ResultsCueImSimOnlyNum.(reg_names{reg}).same.Acc(subj,c)=atanh(nanmean(reg_correl(conds_accurate(:,c)==1,1)));
                        ResultsCueImSimOnlyNum.(reg_names{reg}).same.NonAcc(subj,c)=atanh(nanmean(reg_correl(conds_non_accurate(:,c)==1,1)));
                        ResultsCueImSim.(reg_names{reg}).same.all_items{subj+1,c+1}=atanh(nanmean(reg_correl(conds(:,c)==1,1)));
                        ResultsCueImSim.(reg_names{reg}).same.Acc{subj+1,c+1}=atanh(nanmean(reg_correl(conds_accurate(:,c)==1,1)));
                        ResultsCueImSim.(reg_names{reg}).same.NonAcc{subj+1,c+1}=atanh(nanmean(reg_correl(conds_non_accurate(:,c)==1,1)));
                        
                        %same-others
                        ResultsCueImSimOnlyNum.(reg_names{reg}).same_others.all_items(subj,c)=atanh(nanmean(reg_correl(conds(:,c)==1,3)));
                        ResultsCueImSimOnlyNum.(reg_names{reg}).same_others.Acc(subj,c)=atanh(nanmean(reg_correl(conds_accurate(:,c)==1,3)));
                        ResultsCueImSimOnlyNum.(reg_names{reg}).same_others.NonAcc(subj,c)=atanh(nanmean(reg_correl(conds_non_accurate(:,c)==1,3)));
                        ResultsCueImSim.(reg_names{reg}).same_others.all_items{subj+1,c+1}=atanh(nanmean(reg_correl(conds(:,c)==1,3)));
                        ResultsCueImSim.(reg_names{reg}).same_others.Acc{subj+1,c+1}=atanh(nanmean(reg_correl(conds_accurate(:,c)==1,3)));
                        ResultsCueImSim.(reg_names{reg}).same_others.NonAcc{subj+1,c+1}=atanh(nanmean(reg_correl(conds_non_accurate(:,c)==1,3)));
                    end
                end %ends the if for > 10 voxels - to run the analysis
            end %ends the if for > 10 voxels - to select voxels
        end %EB subject conditional
    end %ends the regions loop
    
end %ends the subjects loop

save(results_fname,'ResultsCueImSim', 'ResultsCueImSimItemsCount', 'ResultsCueImSimOnlyNum');
end






