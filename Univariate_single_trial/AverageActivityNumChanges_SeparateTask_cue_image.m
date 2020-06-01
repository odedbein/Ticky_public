function [ResultsAverageActivity, ResultsAverageActivityItemsCount, ResultsAverageActivityOnlyNum,ActiveVoxels]=AverageActivityNumChanges_SeparateTask_cue_image(engram)

warning('off','all')

if engram
    mydir='/data/Bein';
else
    mydir='/Volumes/data/Bein';
end
proj_dir=fullfile(mydir,'TickyReanalysis');
rmpath('/Volumes/data/Bein/fMRI_course/AnalysisScripts');
results_fname=sprintf('Univar_singleTrials_withVoxelSelection.mat');
results_fname=fullfile(proj_dir,'results',results_fname);
vox_sel={'MedianSplit','High2Thirds'};

%THIS HAS ALL OF THEM: subjects={'AB'; 'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD'; 'JG'; 'JM'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};
%subjects={'JA';'JG';'JW';'YE'; 'AR'}; - the ones that I editted their
%behavioral files
%LD and AK - only 216 trials

subjects={'AB';'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD';'JG'; 'JM'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};

subj_dir=fullfile(proj_dir,'SubData');

trial_types={...%R - changes in room, I - changes in items
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
num_comparisons=10;%0-4 changes+1

reg_names={'hipp',...%1
            'lhipp',...%2
            'rhipp',...%3
            'ant_hipp',...%4
            'mid_hipp',...%5
            'post_hipp',...%6
            'ant_lhipp',...%7
            'mid_lhipp',...%8
            'post_lhipp',...%9
            'ant_rhipp',...%10
            'mid_rhipp',...%11
            'post_rhipp',...%12
            'CA1',...%13
            'CA23DG',...%14
            'subiculum',...%15
            'lCA1',...%16
            'lCA23DG',...%17
            'lsubiculum',...%18
            'rCA1',...%19
            'rCA23DG',...%20
            'rsubiculum'...%21
            'lCA1_ant',...%1
            'lCA1_mid',...%2
            'lCA1_post',...%3
            'rCA1_ant',...%4
            'rCA1_mid',...%5
            'rCA1_post',...%6
            'lCA23DG_ant',...%7
            'lCA23DG_mid',...%8
            'lCA23DG_post',...%9
            'rCA23DG_ant',...%10
            'rCA23DG_mid',...%11
            'rCA23DG_post',...%12
            'lSub_ant',...%13
            'lSub_mid',...%14
            'lSub_post',...%15
            'rSub_ant',...%16
            'rSub_mid',...%17
            'rSub_post',...%18
            'Para',...%22
            'Peri',...%23
            'Ent',...%24
            'rPara',...%25
            'lPara',...%26
            'rPeri',...%27
            'lPeri',...%28
            'rEnt',...%29
            'lEnt',...%30
            'rPara_thirds',...%1
            'lPara_thirds',...%2
            'rPeri_thirds',...%3
            'lPeri_thirds',...%4
            'rMid_ParahippGyrus',...%5
            'lMid_ParahippGyrus',...%6
            'vta_epi','mean_SN_thr60', 'mean_SN_R_thr60','mean_SN_L_thr60','mean_VTA_thr60','mean_VTA_R_thr60','mean_VTA_L_thr60',...
            'mean_SN_thr80', 'mean_SN_R_thr80','mean_SN_L_thr80','mean_VTA_thr80','mean_VTA_R_thr80','mean_VTA_L_thr80',...
            'mean_VTA_thr95','mean_VTA_R_thr95','mean_VTA_L_thr95',...
            'rputamen','lputamen','putamen','rcaudate','lcaudate','caudate','racc','lacc','acc',...
            };
        
        
MTL_regs={'Para',...%22
            'Peri',...%23
            'Ent',...%24
            'rPara',...%25
            'lPara',...%26
            'rPeri',...%27
            'lPeri',...%28
            'rEnt',...%29
            'lEnt',...%30
            };
        
thirds_regs={
            'lCA1_ant',...%1
            'lCA1_mid',...%2
            'lCA1_post',...%3
            'rCA1_ant',...%4
            'rCA1_mid',...%5
            'rCA1_post',...%6
            'lCA23DG_ant',...%7
            'lCA23DG_mid',...%8
            'lCA23DG_post',...%9
            'rCA23DG_ant',...%10
            'rCA23DG_mid',...%11
            'rCA23DG_post',...%12
            'lSub_ant',...%13
            'lSub_mid',...%14
            'lSub_post',...%15
            'rSub_ant',...%16
            'rSub_mid',...%17
            'rSub_post',...%18
            };
        
thirds_MTL_regs={ 
    'rPara_thirds',...%1
    'lPara_thirds',...%2
    'rPeri_thirds',...%3
    'lPeri_thirds',...%4
    'rMid_ParahippGyrus',...%5
    'lMid_ParahippGyrus',...%6
            };        
        
midbrain_regs={'vta_epi','mean_SN_thr60', 'mean_SN_R_thr60','mean_SN_L_thr60','mean_VTA_thr60','mean_VTA_R_thr60','mean_VTA_L_thr60',...
            'mean_SN_thr80', 'mean_SN_R_thr80','mean_SN_L_thr80','mean_VTA_thr80','mean_VTA_R_thr80','mean_VTA_L_thr80',...
            'mean_VTA_thr95','mean_VTA_R_thr95','mean_VTA_L_thr95',...
            };
        
striatum_regs={'rputamen','lputamen','putamen','rcaudate','lcaudate','caudate','racc','lacc','acc',...
             };
                 
ResultsAverageActivity={};
ResultsAverageActivityOnlyNum={};
ResultsAverageActivityItemsCount={};
ActiveVoxels={};
%prepare the header:
for reg=1:numel(reg_names)
    ResultsAverageActivity.(reg_names{reg}).image.all_items{1,1}='subjects';
    ResultsAverageActivity.(reg_names{reg}).image.all_items(1,(2:num_comparisons+1))={'lay: 0-changes','lay: 1-change','lay: 2-changes','lay: 3-changes','lay: 4-changes','item: 0-changes','item: 1-change','item: 2-changes','item: 3-changes','item: 4-changes'};
    ResultsAverageActivity.(reg_names{reg}).image.Acc=ResultsAverageActivity.(reg_names{reg}).image.all_items;
    ResultsAverageActivity.(reg_names{reg}).image.NonAcc=ResultsAverageActivity.(reg_names{reg}).image.all_items;
    
    ResultsAverageActivity.(reg_names{reg}).image.MedianSplit=ResultsAverageActivity.(reg_names{reg}).image;
    ResultsAverageActivity.(reg_names{reg}).image.High2Thirds=ResultsAverageActivity.(reg_names{reg}).image;
    ResultsAverageActivity.(reg_names{reg}).cue=ResultsAverageActivity.(reg_names{reg}).image;
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
    load(fullfile(subj_dir,char(subjects(subj)),'data','reg_dataSubRegionsThirds.mat'),'reg_data');
    Thirds_hipp_data=reg_data;
    load(fullfile(subj_dir,char(subjects(subj)),'data','reg_dataMTL.mat'),'reg_data');
    MTL_data=reg_data;
    load(fullfile(subj_dir,char(subjects(subj)),'data','reg_dataMTLthirds.mat'),'reg_data');
    Thirds_MTL_data=reg_data;
    load(fullfile(subj_dir,char(subjects(subj)),'data','reg_data_midbrain.mat'),'reg_data');
    midbrain_data=reg_data;
    load(fullfile(subj_dir,char(subjects(subj)),'data','reg_data_striatum.mat'),'reg_data');
    striatum_data=reg_data;
    load(fullfile(subj_dir,char(subjects(subj)),'data','reg_data.mat'),'reg_data');
    %average on all voxels to get the mean activity in this region:
    for mreg=1:numel(thirds_regs)
        reg_data.(thirds_regs{mreg})=Thirds_hipp_data.(thirds_regs{mreg});
    end
    for mreg=1:numel(MTL_regs)
        reg_data.(MTL_regs{mreg})=MTL_data.(MTL_regs{mreg});
    end
    
    for mreg=1:numel(thirds_MTL_regs)
        reg_data.(thirds_MTL_regs{mreg})=Thirds_MTL_data.(thirds_MTL_regs{mreg});
    end
    
    for mreg=1:numel(midbrain_regs)
        reg_data.(midbrain_regs{mreg})=midbrain_data.(midbrain_regs{mreg});
    end
    
    for mreg=1:numel(striatum_regs)
        reg_data.(striatum_regs{mreg})=striatum_data.(striatum_regs{mreg});
    end
    
    reg_data_in_voxels=reg_data; %just to keep track of the number of voxels
    
    for reg=1:numel(reg_names)
        reg_data.(reg_names{reg}).cue=mean(reg_data.(reg_names{reg}).cue);
        reg_data.(reg_names{reg}).Image=mean(reg_data.(reg_names{reg}).Image);
    end
    
    
    %load the behavioral data, to know which items go to which condition:
    behav_filename=fullfile(subj_dir,char(subjects(subj)),sprintf('output_subject_%s.txt',subjects{subj}));
    if strcmp(char(subjects(subj)),'AR')
        subj_behavior=textscan(fopen(behav_filename),'%d %d %d %.1f %s %s %d %s %.1f %s');
    else
        subj_behavior=textscan(fopen(behav_filename),'%d %d %d %.1f %s %d %s %.1f %s %s %.3f %.3f %d');
    end
    
    %gett the accuracy and the trial type
    acc=subj_behavior{3}(1:num_trials);
    type=subj_behavior{5}(1:num_trials);
    if strcmp(char(subjects(subj)),'AR')
        task=subj_behavior{8}(1:num_trials);
    else
        task=subj_behavior{7}(1:num_trials);
    end
    
    if strcmp(char(subjects(subj)),'LD') %need to cut the first two sessions - brain data only from session 3.
        acc=subj_behavior{3}(((270/10*2)+1):end);
        type=subj_behavior{5}(((270/10*2)+1):end);
        task=subj_behavior{7}(((270/10*2)+1):end);
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
    ResultsAverageActivityItemsCount.all_items{subj,1}=subjects{subj};
    ResultsAverageActivityItemsCount.all_items(subj,(2:num_comparisons+1))=num2cell(sum(conds));
    
    %prepare the matrix with only accuracte responses
    conds_accurate=conds;
    conds_accurate(~acc,:)=0;
    %save the number of items per each subject
    ResultsAverageActivityItemsCount.Acc{subj,1}=subjects{subj};
    ResultsAverageActivityItemsCount.Acc(subj,(2:num_comparisons+1))=num2cell(sum(conds_accurate));
    
    %prepare the matrix with only NON accuracte responses
    conds_non_accurate=conds;
    conds_non_accurate(acc==1,:)=0;
    %save the number of items per each subject
    ResultsAverageActivityItemsCount.NonAcc{subj,1}=subjects{subj};
    ResultsAverageActivityItemsCount.NonAcc(subj,(2:num_comparisons+1))=num2cell(sum(conds_non_accurate));
    
    %now get the data for each region and calculate the average activity per condition:
    for reg=1:numel(reg_names)
        %put the subj_name in the reuslts structure:
        ResultsAverageActivity.(reg_names{reg}).cue.all_items(subj+1,1)=subjects(subj);
        ResultsAverageActivity.(reg_names{reg}).image.all_items(subj+1,1)=subjects(subj);
        ResultsAverageActivity.(reg_names{reg}).Acc(subj+1,1)=subjects(subj);
        ResultsAverageActivity.(reg_names{reg}).NonAcc(subj+1,1)=subjects(subj);
        ResultsAverageActivityVoxelsCount.(reg_names{reg}).all_items(subj,1)=subjects(subj);
        ResultsAverageActivityVoxelsCount.(reg_names{reg}).all_items(subj,2)=num2cell(size(reg_data_in_voxels.(reg_names{reg}).cue,1));
        
        cue_data1=reg_data.(reg_names{reg}).cue;
        
        %check the data for no nans:
        if isempty(cue_data1)
            fprintf(sprintf('region %s has no cue data \n',reg_names{reg}))
        elseif ~isempty(find(isnan(cue_data1)))
            fprintf(sprintf('cue data has nans region %s \n',reg_names{reg}))
        elseif size(reg_data_in_voxels.(reg_names{reg}).cue,1)<10
            fprintf(sprintf('cue data has less than 10 voxels region %s \n',reg_names{reg}))
        else %reginio is fine
%             if strcmp(char(subjects(subj)),'EB') && reg>51
%                 fprintf(sprintf('analyzing reg %s \n',reg_names{reg1}))
%             end
            %%compute cue activity per task:
            for c=1:size(conds,2) %for each condition
                %compute activation for all_items:
                curr_trials=(conds(:,c)==1);
                AvReg1=nanmean(cue_data1(curr_trials));
                ResultsAverageActivityOnlyNum.(reg_names{reg}).cue.all_items(subj,c)=AvReg1;
                ResultsAverageActivity.(reg_names{reg}).cue.all_items{subj+1,c+1}=AvReg1;
                
                %compute activation for accurate responses:
                curr_trials=(conds_accurate(:,c)==1);
                AvReg1=nanmean(cue_data1(curr_trials));
                ResultsAverageActivityOnlyNum.(reg_names{reg}).cue.Acc(subj,c)=AvReg1;
                ResultsAverageActivity.(reg_names{reg}).cue.Acc{subj+1,c}=AvReg1;
                
                %compute activation for non-accurate responses:
                curr_trials=(conds_non_accurate(:,c)==1);
                if length(find(curr_trials))>=10
                    AvReg1=nanmean(cue_data1(curr_trials));
                    ResultsAverageActivityOnlyNum.(reg_names{reg}).cue.NonAcc(subj,c)=AvReg1;
                    ResultsAverageActivity.(reg_names{reg}).cue.NonAcc{subj+1,c}=AvReg1;
                else
                    ResultsAverageActivityOnlyNum.(reg_names{reg}).cue.NonAcc(subj,c)=nan;
                    ResultsAverageActivity.(reg_names{reg}).cue.NonAcc{subj+1,c+1}=nan;
                end
            end
            
        end
        
        %%compute images activity:
        image_data1=reg_data.(reg_names{reg}).Image; %that's averaged on voxels
        image_data_vox=reg_data_in_voxels.(reg_names{reg}).Image; %that's per voxel, for voxel selection
        cue_data_vox=reg_data_in_voxels.(reg_names{reg}).cue; %that's per voxel, for voxel selection
        image_data_av_cond=[];
        %check the data for no nans:
        if isempty(image_data1)
            fprintf(sprintf('region %s has no image data \n',reg_names{reg}))
        elseif ~isempty(find(isnan(image_data1)))
            fprintf(sprintf('image data has nans region %s \n',reg_names{reg}))
        elseif size(reg_data_in_voxels.(reg_names{reg}).Image,1)<10
            fprintf(sprintf('image data has less than 10 voxels region %s \n',reg_names{reg}))
        else %region is fine - compute
            %%compute activity based on conditions:
            for c=1:size(conds,2) %for each condition
%                 if strcmp(char(subjects(subj)),'EB') && reg>51
%                     fprintf(sprintf('analyzing reg %s \n',reg_names{reg1}))
%                 end
                %compute activation for all_items:
                curr_trials=(conds(:,c)==1);
                AvReg1=nanmean(image_data1(curr_trials));
                ResultsAverageActivityOnlyNum.(reg_names{reg}).image.all_items(subj,c)=AvReg1;
                ResultsAverageActivity.(reg_names{reg}).image.all_items{subj+1,c+1}=AvReg1;
                image_data_av_cond=[image_data_av_cond nanmean(image_data_vox(:,curr_trials),2)];
                %compute activation for accurate responses:
                curr_trials=(conds_accurate(:,c)==1);
                AvReg1=nanmean(image_data1(curr_trials));
                ResultsAverageActivityOnlyNum.(reg_names{reg}).image.Acc(subj,c)=AvReg1;
                ResultsAverageActivity.(reg_names{reg}).image.Acc{subj+1,c+1}=AvReg1;
                
                %compute activation for non-accurate responses:
                curr_trials=(conds_non_accurate(:,c)==1);
                if length(find(curr_trials))>=10
                    AvReg1=nanmean(image_data1(curr_trials));
                    ResultsAverageActivityOnlyNum.(reg_names{reg}).image.NonAcc(subj,c)=AvReg1;
                    ResultsAverageActivity.(reg_names{reg}).image.NonAcc{subj+1,c+1}=AvReg1;
                else
                    ResultsAverageActivityOnlyNum.(reg_names{reg}).image.NonAcc(subj,c)=nan;
                    ResultsAverageActivity.(reg_names{reg}).image.NonAcc{subj+1,c+1}=nan;
                end
            end
            
            %% compute selective voxels:
            
            %compute average activation across both tasks and level of
            %changes:
            
            curr_reg=mean(image_data_av_cond,2);
            %median split:
            med_split=find(curr_reg>median(curr_reg));
            ActiveVoxels.(reg_names{reg}).MedianSplit.(subjects{subj})=med_split;
            
            %thirds split:
            sort_curr_reg=sort(curr_reg);
            curr_reg_size=length(sort_curr_reg);
            low_third_thresh=sort_curr_reg(floor(curr_reg_size/3));
            high_third_thresh=sort_curr_reg(floor(curr_reg_size*2/3));
            HighThird=find(curr_reg>high_third_thresh);
            High2Thirds=find(curr_reg>low_third_thresh);
            ActiveVoxels.(reg_names{reg}).HighThird.(subjects{subj})=HighThird;
            ActiveVoxels.(reg_names{reg}).High2Thirds.(subjects{subj})=High2Thirds;
            
            %t>0.5:
            more_than1=find(curr_reg>0.5);
            ActiveVoxels.(reg_names{reg}).TstatMoreThan05.(subjects{subj})=more_than1;
            
            %save number of voxels:
            ResultsAverageActivityVoxelsCount.(reg_names{reg}).MedianSplit(subj,1)=subjects(subj);
            ResultsAverageActivityVoxelsCount.(reg_names{reg}).MedianSplit(subj,2)=num2cell(length(med_split));
            ResultsAverageActivityVoxelsCount.(reg_names{reg}).TstatMoreThan05(subj,1)=subjects(subj);
            ResultsAverageActivityVoxelsCount.(reg_names{reg}).TstatMoreThan05(subj,2)=num2cell(length(more_than1));
            ResultsAverageActivityVoxelsCount.(reg_names{reg}).HighThird(subj,1)=subjects(subj);
            ResultsAverageActivityVoxelsCount.(reg_names{reg}).HighThird(subj,2)=num2cell(length(HighThird));
            ResultsAverageActivityVoxelsCount.(reg_names{reg}).High2Thirds(subj,1)=subjects(subj);
            ResultsAverageActivityVoxelsCount.(reg_names{reg}).High2Thirds(subj,2)=num2cell(length(High2Thirds));
            
            
        for curr_sel=1:2 %median, then 2highthirds
            curr_vox=ActiveVoxels.(reg_names{reg}).(vox_sel{curr_sel}).(subjects{subj});
            image_data1=nanmean(image_data_vox(curr_vox,:)); %select voxels, then average on voxels
            cue_data1=nanmean(cue_data_vox(curr_vox,:));  
            for c=1:size(conds,2) %for each condition
%                 if strcmp(char(subjects(subj)),'EB') && reg>51
%                     fprintf(sprintf('analyzing reg %s \n',reg_names{reg1}))
%                 end
                %compute activation for all_items:
                curr_trials=(conds(:,c)==1);
                AvReg1=nanmean(image_data1(curr_trials));
                ResultsAverageActivityOnlyNum.(reg_names{reg}).image.(vox_sel{curr_sel}).all_items(subj,c)=AvReg1;
                ResultsAverageActivity.(reg_names{reg}).image.(vox_sel{curr_sel}).all_items{subj+1,c+1}=AvReg1;
                %cue
                AvReg1=nanmean(cue_data1(curr_trials));
                ResultsAverageActivityOnlyNum.(reg_names{reg}).cue.(vox_sel{curr_sel}).all_items(subj,c)=AvReg1;
                ResultsAverageActivity.(reg_names{reg}).cue.(vox_sel{curr_sel}).all_items{subj+1,c+1}=AvReg1;
                
                %compute activation for accurate responses:
                curr_trials=(conds_accurate(:,c)==1);
                AvReg1=nanmean(image_data1(curr_trials));
                ResultsAverageActivityOnlyNum.(reg_names{reg}).image.(vox_sel{curr_sel}).Acc(subj,c)=AvReg1;
                ResultsAverageActivity.(reg_names{reg}).image.(vox_sel{curr_sel}).Acc{subj+1,c+1}=AvReg1;
                %cue
                AvReg1=nanmean(cue_data1(curr_trials));
                ResultsAverageActivityOnlyNum.(reg_names{reg}).cue.(vox_sel{curr_sel}).Acc(subj,c)=AvReg1;
                ResultsAverageActivity.(reg_names{reg}).cue.(vox_sel{curr_sel}).Acc{subj+1,c+1}=AvReg1;
                
                %compute activation for non-accurate responses:
                curr_trials=(conds_non_accurate(:,c)==1);
                if length(find(curr_trials))>=10
                    AvReg1=nanmean(image_data1(curr_trials));
                    ResultsAverageActivityOnlyNum.(reg_names{reg}).image.(vox_sel{curr_sel}).NonAcc(subj,c)=AvReg1;
                    ResultsAverageActivity.(reg_names{reg}).image.(vox_sel{curr_sel}).NonAcc{subj+1,c+1}=AvReg1;
                    
                    AvReg1=nanmean(cue_data1(curr_trials));
                    ResultsAverageActivityOnlyNum.(reg_names{reg}).cue.(vox_sel{curr_sel}).NonAcc(subj,c)=AvReg1;
                    ResultsAverageActivity.(reg_names{reg}).cue.(vox_sel{curr_sel}).NonAcc{subj+1,c+1}=AvReg1;
                else
                    ResultsAverageActivityOnlyNum.(reg_names{reg}).image.(vox_sel{curr_sel}).NonAcc(subj,c)=nan;
                    ResultsAverageActivity.(reg_names{reg}).image.(vox_sel{curr_sel}).NonAcc{subj+1,c+1}=nan;
                    
                    ResultsAverageActivityOnlyNum.(reg_names{reg}).cue.(vox_sel{curr_sel}).NonAcc(subj,c)=nan;
                    ResultsAverageActivity.(reg_names{reg}).cue.(vox_sel{curr_sel}).NonAcc{subj+1,c+1}=nan;
                end
            end
        end
            
        end
        
        
    end %ends the reg1 loop
    
end %ends the subjects loop

save(results_fname,'ResultsAverageActivity', 'ResultsAverageActivityItemsCount', 'ResultsAverageActivityOnlyNum','ActiveVoxels');

end
        
        
        
        
    