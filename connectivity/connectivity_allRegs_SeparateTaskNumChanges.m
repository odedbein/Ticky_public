function [ResultsBetaSerConnectivity, ResultsBetaSerConnectivityItemsCount, ResultsBetaSerConnectivityOnlyNum]=connectivity_allRegs_SeparateTaskNumChanges()
%this script does not save the output - if you want that, add it
warning('off','all')
rmpath('/Volumes/Oded/Bein/fMRI_course/AnalysisScripts');

[~, hostname]=system('hostname');
if strcmp(hostname(1:6),'joanna')%the hostname command gives 1X7 char output, we only need the first 6.
    proj_dir='/Volumes/Oded/Bein/TickyReanalysis';
else
    proj_dir='/Volumes/davachilab/Bein/TickyReanalysis';
end

FisherTrans=1;

%LD - only 216 trials
%AK - 189 tials - eighth run was running the same block again, that run is excluded
subjects={'AB'; 'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD';'JG'; 'JM'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};

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
num_comparisons=10;%` 0-4*layout/item

reg_names={'CA1',...
            'CA23DG',...
            'lCA1',...
            'lCA23DG',...
            'rCA1',...
            'rCA23DG',...
            'Ent',...
            'rEnt',...
            'lEnt'...
            };
        
        
MTL_regs={'Ent',...
            'rEnt',...
            'lEnt'...
            };
        
ResultsBetaSerConnectivity={};
ResultsBetaSerConnectivityOnlyNum={};
ResultsBetaSerConnectivityItemsCount={};

%prepare the header:
for reg1=1:numel(reg_names)
        for reg2=reg1+1:numel(reg_names)
            ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).all_items{1,1}='subjects';
            ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).all_items(1,(2:num_comparisons+3))={'lay_cue','item_cue','lay: 0-changes','lay: 1-change','lay: 2-changes','lay: 3-changes','lay: 4-changes','item: 0-changes','item: 1-change','item: 2-changes','item: 3-changes','item: 4-changes'};
            ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).Acc{1,1}='subjects';
            ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).Acc(1,(2:num_comparisons+3))={'lay_cue','item_cue','lay: 0-changes','lay: 1-change','lay: 2-changes','lay: 3-changes','lay: 4-changes','item: 0-changes','item: 1-change','item: 2-changes','item: 3-changes','item: 4-changes'};
            ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).NonAcc{1,1}='subjects';
            ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).NonAcc(1,(2:num_comparisons+3))={'lay_cue','item_cue','lay: 0-changes','lay: 1-change','lay: 2-changes','lay: 3-changes','lay: 4-changes','item: 0-changes','item: 1-change','item: 2-changes','item: 3-changes','item: 4-changes'};
        end
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
    load(fullfile(subj_dir,char(subjects(subj)),'data','reg_dataMTL.mat'),'reg_data');
    MTL_data=reg_data;
    load(fullfile(subj_dir,char(subjects(subj)),'data','reg_data.mat'),'reg_data');
    %average on all voxels to get the mean activity in this region:
    for mreg=1:numel(MTL_regs)
        reg_data.(MTL_regs{mreg})=MTL_data.(MTL_regs{mreg});
    end
    
    reg_data_num_voxels=reg_data; %just to keep track of the number of voxels
    
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
    
    %get the accuracy and the trial type
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
    ResultsBetaSerConnectivityItemsCount.all_items{subj,1}=subjects{subj};
    ResultsBetaSerConnectivityItemsCount.all_items(subj,(2:num_comparisons+1))=num2cell(sum(conds));
    
    %prepare the matrix with only accuracte responses
    conds_accurate=conds;
    conds_accurate(~acc,:)=0;
    %save the number of items per each subject
    ResultsBetaSerConnectivityItemsCount.Acc{subj,1}=subjects{subj};
    ResultsBetaSerConnectivityItemsCount.Acc(subj,(2:num_comparisons+1))=num2cell(sum(conds_accurate));
    
    %prepare the matrix with only NON accuracte responses
    conds_non_accurate=conds;
    conds_non_accurate(acc==1,:)=0;
    %save the number of items per each subject
    ResultsBetaSerConnectivityItemsCount.NonAcc{subj,1}=subjects{subj};
    ResultsBetaSerConnectivityItemsCount.NonAcc(subj,(2:num_comparisons+1))=num2cell(sum(conds_non_accurate));
    
    %now get the data for each region and calculate the correlations between regions:
    %calculate
    for reg1=1:numel(reg_names)-1
        if size(reg_data_num_voxels.(reg_names{reg1}).cue,1)<10
            fprintf(sprintf('less than 10 voxels region %s \n',reg_names{reg1}))
        else
            
            for reg2=reg1+1:numel(reg_names)
                if size(reg_data_num_voxels.(reg_names{reg2}).cue,1)<10
                    fprintf(sprintf('less than 10 voxels region %s \n',reg_names{reg2}))
                else
                    %put the subj_name in the reuslts structure:
                    ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).all_items(subj+1,1)=subjects(subj);
                    ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).Acc(subj+1,1)=subjects(subj);
                    ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).NonAcc(subj+1,1)=subjects(subj);
                    
                    cue_data1=reg_data.(reg_names{reg1}).cue;
                    cue_data2=reg_data.(reg_names{reg2}).cue;
                    
                    %check the data for no nans:
                    if ~isempty(find(isnan(cue_data1))) || ~isempty(find(isnan(cue_data2)))
                        fprintf('cue data has nans\n')
                    end
                    
                    %%compute cue connectivity in each task:
                    for t=0:1
                        % all_items:
                        BseriesReg1=cue_data1(~task_binary==t); %stupid hack to have first layout, then item
                        BseriesReg2=cue_data2(~task_binary==t);
                        corr=corrcoef(BseriesReg1,BseriesReg2);
                        if FisherTrans
                            reg_correl=atanh(corr(2,1));
                        else
                            reg_correl=corr(2,1);
                        end
                        
                        ResultsBetaSerConnectivityOnlyNum.([reg_names{reg1} '_' reg_names{reg2}]).all_items(subj,t+1)=reg_correl;
                        ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).all_items{subj+1,t+1+1}=reg_correl;
                        % accurate responses:
                        BseriesReg1=cue_data1(acc==1);
                        BseriesReg2=cue_data2(acc==1);
                        corr=corrcoef(BseriesReg1,BseriesReg2);
                        if FisherTrans
                            reg_correl=atanh(corr(2,1));
                        else
                            reg_correl=corr(2,1);
                        end
                        ResultsBetaSerConnectivityOnlyNum.([reg_names{reg1} '_' reg_names{reg2}]).Acc(subj,t+1)=reg_correl;
                        ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).Acc{subj+1,t+1+1}=reg_correl;
                        % inaccurate responses:
                        BseriesReg1=cue_data1(acc==0);
                        BseriesReg2=cue_data2(acc==0);
                        corr=corrcoef(BseriesReg1,BseriesReg2);
                        if FisherTrans
                            reg_correl=atanh(corr(2,1));
                        else
                            reg_correl=corr(2,1);
                        end
                        ResultsBetaSerConnectivityOnlyNum.([reg_names{reg1} '_' reg_names{reg2}]).NonAcc(subj,t+1)=reg_correl;
                        ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).nonAcc{subj+1,t+1+1}=reg_correl;
                    end
                    
                    %%compute images connectivity:
                    image_data1=reg_data.(reg_names{reg1}).Image;
                    image_data2=reg_data.(reg_names{reg2}).Image;
                    
                    %check the data for no nans:
                    if ~isempty(find(isnan(image_data1))) || ~isempty(find(isnan(image_data2)))
                        fprintf('cue data has nans\n')
                    end
                    
                    %%compute connectivity based on conditions:
                    for c=1:size(conds,2) %for each condition
                        %compute corelations for all_items:
                        curr_trials=(conds(:,c)==1);
                        BseriesReg1=image_data1(curr_trials);
                        BseriesReg2=image_data2(curr_trials);
                        corr=corrcoef(BseriesReg1,BseriesReg2);
                        if FisherTrans
                            reg_correl=atanh(corr(2,1));
                        else
                            reg_correl=corr(2,1);
                        end
                        ResultsBetaSerConnectivityOnlyNum.([reg_names{reg1} '_' reg_names{reg2}]).all_items(subj,c+2)=reg_correl;
                        ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).all_items{subj+1,c+3}=reg_correl;
                        
                        %compute corelations for accurate responses:
                        curr_trials=(conds_accurate(:,c)==1);
                        BseriesReg1=image_data1(curr_trials);
                        BseriesReg2=image_data2(curr_trials);
                        corr=corrcoef(BseriesReg1,BseriesReg2);
                        if FisherTrans
                            reg_correl=atanh(corr(2,1));
                        else
                            reg_correl=corr(2,1);
                        end
                        ResultsBetaSerConnectivityOnlyNum.([reg_names{reg1} '_' reg_names{reg2}]).Acc(subj,c+2)=reg_correl;
                        ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).Acc{subj+1,c+3}=reg_correl;
                        
                        %compute corelations for non-accurate responses:
                        curr_trials=(conds_non_accurate(:,c)==1);
                        if length(find(curr_trials))>=10
                            BseriesReg1=image_data1(curr_trials);
                            BseriesReg2=image_data2(curr_trials);
                            corr=corrcoef(BseriesReg1,BseriesReg2);
                            if FisherTrans
                                reg_correl=atanh(corr(2,1));
                            else
                                reg_correl=corr(2,1);
                            end
                            ResultsBetaSerConnectivityOnlyNum.([reg_names{reg1} '_' reg_names{reg2}]).NonAcc(subj,c+2)=reg_correl;
                            ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).NonAcc{subj+1,c+3}=reg_correl;
                        else
                            ResultsBetaSerConnectivityOnlyNum.([reg_names{reg1} '_' reg_names{reg2}]).NonAcc(subj,c+2)=nan;
                            ResultsBetaSerConnectivity.([reg_names{reg1} '_' reg_names{reg2}]).NonAcc{subj+1,c+3}=nan;
                        end
                    end
                end%%ends the if on whether the reg2 has less then 10 voxels
            end%ends the reg2 loop
        end %ends the if on whether  reg1 has less then 10 voxels
    end %ends the reg1 loop
    
end %ends the subjects loop

        
        
        
        
    