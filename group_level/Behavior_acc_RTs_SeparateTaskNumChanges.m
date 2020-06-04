function [Results,ResultsOnlyNum]=Behavior_acc_RTs_SeparateTaskNumChanges()

warning('off','all')
rmpath('/Volumes/Oded/Bein/fMRI_course/AnalysisScripts');

engram=0;
if engram
    proj_dir='/data/Bein/TickyReanalysis';
else
    proj_dir='/Volumes/data/Bein/TickyReanalysis';
end

%subjects={'JA';'JG';'JW';'YE'; 'AR'}; - the ones that I editted their behavioral files
%LD: 216 trials, AK: 189 trials
%AD - not enough entorhinal voxels - excluded from ms.

%subjects={'AB'; 'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD';'JG'; 'JM'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};
subjects={'AB'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD';'JG'; 'JM'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};

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

Results={};

%prepare the header:


Results.Acc={'subjects','lay: 0-changes','lay: 1-change','lay: 2-changes','lay: 3-changes','lay: 4-changes','item: 0-changes','item: 1-change','item: 2-changes','item: 3-changes','item: 4-changes'};
Results.all_items=Results.Acc;
Results.RT.all_items=Results.Acc;
Results.RT.Acc=Results.Acc;


for subj=1:numel(subjects)
    if strcmp(char(subjects(subj)),'LD')
        num_trials=216;
    elseif strcmp(char(subjects(subj)),'AK')
        num_trials=189;
    else
        num_trials=270;
    end
    fprintf('analyzing subj %s\n',char(subjects(subj)));
    
    
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
        RTs=subj_behavior{7}(1:num_trials);
    else
        task=subj_behavior{7}(1:num_trials);
        RTs=subj_behavior{6}(1:num_trials);
    end
    
    if strcmp(char(subjects(subj)),'LD') %need to cut the first two sessions - brain data only from session 3.
        acc=subj_behavior{3}(((270/10*2)+1):end);
        type=subj_behavior{5}(((270/10*2)+1):end);
        task=subj_behavior{7}(((270/10*2)+1):end);
        RTs=subj_behavior{6}(((270/10*2)+1):end);
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
    Results.all_items{subj+1,1}=subjects{subj};
    Results.all_items(subj+1,(2:num_comparisons+1))=num2cell(sum(conds));
    
    %prepare the matrix with only accuracte responses
    conds_accurate=conds;
    conds_accurate(~acc,:)=0;
    %save the number of items per each subject
    Results.Acc{subj+1,1}=subjects{subj};
    Results.Acc(subj+1,(2:num_comparisons+1))=num2cell((sum(conds_accurate))./sum(conds));
    
    %now RTs for each condition
    Results.RT.all_items{subj+1,1}=subjects{subj};
    Results.RT.Acc{subj+1,1}=subjects{subj};
    for c=1:size(conds,2) %for each condition
        %compute corelations for all_items:
        curr_trials=(conds(:,c)==1);
        Results.RT.all_items{subj+1,c+1}=nanmean(RTs(curr_trials));
        
        curr_trials=(conds_accurate(:,c)==1);
        Results.RT.Acc{subj+1,c+1}=nanmean(RTs(curr_trials));
    end
end %ends the subjects loop


%%now make the onlyNumber struct
ResultsOnlyNum.all_items=cell2mat(Results.all_items(2:end,2:end));
ResultsOnlyNum.Acc=cell2mat(Results.Acc(2:end,2:end));
ResultsOnlyNum.RT.all_items=cell2mat(Results.RT.all_items(2:end,2:end));
ResultsOnlyNum.RT.Acc=cell2mat(Results.RT.Acc(2:end,2:end));


end

        
        
        
        
    