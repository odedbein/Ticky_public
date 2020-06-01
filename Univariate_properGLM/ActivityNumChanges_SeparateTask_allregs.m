function [ResultsAverageActivity, ResultsAverageActivityVoxelsCount, ResultsAverageActivityOnlyNum,ActiveVoxels]=ActivityNumChanges_SeparateTask_allregs(engram)

warning('off','all')
rmpath('/Volumes/Oded/Bein/fMRI_course/AnalysisScripts');


if engram
    mydir='/data/Bein';
else
    mydir='/Volumes/data/Bein';
end
proj_dir=fullfile(mydir,'TickyReanalysis');
rmpath('/Volumes/data/Bein/fMRI_course/AnalysisScripts');
results_fname=sprintf('Univar_properGLM_withVoxelSelection.mat');
results_fname=fullfile(proj_dir,'results',results_fname);

%THIS HAS ALL OF THEM: subjects={'AB'; 'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD'; 'JG'; 'JM'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};
%LD and AK - only 216/189 trials

subjects={'AB';'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD';'JG'; 'JM'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};

subj_dir=fullfile(proj_dir,'SubData');
region_mat_dir='regions_mat'; %a directory to store the matfiles of each region
trial_types={...%R - changes in room, I - changes in items
    'R0I0';%1
    'R0I1';%2
    'R0I2';%3
    'R1I0';%4
    'R1I1';%5
    'R1I2';%6
    'R2I0';%7
    'R2I1';%8
    'R2I2';%9
    '1-changes';
    '2-changes';
    '3-changes'...
    };
nCope=12; %number of total copes - each one will have: tstat1=layout task, tstat2=item task
num_comparisons=nCope*2;%number of copes, in each task
%here all of of them
reg_names={'lCA1',...
    'lCA23DG',...
    'rCA1',...
    'rCA23DG',...
    'rEnt',...
    'lEnt'...
    };

ResultsAverageActivity={};
ResultsAverageActivityOnlyNum={};
ResultsAverageActivityVoxelsCount={};
ActiveVoxels={};

%prepare the header:
for reg=1:numel(reg_names)
    ResultsAverageActivity.(reg_names{reg}).all_items.lay{1,1}='subjects';
    ResultsAverageActivity.(reg_names{reg}).all_items.lay(1,(2:nCope+1))=trial_types';
    ResultsAverageActivity.(reg_names{reg}).all_items.item{1,1}='subjects';
    ResultsAverageActivity.(reg_names{reg}).all_items.item(1,(2:nCope+1))=trial_types';
end

for subj=1:numel(subjects)
    
    load(fullfile(subj_dir,char(subjects(subj)),'data','Univariate_reg_data.mat'),'reg_data');
    %average on all voxels to get the mean activity in this region:
    
    reg_data_in_voxels=reg_data; %to keep track of the number of voxels, and compute selectivity:
    
    for reg=1:numel(reg_names)
        reg_data.(reg_names{reg}).lay=nanmean(reg_data.(reg_names{reg}).lay);
        reg_data.(reg_names{reg}).item=nanmean(reg_data.(reg_names{reg}).item);
    end
    
    %get the data for each region and calculate the average activity per condition:
    for reg=1:numel(reg_names)
        %put the subj_name in the reuslts structure:
        ResultsAverageActivity.(reg_names{reg}).all_items.lay(subj+1,1)=subjects(subj);
        ResultsAverageActivity.(reg_names{reg}).all_items.item(subj+1,1)=subjects(subj);
        ResultsAverageActivityVoxelsCount.(reg_names{reg}).all_items(subj,1)=subjects(subj);
        ResultsAverageActivityVoxelsCount.(reg_names{reg}).all_items(subj,2)=num2cell(size(reg_data_in_voxels.(reg_names{reg}).lay,1));
        
        %check the data for no nans:
        if isempty(reg_data.(reg_names{reg}).lay)
            fprintf(sprintf('region %s has no lay data \n',reg_names{reg}))
        elseif size(reg_data_in_voxels.(reg_names{reg}).lay,1)<10
            fprintf(sprintf('data has less than 10 voxels region %s \n',reg_names{reg}))
        else %region is fine - put the correct cope data in the correct field
            %%compute activity:
            AvReg1=reg_data.(reg_names{reg}).lay;
            ResultsAverageActivityOnlyNum.(reg_names{reg}).all_items.lay(subj,:)=AvReg1;
            ResultsAverageActivity.(reg_names{reg}).all_items.lay(subj+1,2:end)=num2cell(AvReg1);
            AvReg1=reg_data.(reg_names{reg}).item;
            ResultsAverageActivityOnlyNum.(reg_names{reg}).all_items.item(subj,:)=AvReg1;
            ResultsAverageActivity.(reg_names{reg}).all_items.item(subj+1,2:end)=num2cell(AvReg1);
            %% compute selective voxels:
            %compute average activation across both tasks and level of
            %changes:
            curr_reg=reg_data_in_voxels.(reg_names{reg}).lay(:,[1 9:12]);
            %there are 9 trial types. Activation in them is in the first 9 columns
            %- the other colums are specific contrasts. Use only first 9.
            curr_reg=[curr_reg reg_data_in_voxels.(reg_names{reg}).item(:,[1 9:12])];
            
            curr_reg=mean(curr_reg,2);
            
            %thirds split:
            sort_curr_reg=sort(curr_reg);
            curr_reg_size=length(sort_curr_reg);
            low_third_thresh=sort_curr_reg(floor(curr_reg_size/3));
            High2Thirds=find(curr_reg>low_third_thresh);
            ActiveVoxels.(reg_names{reg}).Av5cond.High2Thirds.(subjects{subj})=High2Thirds;
            
            %save number of voxels:
            ResultsAverageActivityVoxelsCount.(reg_names{reg}).High2Thirds(subj,1)=subjects(subj);
            ResultsAverageActivityVoxelsCount.(reg_names{reg}).High2Thirds(subj,2)=num2cell(length(High2Thirds));
        
        end
        
    end %ends the reg loop
    
end %ends the subjects loop
save(results_fname,'ResultsAverageActivity', 'ResultsAverageActivityVoxelsCount', 'ResultsAverageActivityOnlyNum','ActiveVoxels');
end





