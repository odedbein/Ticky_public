function get_Im_data()

[~, hostname]=system('hostname');
if strcmp(hostname(1:6),'joanna')%the hostname command gives 1X7 char output, we only need the first 6.
    proj_dir='/Volumes/Oded/Bein/TickyReanalysis';
    % software and data paths
    warning('off','all')
    addpath(genpath('/Volumes/Oded/Bein/AnalysisMethods/General_scripts'))
    rmpath('/Volumes/Oded/Bein/fMRI_course/AnalysisMethods/AnalysisScripts');

else
    proj_dir='/Volumes/davachilab/Bein/TickyReanalysis';
    % software and data paths
    warning('off','all')
    addpath(genpath('/Volumes/davachilab/Bein/AnalysisMethods/General_scripts'))
    rmpath('/Volumes/davachilab/Bein/fMRI_course/AnalysisMethods/AnalysisScripts');
end

%THIS HAS ALL OF THEM: 
subjects={'AB'; 'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD'; 'JG'; 'JM'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};

subj_dir=fullfile(proj_dir,'SubData');
region_mat_dir='regions_mat'; %a directory to store the matfiles of each region
roi_dir='ROIs';

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

analysis_dir='UnivariateGLM/AllTrials/SubjLevelTasks.gfeat';
tmap_dir='stats';
nCope=12; %number of total copes - each one will have: tstat1=layout task, tstat2=item task
%here all of of them
reg_names={'lCA1',...
    'lCA23DG',...
    'rCA1',...
    'rCA23DG',...
    'rEnt',...
    'lEnt'...
    };

%this loop grabs for each subject the t-map, and greate a matrix file for
%each region

%that's for hipp:
for subj=1:numel(subjects)
            
    fprintf('creating data structure for subj %s\n',char(subjects(subj)));
    reg_data={};
    all_reg={};
    
    %% start with the hippocampus
    %upolad the sub-regions file:
    %R hipp: 1 is CA1, 3 is CA2/3/dentate, 2 is subiculum
    %L hipp: 4 is CA1, 6 is CA2/3/dentate, 5 is subiculum
    fileName=fullfile(subj_dir,char(subjects(subj)),roi_dir,'Hip15_oblique_al');
   
    if ~exist([fileName '.nii'],'file')
       %disp(['unzipping ' fileName])
       unix(['gunzip ' fileName '.nii.gz']);
    else
        disp([fileName ' doesn''t exist'])
    end
    
    %a stupid hack because something in the files were not that good - I
    %read them with spm_vol
    all_sub_reg=spm_read_vols(spm_vol(sprintf('%s.nii',fileName)));
    
    % zip up nifti file - we took the data, so can zip again
    unix(['gzip -f ' fileName '.nii']);
   
    %put them all in one structure, easier for later
    all_reg{1}=find(all_sub_reg==4);
    all_reg{2}=find(all_sub_reg==6);
    all_reg{4}=find(all_sub_reg==1);
    all_reg{5}=find(all_sub_reg==3);
    
    %% now the rest
    for roi=(numel(all_reg)+1):numel(reg_names)
        
        %upolad the regions file:
        fileName=fullfile(subj_dir,char(subjects(subj)),roi_dir,reg_names{roi});
        
        if ~exist([fileName '.nii'],'file')
            %disp(['unzipping ' fileName])
            unix(['gunzip ' fileName '.nii.gz']);
        else
            disp([fileName ' doesn''t exist'])
        end
        
        %a stupid hack because something in the files were not that good - I
        %read them with spm_vol
        curr_reg=spm_read_vols(spm_vol(sprintf('%s.nii',fileName)));
        
        % zip up nifti file - we took the data, so can zip again
        unix(['gzip -f ' fileName '.nii']);
        
        %put them all in one structure, easier for later
        all_reg{roi}=find(curr_reg);
    end
    
    %% for each contrast, retrieve the data
    for i=1:nCope
        %for each cope, retrieve the tmap in the layout task:
        cope_dir=['cope' num2str(i) '.feat'];
        fileName=fullfile(subj_dir,char(subjects(subj)),analysis_dir,cope_dir,tmap_dir,'tstat1');
        if ~exist([fileName '.nii'],'file')
           % disp(['unzipping ' fileName])
            unix(['gunzip ' fileName '.nii.gz']);
        end
        
        %upload the map
        data=spm_read_vols(spm_vol(sprintf('%s.nii',fileName)));
        % zip up nifti file - we took the data, so can zip again
        unix(['gzip -f ' fileName '.nii']);
        
        %get the relevant files for each region
        for reg=1:numel(reg_names)
            reg_data.(reg_names{reg}).lay(:,i)=data(all_reg{reg});
        end
        
        %for each cope, retrieve the tmap in the item task:
        cope_dir=['cope' num2str(i) '.feat'];
        fileName=fullfile(subj_dir,char(subjects(subj)),analysis_dir,cope_dir,tmap_dir,'tstat2');
        if ~exist([fileName '.nii'],'file')
           % disp(['unzipping ' fileName])
            unix(['gunzip ' fileName '.nii.gz']);
        end
        
        %upload the map
        data=spm_read_vols(spm_vol(sprintf('%s.nii',fileName)));
        % zip up nifti file - we took the data, so can zip again
        unix(['gzip -f ' fileName '.nii']);
        
        %get the relevant files for each region
        for reg=1:numel(reg_names)
            reg_data.(reg_names{reg}).item(:,i)=data(all_reg{reg});
        end
        
    end
    save(fullfile(subj_dir,char(subjects(subj)),'data','Univariate_reg_data.mat'),'reg_data');
   
end

end

        
   

