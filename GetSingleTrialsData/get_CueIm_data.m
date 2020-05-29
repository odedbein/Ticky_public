function get_CueIm_data()

[~, hostname]=system('hostname');
if strcmp(hostname(1:6),'joanna')%the hostname command gives 1X7 char output, we only need the first 6.
    proj_dir='/Volumes/Oded/Bein/TickyReanalysis';
    % software and data paths
    warning('off','all')
    addpath(genpath('/Volumes/Oded/Bein/General_scripts'))
    rmpath('/Volumes/Oded/Bein/fMRI_course/AnalysisScripts');

else
    proj_dir='/Volumes/davachilab/Bein/TickyReanalysis';
    % software and data paths
    warning('off','all')
    addpath(genpath('/Volumes/davachilab/Bein/General_scripts'))
    rmpath('/Volumes/davachilab/Bein/fMRI_course/AnalysisScripts');
end

%THIS HAS ALL OF THEM: 
subjects={'AB'; 'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD'; 'JG'; 'JM'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};

subj_dir=fullfile(proj_dir,'SubData');
region_mat_dir='regions_mat'; %a directory to store the matfiles of each region
roi_dir='ROIs';

tmaps_dir='single_trial_Tmaps';
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
            };


%regions extracted:
% 1. all hipp (bilateral,r,l) 
% 2. sub-regions (bilateral,r,l) 
% 3. anterior, mid,post (bilateral,r,l) 
% 
%this loop grabs for each subject the t-map, and greate a matrix file for
%each region

%that's for hipp:
for subj=1:numel(subjects)
    num_trials=270;
    if strcmp(char(subjects(subj)),'LD') || strcmp(char(subjects(subj)),'AK')
        num_trials=216; %only had 8 scanning sessions
    end
            
    fprintf('creating data structure for subj %s\n',char(subjects(subj)));
    reg_data={};
    all_reg={};
    %upolad the sub-regions file:
    %R hipp: 1 is CA1, 3 is CA2/3/dentate, 2 is subiculum
    %L hipp: 4 is CA1, 6 is CA2/3/dentate, 5 is subiculum
    fileName=fullfile(subj_dir,char(subjects(subj)),roi_dir,'Hip15_oblique_al');
   
    if ~exist([fileName '.nii'],'file')
       disp(['unzipping ' fileName])
       unix(['gunzip ' fileName '.nii.gz']);
    end
    
    %a stupid hack because something in the files were not that good - I
    %read them with spm_vol
    all_sub_reg=spm_read_vols(spm_vol(sprintf('%s.nii',fileName)));
    
    % zip up nifti file - we took the data, so can zip again
    unix(['gzip -f ' fileName '.nii']);

    %upload the ant/mid/post
    fileName=fullfile(subj_dir,char(subjects(subj)),roi_dir,'Hip15_AntMidPost');
    %bilateral regions:
    %1-ant, 2-mid,3-post
    if ~exist([fileName '.nii'],'file')
       disp(['unzipping ' fileName])
       unix(['gunzip ' fileName '.nii.gz']);
    end
    
    bilateral_hipp_thirds=spm_read_vols(spm_vol(sprintf('%s.nii',fileName)));
    
    % zip up nifti file - we took the data, so can zip again
    unix(['gzip -f ' fileName '.nii']);
    
    %prepare the bilateral hipp ROI:
    bilateral_hipp=bilateral_hipp_thirds;
    bilateral_hipp(bilateral_hipp>0)=1;
    
    %prepare the left/right hipp:
    lhipp=zeros(size(bilateral_hipp));
    lhipp(all_sub_reg>=4)=1;
    rhipp=zeros(size(bilateral_hipp));
    rhipp(all_sub_reg<=3 & all_sub_reg>0)=1;
    
    %prepare the right/left ant/mid/post:
    %R: ant-1,mid-2,post-3
    %L: ant-4,mid-5,post-6
    hipp_thirds= bilateral_hipp_thirds;
    hipp_thirds((hipp_thirds==1) & (lhipp==1))=4;
    hipp_thirds((hipp_thirds==2) & (lhipp==1))=5;
    hipp_thirds((hipp_thirds==3) & (lhipp==1))=6;
    
    %put them all in one structure, easier for later
    all_reg{1}=find(bilateral_hipp);
    all_reg{2}=find(lhipp);
    all_reg{3}=find(rhipp);
    all_reg{4}=find(bilateral_hipp_thirds==1);
    all_reg{5}=find(bilateral_hipp_thirds==2);
    all_reg{6}=find(bilateral_hipp_thirds==3);
    all_reg{7}=find(hipp_thirds==4);
    all_reg{8}=find(hipp_thirds==5);
    all_reg{9}=find(hipp_thirds==6);
    all_reg{10}=find(hipp_thirds==1);
    all_reg{11}=find(hipp_thirds==2);
    all_reg{12}=find(hipp_thirds==3);
    all_reg{13}=find(all_sub_reg==1 | all_sub_reg==4);
    all_reg{14}=find(all_sub_reg==3 | all_sub_reg==6);
    all_reg{15}=find(all_sub_reg==2 | all_sub_reg==5);
    all_reg{16}=find(all_sub_reg==4);
    all_reg{17}=find(all_sub_reg==6);
    all_reg{18}=find(all_sub_reg==5);
    all_reg{19}=find(all_sub_reg==1);
    all_reg{20}=find(all_sub_reg==3);
    all_reg{21}=find(all_sub_reg==2);
    
    for i=1:num_trials
        %for each cue, retrieve the tmap:
        fileName=fullfile(subj_dir,char(subjects(subj)),tmaps_dir,sprintf('tstat_Cue%d',i));
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
            reg_data.(reg_names{reg}).cue(:,i)=data(all_reg{reg});
        end
        
        %for each image, retrieve the tmap:
        fileName=fullfile(subj_dir,char(subjects(subj)),tmaps_dir,sprintf('tstat_Im%d',i));
        if ~exist([fileName '.nii'],'file')
          %  disp(['unzipping ' fileName])
            unix(['gunzip ' fileName '.nii.gz']);
        end
        
        %upload the map
        data=spm_read_vols(spm_vol(sprintf('%s.nii',fileName)));
        % zip up nifti file - we took the data, so can zip again
        unix(['gzip -f ' fileName '.nii']);
        
        %get the relevant files for each region
        for reg=1:numel(reg_names)
            reg_data.(reg_names{reg}).Image(:,i)=data(all_reg{reg});
        end
        
    end
    save(fullfile(subj_dir,char(subjects(subj)),'data','reg_data.mat'),'reg_data');
   
end

%that's for MTL:
%first, create a box nii file that marks R and L - only need to run once -
%I created them
% subj=1;
% fileName=fullfile(subj_dir,char(subjects(subj)),roi_dir,reg_names{1});
% if ~exist([fileName '.nii'],'file')
%        disp(['unzipping ' fileName])
%        unix(['gunzip ' fileName '.nii.gz']);
% end
% all_sub_reg_header=spm_vol(sprintf('%s.nii',fileName));
% all_sub_reg=spm_read_vols(spm_vol(sprintf('%s.nii',fileName)));
% Right_side=[ones(size(all_sub_reg,1)/2,size(all_sub_reg,2),size(all_sub_reg,3)); zeros(size(all_sub_reg,1)/2,size(all_sub_reg,2),size(all_sub_reg,3))];
% all_sub_reg_header.fname=fullfile(proj_dir,'scripts','Rside.nii');
% spm_write_vol(all_sub_reg_header,Right_side);
% unix(['gzip -f ' all_sub_reg_header.fname]);
% 
% Left_side=[zeros(size(all_sub_reg,1)/2,size(all_sub_reg,2),size(all_sub_reg,3));ones(size(all_sub_reg,1)/2,size(all_sub_reg,2),size(all_sub_reg,3))];
% all_sub_reg_header.fname=fullfile(proj_dir,'scripts','Lside.nii');
% spm_write_vol(all_sub_reg_header,Left_side);
% unix(['gzip -f ' all_sub_reg_header.fname]);

%MTL:
reg_names={'Para',...%1
    'Peri',...%2
    'Ent',...%3
    'rPara',...%4
    'lPara',...%5
    'rPeri',...%6
    'lPeri',...%7
    'rEnt',...%8
    'lEnt',...%9
    };

for subj=3:numel(subjects)
    fprintf('creating data structure for subj %s\n',char(subjects(subj)));
    reg_data={};
    all_reg={};
    
    num_trials=270;
    if strcmp(char(subjects(subj)),'LD') || strcmp(char(subjects(subj)),'AK')
        num_trials=216; %only had 8 scanning sessions
    end
        
    %prepare the structure with all regions:
    for reg=1:numel(reg_names)
        
        %upload the region's file
        fileName=fullfile(subj_dir,char(subjects(subj)),roi_dir,reg_names{reg});
        
        if ~exist([fileName '.nii'],'file')
            disp(['unzipping ' fileName])
            unix(['gunzip ' fileName '.nii.gz']);
        end
        
        %a stupid hack because something in the files were not that good - I
        %read them with spm_vol
        curr_reg=spm_read_vols(spm_vol(sprintf('%s.nii',fileName)));
        
        % zip up nifti file - we took the data, so can zip again
        unix(['gzip -f ' fileName '.nii']);
        
        %put them all in one structure, easier for later
        all_reg{reg}=find(curr_reg);
    end
    
    for i=1:num_trials
        %for each cue, retrieve the tmap:
        fileName=fullfile(subj_dir,char(subjects(subj)),tmaps_dir,sprintf('tstat_Cue%d',i));
        if ~exist([fileName '.nii'],'file')
           % disp(['unzipping ' fileName])
            unix(['gunzip ' fileName '.nii.gz']);
        end
        
        %upload the map
        data=spm_read_vols(spm_vol(sprintf('%s.nii',fileName)));
        % zip up nifti file - we took the data, so can zip again
        unix(['gzip -f ' fileName '.nii']);
        
        %get the relevant data for each region
        for reg=1:numel(reg_names)
            reg_data.(reg_names{reg}).cue(:,i)=data(all_reg{reg});
        end
        
        %for each image, retrieve the tmap:
        fileName=fullfile(subj_dir,char(subjects(subj)),tmaps_dir,sprintf('tstat_Im%d',i));
        if ~exist([fileName '.nii'],'file')
          %  disp(['unzipping ' fileName])
            unix(['gunzip ' fileName '.nii.gz']);
        end
        
        %upload the map
        data=spm_read_vols(spm_vol(sprintf('%s.nii',fileName)));
        % zip up nifti file - we took the data, so can zip again
        unix(['gzip -f ' fileName '.nii']);
        
        %get the relevant files for each region
        for reg=1:numel(reg_names)
            reg_data.(reg_names{reg}).Image(:,i)=data(all_reg{reg});
        end
        
    end
    save(fullfile(subj_dir,char(subjects(subj)),'data','reg_dataMTL.mat'),'reg_data');
   
end
    
        
   

