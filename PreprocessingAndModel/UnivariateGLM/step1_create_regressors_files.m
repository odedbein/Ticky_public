function step1_create_regressors_files()
%creates the regressor files for each model, for each subject

proj_dir='/Volumes/Oded/Bein/TickyReanalysis';
subjects={'AB'; 'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD'; 'JG'; 'JM'; 'JR'; 'JW'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};
%subjects={'JA';'JG';'JW';'YE'};
%LD - only 216 trials
%AK - 189 tials - eighth run was running the same block again, excluded
% subjects={'LD'};
% subjects={'AK'};

imported=0;
trials_per_sess=27;
cue_duration=1.5;
im_duration=4;
im_delay=cue_duration+1;%images appeared 1 sec after offset of cue
trial_types={...
    'R0I0';
    'R0I1';
    'R0I2';
    'R1I0';
    'R1I1';
    'R1I2';
    'R2I0';
    'R2I1';
    'R2I2'...
    };
for subj=1:numel(subjects)
    fprintf('creating regressors for subject %s\n',subjects{subj});
    subj_dir=fullfile(proj_dir,'SubData',subjects{subj});
    reg_dir=fullfile(subj_dir,'regressorsUnivariate');
    reg_cor_dir=(fullfile(reg_dir,'CorIncor'));
    reg_all_dir=(fullfile(reg_dir,'All'));
    if exist(reg_dir)
        rmdir(reg_dir,'s');
    end
    if ~exist(reg_cor_dir)
        mkdir (reg_cor_dir);
    end
    if ~exist(reg_all_dir)
        mkdir (reg_all_dir);
    end
    behav_filename=fullfile(subj_dir,sprintf('output_subject_%s.txt',subjects{subj}));
    if strcmp(char(subjects(subj)),'AR')
        subj_behavior=textscan(fopen(behav_filename),'%d %d %d %.1f %s %s %d %s %.1f %s');
    else
        subj_behavior=textscan(fopen(behav_filename),'%d %d %d %.1f %s %d %s %.1f %s %s %.3f %.3f %d');
    end
    
    if strcmp(char(subjects(subj)),'LD')
        num_trials=216;
    elseif strcmp(char(subjects(subj)),'AK')
        num_trials=189;
    else
        num_trials=270;
    end
    
    if ~imported
        %subj_behavior=textscan(fopen(behav_filename),'%d %d %d %.1f %s %d %s %.1f %s %s %.3f %.3f %d');
        if strcmp(subjects{subj},'LD') %need to take session 3-10 and correct the timing of the 9th repetition.
            timing=subj_behavior{4}(2*trials_per_sess+1:end);
            type=subj_behavior{5}(2*trials_per_sess+1:end);
            acc=subj_behavior{3}(((270/10*2)+1):end);
            %change the timing of the items in the 9th session (7th after removing 2), to mach
            %with Katherine's regressors (5 secs delay).
            timing(6*trials_per_sess+1:7*trials_per_sess)=timing(6*trials_per_sess+1:7*trials_per_sess)+5;
        else
            timing=subj_behavior{4}(1:num_trials);
            type=subj_behavior{5}(1:num_trials);
            acc=subj_behavior{3}(1:num_trials);
        end
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %if imported - make sure you import the relevant data manually
%         subj_behavior=eval(sprintf('outputsubject%s',subjects{subj}));
%         num_trials=size(subj_behavior,1);
%         timing=[];
%         for i=1:num_trials
%             timing=[timing; subj_behavior{i,4}];
%         end
%         type={};
%         for i=1:num_trials
%             type=[type; subj_behavior{i,5}];
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    for sess=1:num_trials/trials_per_sess;
        if strcmp(subjects{subj},'AR') %need to correct timing
            timing(((sess-1)*trials_per_sess+1):((sess-1)*trials_per_sess+trials_per_sess))=timing(((sess-1)*trials_per_sess+1):((sess-1)*trials_per_sess+trials_per_sess))-(300*(sess-1));
        end
        %% spearate to corr and incorr
        %create one regressor for correct cues:
        reg_file=sprintf('run%d_CorrectCues.txt',sess);
        fid = fopen(fullfile(reg_cor_dir,reg_file), 'w');
        for c=1:trials_per_sess
            cue=(sess-1)*trials_per_sess+c;
            if acc(cue)==1
                fprintf(fid,'%.1f\t%.1f\t%d\n',timing(cue),cue_duration,1);
            end
        end
        fclose(fid);
        
        %create one regressor for incorrect cues:
        reg_file=sprintf('run%d_IncorrectCues.txt',sess);
        fid = fopen(fullfile(reg_cor_dir,reg_file), 'w');
        ii=0;
        for c=1:trials_per_sess
            cue=(sess-1)*trials_per_sess+c;
            if acc(cue)==0
                fprintf(fid,'%.1f\t%.1f\t%d\n',timing(cue),cue_duration,1);
                ii=ii+1;
            end
        end
        if ~ii %this is an empty file
            fprintf(fid,'%.1f\t%.1f\t%d\n',0,0,0);
        end
        fclose(fid);
        
        %images - accurate are divided to different conditions:
        curr_im=type((sess-1)*trials_per_sess+1:(sess-1)*trials_per_sess+trials_per_sess); %get all trial types in the current session
        
        for tt=1:numel(trial_types)
            ii=0;
            ImType=trial_types{tt};
            reg_file=sprintf('run%d_ImType_%s.txt',sess,ImType);
            fid = fopen(fullfile(reg_cor_dir,reg_file), 'w');
            curr_trials=find(strcmp(curr_im,trial_types{tt}))+(sess-1)*trials_per_sess;
            for curr_im_type=1:length(curr_trials)
                if acc(curr_trials(curr_im_type))==1
                    fprintf(fid,'%.1f\t%.1f\t%d\n',timing(curr_trials(curr_im_type))+im_delay,im_duration,1);
                    ii=ii+1;
                end
            end
            if ~ii %this is an empty file
                fprintf(fid,'%.1f\t%.1f\t%d\n',0,0,0);
            end
            fclose(fid);
        end
        
        %images - all inaccurate in a separated regressor, together:
        reg_file=sprintf('run%d_IcorrectIm.txt',sess);
        fid = fopen(fullfile(reg_cor_dir,reg_file), 'w');
        curr_trials=(sess-1)*trials_per_sess+1:(sess-1)*trials_per_sess+trials_per_sess;
        ii=0;
        for curr_im_type=1:length(curr_trials)
            if acc(curr_trials(curr_im_type))==0
                fprintf(fid,'%.1f\t%.1f\t%d\n',timing(curr_trials(curr_im_type))+im_delay,im_duration,1);
                ii=ii+1;
            end
        end
        if ~ii %this is an empty file
            fprintf(fid,'%.1f\t%.1f\t%d\n',0,0,0);
        end
        fclose(fid);
        
        %% all trials - collapse on cor/incor
        %create one regressor for correct cues:
        reg_file=sprintf('run%d_Cues.txt',sess);
        fid = fopen(fullfile(reg_all_dir,reg_file), 'w');
        for c=1:trials_per_sess
            cue=(sess-1)*trials_per_sess+c;
            fprintf(fid,'%.1f\t%.1f\t%d\n',timing(cue),cue_duration,1);
            
        end
        fclose(fid);
        
        %images - accurate are divided to different conditions:
        curr_im=type((sess-1)*trials_per_sess+1:(sess-1)*trials_per_sess+trials_per_sess); %get all trial types in the current session
        
        for tt=1:numel(trial_types)
            ImType=trial_types{tt};
            reg_file=sprintf('run%d_ImType_%s.txt',sess,ImType);
            fid = fopen(fullfile(reg_all_dir,reg_file), 'w');
            curr_trials=find(strcmp(curr_im,trial_types{tt}))+(sess-1)*trials_per_sess;
            for curr_im_type=1:length(curr_trials)
                fprintf(fid,'%.1f\t%.1f\t%d\n',timing(curr_trials(curr_im_type))+im_delay,im_duration,1);
            end
            fclose(fid);
        end
        
    end %ends all sess for current subj
end %ends all subjs

end
  