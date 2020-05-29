function step4_create_regressors_files()
%creates the regressor files for each model, for each subject

proj_dir='/Volumes/Oded/Bein/TickyReanalysis';
subjects={'AB'; 'AD'; 'AK'; 'AR'; 'AT'; 'BW'; 'CR'; 'DH'; 'DM'; 'EB'; 'JA'; 'JD'; 'JG'; 'JM'; 'JR'; 'JW'; 'JR'; 'JW'; 'KZ'; 'LD'; 'SB'; 'YE'};
%subjects={'JA';'JG';'JW';'YE'};
subjects={'LD'};
imported=1;
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
    reg_dir=fullfile(subj_dir,'regressors');
    if ~exist(reg_dir)
        mkdir (reg_dir);
    end
    behav_filename=fullfile(subj_dir,sprintf('output_subject_%s.txt',subjects{subj}));
    
    if ~imported
        subj_behavior=textscan(fopen(behav_filename),'%d %d %d %.1f %s %d %s %.1f %s %s %.3f %.3f %d');
        if strcmp(subjects{subj},'LD') %need to take session 3-10 and correct the timing of the 9th repetition.
            timing=subj_behavior{4}(2*trials_per_sess+1:end);
            type=subj_behavior{5}(2*trials_per_sess+1:end);
            num_items=length(type);
            %change the timing of the items in the 9th session (7th after removing 2), to mach
            %with Katherine's regressors (5 secs delay).
            timing(6*trials_per_sess+1:7*trials_per_sess)=timing(6*trials_per_sess+1:7*trials_per_sess)+5;
        else
            timing=subj_behavior{4};
            type=subj_behavior{5};
            num_items=size(subj_behavior{1},1);
        end
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %if imported - make sure you import the relevant data manually
        subj_behavior=eval(sprintf('outputsubject%s',subjects{subj}));
        num_items=size(subj_behavior,1);
        timing=[];
        for i=1:num_items
            timing=[timing; subj_behavior{i,4}];
        end
        type={};
        for i=1:num_items
            type=[type; subj_behavior{i,5}];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    for sess=1:num_items/trials_per_sess;
        if strcmp(subjects{subj},'AR') %need to correct timing
            timing(((sess-1)*trials_per_sess+1):((sess-1)*trials_per_sess+trials_per_sess))=timing(((sess-1)*trials_per_sess+1):((sess-1)*trials_per_sess+trials_per_sess))-(300*(sess-1));
        end
        %create regressors files for each trial, in each session
        for i=1:trials_per_sess
            trial_num=(sess-1)*trials_per_sess+i;
            %% cue model:
            %cue regressor
            reg_file=sprintf('Cue%dModel_cue.txt',trial_num);
            fid = fopen(fullfile(reg_dir,reg_file), 'w');
            fprintf(fid,'%.1f\t%.1f\t%d',timing(trial_num),cue_duration,1);
            fclose(fid);
            
            %other cues:
            reg_file=sprintf('Cue%dModel_AllOtherCues.txt',trial_num);
            fid = fopen(fullfile(reg_dir,reg_file), 'w');
            for c=1:trials_per_sess
                cue=(sess-1)*trials_per_sess+c;
                if cue~=trial_num %all other cues but the target cue
                   fprintf(fid,'%.1f\t%.1f\t%d\n',timing(cue),cue_duration,1);
                end
            end
            fclose(fid);
            
            %images - divided to different conditions - we know that there
            %were activity differences between conditions, so should model
            %differently
            curr_im=type((sess-1)*trials_per_sess+1:(sess-1)*trials_per_sess+trials_per_sess);
            for tt=1:numel(trial_types)
                reg_file=sprintf('Cue%dModel_ImType%d.txt',trial_num,tt);
                fid = fopen(fullfile(reg_dir,reg_file), 'w');
                curr_trials=find(strcmp(curr_im,trial_types{tt}))+(sess-1)*trials_per_sess;
                for curr_im_type=1:length(curr_trials)
                    fprintf(fid,'%.1f\t%.1f\t%d\n',timing(curr_trials(curr_im_type))+im_delay,im_duration,1);
                end
                fclose(fid);
            end
            
            %% image model:
            %image regressor
            reg_file=sprintf('Im%dModel_image.txt',trial_num);
            fid = fopen(fullfile(reg_dir,reg_file), 'w');
            fprintf(fid,'%.1f\t%.1f\t%d',timing(trial_num)+im_delay,im_duration,1);
            fclose(fid);
            
            %all_cues:
            reg_file=sprintf('Im%dModel_AllCues.txt',trial_num);
            fid = fopen(fullfile(reg_dir,reg_file), 'w');
            for c=1:trials_per_sess
                cue=(sess-1)*trials_per_sess+c;
                fprintf(fid,'%.1f\t%.1f\t%d\n',timing(cue),cue_duration,1);
            end
            fclose(fid);
            
            %images - divided to different conditions - we know that there
            %were activity differences between conditions, so should model
            %differently
            curr_im=type((sess-1)*trials_per_sess+1:(sess-1)*trials_per_sess+trials_per_sess);
            for tt=1:numel(trial_types)
                reg_file=sprintf('Im%dModel_ImType%d.txt',trial_num,tt);
                fid = fopen(fullfile(reg_dir,reg_file), 'w');
                curr_trials=find(strcmp(curr_im,trial_types{tt}))+(sess-1)*trials_per_sess;
                for curr_im_type=1:length(curr_trials)
                    if curr_trials(curr_im_type)~=trial_num %if it's not the current image
                        fprintf(fid,'%.1f\t%.1f\t%d\n',timing(curr_trials(curr_im_type))+im_delay,im_duration,1);
                    end
                end
                fclose(fid);
            end
            
        end
    end %ends all sess for current subj
end %ends all subjs

end
  