
behav_filename=fullfile(subj_dir,subjects{subj},sprintf('output_subject_%s1.txt',subjects{subj}));
fid=fopen(behav_filename,'w');
data=outputsubjectAR;
%fprintf(fid,'%d %d %d %.1f %s %d %s %.1f %s %s %.3f %.3f %d\n',

%for AR:
for r=1:270
fprintf(fid,'%d\t%d\t%d\t%.1f\t%s\t%s\t%d\t%s\t%.1f\t%s\n',...
        data{r,1},data{r,2},data{r,3},data{r,4},data{r,5}(2:end),data{r,6}(2:end),data{r,7},data{r,8}(2:end),data{r,9},data{r,10}(2:end));
end
fclose(fid);

%for AR to read:subj_behavior=textscan(fopen(behav_filename),'%d %d %d %.1f %s %s %d %s %.1f %s');

%for the rest
subjects={'JA';'JG';'JW';'YE'};
for subj=1:numel(subjects)
    subj_dir=fullfile(proj_dir,'SubData',subjects{subj});
    behav_filename=fullfile(subj_dir,sprintf('output_subject_%s1.txt',subjects{subj}));
    fid=fopen(behav_filename,'w');
    data=eval(sprintf('outputsubject%s',subjects{subj}));
    for r=1:270
    fprintf(fid,'%d\t%d\t%d\t%.1f\t%s\t%d\t%s\t%.1f\t%d\t%s\t%.3f\t%.3f\t%d\n',...
            data{r,1},data{r,2},data{r,3},data{r,4},data{r,5},data{r,6},data{r,7},data{r,8},data{r,9},data{r,10},data{r,11},data{r,12},data{r,13});
    end
    fclose(fid);
end

%for JW
subjects={'JW'};
for subj=1:numel(subjects)
    subj_dir=fullfile(proj_dir,'SubData',subjects{subj});
    behav_filename=fullfile(subj_dir,sprintf('output_subject_%s1.txt',subjects{subj}));
    fid=fopen(behav_filename,'w');
    data=eval(sprintf('outputsubject%s',subjects{subj}));
    for r=1:270
    fprintf(fid,'%d\t%d\t%d\t%.1f\t%s\t%d\t%s\t%.1f\t%s\t%s\t%.3f\t%.3f\t%d\n',...
            data{r,1},data{r,2},data{r,3},data{r,4},data{r,5}(2:end),data{r,6},data{r,7}(2:end),data{r,8},data{r,9}(2:end),data{r,10}(2:end),data{r,11},data{r,12},data{r,13});
    end
    fclose(fid);
end