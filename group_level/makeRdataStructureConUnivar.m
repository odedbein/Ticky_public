function makeRdataStructureConUnivar(ResultsBetaSerConnectivityOnlyNum,ResultsAverageActivityOnlyNum,reg,engram,acc)
%created by Oded Bein
%assumes that connectivity analyses and activity analyses were
%run. Takes the data and write it to a txt file that will be used to R
%analysis

%this script assumes that the connectivity and activity analyses were on
%separated tasks, then colapsed - it takes the average.

%input:
%ResultsBetaSerConnectivityOnlyNum:     connectivity data structure
%ResultsAverageActivityOnlyNum:         mean unvariate activity data per subject
%reg:                                   regions (a string - should be reg1_reg2 - as appears in the connectivity data structure
%acc:                                   [0/1] - only accurate responses(1), or all responses(0)

%data structures to load:
%load('/Volumes/data/Bein/TickyReanalysis/results/data_in_paper/All_regs_connectivityAKcorrectFisherTransformed.mat')
%load('/Volumes/data/Bein/TickyReanalysis/results/data_in_paper/Univar_cue_image_single_trials.mat')
%-----------------------------------------------------------

if engram
    mydir='/data/Bein';
else
    mydir='/Volumes/data/Bein';
end

rmpath('/Volumes/data/Bein/fMRI_course/AnalysisScripts');

sing_trials=1;
if sing_trials %univariate was based on single trials
    sing_ttl='single_trials';
else
    sing_ttl='ProperGLM';
end

if acc
    acc_fname='onlyAcc';
else
    acc_fname='allItems';
end
proj_dir=fullfile(mydir,'TickyReanalysis');
results_fname=sprintf('%s_%sExcAD_AKcorrect_Univar_%s_withAccRT.txt',reg,acc_fname,sing_ttl);
results_dir=fullfile(proj_dir,'results','connectivity');
results_fname=fullfile(results_dir,results_fname);

Exc_AD=1;
Exc_AK=0;

%get the data:
reg1=reg(1:strfind(reg,'_')-1);
reg2=reg(strfind(reg,'_')+1:end);
if acc
    con_data=ResultsBetaSerConnectivityOnlyNum.(reg).Acc;
else
    con_data=ResultsBetaSerConnectivityOnlyNum.(reg).all_items;
end

if ~sing_trials
    if acc
        task_data1=ResultsAverageActivityOnlyNum.(reg1).Acc;
        task_data2=ResultsAverageActivityOnlyNum.(reg2).Acc;
    else
        task_data1=ResultsAverageActivityOnlyNum.(reg1).all_items;
        task_data2=ResultsAverageActivityOnlyNum.(reg2).all_items;
    end
    %construct the data per task from subj' contrasts:
    univariate_data1=[task_data1.lay(:,[1 10:12 9]) task_data1.item(:,[1 10:12 9])];
    univariate_data2=[task_data2.lay(:,[1 10:12 9]) task_data2.item(:,[1 10:12 9])];
else %activation if from single trials:
    if acc
        univariate_data1=ResultsAverageActivityOnlyNum.(reg1).image.Acc;
        univariate_data2=ResultsAverageActivityOnlyNum.(reg2).image.Acc;
    else
        univariate_data1=ResultsAverageActivityOnlyNum.(reg1).image.all_items;
        univariate_data2=ResultsAverageActivityOnlyNum.(reg2).image.all_items;
    end

end

%TAKE OUT THE CUE DATA:
con_data=con_data(:,3:end);
if Exc_AD
    con_data=con_data([1 3:end],:);
    univariate_data1=univariate_data1([1 3:end],:);
    univariate_data2=univariate_data2([1 3:end],:);
end


%set up some fixed columns:
n=size(con_data,1);
nCond=5;
S=repmat([1:n]',nCond,1);
linearModel=[];
for c=1:nCond
    linearModel=[linearModel;ones(n,1)*c];
end
all_or_none=[ones(n,1)*(-1);ones(n*(nCond-1),1)*1/(nCond-1)];


%accuracy model - 0 and 4 vs. 1-3:
acc_mdl=[ones(n,1)*0.5;ones(n*3,1)*-(1/3);ones(n,1)*0.5];
%RT model - 
RT_mdl=[ones(n*3,1)*(1/3);ones(n*2,1)*-(0.5)];



%start building the matrix to be written to a file:
M=[S linearModel all_or_none acc_mdl RT_mdl];

%get the connectivity data: average across tasks
lay=1:5;
item=6:10;
lay_data=con_data(:,lay);
item_data=con_data(:,item);
data=zeros([size(lay_data),2]);
data(:,:,1)=lay_data;
data(:,:,2)=item_data;
data=mean(data,3);
M=[M reshape(data,size(data,1)*size(data,2),1)];

%get the activity data for reg1: average across tasks
if sing_trials
    lay_data=univariate_data1(:,lay);
    item_data=univariate_data1(:,item);
else
    lay_data=univariate_data1.lay(:,[1 9:12]); %0-changes,4-change, contrasts of 1-3 changes
    item_data=univariate_data1.item(:,[1 9:12]);
end
data=zeros([size(lay_data),2]);
data(:,:,1)=lay_data;
data(:,:,2)=item_data;
data=mean(data,3);
M=[M reshape(data,size(data,1)*size(data,2),1)];

%get the activity data for reg2: average across tasks
if sing_trials
    lay_data=univariate_data2(:,lay);
    item_data=univariate_data2(:,item);
else
    lay_data=univariate_data2.lay(:,[1 9:12]); %0-changes,4-change, contrasts of 1-3 changes
    item_data=univariate_data2.item(:,[1 9:12]);
end

data=zeros([size(lay_data),2]);
data(:,:,1)=lay_data;
data(:,:,2)=item_data;
data=mean(data,3);
M=[M reshape(data,size(data,1)*size(data,2),1)];

%write it up:
fid = fopen(results_fname, 'w');
%set up the header:
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
            'subject','linearModel','all_or_none','acc_model','RT_model','BseriesCor',sprintf('%s_activity',reg1),sprintf('%s_activity',reg2));
dlmwrite(results_fname,M,'-append','delimiter','\t','precision','%.4f');        
end