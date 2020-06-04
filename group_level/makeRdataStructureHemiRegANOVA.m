function makeRdataStructureHemiRegANOVA(ResultsBetaSerConnectivityOnlyNum,acc)
%created by Oded Bein
%assumes that connectivity analyses and single-trial activity analyses were
%run. Takes the data and write it to a txt file that will be used to R
%analysis

%this script assumes that the connectivity and activity analyses were on
%separated tasks, then colapsed - it takes the average.

%input:
%ResultsBetaSerConnectivityOnlyNum:     connectivity data structure
%ResultsAverageActivityOnlyNum:         mean unvariate activity data per subject
%reg:                                   regions (a string - should be reg1_reg2 - as appears in the connectivity data structure
%acc:                                   [0/1] - only accurate responses(1), or all responses(0)

%-----------------------------------------------------------
proj_dir='/Volumes/data/Bein/TickyReanalysis';

Exc_AD=1;
%Exc_AK=0;
if acc
    acc_fname='onlyAcc';
else
    acc_fname='allItems';
end

data_type='ProperGLM_AvMatchMis_High2Thirds';

fname=sprintf('LeftRightByEntCA3_%s_%s_ExcAD_AKcorrect_noUnivar.txt',data_type,acc_fname);
results_dir=fullfile(proj_dir,'results','connectivity');

%header={'subject','linearModel','all_or_none','BseriesCor',sprintf('%s_activity',reg1),sprintf('%s_activity',reg2)};
%regions
reg_ttl={'reg1lh','reg2lh','reg1rh','reg2rh'};
%get the data:
regs.reg1lh='lCA1_lEnt';
regs.reg2lh='lCA1_lCA23DG';
regs.reg1rh='rCA1_rEnt';
regs.reg2rh='rCA1_rCA23DG';

if acc
    for reg=1:numel(reg_ttl)
        con_data.(reg_ttl{reg})=ResultsBetaSerConnectivityOnlyNum.(regs.(reg_ttl{reg})).Acc;
    end
else
    for reg=1:numel(reg_ttl)
        con_data.(reg_ttl{reg})=ResultsBetaSerConnectivityOnlyNum.(regs.(reg_ttl{reg})).all_items;
    end
end

if Exc_AD
    for reg=1:numel(reg_ttl)
        con_data.(reg_ttl{reg})=con_data.(reg_ttl{reg})([1 3:end],:);
    end
end


%set up some fixed columns:
n=size(con_data.(reg_ttl{reg}),1);
nCond=5;
S=repmat([1:n]',nCond,1);
linearModel=[];
for c=1:nCond
    linearModel=[linearModel;ones(n,1)*c];
end
all_or_none=[ones(n,1)*(-1);ones(n*(nCond-1),1)*1/(nCond-1)];

%start building the matrix to be written to a file:
M=[S linearModel all_or_none];
M=repmat(M,4,1);

%add the hemisphere/region regressors
reg_regress=repmat([ones(length(linearModel),1); ones(length(linearModel),1)*2],2,1);
hem_regress=[ones(length(linearModel)*2,1); ones(length(linearModel)*2,1)*2];
M=[M hem_regress reg_regress];

%get the connectivity data: average across tasks
Mcon=[];
lay=3:7;
item=8:12;

for reg=1:numel(reg_ttl)
    curr_con_data=con_data.(reg_ttl{reg});
    lay_data=curr_con_data(:,lay);
    item_data=curr_con_data(:,item);
    data=zeros([size(lay_data),2]);
    data(:,:,1)=lay_data;
    data(:,:,2)=item_data;
    data=mean(data,3);
    Mcon=[Mcon; reshape(data,size(data,1)*size(data,2),1)];
end

M=[M Mcon];

%write it up:
filename=fullfile(results_dir,fname);
fid = fopen(filename, 'w');
%set up the header:
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n', ...
            'subject','linearModel','all_or_none','hemisphere','reg','BseriesCor');
dlmwrite(filename,M,'-append','delimiter','\t','precision','%.4f');        
end