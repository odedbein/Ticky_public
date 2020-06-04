function [data, curr_data]=plot_Behavior_SeparateTask_NumChanges(ResultsOnlyNum,acc,Rstruct,closePrev)
%% I USED THIS FILE INITIALLY,BUT REPORTED STATS ARE NOW BASED ON R SCRIPT. THIS FILE WAS USED TO PRODUCE THE FILES FOR THE R ANALYSIS
%subj AD (subj #2) had only 12 voxels in the left Ent, for these analyses,
%may want to exclude him
contrasts=0;
Exc_AD=0; %i removed from the original analysis.
Exc_AK=0;
if closePrev
    close all
end

%% folders etc, for the Rfile:
engram=0;
if engram
    mydir='/data/Bein';
else
    mydir='/Volumes/data/Bein';
end
proj_dir=fullfile(mydir,'TickyReanalysis');
results_fname=sprintf('Behavior_ExcAD_AKcorrect.txt');
results_dir=fullfile(proj_dir,'results','behavior');
results_fname=fullfile(results_dir,results_fname);
%% analyze accuracy data:

curr_data=ResultsOnlyNum.Acc;

if Exc_AD && Exc_AK
    curr_data=curr_data([1 4:end],:);
    display ('AD and AK are excluded from analysis')
elseif Exc_AD && ~Exc_AK
    curr_data=curr_data([1 3:end],:);
    display ('AD is excluded from analysis')
elseif ~Exc_AD && Exc_AK
    curr_data=curr_data([1:2 4:end],:);
    display ('AK is excluded from analysis')
end

%% make behavioral data R sturcture:
if Rstruct
    %set up some fixed columns:
    n=size(curr_data,1);
    nCond=5;
    Scol=repmat([1:n]',nCond,1);
    linearModel=[];
    for c=1:nCond
        linearModel=[linearModel;ones(n,1)*c];
    end
end

averageCon=mean(curr_data);

lay=1:5;
item=6:10;
%plot 0-4 changes:
%compute within-subject SEM for both tasks:
n=size(curr_data,1);
subjs_av=mean(curr_data(:,lay),2);
withinSEM=curr_data(:,lay)-repmat(subjs_av,[1,5]);
SEM=abs(std(withinSEM)/sqrt(n));
subjs_av=mean(curr_data(:,item),2);
withinSEM=curr_data(:,item)-repmat(subjs_av,[1,5]);
SEM=[SEM abs(std(withinSEM)/sqrt(n))];
%set up the graph:
greyscaleColors=repmat(linspace(1,0,5)',2,3);
%Xticks=ResultsBetaSerConnectivity.hipp_Para.Acc(1,3:7);
Xticks={'lay: 0','lay: 1','lay: 2','lay: 3','lay: 4','item: 0','item: 1','item: 2','item: 3','item: 4'};
bar(1:10,averageCon);
set(gca,'XTickLabel',Xticks)
ylabel('Accuracy rates','Fontsize',16)
xlabel('Total number of Changes','Fontsize',14);

hold on
for i=1:10
    bar(i,averageCon(i),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i),SEM(i),'k.');
end
hold off
xlim([0,11]);
%stats:

if contrasts
    %%fitting a linear contrast:
    display(sprintf('linear contrast layout task: \n'));
    lin=[-2 -1 0 1 2];
    subj_con=curr_data(:,lay).*repmat(lin,size(curr_data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    display(sprintf('linear contrast item task: \n'));
    lin=[-2 -1 0 1 2];
    subj_con=curr_data(:,item).*repmat(lin,size(curr_data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    %%fitting a quadratic contrast:
    display(sprintf('quadratic contrast layout task: \n'));
    con=[2 -1 -2 1 2];
    con=con-mean(con);
    subj_con=curr_data(:,lay).*repmat(con,size(curr_data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    display(sprintf('quadratic contrast item task: \n'));
    con=[2 -1 -2 1 2];
    con=con-mean(con);
    subj_con=curr_data(:,item).*repmat(con,size(curr_data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    %%fitting an all-or-none contrast:
    display(sprintf('all-or-none contrast layout task: \n'));
    lin=[-1 0.25 0.25 0.25 0.25];
    subj_con=curr_data(:,lay).*repmat(lin,size(curr_data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    display(sprintf('all-or-none contrast item task: \n'));
    lin=[-1 0.25 0.25 0.25 0.25];
    subj_con=curr_data(:,item).*repmat(lin,size(curr_data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    display(sprintf('0-4 ttest layout: \n'));
    [h,p,ci,stats] = ttest(curr_data(:,lay(1)),curr_data(:,lay(5)));
    p
    
    display(sprintf('0-4 ttest item: \n'));
    [h,p,ci,stats] = ttest(curr_data(:,item(1)),curr_data(:,item(5)));
    p
    
end
%% ANOVA lay_item
%run the ANOVA:
lay_data=curr_data(:,lay);
item_data=curr_data(:,item);
n=size(lay_data,1);
Y=reshape(lay_data,size(lay_data,1)*size(lay_data,2),1);
Y=[Y; reshape(item_data,size(item_data,1)*size(item_data,2),1)];
S=repmat([1:n]',10,1);
F1=[ones(n*5,1);ones(n*5,1)*2];%task
F2=repmat([ones(n,1);ones(n,1)*2;ones(n,1)*3;ones(n,1)*4;ones(n,1)*5],2,1);%number of changes
display(sprintf('ANOVA Accuracy rates layout vs. item task\n'));
stats = rm_anova2(Y,S,F1,F2,{'task:(1)layout/(2)item','#changes'})
My_alpha=.05;
[P,MSEAB,F3] = RMAOV2_mod([Y,F1,F2,S],My_alpha,1);

data=zeros([size(lay_data),2]);
data(:,:,1)=lay_data;
data(:,:,2)=item_data;
data=mean(data,3);

if Rstruct
    M=[Scol linearModel reshape(data,[size(data,1)*size(data,2),1])];
end

%do the analysis average tasks:
if contrasts
    %%fitting a linear contrast:
    
    display(sprintf('linear contrast average tasks: \n'));
    lin=[-2 -1 0 1 2];
    subj_con=data.*repmat(lin,size(data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    display(sprintf('all-or-none contrast average tasks: \n'));
    lin=[-1 0.25 0.25 0.25 0.25];
    subj_con=data.*repmat(lin,size(data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    %%fitting a quadratic contrast:
    display(sprintf('quadratic contrast average tasks: \n'));
    con=[2 -1 -2 1 2];
    con=con-mean(con);
    subj_con=data.*repmat(con,size(data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    display(sprintf('0-4 ttest: \n'));
    [h,p,ci,stats] = ttest(data(:,1),data(:,5));
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
end

%plot it:
%compute within-subject SEM:
n=size(data,1);
subjs_av=mean(data,2);
withinSEM=data-repmat(subjs_av,[1,5]);
SEM=abs(std(withinSEM)/sqrt(n));
%set up the graph:
greyscaleColors=repmat(linspace(1,0,5)',1,3);
averageCon=mean(data);
Xticks={'0','1','2','3','4'};
figure;
bar(1:5,averageCon);
set(gca,'XTickLabel',Xticks)
ylabel('Accuracy','Fontsize',16)
xlabel('Total number of Changes','Fontsize',16);


title('Accuracy (av tasks): only accurate resp','Fontsize',20);

hold on
for i=1:5
    bar(i,averageCon(i),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i),SEM(i),'k.');
end
hold off


% ttests:
disp('ACCURACY DATA:');
for i=1:4
    for ii=(i+1):5
        fprintf('%d-%d ttest: \n', (i-1), (ii-1))
        [h,p,ci,stats] = ttest(data(:,i),data(:,ii));
        p
    end
end


%% analyze RTs data in seconds:

if acc
    curr_data=ResultsOnlyNum.RT.Acc;
else
    curr_data=ResultsOnlyNum.RT.all_items;
end

if Exc_AD && Exc_AK
    curr_data=curr_data([1 4:end],:);
    display ('AD and AK are excluded from analysis')
elseif Exc_AD && ~Exc_AK
    curr_data=curr_data([1 3:end],:);
    display ('AD is excluded from analysis')
elseif ~Exc_AD && Exc_AK
    curr_data=curr_data([1:2 4:end],:);
    display ('AK is excluded from analysis')
end

averageCon=mean(curr_data);

lay=1:5;
item=6:10;
%plot 0-4 changes:
%compute within-subject SEM for both tasks:
n=size(curr_data,1);
subjs_av=mean(curr_data(:,lay),2);
withinSEM=curr_data(:,lay)-repmat(subjs_av,[1,5]);
SEM=abs(std(withinSEM)/sqrt(n));
subjs_av=mean(curr_data(:,item),2);
withinSEM=curr_data(:,item)-repmat(subjs_av,[1,5]);
SEM=[SEM abs(std(withinSEM)/sqrt(n))];
%set up the graph:
greyscaleColors=repmat(linspace(1,0,5)',2,3);
%Xticks=ResultsBetaSerConnectivity.hipp_Para.Acc(1,3:7);
Xticks={'lay: 0','lay: 1','lay: 2','lay: 3','lay: 4','item: 0','item: 1','item: 2','item: 3','item: 4'};
figure;
bar(1:10,averageCon);
set(gca,'XTickLabel',Xticks)
ylabel('RTs(s)','Fontsize',16)
xlabel('Total number of Changes','Fontsize',14);

hold on
for i=1:10
    bar(i,averageCon(i),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i),SEM(i),'k.');
end
hold off
xlim([0,11]);
%stats:
if contrasts
    %%fitting a linear contrast:
    display(sprintf('linear contrast layout task: \n'));
    lin=[-2 -1 0 1 2];
    subj_con=curr_data(:,lay).*repmat(lin,size(curr_data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    display(sprintf('linear contrast item task: \n'));
    lin=[-2 -1 0 1 2];
    subj_con=curr_data(:,item).*repmat(lin,size(curr_data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    %%fitting a quadratic contrast:
    display(sprintf('quadratic contrast layout task: \n'));
    con=[2 -1 -2 1 2];
    con=con-mean(con);
    subj_con=curr_data(:,lay).*repmat(con,size(curr_data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    display(sprintf('quadratic contrast item task: \n'));
    con=[2 -1 -2 1 2];
    con=con-mean(con);
    subj_con=curr_data(:,item).*repmat(con,size(curr_data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    %%fitting an all-or-none contrast:
    display(sprintf('all-or-none contrast layout task: \n'));
    lin=[-1 0.25 0.25 0.25 0.25];
    subj_con=curr_data(:,lay).*repmat(lin,size(curr_data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    display(sprintf('all-or-none contrast item task: \n'));
    lin=[-1 0.25 0.25 0.25 0.25];
    subj_con=curr_data(:,item).*repmat(lin,size(curr_data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    display(sprintf('0-4 ttest layout: \n'));
    [h,p,ci,stats] = ttest(curr_data(:,lay(1)),curr_data(:,lay(5)));
    p
    
    display(sprintf('0-4 ttest item: \n'));
    [h,p,ci,stats] = ttest(curr_data(:,item(1)),curr_data(:,item(5)));
    p
end

%% ANOVA lay_item
%run the ANOVA:
lay_data=curr_data(:,lay);
item_data=curr_data(:,item);
n=size(lay_data,1);
Y=reshape(lay_data,size(lay_data,1)*size(lay_data,2),1);
Y=[Y; reshape(item_data,size(item_data,1)*size(item_data,2),1)];
S=repmat([1:n]',10,1);
F1=[ones(n*5,1);ones(n*5,1)*2];%task
F2=repmat([ones(n,1);ones(n,1)*2;ones(n,1)*3;ones(n,1)*4;ones(n,1)*5],2,1);%number of changes
display(sprintf('ANOVA RTs layout vs. item task\n'));
stats = rm_anova2(Y,S,F1,F2,{'task:(1)layout/(2)item','#changes'})
My_alpha=.05;
[P,MSEAB,F3] = RMAOV2_mod([Y,F1,F2,S],My_alpha,1);

data=zeros([size(lay_data),2]);
data(:,:,1)=lay_data;
data(:,:,2)=item_data;
data=mean(data,3);

%do the analysis average tasks:
if contrasts
    %%fitting a linear contrast:
    display(sprintf('linear contrast average tasks: \n'));
    lin=[-2 -1 0 1 2];
    subj_con=data.*repmat(lin,size(data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    display(sprintf('all-or-none contrast average tasks: \n'));
    lin=[-1 0.25 0.25 0.25 0.25];
    subj_con=data.*repmat(lin,size(data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    %%fitting a quadratic contrast:
    display(sprintf('quadratic contrast average tasks: \n'));
    con=[2 -1 -2 1 2];
    con=con-mean(con);
    subj_con=data.*repmat(con,size(data,1),1);
    subj_con=sum(subj_con,2);
    [h,p,ci,stats]=ttest(subj_con);
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
    
    display(sprintf('0-4 ttest: \n'));
    [h,p,ci,stats] = ttest(data(:,1),data(:,5));
    stats
    p
    CohenD=mean(subj_con)/std(subj_con)
end

%plot it:
%compute within-subject SEM:
n=size(data,1);
subjs_av=mean(data,2);
withinSEM=data-repmat(subjs_av,[1,5]);
SEM=abs(std(withinSEM)/sqrt(n));
%set up the graph:
greyscaleColors=repmat(linspace(1,0,5)',1,3);
averageCon=mean(data);
Xticks={'0','1','2','3','4'};
figure;
bar(1:5,averageCon);
set(gca,'XTickLabel',Xticks)
ylabel('RTs(s)','Fontsize',16)
xlabel('Total number of Changes','Fontsize',16);

if acc
    title('RTs (av tasks): only accurate resp','Fontsize',20);
else
    title('RTs (av tasks): all items','Fontsize',20);
end

hold on
for i=1:5
    bar(i,averageCon(i),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i),SEM(i),'k.');
end
hold off


% ttests:
disp('RT DATA:');
for i=1:4
    for ii=(i+1):5
        fprintf('%d-%d ttest: \n', (i-1), (ii-1))
        [h,p,ci,stats] = ttest(data(:,i),data(:,ii));
        p
    end
end

%% plot both accuracy and RT data averaged on tasks for figure 1:

%ACCURACY
curr_data=ResultsOnlyNum.Acc;

if Exc_AD && Exc_AK
    curr_data=curr_data([1 4:end],:);
    display ('AD and AK are excluded from analysis')
elseif Exc_AD && ~Exc_AK
    curr_data=curr_data([1 3:end],:);
    display ('AD is excluded from analysis')
elseif ~Exc_AD && Exc_AK
    curr_data=curr_data([1:2 4:end],:);
    display ('AK is excluded from analysis')
end

lay=1:5;
item=6:10;
lay_data=curr_data(:,lay);
item_data=curr_data(:,item);
data=zeros([size(lay_data),2]);
data(:,:,1)=lay_data;
data(:,:,2)=item_data;
data=mean(data,3);

%plot it:
%compute within-subject SEM:
n=size(data,1);
subjs_av=mean(data,2);
withinSEM=data-repmat(subjs_av,[1,5]);
SEM=abs(std(withinSEM)/sqrt(n));
%set up the graph:
greyscaleColors=repmat(linspace(0.8,0.2,5)',1,3);
averageCon=mean(data);
Xticks={'0','1','2','3','4'};
h=figure;
set(h,'PaperUnits','centimeters');
set(h,'PaperSize',[5 4]);

subplot(2,1,1);
bar(1:5,averageCon);
set(gca,'XTickLabel',Xticks)
ylabel('accuracy','Fontsize',16)
xlabel('changes','Fontsize',16);
%title('Accuracy (av tasks): only accurate resp','Fontsize',20);

hold on
for i=1:5
    bar(i,averageCon(i),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i),SEM(i),'k.');
end
ylim([0.5 1]);
hold off

%RT
curr_data=ResultsOnlyNum.RT.Acc/1000;

if Exc_AD && Exc_AK
    curr_data=curr_data([1 4:end],:);
    display ('AD and AK are excluded from analysis')
elseif Exc_AD && ~Exc_AK
    curr_data=curr_data([1 3:end],:);
    display ('AD is excluded from analysis')
elseif ~Exc_AD && Exc_AK
    curr_data=curr_data([1:2 4:end],:);
    display ('AK is excluded from analysis')
end

lay_data=curr_data(:,lay);
item_data=curr_data(:,item);
data=zeros([size(lay_data),2]);
data(:,:,1)=lay_data;
data(:,:,2)=item_data;
data=mean(data,3);
subplot(2,1,2);
n=size(data,1);
subjs_av=mean(data,2);
withinSEM=data-repmat(subjs_av,[1,5]);
SEM=abs(std(withinSEM)/sqrt(n));
%set up the graph:
greyscaleColors=repmat(linspace(0.8,0.2,5)',1,3);
averageCon=mean(data);
Xticks={'0','1','2','3','4'};

bar(1:5,averageCon);
set(gca,'XTickLabel',Xticks)
ylabel('RT(sec)','Fontsize',16)
xlabel('changes','Fontsize',16);

hold on
for i=1:5
    bar(i,averageCon(i),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i),SEM(i),'k.');
end
ylim([2 2.5]);
hold off
set(h,'PaperUnits','centimeters');
set(h,'PaperSize',[5 4]);

if Rstruct
    M=[M reshape(data,[size(data,1)*size(data,2),1])];
    
    %add the RTs for all items:
    curr_data=ResultsOnlyNum.RT.all_items/1000;
    lay_data=curr_data(:,lay);
    item_data=curr_data(:,item);
    data=zeros([size(lay_data),2]);
    data(:,:,1)=lay_data;
    data(:,:,2)=item_data;
    data=mean(data,3);
    
    M=[M reshape(data,[size(data,1)*size(data,2),1])];
    
    
    %write it up:
    fid = fopen(results_fname, 'w');
    %set up the header:
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\n', ...
        'subject','linearModel','accuracy_rates','RT_accurate','RT_all_items');
    dlmwrite(results_fname,M,'-append','delimiter','\t','precision','%.4f');
    
end




end