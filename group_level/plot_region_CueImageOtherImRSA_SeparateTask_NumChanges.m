function [data, curr_data]=plot_region_CueImageOtherImRSA_SeparateTask_NumChanges(ResultsCueImSimOnlyNum,reg,acc,same_others,closePrev)
%subj AD (subj #2) had only 12 voxels in the left Ent,removed from analysis

%same/same_others: same is only the cue with the image on that trial.
%same_others - other trials subtracted to control for specificity -
%reported in the supplementary.

%the data currently reported in the supplementary
%without voxel selection:
%load('/Volumes/data/Bein/TickyReanalysis/results/rsa_PE/CueImageVsOtherRoomsExcludeSameRoom_Ftrans.mat')

%reported in the appeal - with voxel selection:
% load('/Volumes/data/Bein/TickyReanalysis/results/data_in_paper/CueImageVsOtherRoomsExcludeSameRoom_Ftrans_ProperGLMAv5CondHigh2Thirds.mat')
Exc_AD=1;
if closePrev
    close all
end

if same_others
    
    if acc==1
        curr_data=ResultsCueImSimOnlyNum.(reg).same_others.Acc;
    elseif acc==3
        curr_data=ResultsCueImSimOnlyNum.(reg).same_others.NonAcc;
    else
        curr_data=ResultsCueImSimOnlyNum.(reg).same_others.all_items;
    end
    
else
    if acc==1
        curr_data=ResultsCueImSimOnlyNum.(reg).same.Acc;
    elseif acc==3
        curr_data=ResultsCueImSimOnlyNum.(reg).same.NonAcc;
    else
        curr_data=ResultsCueImSimOnlyNum.(reg).same.all_items;
    end
end

    
if Exc_AD
    curr_data=curr_data([1 3:end],:);
end

averageCon=nanmean(curr_data);

figure;
lay=1:5;
item=6:10;
%plot 0-4 changes:
%compute within-subject SEM for both tasks:
n=size(curr_data,1);
subjs_av=nanmean(curr_data(:,lay),2);
withinSEM=curr_data(:,lay)-repmat(subjs_av,[1,5]);
SEM=abs(nanstd(withinSEM)/sqrt(n));
subjs_av=nanmean(curr_data(:,item),2);
withinSEM=curr_data(:,item)-repmat(subjs_av,[1,5]);
SEM=[SEM abs(nanstd(withinSEM)/sqrt(n))];
%set up the graph:
greyscaleColors=repmat(linspace(1,0,5)',2,3);
%Xticks=ResultsBetaSerConnectivity.hipp_Para.Acc(1,3:7);
Xticks={'lay: 0','lay: 1','lay: 2','lay: 3','lay: 4','item: 0','item: 1','item: 2','item: 3','item: 4'};
bar(1:10,averageCon);
set(gca,'XTickLabel',Xticks)
ylabel('Cue-Image similarity (Fisher-transformed)','Fontsize',16)
xlabel('Total number of Changes','Fontsize',14);

if acc==1
    title(sprintf('%s Cue-Image similarity: only accurate resp',reg),'Fontsize',20);
elseif acc==3
    title(sprintf('%s Cue-Image similarity: only inaccurate resp',reg),'Fontsize',20);

else
    title(sprintf('%s Cue-Image similarity: all items',reg),'Fontsize',20);
end

hold on
for i=1:10
    bar(i,averageCon(i),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i),SEM(i),'k.');
end
hold off
ylim([0.28, 0.35]);
%ylim([0.30, 0.33]);
xlim([0,11]);
%stats:

%%fitting a linear contrast:
fprintf('linear contrast layout task: \n')
lin=[-2 -1 0 1 2];
subj_con=curr_data(:,lay).*repmat(lin,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p

fprintf('linear contrast item task: \n')
lin=[-2 -1 0 1 2];
subj_con=curr_data(:,item).*repmat(lin,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p


%%fitting a quadratic contrast:
fprintf('quadratic contrast layout task: \n')
con=[2 -1 -2 1 2];
con=con-nanmean(con);
subj_con=curr_data(:,lay).*repmat(con,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p

fprintf('quadratic contrast item task: \n')
con=[2 -1 -2 1 2];
con=con-nanmean(con);
subj_con=curr_data(:,item).*repmat(con,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p

%%fitting an all-or-none contrast:
fprintf('all-or-none contrast layout task: \n')
lin=[-1 0.25 0.25 0.25 0.25];
subj_con=curr_data(:,lay).*repmat(lin,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p

fprintf('all-or-none contrast item task: \n')
lin=[-1 0.25 0.25 0.25 0.25];
subj_con=curr_data(:,item).*repmat(lin,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p

fprintf('0-1 ttest layout: \n')
[h,p,ci,stats] = ttest(curr_data(:,lay(1)),curr_data(:,lay(2)));
p
fprintf('0-2 ttest layout: \n')
[h,p,ci,stats] = ttest(curr_data(:,lay(1)),curr_data(:,lay(3)));
p

fprintf('0-3 ttest layout: \n')
[h,p,ci,stats] = ttest(curr_data(:,lay(1)),curr_data(:,lay(4)));
p

fprintf('0-4 ttest layout: \n')
[h,p,ci,stats] = ttest(curr_data(:,lay(1)),curr_data(:,lay(5)));
p

fprintf('0-1 ttest item: \n')
[h,p,ci,stats] = ttest(curr_data(:,item(1)),curr_data(:,item(2)));
p
fprintf('0-2 ttest item: \n')
[h,p,ci,stats] = ttest(curr_data(:,item(1)),curr_data(:,item(3)));
p
fprintf('0-3 ttest item: \n')
[h,p,ci,stats] = ttest(curr_data(:,item(1)),curr_data(:,item(4)));
p
fprintf('0-4 ttest item: \n')
[h,p,ci,stats] = ttest(curr_data(:,item(1)),curr_data(:,item(5)));
p

%% ANOVA lay_item
%run the ANOVA:
lay_data=curr_data(:,lay);
item_data=curr_data(:,item);
n=size(lay_data,1);
Y=reshape(lay_data,size(lay_data,1)*size(lay_data,2),1);
Y=[Y; reshape(item_data,size(item_data,1)*size(item_data,2),1)];
S=repmat([1:n]',10,1);
F1=[ones(n*5,1);zeros(n*5,1)];%task
F2=repmat([ones(n,1);ones(n,1)*2;ones(n,1)*3;ones(n,1)*4;ones(n,1)*5],2,1);%number of changes
fprintf('ANOVA Cue-Image similarity %s layout vs. item task\n',reg)
stats = rm_anova2(Y,S,F1,F2,{'task:(1)layout/(2)item','#changes'})

%do the analysis average tasks:
%%fitting a linear contrast:
data=zeros([size(lay_data),2]);
data(:,:,1)=lay_data;
data(:,:,2)=item_data;
data=nanmean(data,3);


%fprintf('linear contrast average tasks: \n')
lin=[-.5 -.25 0 .25 .5];
subj_con=data.*repmat(lin,size(data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
fprintf('linear contrast average tasks: M = %.3f, SD = %.3f\n',mean(subj_con),std(subj_con))
ci
stats
p
CohenD=nanmean(subj_con)/nanstd(subj_con)

%fprintf('all-or-none contrast average tasks: \n')
lin=[-1 0.25 0.25 0.25 0.25];
subj_con=data.*repmat(lin,size(data,1),1);
subj_con=sum(subj_con,2);
fprintf('all-or-none contrast average tasks: M = %.3f, SD = %.3f\n',mean(subj_con),std(subj_con))
[h,p,ci,stats]=ttest(subj_con);
ci
n_pos=sum(subj_con < 0);
CohenD=nanmean(subj_con)/std(subj_con)
fprintf('all-or-none contrast average tasks: t = %.4f, p = %.4f \n N negative: %d, %.2f\n',stats.tstat,p,n_pos,n_pos/19)

%%fitting a quadratic contrast:
fprintf('quadratic contrast average tasks: \n')
con=[2 -1 -2 -1 2];
con=con-nanmean(con);
subj_con=data.*repmat(con,size(data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p


fprintf('0-1 ttest average tasks: \n')
[h,p,ci,stats] = ttest(data(:,1),data(:,2));
p
fprintf('0-2 ttest average tasks: \n')
[h,p,ci,stats] = ttest(data(:,1),data(:,3));
stats
p
fprintf('0-3 ttest average tasks: \n')
[h,p,ci,stats] = ttest(data(:,1),data(:,4));
stats
p

fprintf('0-4 ttest average tasks: \n')
[h,p,ci,stats] = ttest(data(:,1),data(:,5));
p

%plot it:
n=size(data,1);
subjs_av=nanmean(data,2);
%withinSEM=data-repmat(subjs_av,[1,5]);
SEM=abs(nanstd(data)/sqrt(n));
%set up the graph:
greyscaleColors=repmat(linspace(1,0,5)',1,3);
averageCon=nanmean(data);
Xticks={'0','1','2','3','4'};
figure;
bar(1:5,zeros(1,5),'FaceColor', 'none');
set(gca,'XTickLabel',Xticks)
ylabel('Cue-Image similarity','Fontsize',16)
xlabel('Total number of Changes','Fontsize',16);

if acc==1
    title(sprintf('%s Cue-Image similarity (av tasks): only accurate resp',reg),'Fontsize',20);
elseif acc==3
    title(sprintf('%s Cue-Image similarity (av tasks): only inaccurate resp',reg),'Fontsize',20);
else
    title(sprintf('%s Cue-Image similarity (av tasks): all items',reg),'Fontsize',20);
end

hold on
for i=1:5
    bar(i,averageCon(i),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i),SEM(i),'k');
end
%ylim([0.28, 0.35]);
ylim([0.3, 0.34]);
hold off

end