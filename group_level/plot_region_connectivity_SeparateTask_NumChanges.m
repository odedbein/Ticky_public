function [data, curr_data]=plot_region_connectivity_SeparateTask_NumChanges(ResultsBetaSerConnectivityOnlyNum,reg,acc,closePrev)
%subj AD (subj #2) had only 12 voxels in the left Ent, removed from
%analyses.

%results in NatComms published paper:
%load('/Volumes/data/Bein/TickyReanalysis/results/data_in_paper/All_regs_connectivityAKcorrectFisherTransformed.mat')
Exc_AD=1;

if closePrev
    close all
end

if acc==1
    curr_data=ResultsBetaSerConnectivityOnlyNum.(reg).Acc;
elseif acc==3
    curr_data=ResultsBetaSerConnectivityOnlyNum.(reg).NonAcc;
else
    curr_data=ResultsBetaSerConnectivityOnlyNum.(reg).all_items;
end

if Exc_AD
    curr_data=curr_data([1 3:end],:);
    display ('AD is excluded from analysis')
end

averageCon=nanmean(curr_data);

lay=3:7;
item=8:12;
%plot 0-4 changes:
%compute within-subject SEM for both tasks:
n=size(curr_data,1);
subjs_av=nanmean(curr_data(:,lay),2);
withinSEM=curr_data(:,lay)-repmat(subjs_av,[1,5]);
SEM=abs(std(withinSEM)/sqrt(n));
subjs_av=nanmean(curr_data(:,item),2);
withinSEM=curr_data(:,item)-repmat(subjs_av,[1,5]);
SEM=[SEM abs(std(withinSEM)/sqrt(n))];
%set up the graph:
greyscaleColors=repmat(linspace(1,0,5)',2,3);
%Xticks=ResultsBetaSerConnectivity.hipp_Para.Acc(1,3:7);
Xticks={'lay: 0','lay: 1','lay: 2','lay: 3','lay: 4','item: 0','item: 1','item: 2','item: 3','item: 4'};
bar(1:10,zeros(1,10),'FaceColor','none');
set(gca,'XTickLabel',Xticks)
ylabel('Beta Series Correlation','Fontsize',16)
xlabel('Total number of Changes','Fontsize',14);

if acc==1
    title(sprintf('%s connectivity: only accurate resp',reg),'Fontsize',20);
elseif acc==3
    title(sprintf('%s connectivity: only inaccurate resp',reg),'Fontsize',20);
else
    title(sprintf('%s connectivity: all items',reg),'Fontsize',20);
end

hold on
for i=1:10
    if (acc==3) && (i==1 || i==5 || i==6 || i==10)
    else
    bar(i,averageCon(i+2),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i+2),SEM(i),'k.');
    end
end
xlim([0,11]);
%ylim([-0.15,.25]);
hold off
%stats:

%%fitting a linear contrast:
fprintf('linear contrast layout task: \n')
lin=[-2 -1 0 1 2];
subj_con=curr_data(:,lay).*repmat(lin,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
stats
p
CohenD=nanmean(subj_con)/std(subj_con)

fprintf('linear contrast item task: \n')
lin=[-2 -1 0 1 2];
subj_con=curr_data(:,item).*repmat(lin,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
stats
p
CohenD=nanmean(subj_con)/std(subj_con)

%without 4-changes condition
% fprintf('linear contrast item task: \n')
% lin=[-2 -1 1 2];
% subj_con=curr_data(:,item(1:4)).*repmat(lin,size(curr_data,1),1);
% subj_con=sum(subj_con,2);
% [h,p,ci,stats]=ttest(subj_con);
% stats
% p
% CohenD=mean(subj_con)/std(subj_con)

%%fitting a quadratic contrast:
fprintf('quadratic contrast layout task: \n')
con=[2 -1 -2 1 2];
con=con-nanmean(con);
subj_con=curr_data(:,lay).*repmat(con,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p
CohenD=nanmean(subj_con)/std(subj_con)

fprintf('quadratic contrast item task: \n')
con=[2 -1 -2 1 2];
con=con-nanmean(con);
subj_con=curr_data(:,item).*repmat(con,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p
CohenD=nanmean(subj_con)/std(subj_con)

%%fitting an all-or-none contrast:
fprintf('all-or-none contrast layout task: \n')
lin=[-1 0.25 0.25 0.25 0.25];
subj_con=curr_data(:,lay).*repmat(lin,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
stats
p
CohenD=nanmean(subj_con)/std(subj_con)

fprintf('all-or-none contrast item task: \n')
lin=[-1 0.25 0.25 0.25 0.25];
subj_con=curr_data(:,item).*repmat(lin,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
stats
p
CohenD=nanmean(subj_con)/std(subj_con)


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
if acc ~= 3
    fprintf('ANOVA connectivity %s layout vs. item task\n',reg)
    stats = rm_anova2(Y,S,F1,F2,{'task:(1)layout/(2)item','#changes'})
end

%do the analysis average tasks:
%%fitting a linear contrast:
data=zeros([size(lay_data),2]);
data(:,:,1)=lay_data;
data(:,:,2)=item_data;
data=nanmean(data,3);
curr_data = data; %just to save it
%one way ANOVA:
n=size(data,1);
Y=reshape(data,size(data,1)*size(data,2),1);
S=repmat([1:n]',5,1);
F1=[ones(n,1);ones(n,1)*2;ones(n,1)*3;ones(n,1)*4;ones(n,1)*5];%number of changes
X=[Y F1 S];

%run the anova:
alpha=.05;
showtable=1;
P = RMAOV1_mod_oded(X,alpha,showtable);

%fprintf('linear contrast average tasks: \n')
lin=[-2 -1 0 1 2];
subj_con=data.*repmat(lin,size(data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
n_pos=sum(subj_con > 0);
CohenD=nanmean(subj_con)/std(subj_con)
fprintf('linear contrast average task: p = %.3f \n N positive: %d, %.2f\n',p,n_pos,n_pos/19)

%fprintf('all-or-none contrast average tasks: \n')
lin=[-1 0.25 0.25 0.25 0.25];
subj_con=data.*repmat(lin,size(data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
n_pos=sum(subj_con > 0);
CohenD=nanmean(subj_con)/std(subj_con)
fprintf('all-or-none contrast average tasks: p = %.3f \n N positive: %d, %.2f\n',p,n_pos,n_pos/19)

%%fitting a quadratic contrast:
con=[2 -1 -2 1 2];
con=con-nanmean(con);
subj_con=data.*repmat(con,size(data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);

CohenD=nanmean(subj_con)/std(subj_con)
n_pos=sum(subj_con > 0);
fprintf('quadratic contrast average tasks: p = %.3f \n N positive: %d, %.2f\n',p,n_pos,n_pos/19)

for i=2:5
    fprintf('0-%d ttest: \n', (i-1))
    [h,p,ci,stats] = ttest(data(:,1),data(:,i));
    p
end


fprintf('1-4 ttest, average tasks: \n');
[h,p,ci,stats] = ttest(data(:,2),data(:,5));
p

fprintf('1-3 ttest, average tasks: \n');
[h,p,ci,stats] = ttest(data(:,2),data(:,4));
p

fprintf('1-2 ttest, average tasks: \n');
[h,p,ci,stats] = ttest(data(:,2),data(:,3));
p

fprintf('2-4 ttest, average tasks: \n');
[h,p,ci,stats] = ttest(data(:,3),data(:,5));
p

fprintf('2-3 ttest, average tasks: \n');
[h,p,ci,stats] = ttest(data(:,3),data(:,4));
p

fprintf('3-4 ttest, average tasks: \n');
[h,p,ci,stats] = ttest(data(:,4),data(:,5));
p

%plot it:
%compute within-subject SEM:
n=size(data,1);
subjs_av=nanmean(data,2);
withinSEM=data-repmat(subjs_av,[1,5]);
SEM=abs(std(withinSEM)/sqrt(n));
%set up the graph:
greyscaleColors=repmat(linspace(1,0,5)',1,3);
averageCon=nanmean(data);
Xticks={'0','1','2','3','4'};
figure;
bar(1:5,zeros(1,5));
set(gca,'XTickLabel',Xticks)
ylabel('Beta Series Correlation','Fontsize',16)
xlabel('Total number of Changes','Fontsize',16);

if acc==1
    title(sprintf('%s connectivity (av tasks): only accurate resp',reg),'Fontsize',20);
elseif acc==3
     title(sprintf('%s connectivity (av tasks): only inaccurate resp',reg),'Fontsize',20);
else
    title(sprintf('%s connectivity (av tasks): all items',reg),'Fontsize',20);
end

hold on
for i=1:5
    if (acc==3) && (i==1 || i==5)
    else
    bar(i,averageCon(i),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i),SEM(i),'k');
    end
end
ylim([-0.15,.2]);
hold off

%plot the match mismatch:
%compute within-subject SEM:
data=[data(:,1) nanmean(data(:,2:end),2)];
n=size(data,1);
subjs_av=nanmean(data,2);
withinSEM=data-repmat(subjs_av,[1,2]);
SEM=abs(std(withinSEM)/sqrt(n));
%set up the graph:
greyscaleColors=repmat([1,0]',1,3);
averageCon=nanmean(data);

figure;
bar(1:2,averageCon);

hold on
for i=1:2
    bar(i,averageCon(i),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i),SEM(i),'k.');
end
hold off




end