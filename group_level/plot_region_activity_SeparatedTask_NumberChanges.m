function [data, curr_data]=plot_region_activity_SeparatedTask_NumberChanges(ResultsAverageActivityOnlyNum,reg,acc,cue_im,new_fig)

%load('/Volumes/data/Bein/TickyReanalysis/results/univariate/Univar_cue_image_single_trials.mat')

if new_fig
    figure;
end
Exc_AD=1; % not enough entorhinal voxels - excluded from analysis

%cue_im: 1:image, 0:cue
if cue_im
    if acc
        curr_data=ResultsAverageActivityOnlyNum.(reg).image.Acc;
    else
        curr_data=ResultsAverageActivityOnlyNum.(reg).image.all_items;
    end
else
    if acc
        curr_data=ResultsAverageActivityOnlyNum.(reg).cue.Acc;
    else
        curr_data=ResultsAverageActivityOnlyNum.(reg).cue.all_items;
    end
end



if Exc_AD
    curr_data=curr_data([1 3:end],:);
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
ylabel('t parameter (beta estimates)','Fontsize',16);
xlabel('Total number of Changes','Fontsize',14);

if acc
    title(sprintf('%s activity: only accurate resp',reg),'Fontsize',20);
else
    title(sprintf('%s activity: all items',reg),'Fontsize',20);
end

hold on
for i=1:10
    bar(i,averageCon(i),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i),SEM(i),'k');
end
hold off
xlim([0,11]);
%stats:

%%fitting a linear contrast:
display(sprintf('linear contrast layout task: \n'));
lin=[-2 -1 0 1 2];
subj_con=curr_data(:,lay).*repmat(lin,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p

display(sprintf('linear contrast item task: \n'));
lin=[-2 -1 0 1 2];
subj_con=curr_data(:,item).*repmat(lin,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p


%%fitting a quadratic contrast:
display(sprintf('quadratic contrast layout task: \n'));
con=[2 -1 -2 1 2];
con=con-mean(con);
subj_con=curr_data(:,lay).*repmat(con,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p

display(sprintf('quadratic contrast item task: \n'));
con=[2 -1 -2 1 2];
con=con-mean(con);
subj_con=curr_data(:,item).*repmat(con,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p

%%fitting an all-or-none contrast:
display(sprintf('all-or-none contrast layout task: \n'));
lin=[-1 0.25 0.25 0.25 0.25];
subj_con=curr_data(:,lay).*repmat(lin,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p

display(sprintf('all-or-none contrast item task: \n'));
lin=[-1 0.25 0.25 0.25 0.25];
subj_con=curr_data(:,item).*repmat(lin,size(curr_data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p


display(sprintf('0-4 ttest layout: \n'));
[h,p,ci,stats] = ttest(curr_data(:,lay(1)),curr_data(:,lay(5)));
p

display(sprintf('0-4 ttest item: \n'));
[h,p,ci,stats] = ttest(curr_data(:,item(1)),curr_data(:,item(5)));
p


%% run the ANOVA for two tasks:
lay_data=curr_data(:,lay);
item_data=curr_data(:,item);
n=size(lay_data,1);
Y=reshape(lay_data,size(lay_data,1)*size(lay_data,2),1);
Y=[Y; reshape(item_data,size(item_data,1)*size(item_data,2),1)];
S=repmat([1:n]',10,1);
F1=[ones(n*5,1);zeros(n*5,1)];%task
F2=repmat([ones(n,1);ones(n,1)*2;ones(n,1)*3;ones(n,1)*4;ones(n,1)*5],2,1);%number of changes
display('ANOVA activity %s layout vs. item task\n',reg);
stats = rm_anova2(Y,S,F1,F2,{'task:(1)layout/(2)item','#changes'})

%average tasks:
data=zeros([size(lay_data),2]);
data(:,:,1)=lay_data;
data(:,:,2)=item_data;
data=mean(data,3);

%%fitting a match > mismatch contrast:
disp('match < mismatch contrast average tasks: \n');
lin=[-1 0.25 0.25 0.25 0.25];
subj_con=data.*repmat(lin,size(data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p
CohenD=nanmean(subj_con)/nanstd(subj_con)


%%fitting a linear contrast:
disp('linear contrast average tasks: \n');
lin=[-2 -1 0 1 2];
subj_con=data.*repmat(lin,size(data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
ci
stats
p
CohenD=nanmean(subj_con)/nanstd(subj_con)


%%fitting a quadratic contrast:
disp('quadratic contrast average tasks: \n');
con=[2 -1 -2 1 2];
con=con-mean(con);
subj_con=data.*repmat(con,size(data,1),1);
subj_con=sum(subj_con,2);
[h,p,ci,stats]=ttest(subj_con);
p


for i=1:4
    for ii=(i+1):5
    
    [h,p,ci,stats] = ttest(data(:,i),data(:,ii));
    fprintf('%d-%d ttest: p = %.3f, t = %.3f \n',(i-1),(ii-1),p,stats.tstat)
    end
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
ylabel('Beta Series Correlation','Fontsize',16)
xlabel('Total number of Changes','Fontsize',16);

if acc
    title(sprintf('%s activity (av tasks): only accurate resp',reg),'Fontsize',20);
else
    title(sprintf('%s activity (av tasks): all items',reg),'Fontsize',20);
end

hold on
for i=1:5
    bar(i,averageCon(i),'FaceColor',greyscaleColors(i,:));
    errorbar(i,averageCon(i),SEM(i),'k');
end
hold off

end