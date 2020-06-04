function [con_data, sim_data]=Con_CueMatchImVsOtherRoomsPredictabilityRSA_AcrossSubCorrel(ResultsCueImSimOnlyNum,ResultsBetaSerConnectivityOnlyNum,reg_con,reg_sim,acc,closePrev,onlyConContrasts)

%to select voxels for the similarity analysis:
%currently reported in the appeal:
%load('/Volumes/data/Bein/TickyReanalysis/results/data_in_paper/CueMatchImPredCompOtherRoomsSelVox_ProperGLM_Av5cond_High2Thirds.mat')
%load('/Volumes/data/Bein/TickyReanalysis/results/data_in_paper/All_regs_connectivityAKcorrectFisherTransformed.mat')

Exc_AD=1; %AD is excluded - no entorhinal voxels
type_correl='Pearson';
if closePrev
    close all
end
plot_contorls=0;

% in the paper we report all items, but accuracy was checked in response to
% reviewer
if acc==1
    sim_data=ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_sim).AccCueOnly;
    con_data=ResultsBetaSerConnectivityOnlyNum.(reg_con).Acc(:,3:end); %don't take cue data
elseif acc==2
    sim_data=ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_sim).AccCueAnd0Changes;
    con_data=ResultsBetaSerConnectivityOnlyNum.(reg_con).Acc(:,3:end); %don't take cue data
elseif acc==3
    sim_data=ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_sim).NonAcc;
    con_data=ResultsBetaSerConnectivityOnlyNum.(reg_con).Acc(:,3:end); %don't take cue data
else
    sim_data=ResultsCueImSimOnlyNum.AvPerTaskNo0Changes.(reg_sim).all_items;
    con_data=ResultsBetaSerConnectivityOnlyNum.(reg_con).all_items(:,3:end); %don't take cue data
end


if Exc_AD
    display ('AD is excluded from analysis')
    sim_data=sim_data([1 3:end],:);
    con_data=con_data([1 3:end],:);
end


lay=1:5;
item=6:10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the similarity measure:
%for each task:
color=[0.5    0.5    0.5
       1      1      1
       0.5    0.5    0.5
       1      1      1];

%first, plot whether there is a differnece between same-different scene in
%predictions in that region:
curr_data=sim_data(:,[1:2 4:5]);%
num_cond=size(curr_data,2);
num_subj=size(curr_data,1);
meanSim=mean(curr_data);
subjAv=mean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=std(withinEr)/sqrt(size(curr_data,1));
f=figure;
set(f,'name','predictability per task: same vs. diff. room','numbertitle','off');

bar(1:num_cond,meanSim,'FaceColor',[1,1,1]);
hold on
for i=1:length(meanSim) %
    bar(i,meanSim(i),'FaceColor',color(i,:));
    errorbar(i,meanSim(i),SEM(i),'k');
end
ylim([-0.03 0.05]);
Xticks={'lay: same room','lay: other rooms','item: same room','item: other rooms'};
set(gca,'XTickLabel',Xticks)
set(gca, 'FontSize', 14)
title(reg_sim)

% Plot subject data
subjDataX = [ones(1,num_subj), ones(1,num_subj)*2,ones(1,num_subj)*3,ones(1,num_subj)*4];
scatter(subjDataX,[curr_data(:,1);curr_data(:,2);curr_data(:,3);curr_data(:,4)]', 22, 'ko', 'filled')
% Plot lines connecting rPlus/rMinus and assocPlus/assocMinus
for ii = 1:num_subj
    plot(1:4, [curr_data(ii,1), curr_data(ii,2),curr_data(ii,3),curr_data(ii,4)])
end
hold off

%do the t-test:
fprintf('lay: same/other t-tests \n');
[h,p,ci,stats]=ttest(curr_data(:,1),curr_data(:,2))
fprintf('lay - mean same: %.3f, SD same: %.3f,\n', mean(curr_data(:,1)), std(curr_data(:,1)));
fprintf('lay - mean other: %.4f, SD other: %.3f,\n', mean(curr_data(:,2)), std(curr_data(:,2)));

fprintf('item: same/other t-tests \n');
[h,p,ci,stats]=ttest(curr_data(:,3),curr_data(:,4))
fprintf('item - mean same: %.3f, SD same: %.3f,\n', mean(curr_data(:,3)), std(curr_data(:,3)));
fprintf('item - mean other: %.4f, SD other: %.3f,\n', mean(curr_data(:,4)), std(curr_data(:,4)));

%% now average similarity across tasks:
color=[0.5    0.5    0.5
    1      1      1];
%first, plot whether there is a differnece between same-different scene in
%predictions in that region:
av_sim_data=[];
for cols=1:3
    av_sim_data=[av_sim_data mean([sim_data(:,cols),sim_data(:,(cols+3))],2)];
end
curr_data=av_sim_data(:,1:2);
num_cond=size(curr_data,2);
num_subj=size(curr_data,1);
meanSim=mean(curr_data);
subjAv=mean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=std(withinEr)/sqrt(size(curr_data,1));
f=figure;
set(f,'name','predictability av. tasks: same vs. diff. room','numbertitle','off');

bar(1:num_cond,meanSim,'FaceColor',[1,1,1]);
hold on
for i=1:length(meanSim) %
    bar(i,meanSim(i),'FaceColor',color(i,:));
    errorbar(i,meanSim(i),SEM(i),'k');
end
ylim([-0.03 0.05]);
Xticks={'same room','other rooms'};
set(gca,'XTickLabel',Xticks)
set(gca, 'FontSize', 14)
title(reg_sim)

% Plot subject data
subjDataX = [ones(1,num_subj), ones(1,num_subj)*2];
scatter(subjDataX,[curr_data(:,1);curr_data(:,2)]', 22, 'ko', 'filled')
% Plot lines connecting rPlus/rMinus and assocPlus/assocMinus
for ii = 1:num_subj
    plot([1,2], [curr_data(ii,1), curr_data(ii,2)])
end
hold off

%do the t-test:
fprintf('average on tasks: same/other t-tests \n');
[h,p,ci,stats]=ttest(curr_data(:,1),curr_data(:,2))
fprintf('mean same: %.3f, SD same: %.3f,\n', mean(curr_data(:,1)), std(curr_data(:,1)));
fprintf('mean other: %.4f, SD other: %.3f,\n', mean(curr_data(:,2)), std(curr_data(:,2)));

%plot the difference:
curr_data=av_sim_data(:,3);%
num_cond=size(curr_data,2);
num_subj=size(curr_data,1);
meanSim=mean(curr_data);
%subjAv=mean(curr_data,2);
%withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=std(curr_data)/sqrt(size(curr_data,1));
f=figure;
set(f,'name','predictability av. tasks: same vs. diff. room','numbertitle','off');

bar(1:num_cond,meanSim,'FaceColor',[1,1,1]);
hold on
for i=1:length(meanSim) %
    bar(i,meanSim(i),'FaceColor',color(i,:));
    errorbar(i,meanSim(i),SEM(i),'k');
end
ylim([[-0.015, 0.021]]);
Xticks={'same minus other room'};
set(gca,'XTickLabel',Xticks)
set(gca, 'FontSize', 14)
title(reg_sim)
num_pos=sum(curr_data > 0);

% Plot subject data
subjDataX = ones(1,num_subj);
scatter(subjDataX,curr_data', 26, 'ko', 'filled')
%effect size of the difference
CohenD=nanmean(curr_data)/nanstd(curr_data);
fprintf('Cohen''s D: same/other: %.02f \n',CohenD);
fprintf('same-other N positive: %.3f, %.4f,\n', num_pos, num_pos/19);

%% do the correlation with connectivity:
% In the paper I do average across tasks - these are just checks to see diff tasks/connectivity contrast look similar...
f=figure;
%compute the linear score:
%lay task:
lin=[-2 -1 0 1 2];
subj_con=con_data(:,lay).*repmat(lin,size(con_data,1),1);
con=sum(subj_con,2);
sim=sim_data(:,3);
[R,P] = corr(sim,con,'type',type_correl);
fprintf(sprintf('lay: con linear score with average predictability: \n %s: %0.3f \n p: %0.3f \n',type_correl,R,P));
subplot(2,2,1);
hold on
scatter(sim,con,'black')
coef_fit = polyfit(sim,con,1);
y_fit = polyval(coef_fit,xlim);
plot(xlim,y_fit,'black');

title('lay: con linear score with average predictability','Fontsize',14);
xlabel('worse prediction   better prediction','Fontsize',12);
ylabel('less lin increase        more lin increase','Fontsize',12);
hold off

%item task:
lin=[-2 -1 0 1 2];
subj_con=con_data(:,item).*repmat(lin,size(con_data,1),1);
con=sum(subj_con,2);
sim=sim_data(:,6);
[R,P] = corr(sim,con,'type',type_correl);
fprintf(sprintf('item: con linear score with average predictability: \n %s: %0.3f \n p: %0.3f \n',type_correl,R,P));
subplot(2,2,2);
hold on
scatter(sim,con,'black')
coef_fit = polyfit(sim,con,1);
y_fit = polyval(coef_fit,xlim);
plot(xlim,y_fit,'black');

title('item: con linear score with average predictability','Fontsize',14);
xlabel('worse prediction   better prediction','Fontsize',12);
ylabel('less lin increase        more lin increase','Fontsize',12);
hold off


%compute the match-mismatch score:
%lay task:
lin=[-1 0.25 0.25 0.25 0.25];
subj_con=con_data(:,lay).*repmat(lin,size(con_data,1),1);
con=sum(subj_con,2);
sim=sim_data(:,3);
[R,P] = corr(sim,con,'type',type_correl);
fprintf(sprintf('lay: con match/mismatch score with average predictability: \n %s: %0.3f \n p: %0.3f \n',type_correl,R,P));
subplot(2,2,3);
hold on
scatter(sim,con,'black')
coef_fit = polyfit(sim,con,1);
y_fit = polyval(coef_fit,xlim);
plot(xlim,y_fit,'black');

title('lay: con match/mismatch score with average predictability','Fontsize',14);
xlabel('worse prediction   better prediction','Fontsize',12);
ylabel('less match/mismatch increase        more match/mismatch increase','Fontsize',12);
hold off


sim(17)=[];con(17)=[];
[R,P] = corr(sim,con,'type',type_correl);
fprintf(sprintf('lay: w/o outlier: \n %s: %0.3f \n p: %0.3f \n',type_correl,R,P));
[b,stats] = robustfit(sim,con);
fprintf('robust reg: b %0.3f p: %0.3f \n',b(2),stats.p(2));

%item task:
lin=[-1 0.25 0.25 0.25 0.25];
subj_con=con_data(:,item).*repmat(lin,size(con_data,1),1);
con=sum(subj_con,2);
sim=sim_data(:,6);
[R,P] = corr(sim,con,'type',type_correl);
fprintf(sprintf('item: con match/mismatch score with average predictability: \n %s: %0.3f \n p: %0.3f \n',type_correl,R,P));
subplot(2,2,4);
hold on
scatter(sim,con,'black')
coef_fit = polyfit(sim,con,1);
y_fit = polyval(coef_fit,xlim);
plot(xlim,y_fit,'black');

title('item: con match/mismatch score with average predictability','Fontsize',14);
xlabel('worse prediction   better prediction','Fontsize',12);
ylabel('less match/mismatch increase        more match/mismatch increase','Fontsize',12);
hold off

%% average on both tasks

con_data=zeros([size(con_data_bothTasks(:,lay)),2]);
con_data(:,:,1)=con_data_bothTasks(:,lay);
con_data(:,:,2)=con_data_bothTasks(:,item);
con_data=mean(con_data,3);

f=figure;
%compute the linear score - just to check that looks similar, match>mismatch below contrast is reported:
lin=[-2 -1 0 1 2];
subj_con=con_data.*repmat(lin,size(con_data,1),1);
con=sum(subj_con,2);
sim=mean([sim_data(:,3),sim_data(:,6)],2);
[R,P] = corr(sim,con,'type',type_correl); %,'Tail','right'
fprintf(sprintf('av tasks: con linear score with average predictability (two-tailed): \n %s: %0.3f \n p: %0.3f \n',type_correl,R,P));
subplot(1,2,1);
hold on
scatter(sim,con,'black')
coef_fit = polyfit(sim,con,1);
y_fit = polyval(coef_fit,xlim);
plot(xlim,y_fit,'black');

title('av tasks: con linear score with average predictability','Fontsize',14);
xlabel('worse prediction   better prediction','Fontsize',12);
ylabel('less lin increase        more lin increase','Fontsize',12);
hold off

%THIS IS REPORTED: compute the match-mismatch score:
lin=[-1 0.25 0.25 0.25 0.25];
subj_con=con_data.*repmat(lin,size(con_data,1),1);
con=sum(subj_con,2);
sim=mean([sim_data(:,3),sim_data(:,6)],2);
[R,P] = corr(sim,con,'type',type_correl); %,'Tail','right'
%[~,~,RL,RU] = corrcoef(sim,con);
fprintf('av tasks: con match/mismatch score with average predictability (two-tailed): \n %s: %0.3f \n p: %0.3f \n',type_correl,R,P);
[b,stats] = robustfit(sim,con);
fprintf('robust reg: b %0.3f p: %0.3f \n',b(2),stats.p(2));
subplot(1,2,2);
hold on
scatter(sim,con,'black')
coef_fit = polyfit(sim,con,1);
y_fit = polyval(coef_fit,xlim);
plot(xlim,y_fit,'black');

title('av tasks: con match/mismatch score with average predictability','Fontsize',14);
xlabel('worse prediction   better prediction','Fontsize',12);
ylabel('less match/mismatch increase        more match/mismatch increase','Fontsize',12);
hold off


%% plot only the match/mismatch:
figure;
hold on
scatter(sim,con,'black')
coef_fit = polyfit(sim,con,1);
y_fit = polyval(coef_fit,xlim);
plot(xlim,y_fit,'black');
hold off


end