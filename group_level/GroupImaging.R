#setwd('~')
# I used this code to run the analyses and to plot - some packages need to be loaded/commented out depands on what you want to do
#(sorry, bad coding... just began my R journey with this)
# clear workspace
rm(list = ls())

#load some libraries
library(dplyr)
library(tidyverse)
#library(Matrix)
library(lme4)
library(readxl)
library(effsize)
library(ggplot2)
library(data.table)
library(knitr)
library(reshape)
library(standardize)
library(ggeffects)
library(sjPlot)
library(sjmisc)
library(sjstats)
#library(MASS) #for the robust regression
library(psychometric) #for CI of Pearson's r
library(Hmisc) #for stats on corr
library(ggpubr)
library(gridExtra)
library(RColorBrewer)

############# definitons, load data #######
exp_name <- "TickyReanalysis"
proj_dir <- sprintf("/Volumes/data/Bein/%s", exp_name)
analysis_scripts_dir <- sprintf("%s/scripts/Rscritps",proj_dir)
data_dir <- sprintf("%s/results/connectivity",proj_dir)
behav_data_dir <- sprintf("%s/results/behavior",proj_dir)
setwd(proj_dir)

##blue colors:
blue_cols <- c(rgb(210/255,224/255,244/255),
               rgb(162/255,193/255,232/255),
               rgb(110/255,159/255,220/255),
               rgb(36/255,71/255,113/255),
               rgb(23/255,49/255,80/255))

green_cols<- c(rgb(147/255,227/255,185/255),
               rgb(93/255,211/255,149/255),
               rgb(23/255,206/255,109/255),
               rgb(15/255,163/255,84/255),
               rgb(5/255,103/255,49/255))

stdcoeff <- function (MOD) {
  b <- summary(MOD)$coef[-1, 1]
  sx <- sd(MOD$model[-1])
  sy <- sd(MOD$model[1])
  beta <- b * sx/sy
  return(beta)
  }

############## Connectivity main analysis ####################

## threeway analysis by region,hemisphere,and #changes (reported in the Results)
data <- read_excel(paste(data_dir, "connectivity_3way_ANOVA_hemisphere_reg_changes.xlsx", sep="/"), sheet = "all_reg")

data$num_changes <- factor(data$linearModel)
data$subID <- factor(data$subject)
#factoring the linear model gives us the number of changes
results.threeway_aov=aov(BseriesCor ~ hemisphere*reg*num_changes + Error(factor(subject) / (hemisphere*reg*num_changes)), data = data)
summary(results.threeway_aov)
eta_sq(results.threeway_aov, partial = TRUE)

#within left hemisphere (reported in the Results):
curr_data <- data %>%
  filter(hemisphere == "left")
results.twoway_aov=aov(BseriesCor ~ reg*num_changes + Error(factor(subject) / (reg*num_changes)), data = curr_data)
summary(results.twoway_aov)
eta_sq(results.twoway_aov, partial = TRUE)

#within each hemisphere (reported in the Results):
curr_data <- data %>%
  filter(hemisphere == "right")
results.twoway_aov=aov(BseriesCor ~ reg*factor(num_changes) + Error(factor(subject) / (reg*factor(num_changes))), data = curr_data)
summary(results.twoway_aov)

#one-way ANOVA within each region (reported in the Results):
curr_data <- data %>%
  filter(hemisphere == "left") %>%
    filter(reg == "ca3")
results.oneway_aov=aov(BseriesCor ~ factor(num_changes) + Error(factor(subject) / factor(num_changes)), data = curr_data)
summary(results.oneway_aov)
eta_sq(results.oneway_aov, partial = TRUE)

#one-way ANOVA within each region (reported in the Results):
curr_data <- data %>%
  filter(hemisphere == "left") %>%
  filter(reg == "ent")
results.oneway_aov=aov(BseriesCor ~ factor(num_changes) + Error(factor(subject) / factor(num_changes)), data = curr_data)
summary(results.oneway_aov)
eta_sq(results.oneway_aov, partial = TRUE)

############## lCA1-lEnt connectivity ####################

#plots for Figure 2:
all_subjs_filename <- sprintf("%s/lCA1_lEnt_allItemsExcAD_AKcorrect.txt", data_dir )
data <- read.table(all_subjs_filename, header=T)
data$fac_linearModel <- factor(data$linearModel) 

data_sum <- data %>%
  group_by(fac_linearModel) %>%
  summarise(mean_cor = mean(BseriesCor),
            SEM_cor = sd(BseriesCor)/sqrt(n()),
            n_data=n())

data_con_temp <- select(data,c(-linearModel,-lCA1_activity,-lEnt_activity)) %>%
  mutate(mm_con = all_or_none*BseriesCor)

data_con  <- data_con_temp %>%
  group_by(subject) %>%
  summarise(subj_con = sum(mm_con))

#position = position_nudge(x = .1)
## plot (with points):
# plt_ent <- ggplot(data_sum, aes(x = fac_linearModel, y = mean_cor, fill=fac_linearModel)) +
#   geom_col(color = "black") +
#   geom_point(data = data, aes(x = fac_linearModel, y = BseriesCor,fill=fac_linearModel), alpha = .4,
#              position = position_jitterdodge(jitter.width = .7,dodge.width = -.1)) +
#   geom_errorbar(aes(ymin = mean_cor - SEM_cor, ymax = mean_cor + SEM_cor),width=0,position = position_nudge(x = .2)) +
#   scale_fill_manual(values = green_cols) +
#   scale_color_manual(values = "black") +
#   labs(x = '# of changes', y = 'Functional connectivity') + 
#   theme_classic() + theme(legend.position = "none")

## plot (no points):
plt_ent <- ggplot(data_sum, aes(x = fac_linearModel, y = mean_cor, fill=fac_linearModel)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = mean_cor - SEM_cor, ymax = mean_cor + SEM_cor),width=0) +
  scale_fill_manual(values = green_cols) +
  scale_color_manual(values = "black") +
  labs(x = '# of changes', y = 'Functional connectivity') + 
  theme_classic() + theme(legend.position = "none")


plt_ent_subj_con <- ggplot(data_con, aes(x = 1, y = subj_con, color = "1")) + 
  geom_hline(yintercept=0, colour="black", linetype = "dotted") +
  xlim(0.5,1.5) +
  geom_jitter(width = .08,size = 2) + 
  scale_color_manual(values = green_cols[4]) +
  labs(y = 'match < mismatch') + 
  theme_classic() + theme(legend.position = "none",
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank())

#mixed effects linear models (Reported in Results):
## analysis with only a random intercept:
data$all_or_none[data$all_or_none==-1]=0
data$all_or_none[data$all_or_none==0.25]=1
#fit a linear model
results.linearModel=lmer(BseriesCor~linearModel+(1|as.factor(data$subject)),data,REML=F)
summary(results.linearModel)


#fit the all-or-none model
results.all_or_none=lmer(BseriesCor~as.factor(all_or_none)+(1|as.factor(data$subject)),data,REML=F)
summary(results.all_or_none)

#compare models
anova(results.linearModel,results.all_or_none)
anova(results.all_or_none,results.linearModel)

#a model with both
results.all_or_noneAndLinear=lmer(BseriesCor~as.factor(all_or_none)+linearModel+(1|as.factor(data$subject)),data,REML=F)
summary(results.all_or_noneAndLinear)

#compare models
anova(results.all_or_noneAndLinear,results.linearModel)
anova(results.all_or_noneAndLinear,results.all_or_none)


############## control for univariate ############
data_dir <- sprintf("%s/results/connectivity",proj_dir)
#mixed-effect models controling for univar (Reported in Results): 
#The script that creates this file is: makeRdataStructureConUnivar.m
all_subjs_filename <- sprintf("%s/lCA1_lEnt_allItemsExcAD_AKcorrect_Univar_single_trials_withAccRT.txt", data_dir )
data <- read.table(all_subjs_filename, header=T)

#first, see that the all_or_none works alone:
null_mdl=lmer(BseriesCor~1+(1|as.factor(data$subject)),data,REML=F)
all_or_noneModel=lmer(BseriesCor~all_or_none+(1|as.factor(data$subject)),data,REML=F)
summary(all_or_noneModel)
anova(null_mdl,all_or_noneModel)

#fit a null model with only the activity:
null_mdl=lmer(BseriesCor~scale(lCA1_activity)*scale(lEnt_activity)+(1|as.factor(data$subject)),data,REML=F)

#fit a linear model - to see that this works as well
linearModel=lmer(BseriesCor~linearModel+scale(lCA1_activity)*scale(lEnt_activity)+(1|as.factor(data$subject)),data,REML=F)
summary(linearModel)

#fit all or none model - this is the main one, REPORTED IN RESULTS:
all_or_noneModel=lmer(BseriesCor~all_or_none+scale(lCA1_activity)*scale(lEnt_activity)+(1|as.factor(data$subject)),data,REML=F)
summary(all_or_noneModel)
anova(null_mdl,all_or_noneModel)
comp=anova(null_mdl,all_or_noneModel)
AIC_dif=comp$AIC[1]-comp$AIC[2]
BIC_dif=comp$BIC[1]-comp$BIC[2]
comp
AIC_dif
BIC_dif


###### ACC/RT actual data per participant ######
#mixed-effect models controling for univar (Reported in Supplementary Information Note 2): 
all_subjs_filename <- sprintf("%s/lCA1_lEnt_allItemsExcAD_AKcorrect_Univar_single_trials_withAccRT.txt", data_dir )

data <- read.table(all_subjs_filename, header=T)
behav_filename <- sprintf("%s/Behavior_ExcAD_AKcorrect.txt", behav_data_dir )
bdata <- read.table(behav_filename, header=T)
#add some cols from bdata to data
data$acc=bdata$accuracy_rates
data$RT_acc=bdata$RT_accurate

#first, see that the all_or_none works alone:
null_mdl=lmer(BseriesCor~1+(1|as.factor(data$subject)),data,REML=F)
all_or_noneModel=lmer(BseriesCor~all_or_none+(1|as.factor(data$subject)),data,REML=F)
summary(all_or_noneModel)
Lin_mdl=lmer(BseriesCor~linearModel+(1|as.factor(data$subject)),data,REML=F)
summary(Lin_mdl)
onlyAccRTacc=lmer(BseriesCor~acc+RT_acc+(1|as.factor(data$subject)),data,REML=F)
summary(onlyAccRTacc)
onlyAcc=lmer(BseriesCor~acc+(1|as.factor(data$subject)),data,REML=F)
summary(onlyAcc)
onlyRTacc=lmer(BseriesCor~RT_acc+(1|as.factor(data$subject)),data,REML=F)
summary(onlyRTacc)

#comparison for the model with both acc and RT: THIS IS REPORTED
full_mdl=lmer(BseriesCor~all_or_none+acc+RT_acc+(1|as.factor(data$subject)),data,REML=F)
anova(all_or_noneModel,full_mdl)
anova(full_mdl,onlyAccRTacc)

# these are just checks - not reported
full_mdl=lmer(BseriesCor~linearModel+acc+RT_acc+(1|as.factor(data$subject)),data,REML=F)
anova(Lin_mdl,full_mdl)
anova(full_mdl,onlyAccRTacc)
#comparison for the model with only Acc:
full_mdl=lmer(BseriesCor~all_or_none+acc+(1|as.factor(data$subject)),data,REML=F)
anova(all_or_noneModel,full_mdl)
anova(full_mdl,onlyAcc)

#comparison for the model with only RT:
full_mdl=lmer(BseriesCor~all_or_none+RT_acc+(1|as.factor(data$subject)),data,REML=F)
anova(all_or_noneModel,full_mdl)
anova(full_mdl,onlyRTacc)


############## lCA1-lCA23DG connectivity ####################

#plots for Figure 2:
all_subjs_filename <- sprintf("%s/lCA1_lCA23DG_allItemsExcAD_AKcorrect_Univar_single_trials_withAccRT.txt", data_dir )
data <- read.table(all_subjs_filename, header=T)

data$fac_linearModel <- factor(data$linearModel) 

data_sum <- data %>%
  group_by(fac_linearModel) %>%
  summarise(mean_cor = mean(BseriesCor),
            SEM_cor = sd(BseriesCor)/sqrt(n()),
            n_data=n())
data_con_temp <- select(data,c(-all_or_none,-lCA1_activity,-lCA23DG_activity)) %>%
  mutate(linearModel = linearModel-3) %>%
  mutate(subj_con = linearModel*BseriesCor)

data_con  <- data_con_temp %>%
  group_by(subject) %>%
  summarise(subj_con = sum(subj_con))

## plot with data points:
# plt_ca23 <- ggplot(data_sum, aes(x = fac_linearModel, y = mean_cor, fill=fac_linearModel)) +
#   geom_col(color = "black") +
#   geom_point(data = data, aes(x = fac_linearModel, y = BseriesCor,fill=fac_linearModel), alpha = .4,
#              position = position_jitterdodge(jitter.width = .7,dodge.width = -.1)) +
#   geom_errorbar(aes(ymin = mean_cor - SEM_cor, ymax = mean_cor + SEM_cor),width=0,position = position_nudge(x = .2)) +
#   scale_fill_manual(values = blue_cols) +
#   scale_color_manual(values = "black") +
#   labs(x = '# of changes', y = 'Functional connectivity') + 
#   theme_classic() + theme(legend.position = "none")

## plot without data points:
plt_ca23 <- ggplot(data_sum, aes(x = fac_linearModel, y = mean_cor, fill=fac_linearModel)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = mean_cor - SEM_cor, ymax = mean_cor + SEM_cor),width=0) +
  scale_fill_manual(values = blue_cols) +
  scale_color_manual(values = "black") +
  labs(x = '# of changes', y = 'Functional connectivity') +
  theme_classic() + theme(legend.position = "none")



## plot the contrast:
plt_ca23_subj_con <- ggplot(data_con, aes(x = 1, y = subj_con, color = "1")) + 
  geom_hline(yintercept=0, colour="black", linetype = "dotted") +
  xlim(0.5,1.5) +
  geom_jitter(width = .08,size = 2) + 
  scale_color_manual(values = blue_cols[3]) +
  labs(y = 'linear trend') + 
  theme_classic() + theme(legend.position = "none",
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank())

# this arranges all plots:
grid.arrange(plt_ca23_subj_con,plt_ca23,plt_ent,plt_ent_subj_con,nrow=1)

#mixed effects linear models (Reported in Results):
#basic model:
results.nullModel=lmer(BseriesCor~1|as.factor(data$subject),data,REML=F)
#fit a linear model
results.linearModel=lmer(BseriesCor~linearModel+(1|as.factor(data$subject)),data,REML=F)
summary(results.linearModel)
anova(results.linearModel,results.nullModel)

#fit the all-or-none model
results.all_or_none=lmer(BseriesCor~as.factor(all_or_none)+(1|as.factor(data$subject)),data,REML=F)
summary(results.all_or_none)

#compare models
anova(results.all_or_none,results.linearModel)

#a model with both
results.all_or_noneAndLinear=lmer(BseriesCor~as.factor(all_or_none)+linearModel+(1|as.factor(data$subject)),data,REML=F)
summary(results.all_or_noneAndLinear)

#compare models
anova(results.all_or_noneAndLinear,results.linearModel)
anova(results.all_or_noneAndLinear,results.all_or_none)

############## control for univariate ############
# reported in Results
all_subjs_filename <- sprintf("%s/lCA1_lCA23DG_allItemsExcAD_AKcorrect_Univar_single_trials_withAccRT.txt", data_dir )

data <- read.table(all_subjs_filename, header=T)
#first, see that the linear works alone:
null_mdl=lmer(BseriesCor~1+(1|as.factor(data$subject)),data,REML=F)
linearModel=lmer(BseriesCor~linearModel+(1|as.factor(data$subject)),data,REML=F)
summary(all_or_noneModel)
anova(null_mdl,all_or_noneModel)

#fit a linear model - full
linearModel=lmer(BseriesCor~linearModel+scale(lCA1_activity)*scale(lCA23DG_activity)+(1|as.factor(data$subject)),data,REML=F)
summary(linearModel)

#fit null model - just activation
null_mdl=lmer(BseriesCor~scale(lCA1_activity)*scale(lCA23DG_activity)+(1|as.factor(data$subject)),data,REML=F)
# just to check - this is the all_none_model (not reported_
all_or_noneModel=lmer(BseriesCor~all_or_none+scale(lCA1_activity)*scale(lCA23DG_activity)+(1|as.factor(data$subject)),data,REML=F)
summary(all_or_noneModel)

# THIS IS REPORTED in the RESULTS:
anova(null_mdl,linearModel)
comp=anova(null_mdl,linearModel)
AIC_dif=comp$AIC[1]-comp$AIC[2]
BIC_dif=comp$BIC[1]-comp$BIC[2]
comp
AIC_dif
BIC_dif

###### ACC/RT actual data per participant ######
# THIS IS REPORTED IN SUPPLEMENTARY NOTE 2
all_subjs_filename <- sprintf("%s/lCA1_lCA23DG_allItemsExcAD_AKcorrect_Univar_single_trials_withAccRT.txt", data_dir )
data <- read.table(all_subjs_filename, header=T)
behav_filename <- sprintf("%s/Behavior_ExcAD_AKcorrect.txt", behav_data_dir )
bdata <- read.table(behav_filename, header=T)
#add some cols from bdata to data
data$acc=bdata$accuracy_rates
data$RT_acc=bdata$RT_accurate
data$RT_allitems=bdata$RT_all_items

#first, see that the linear works alone:
null_mdl=lmer(BseriesCor~1+(1|as.factor(data$subject)),data,REML=F)
linear_mdl=lmer(BseriesCor~linearModel+(1|as.factor(data$subject)),data,REML=F)
summary(linear_mdl)
onlyAccRTacc=lmer(BseriesCor~acc+RT_acc+(1|as.factor(data$subject)),data,REML=F)
summary(onlyAccRTacc)
full_mdl=lmer(BseriesCor~linearModel+acc+RT_acc+(1|as.factor(data$subject)),data,REML=F)
anova(linear_mdl,full_mdl)
anova(full_mdl,onlyAccRTacc)
anova(null_mdl,onlyAccRTacc)

onlyAccRTall=lmer(BseriesCor~acc+RT_allitems+(1|as.factor(data$subject)),data,REML=F)
full_mdl=lmer(BseriesCor~linearModel+acc+RT_allitems+(1|as.factor(data$subject)),data,REML=F)
summary(onlyAccRTall)
anova(full_mdl,onlyAccRTall)

#comparison for the model with only Acc:
full_mdl=lmer(BseriesCor~linearModel+acc+(1|as.factor(data$subject)),data,REML=F)
anova(linear_mdl,full_mdl)
anova(full_mdl,onlyAcc)

#comparison for the model with only RT:
full_mdl=lmer(BseriesCor~linearModel+RT_acc+(1|as.factor(data$subject)),data,REML=F)
anova(linear_mdl,full_mdl)
anova(full_mdl,onlyRTacc)

###################  PE_RSA  #######################################
## average per participant
data_dir <- sprintf("%s/results/filesForR",proj_dir)
all_subjs_filename <- sprintf("%s/Ticky_dataR.xlsx", data_dir )


## this controls for univariate with the match mismatch contrast:
#this uploads the univar data based on single trials
get_data <- read_excel(all_subjs_filename, sheet = "lCA1PE_matmis_Strials_uni2h3rd")

#plots the match-mismatch contrast (FIGURE 3):
plt_PE_subj_con <- ggplot(get_data, aes(x = 1, y = -sim, color = "1")) + 
  geom_hline(yintercept=0, colour="black") +
  xlim(0.5,1.5) +
  geom_jitter(width = .205,size = 3) + 
  scale_color_manual(values = green_cols[4]) +
  labs(y = 'match < mismatch') + 
  theme_classic() + theme(legend.position = "none",
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank())

# CONTROL FOR UNIVARIATE - REPORTED IN SUPP
#test the all_none_contrast
lm_mdl <- lm(sim ~ 1, data = get_data)  # lm reg model
summary(lm_mdl)

lm_univar_mdl <- lm(sim ~ 1+scale(univar_image_match_mismatch), data = get_data)  # lm reg model
summary(lm_univar_mdl)

lm_univar_mdl <- lm(sim ~ 1+scale(univar_cue_match_mismatch), data = get_data)  # lm reg model
summary(lm_univar_mdl)

lm_univar_mdl <- lm(sim ~ 1+scale(univar_cue_match_mismatch)+scale(univar_image_match_mismatch), data = get_data)  # lm reg model
summary(lm_univar_mdl)
confint(lm_univar_mdl)


## this controls for univariate all conditions, mixed-level model, and model comparisons
#this is like I treat the connectivity - but i didn't report this, I reported just the lm above, this works as well though:
data_dir <- sprintf("%s/results/filesForR",proj_dir)
##single trials
all_subjs_filename <- sprintf("%s/PE_lCA1_allItemsExcAD_AKcorrect_Univar_single_trials_2highThirds.txt", data_dir )
data <- read.table(all_subjs_filename, header=T)
#first, see that the all_or_none works alone: great, it does
null_mdl=lmer(PE_same~1+(1|as.factor(data$subject)),data,REML=F)
all_or_noneModel=lmer(PE_same~all_or_none+(1|as.factor(data$subject)),data,REML=F)
summary(all_or_noneModel)
anova(null_mdl,all_or_noneModel)

#fit a null model with only the activity:
null_mdl=lmer(PE_same~scale(image_activity)*scale(cue_activity)+(1|as.factor(data$subject)),data,REML=F)

#fit a linear model
linearModel=lmer(PE_same~linearModel+scale(image_activity)+scale(cue_activity)+(1|as.factor(data$subject)),data,REML=F)
summary(linearModel)

#fit all or none model
all_or_noneModel=lmer(PE_same~all_or_none+scale(image_activity)*scale(cue_activity)+(1|as.factor(data$subject)),data,REML=F)
summary(all_or_noneModel)
anova(null_mdl,all_or_noneModel)
comp=anova(null_mdl,linearModel)
AIC_dif=comp$AIC[1]-comp$AIC[2]
BIC_dif=comp$BIC[1]-comp$BIC[2]
comp
AIC_dif
BIC_dif


########### Prediction strength #########
#THIS IS REPORTED: control for univariate activity during the cue/image:
data_dir <- sprintf("%s/results/filesForR",proj_dir)
all_subjs_filename <- sprintf("%s/CA1_prediction_CA1_ent_con.xlsx", data_dir )

get_data <- read_excel(all_subjs_filename, sheet = "Pred2HThirds_PropGLM_conALL_v")
#get_data <- read_excel(all_subjs_filename, sheet = "PredMedSplit_PropGLM_conALL_v")
get_data$sc_cue_act=scale(get_data$cue_act)
get_data$sc_intact_image_act_same=scale(get_data$intact_image_act_same)
get_data$sc_intact_image_act_diff=scale(get_data$intact_image_act_diff)
null_mdl <- lm(sim ~ 1, data = get_data)  # lm reg model
summary(null_mdl)

#this is fine - THIS IS REPORTED:
lm_mdl <- lm(sim ~ 1+sc_cue_act+sc_intact_image_act_same, data = get_data)  # lm reg model
confint(lm_mdl)
summary(lm_mdl)

########## correlation prediction strength connectivity #######
## note that here - bc we're not estimating the intercept - no need to scale. I checked,
##scaling doesn't change anything
data_dir <- sprintf("%s/results/filesForR",proj_dir)
all_subjs_filename <- sprintf("%s/CA1_prediction_CA1_ent_con.xlsx", data_dir )

get_data <- read_excel(all_subjs_filename, sheet = "Pred2HThirds_PropGLM_conALL_v")
rcorr(get_data$sim,get_data$con)
zr <- .5*log((1+.4)/(1-.4))
n <- 19
r_se <- 1/sqrt(n-3)
curr_alpha <- .1 #to match one-tailed t-test
zcrit <- qnorm((1-curr_alpha/2), mean=0, sd=1)
CI_lower = zr - zcrit * r_se
CI_upper = zr + zcrit * r_se
#this is still in z, now we need to convert back to r:
CI_lower_r=tanh(CI_lower)
CI_upper_r=tanh(CI_upper)

#and - all of the above is identical to:
CIr(.40,19,.90)

# linear models:
lm_mdl <- lm(con ~ sim, data = get_data)  # lm reg model
summary(lm_mdl)

#THIS IS REPORTED (standardize): control for univariate activation for the similarity measure and for connectivity (as above)
lm_mdl <- lm(con ~ sim+cue_act+intact_image_act_same+ca1_match_mis_univar_single_trials+ent_match_mis_univar_single_trials, data = get_data)  # lm reg model
summary(lm_mdl)
a <- stdcoeff(lm_mdl)
confint(lm_mdl, level = .90)
#standardize for reporting standardize betas
sc_lm_mdl <- lm(scale(con) ~ scale(sim)+scale(cue_act)+scale(intact_image_act_same)+scale(ca1_match_mis_univar_single_trials)+scale(ent_match_mis_univar_single_trials), data = get_data)  # lm reg model
summary(sc_lm_mdl)
confint(sc_lm_mdl, level = .90)


#robust regression:
rlm_mdl <- rlm(con ~ sim+cue_act+intact_image_act_same+ca1_match_mis_univar_single_trials+ent_match_mis_univar_single_trials, data = get_data)  #  reg model
summary(rlm_mdl)
sc_rlm_mdl <- rlm(scale(con) ~ scale(sim)+scale(cue_act)+scale(intact_image_act_same)+scale(ca1_match_mis_univar_single_trials)+scale(ent_match_mis_univar_single_trials), data = get_data)  
summary(sc_rlm_mdl)
confint.default(object = sc_rlm_mdl, parm = "scale(sim)", level = 0.90)
t.stats=coef(summary(rlm_mdl))[2,'t value']
p.value = 2*pt(abs(t.stats), 13,lower.tail=FALSE)
p.value

######### Univariate ###########
# this is for supplementary figure 1, stats are also reported in the main text in the Results
all_spectral<-brewer.pal(9,"Greys")
myspec <-all_spectral[c(1,3,5,7,9)]
## grab it from the connectivity:
all_subjs_filename <- sprintf("%s/lCA1_lCA23DG_allItemsExcAD_AKcorrect_Univar_single_trials_withAccRT.txt", data_dir )
data <- read.table(all_subjs_filename, header=T)

data$fac_linearModel <- factor(data$linearModel) 
#get ent:
all_subjs_filename <- sprintf("%s/lCA1_lEnt_allItemsExcAD_AKcorrect_Univar_single_trials_withAccRT.txt", data_dir )
data_ent <- read.table(all_subjs_filename, header=T)
data$lEnt_activity <- data_ent$lEnt_activity

data_sum <- data %>%
  group_by(fac_linearModel) %>%
  summarise(mean_d1 = mean(lCA1_activity),
            SEM_d1 = sd(lCA1_activity)/sqrt(n()),
            mean_d23 = mean(lCA23DG_activity),
            SEM_d23 = sd(lCA23DG_activity)/sqrt(n()),
            mean_dEnt = mean(lEnt_activity),
            SEM_dEnt = sd(lEnt_activity)/sqrt(n()),
            n_data=n())

data_con_temp <- select(data,c(-all_or_none,-BseriesCor,-lCA23DG_activity)) %>%
  mutate(linearModel = linearModel-3) %>%
  mutate(subj_con = linearModel*lCA1_activity)

data_con  <- data_con_temp %>%
  group_by(subject) %>%
  summarise(subj_con = sum(subj_con))

t.test(data_con$subj_con,mu=0)

## plot without data points:
plt_b1 <- ggplot(data_sum, aes(x = fac_linearModel, y = mean_d1, fill=fac_linearModel)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = mean_d1 - SEM_d1, ymax = mean_d1 + SEM_d1),width=0) +
  scale_fill_manual(values = myspec) +
  scale_color_manual(values = "black") +
  labs(x = '# of changes', y = 'Univariate activation (t)') +
  theme_classic() + theme(legend.position = "none")

## plot the contrast:
plt_subj_con <- ggplot(data_con, aes(x = 1, y = subj_con, color = "1")) + 
  geom_hline(yintercept=0, colour="black", linetype = "dotted") +
  xlim(0.5,1.5) +
  geom_jitter(width = .1,size = 2) + 
  scale_color_manual(values = myspec[3]) +
  labs(y = 'linear increase') + 
  theme_classic() + theme(legend.position = "none",
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank())

plt_b2 <- ggplot(data_sum, aes(x = fac_linearModel, y = mean_d23, fill=fac_linearModel)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = mean_d23 - SEM_d23, ymax = mean_d23 + SEM_d23),width=0) +
  scale_fill_manual(values = blue_cols) +
  scale_color_manual(values = "black") +
  labs(x = '# of changes', y = 'Univariate activation (t)') +
  theme_classic() + theme(legend.position = "none")

plt_b3 <- ggplot(data_sum, aes(x = fac_linearModel, y = mean_dEnt, fill=fac_linearModel)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = mean_dEnt - SEM_dEnt, ymax = mean_dEnt + SEM_d23),width=0) +
  scale_fill_manual(values = green_cols) +
  scale_color_manual(values = "black") +
  labs(x = '# of changes', y = 'Univariate activation (t)') +
  theme_classic() + theme(legend.position = "none")

grid.arrange(plt_subj_con,plt_b1,plt_b2,plt_b3,nrow=1)
