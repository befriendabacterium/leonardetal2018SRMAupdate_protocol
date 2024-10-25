# PREAMBLE ------------------------------------------------------------

#N.B. This code is intended only to demonstrate the core meta-analytical techniques to be used in the new analysis and their visualisation, and is not intended to be definitive
# Whilst the new review will focus on comparisons of estimate cases of infection in people exposed to the lowest and higher levels of polluted seawater, this demo uses the data from the Leonard et al (2018) meta-analysis, which focuses on comparisons of those not exposed and those exposed to seawater per se
# As such, many of the pre-processing steps are solely for the purpose of wrangling the already part-processed Leonard et al 2018 data into a format suitable for demonstrating the analysis and visualisation we intend

#install.packages('ggplot2')
library(ggplot2)
#install.packages('dplyr')
library(dplyr)
#install.packages('forcats')
library(forcats)
#install.packages('rockchalk')
library(rockchalk)
#install.packages('effectsize')
library(effectsize)
#install.packages('patchwork')
library(patchwork)

# PRE-PROCESSING ------------------------------------------------------------

## READ IN DATA ------------------------------------------------------------

# this is just the dataset that was deposited in Open Research Exeter when the last manuscript was published (https://ore.exeter.ac.uk/repository/handle/10871/22805)
leonardetal2018_metaanalysis_df<-read.csv('data/metaanalysis dataset.csv')

## MANUALLY MAKE SOME MINOR CORRECTIONS/EDITS ------------------------------------------------------------

#for Calderon 1982 the number exposed (2)/unexposed (11) wasn't recorded in the dataset, so we add it manually (for forest plot)
leonardetal2018_metaanalysis_df$numberexposed[leonardetal2018_metaanalysis_df$studyid=='Calderon 1982']<-2
leonardetal2018_metaanalysis_df$numberofunexposed[leonardetal2018_metaanalysis_df$studyid=='Calderon 1982']<-11

leonardetal2018_metaanalysis_df[which(colnames(leonardetal2018_metaanalysis_df)=='numberofunexposedcaes')]

#correct typo on one of column names
leonardetal2018_metaanalysis_df<-rename(leonardetal2018_metaanalysis_df, 'numberofunexposedcases'='numberofunexposedcaes')

#rename 'Enteric' health outcome category to 'Gastrointestinal' as is more commonly used term
leonardetal2018_metaanalysis_df$healthoutcomecategory[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Enteric']<-'Gastrointestinal'

## DISTILL DATASET TO THE STUDIES AND EXPOSURE-OUTCOME COMBINATIONS INCLUDED IN THE LEONARD ET AL 2018 META-ANALYSIS  ------------------------------------------------------------

#check number of studies in the meta-analysis dataset (should be 19, as reported in the paper)
length(unique(leonardetal2018_metaanalysis_df$studyid))

#it's actually 27 - but this is because several studies were excluded 
# Harrington 1993 is excluded because no odds ratios reported or calculable for the bathing/non bathing comparison
# 6 excluded because of microbe specific outcomes analysed separately - Begier 1997 (echovirus); 

#list of excluded studies
excludes<-c('Alexander 1992','Begier 2008','Brown 1987','Charoenca 1995','Dwight 2004','Fewtrell 1994','Fleming 2004','Gammie 1997','Haile 1999','Harder-Lauridsen','Harding 2015','Harrington 1993','Hoque 2002','Ihekweazu 2006','Lepesteur 2006','Morens 1994','Nelson 1997','Reed 2006','Roy 2004','Soraas 2013')

#remove the 7 excluded studies
leonardetal2018_metaanalysis_df<-leonardetal2018_metaanalysis_df[!leonardetal2018_metaanalysis_df$studyid%in%excludes,]

#remove symptoms of dysphagia (removed from original analysis) and where case definitions not reported
leonardetal2018_metaanalysis_df<-leonardetal2018_metaanalysis_df[!leonardetal2018_metaanalysis_df$symptom%in%c('Dysphagia','Ear infection (case definition not reported)','Gastrointestinal infection (case definition not reported)'),]

#remove head immersion from main analysis as in Leonard et al 2019
leonardetal2018_metaanalysis_df<-leonardetal2018_metaanalysis_df[leonardetal2018_metaanalysis_df$exposureanalysis%in%c('any','both'),]

## PROCESSING OF EODDS RATIOS  ------------------------------------------------------------

# The metafor R package requires an 'escalc' dataframe to work with, which is the output of the 'escalc' function for calculating estimates of the association/effect
# However, we're started with the standard dataframe from Leonard et al (2018), which already includes estimates of the association/effect (odds ratios) calculated in STATA
# Therefore, to get this working in metafor, we will use the workaround solution of re-calculating the effect sizes using escalc, which coerces the dataframe to 'escalc' class. 
# However, we will replacing these with the original estimates of the association/effect for this demonstrative analysis/visualisation
# N.B. In the actual review this won't be an issue as we will start afresh

#coerce the dataframe to escalc via re-calculating odds ratios (odds ratios and their variance will be outputted as 'yi' and 'vi') and appending them to the dataframe
leonardetal2018_metaanalysis_df <- metafor::escalc(measure="OR", 
                                 ai=numberofexposedcases, bi=numberofexposednoncases,
                                 ci=numberofunexposedcases, di=numberofunexposednoncases,
                                 #n1i=numberexposed, n2i=numberofunexposed,
                                 data=leonardetal2018_metaanalysis_df, add=T)

#calc standard error of original ORs 
leonardetal2018_metaanalysis_df$or_se<-(leonardetal2018_metaanalysis_df$uor-leonardetal2018_metaanalysis_df$lor)/3.92
#calc standard deviation of original ORs 
leonardetal2018_metaanalysis_df$or_sd<-leonardetal2018_metaanalysis_df$or_se*sqrt(leonardetal2018_metaanalysis_df$studysize)

#replace recalculated effect sizes with the original ones (note these are a mix of adjusted and unadjusted, preferring the former where available
leonardetal2018_metaanalysis_df$yi<-leonardetal2018_metaanalysis_df$logOR
leonardetal2018_metaanalysis_df$vi<-leonardetal2018_metaanalysis_df$logORse^2 #square for variance

#add an effect size ID
leonardetal2018_metaanalysis_df$esid<-1:nrow(leonardetal2018_metaanalysis_df)

#add an observation ID
leonardetal2018_metaanalysis_df$observation_id<-as.integer(unlist(tapply(leonardetal2018_metaanalysis_df$studyid,leonardetal2018_metaanalysis_df$studyid,function(x){1:length(x)})))

## CONVERSION OF RISK RATIOS TO ODDS RATIOS  ------------------------------------------------------------

# Next we're going to convert the Leonard et al 2018 odds ratios and their variances (now held in  'yi' and 'vi', respectively) to risk ratios, because they're easier to interpret

# For this we first define the baserate for each estimate (i.e. case rate in unexpoed population)
leonardetal2018_metaanalysis_df$baserate<-leonardetal2018_metaanalysis_df$numberofunexposedcases/leonardetal2018_metaanalysis_df$numberofunexposed

#calculate mean baserate per health outcome category
meanbaserate_perhealthoutcome<-
       tapply(leonardetal2018_metaanalysis_df$baserate,
       leonardetal2018_metaanalysis_df$healthoutcomecategory,
       mean, na.rm=T)

#temp fill in NAs with mean from all studies
#leonardetal2018_metaanalysis_df$baserate[is.na(leonardetal2018_metaanalysis_df$baserate)]<-mean(leonardetal2018_metaanalysis_df$baserate, na.rm = T)

#convert 'any' infection odds ratio point estimates to risk ratios 
leonardetal2018_metaanalysis_df$yi[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Any']<-metafor::transf.lnortorr(leonardetal2018_metaanalysis_df$yi[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Any'], 
                                                                                       pc=meanbaserate_perhealthoutcome[names(meanbaserate_perhealthoutcome)=='Any'])
#convert 'any' infection odds ratio variance estimates to risk ratios 
leonardetal2018_metaanalysis_df$vi[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Any']<-metafor::transf.lnortorr(leonardetal2018_metaanalysis_df$vi[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Any'], 
                                                                                       pc=meanbaserate_perhealthoutcome[names(meanbaserate_perhealthoutcome)=='Any'])

#convert 'any' infection odds ratio point estimates to risk ratios 
leonardetal2018_metaanalysis_df$yi[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Ear']<-metafor::transf.lnortorr(leonardetal2018_metaanalysis_df$yi[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Ear'], 
                                                                                       pc=meanbaserate_perhealthoutcome[names(meanbaserate_perhealthoutcome)=='Ear'])
#convert 'any' infection odds ratio variance estimates to risk ratios 
leonardetal2018_metaanalysis_df$vi[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Ear']<-metafor::transf.lnortorr(leonardetal2018_metaanalysis_df$vi[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Ear'], 
                                                                                       pc=meanbaserate_perhealthoutcome[names(meanbaserate_perhealthoutcome)=='Ear'])

#convert 'any' infection odds ratio point estimates to risk ratios 
leonardetal2018_metaanalysis_df$yi[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Gastrointestinal']<-metafor::transf.lnortorr(leonardetal2018_metaanalysis_df$yi[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Gastrointestinal'], 
                                                                                       pc=meanbaserate_perhealthoutcome[names(meanbaserate_perhealthoutcome)=='Gastrointestinal'])
#convert 'any' infection odds ratio variance estimates to risk ratios 
leonardetal2018_metaanalysis_df$vi[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Gastrointestinal']<-metafor::transf.lnortorr(leonardetal2018_metaanalysis_df$vi[leonardetal2018_metaanalysis_df$healthoutcomecategory=='Gastrointestinal'], 
                                                                                       pc=meanbaserate_perhealthoutcome[names(meanbaserate_perhealthoutcome)=='Gastrointestinal'])

#remove exposure-outcome combinations for which estimates of the association/effect could not be calculated
leonardetal2018_metaanalysis_df<-leonardetal2018_metaanalysis_df[!is.na(leonardetal2018_metaanalysis_df$yi),]
leonardetal2018_metaanalysis_df<-leonardetal2018_metaanalysis_df[!is.na(leonardetal2018_metaanalysis_df$vi),]

## FINAL CHECKS  ------------------------------------------------------------

#check number of studies again (should be 19 as reported in paper)
length(unique(leonardetal2018_metaanalysis_df$studyid))

#set seed so colour palette generation reproduces previous version in tidy_metadata()
set.seed(123)
#(re)generate global palette of colours for each study using microshades package
study_palette<-sample(unlist(microshades::microshades_cvd_palettes),
                      length(unique(leonardetal2018_metaanalysis_df$studyid)), replace=F)

#assign the study names to the palette
names(study_palette)<-levels(as.factor(leonardetal2018_metaanalysis_df$studyid))

#sample study palette colours based on studyID
leonardetal2018_metaanalysis_df$study_colour<-study_palette[match(leonardetal2018_metaanalysis_df$studyid,names(study_palette))]

# MODEL FOR OVERALL EFFECT (INTERCEPT-ONLY MODEL) -------------------------------------------------------------

#main model
model <- metafor::rma.mv(yi, vi,
                         random = list(~1|studyid, ~1|studyid/observation_id),
                         data = leonardetal2018_metaanalysis_df,
                         control=list(rel.tol=1e-8))

model

exp(cbind(model$b,model$ci.lb,model$ci.ub))

#model with random effect removed for calculating I2 via Jackson apporach (see below)
model_noranf <- metafor::rma.mv(yi, vi,
                                data = leonardetal2018_metaanalysis_df,
                                control=list(rel.tol=1e-8))

#use Jackson approach to calculate I2 - https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
I2<-c(100 * (vcov(model)[1,1] - vcov(model_noranf)[1,1]) / vcov(model)[1,1])


#estimates<-exp(cbind(model$b,model$ci.lb,model$ci.ub)) #for log odds
estimates<-cbind(model$b,model$ci.lb,model$ci.ub)
estimates<-round(estimates, 2)
estimates_rev<-estimates

# ORCHARD PLOT FOR OVERALL EFFECT (INTERCEPT-ONLY MODEL) -------------------------------------------------------------

# code to return the study study_palette in later code
microshades_cvd<-tapply(leonardetal2018_metaanalysis_df$study_colour,leonardetal2018_metaanalysis_df$studyid, unique)
#subset the microshades per study palette to just the studies included in the current model/sub-dataset
cbpl_temp<-microshades_cvd[match(unique(leonardetal2018_metaanalysis_df$studyid), names(microshades_cvd))]

#make orchard plot
orchy_H1a<-orchaRd::orchard_plot(model, 
                             group='studyid', 
                             #mod='symptom_simplified', 
                             #tree.order = neworder_rev,
                             xlab='Risk ratio', 
                             flip=F, colour = T,
                             angle=0, k.pos=5.5, legend.pos = 'none')+
  ggplot2::scale_fill_manual(values = cbpl_temp) +
  ggplot2::scale_colour_manual(values = cbpl_temp)+
  scale_x_discrete(labels = 'Overall') +
  scale_y_continuous(limits = c(-1, 7), breaks=c(0:7)) +
  theme(axis.text.x = ggplot2::element_text(size = 15, colour ="black",
                                            hjust = 0.5,
                                            vjust = 0,
                                            angle=0))+
  theme(axis.text.y = ggplot2::element_text(size = 12, colour ="black",
                                            hjust = 0,
                                            vjust = 0.5,
                                            angle=0))+
  annotate("text", x=1, y=-0.5, label=paste(estimates_rev[,1],' (',estimates_rev[,2],',',estimates_rev[,3],')', sep=''))
  #annotate("text", x=1, y=-3, label="Odds ratio (95% CI)",fontface = 2)

#remove the 'no association/effect' reference line that orchard plot puts at zero, because for the risk ratio it should be at 1 (we'll add another for 0, see below)
orchy_H1a$layers[[2]]<-NULL # 

#add reference lines in manually
orchy_H1a<-orchy_H1a+
           ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "black")+ #add 'no association/effect' reference line at 1
           ggplot2::geom_hline(yintercept = 0, linetype = 1, colour = "black") #add reference line at 0 representing minimum RR (cannot be less than 0)

#check final output
orchy_H1a

#save final output
ggsave('figures/orchard_H1a.tiff', plot=last_plot(), width=8, height=6)

# UPDATED APPROACH, SIMPLIFIED SYMPTOM CATEGORIES -------------------------------------------

#leonardetal2018_metaanalysis_df$symptom[grep('Ear discharge \\(single symptom case definition\\)|Ear infection \\(sensitive case definition\\)|Earache \\(single symptom case definition\\)',leonardetal2018_metaanalysis_df$symptom)]<-'Symptoms of ear ailments (sensitive case definitions)'

#main model_H2
model_H2 <- metafor::rma.mv(yi, vi,
                            mod = ~ 0 + healthoutcomecategory, 
                            random = list(~1|studyid, ~1|studyid/observation_id),
                            data = leonardetal2018_metaanalysis_df,
                            control=list(rel.tol=1e-8))

model_H2
#main model_H2
model_H2_noranf <- metafor::rma.mv(yi, vi,
                                   mod = ~ 0 + healthoutcomecategory, 
                                   data = leonardetal2018_metaanalysis_df,
                                   control=list(rel.tol=1e-8))

#use Jackson approach to calculate I2 - https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
I2<-c(100 * (vcov(model_H2)[1,1] - vcov(model_H2_noranf)[1,1]) / vcov(model_H2)[1,1])

#use Wolfgang's approach to calculate pseudo-R2 of 3 variance components (2 sigmas and tau) - https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-October/000247.html; https://stackoverflow.com/questions/22356450/getting-r-squared-from-a-mixed-effects-multilevel-model_H2-in-metafor; https://gist.github.com/wviechtb/6fbfca40483cb9744384ab4572639169
R2<-100*(max(0,(sum(model$sigma2, model$tau2) -
                  sum(model_H2$sigma2,model_H2$tau2)) / 
               sum(model$sigma2,model$tau2)))

#estimates<-exp(cbind(model_H2$b,model_H2$ci.lb,model_H2$ci.ub)) #for log odds
estimates<-cbind(model_H2$b,model_H2$ci.lb,model_H2$ci.ub)
estimates<-round(estimates, 2)

unique(leonardetal2018_metaanalysis_df$symptom_simplified)
neworder<-c('Any',
            'Ear',
            'Gastrointestinal')

neworder_rev<-neworder

estimates_rev<-estimates[match(neworder_rev,levels(as.factor(leonardetal2018_metaanalysis_df$healthoutcomecategory))),]

#make orchard plot
orchy_H1b<-orchaRd::orchard_plot(model_H2, 
                             group='studyid', 
                             mod='healthoutcomecategory', 
                             tree.order = neworder_rev,
                             xlab='Odds ratio',
                             flip=F,
                             colour=T,
                             angle=0, k.pos=5.5, legend.pos = 'top.right')+
  ggplot2::scale_fill_manual(values = cbpl_temp) +
  ggplot2::scale_colour_manual(values = cbpl_temp)+
  #expand_limits(x=c(1,3))+
  scale_y_continuous(limits = c(-1, 7), breaks=c(0:7)) +
  theme(axis.text.x = ggplot2::element_text(size = 15, colour ="black",
                                            hjust = 0.5,
                                            vjust = 0,
                                            angle=0))+
  theme(axis.text.y = ggplot2::element_text(size = 12, colour ="black",
                                            hjust = 0,
                                            vjust = 0.5,
                                            angle=0))+
  annotate("text", x=1:3, y=-0.5, label=paste(estimates_rev[,1],' (',estimates_rev[,2],',',estimates_rev[,3],')', sep=''))
  #annotate("text", x=3.5, y=-2.5, label="Odds ratio (95% CI)",fontface = 2)

#remove the 'no association/effect' reference line that orchard plot puts at zero, because for the risk ratio it should be at 1 (we'll add another for 0, see below)
orchy_H1b$layers[[2]]<-NULL # 

#add reference lines in manually
orchy_H1b<-orchy_H1b+
  ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "black")+ #add 'no association/effect' reference line at 1
  ggplot2::geom_hline(yintercept = 0, linetype = 1, colour = "black") #add reference line at 0 representing minimum RR (cannot be less than 0)

#check final output
orchy_H1b

#save final output
ggsave('figures/orchard_H1a&H1b.tiff', plot=last_plot(), width=8, height=4)
ggsave('figures/orchard_H1b.tiff', plot=last_plot(), width=8, height=6)

# COMBINED PLOT -----------------------------------------------------------

#combine plots using the '+' operators from the patchwork R package
orchy_H1a+
  orchy_H1b+ggpubr::rremove('ylab')+ggpubr::rremove('y.text')+
  plot_layout(widths = c(1,3))

#save final output
ggsave('figures/orchard_H1a&H1b.tiff', plot=last_plot(), width=8, height=5)