#install.packages('ggplot2')
library(ggplot2)
#install.packages('dplyr')
library(dplyr)
#install.packages('forcats')
library(forcats)
#install.packages('rockchalk')
library(rockchalk)

# READ IN DATA ------------------------------------------------------------

# read in meta-analysis dataset
# this is just the dataset that was deposited in Open Research Exeter when the last manuscript was published (https://ore.exeter.ac.uk/repository/handle/10871/22805)
metadata_2018<-read.csv('data/metaanalysis dataset.csv')

#for Calderon 1982 the number exposed (2)/unexposed (11) wasn't recorded in the dataset, so we add it manually (for forest plot)
metadata_2018$numberexposed[metadata_2018$studyid=='Calderon 1982']<-2
metadata_2018$numberofunexposed[metadata_2018$studyid=='Calderon 1982']<-11

#check number of studies in the meta-analysis dataset (should be 19, as reported in the paper)
length(unique(metadata_2018$studyid))

#it's actually 27 - but this is because several studies were excluded 
# Harrington 1993 is excluded because no odds ratios reported or calculable for the bathing/non bathing comparison
# 6 excluded because of microbe specific outcomes analysed separately - Begier 1997 (echovirus); 

#list of excluded studies
excludes<-c('Alexander 1992','Begier 2008','Brown 1987','Charoenca 1995','Dwight 2004','Fewtrell 1994','Fleming 2004','Gammie 1997','Haile 1999','Harder-Lauridsen','Harding 2015','Harrington 1993','Hoque 2002','Ihekweazu 2006','Lepesteur 2006','Morens 1994','Nelson 1997','Reed 2006','Roy 2004','Soraas 2013')
#remove the 7 excluded studies
metadata_2018<-metadata_2018[!metadata_2018$studyid%in%excludes,]


#calc standard deviation of ORs
metadata_2018$or_se<-(metadata_2018$uor-metadata_2018$lor)/3.92

#calc standard deviation of ORs
metadata_2018$or_sd<-metadata_2018$or_se*sqrt(metadata_2018$studysize)

#add an effect size ID
metadata_2018$esid<-1:nrow(metadata_2018)

#generate an escalc dataframe as a workaround because need escalc class for metafor
metadata_2018 <- metafor::escalc(measure="OR", 
                                 ai=numberofexposedcases, bi=numberofexposednoncases,
                                 ci=numberofunexposedcaes, di=numberofunexposednoncases,
                                 #n1i=numberexposed, n2i=numberofunexposed,
                                 data=metadata_2018, add=T)

#replace generated effect sizes with the old ones (workaround cos need class of escalc) - in real one will only replace those for which we have marginal ORs
metadata_2018$yi<-metadata_2018$logOR
metadata_2018$vi<-metadata_2018$logORse^2 #square for variance

#remove NAs (missing yi or vi)
metadata_2018<-metadata_2018[!is.na(metadata_2018$yi),]
metadata_2018<-metadata_2018[!is.na(metadata_2018$vi),]

#add an observation ID
metadata_2018$observation_id<-as.integer(unlist(tapply(metadata_2018$studyid,metadata_2018$studyid,function(x){1:length(x)})))

#remove symptoms of dysphagia (removed from original analysis) and where case definitions not reported
metadata_2018<-metadata_2018[!metadata_2018$symptom%in%c('Dysphagia','Ear infection (case definition not reported)','Gastrointestinal infection (case definition not reported)'),]

#rename symptoms to align with original
metadata_2018$symptom[metadata_2018$symptom=='Any symptom of infection']<-'Symptoms of any illness'
metadata_2018$symptom[metadata_2018$symptom=='Ear infection (sensitive case definition)']<-'Symptoms of ear ailments (sensitive case definitions)'
metadata_2018$symptom[metadata_2018$symptom=='Earache (single symptom case definition)']<-'Ear ache'
metadata_2018$symptom[metadata_2018$symptom=='Ear discharge (single symptom case definition)']<-'Ear discharge'
metadata_2018$symptom[metadata_2018$symptom=='Gastrointestinal infection (case definition requires at least one of a list of symptoms)']<-'Symptoms of gastrointestinal illness (sensitive case definitions)'
metadata_2018$symptom[metadata_2018$symptom=='Gastrointestinal infection (case definition requires two or more symptoms)']<-'Symptoms of gastrointestinal illness (specific case definitions)'

#make factor
metadata_2018$symptom<-as.factor(metadata_2018$symptom)
#make factor
metadata_2018$symptom_simplified<-as.factor(metadata_2018$symptom)

#create a vector of simplified symptom categories as intended to be used in the update, by merging previous categories based on site of putative infection

#merge any symptoms
metadata_2018$symptom_simplified<-rockchalk::combineLevels(fac = metadata_2018$symptom_simplified, 
                                                           levs = c("Symptoms of any illness"), 
                                                           newLabel="Any")
#merge ear symptoms
metadata_2018$symptom_simplified<-rockchalk::combineLevels(fac = metadata_2018$symptom_simplified, 
                                                           levs = c("Ear ache","Ear discharge","Symptoms of ear ailments (sensitive case definitions)"), 
                                                           newLabel="Ear")

#merge gastro symptoms
metadata_2018$symptom_simplified<-rockchalk::combineLevels(fac = metadata_2018$symptom_simplified, 
                                                           levs = c("Stomach ache","Vomiting","Diarrhoea","Nausea","Symptoms of gastrointestinal illness (sensitive case definitions)","Symptoms of gastrointestinal illness (specific case definitions)"), 
                                                           newLabel="Gastrointestinal")


levels(metadata_2018$symptom_simplified)

#note that the simplified symptom categories are equivalent to the previous 'health outcome categories' that were in the original dataset, I just didn't see them before. So either can be used
tapply(metadata_2018$symptom_simplified,
       metadata_2018$healthoutcomecategory,
       function(x){print(unique(x))})

#check number of studies (should be 19 as reported in paper)
length(unique(metadata_2018$studyid))

# STUDY-LEVEL EFFECTS (FOREST PLOT) -------------------------------------------------------------

metadata_2018_anyexp<-metadata_2018[metadata_2018$exposureanalysis%in%c('any','both'),]

#main model
model <- metafor::rma.mv(yi, vi,
                         random = list(~1|studyid, ~1|studyid/observation_id),
                         data = metadata_2018_anyexp,
                         control=list(rel.tol=1e-8))

model

exp(cbind(model$b,model$ci.lb,model$ci.ub))

#model with random effect removed for calculating I2 via Jackson apporach (see below)
model_noranf <- metafor::rma.mv(yi, vi,
                                data = metadata_2018_anyexp,
                                control=list(rel.tol=1e-8))

#use Jackson approach to calculate I2 - https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
I2<-c(100 * (vcov(model)[1,1] - vcov(model_noranf)[1,1]) / vcov(model)[1,1])

agg <- metafor::aggregate.escalc(metadata_2018_anyexp, cluster=studyid, V=vcov(model, type="obs"), addk=TRUE, var.names=c("yi","vi"))

k <- nrow(agg)

res <- metafor::rma(yi, vi, method="EE", data=agg, digits=3)
res

#FOREST PLOT WITH BACK-TRANSFORMED ORs FOR EACH STUDY AND OVERALL OR

grDevices::tiff('figures/newapproach_forestplot.tiff', res=300, units='in', width=10, height=8)

par(mar = c(5, 4, 0, 0)) # Reduce top margin

sav<-metafor::forest.rma(res,
                     atransf = exp,
                     slab=studyid,
                     efac=c(0,1),
                     xlim=c(-10,6),
                     alim=log(c(0.13,20)),
                     at=log(exp(c(-2,-1,0,1,2,3))),
                     ilab=cbind(paste(round(agg$numberexposed)),paste(round(agg$numberofunexposed))),
                     ilab.xpos=c(-5,-3),
                     header=TRUE, 
                     xlab='Odds ratio',
                     mlab="Overall", 
                     cex=1,
                     col='black',
                     shade='zebra')

text(sav$ilab.xpos[1:2], rep(k,2)+c(1.7), c("Exposed to  \n seawater","Not exposed to  \n seawater"), cex=0.75)
#segments(sav$ilab.xpos[1]-0.22, k+1.75, sav$ilab.xpos[2]+0.13, k+1.75)
text(mean(sav$ilab.xpos[1:2]), k+3, "Average numbers \n in each study arm*",cex=1)

text(sav$xlim[1], -1.75, pos=4, 
     bquote(paste(I^2, " = ", .(round(I2)), "%", "; ",
                  tau^2, " = ", .(metafor::fmtx(model$tau2, digits=2)), "; ",
                  chi^2, " = ", .(metafor::fmtx(model$QE, digits=2)),
                  ", df = ", .(model$k - model$p), ", ",
                  .(metafor::fmtp(model$QEp, digits=2, pname="P", add0=TRUE, equal=TRUE)))))

mtext('*Effects from multiple symptom categories within studies are aggregated', side=1, line=3.5, at=-10, adj=0, cex=1)

dev.off()

# UPDATED APPROACH, SIMPLIFIED SYMPTOM CATEGORIES -------------------------------------------

#metadata_2018_anyexp$symptom[grep('Ear discharge \\(single symptom case definition\\)|Ear infection \\(sensitive case definition\\)|Earache \\(single symptom case definition\\)',metadata_2018_anyexp$symptom)]<-'Symptoms of ear ailments (sensitive case definitions)'

#main model_H2
model_H2 <- metafor::rma.mv(yi, vi,
                            mod = ~ 0 + symptom_simplified, 
                            random = list(~1|studyid, ~1|studyid/observation_id),
                            data = metadata_2018_anyexp,
                            control=list(rel.tol=1e-8))

model_H2
#main model_H2
model_H2_noranf <- metafor::rma.mv(yi, vi,
                                   mod = ~ 0 + symptom_simplified, 
                                   data = metadata_2018_anyexp,
                                   control=list(rel.tol=1e-8))

#use Jackson approach to calculate I2 - https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
I2<-c(100 * (vcov(model_H2)[1,1] - vcov(model_H2_noranf)[1,1]) / vcov(model_H2)[1,1])

#use Wolfgang's approach to calculate pseudo-R2 of 3 variance components (2 sigmas and tau) - https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-October/000247.html; https://stackoverflow.com/questions/22356450/getting-r-squared-from-a-mixed-effects-multilevel-model_H2-in-metafor; https://gist.github.com/wviechtb/6fbfca40483cb9744384ab4572639169
R2<-100*(max(0,(sum(model$sigma2, model$tau2) -
                  sum(model_H2$sigma2,model_H2$tau2)) / 
               sum(model$sigma2,model$tau2)))

estimates<-exp(cbind(model_H2$b,model_H2$ci.lb,model_H2$ci.ub))
estimates<-round(estimates, 2)

unique(metadata_2018_anyexp$symptom_simplified)
neworder<-c('Any',
            'Ear',
            'Gastrointestinal')

neworder_rev<-rev(neworder)

estimates_rev<-estimates[match(neworder_rev,levels(as.factor(metadata_2018_anyexp$symptom_simplified))),]

#neworder_rev<-neworder
#estimates_rev<-estimates

orchy<-orchaRd::orchard_plot(model_H2, 
                             group='studyid', 
                             mod='symptom_simplified', 
                             tree.order = neworder_rev,
                             xlab='Odds ratio', 
                             angle=0, k.pos=3, legend.pos = 'bottom.left')+
  expand_limits(x=c(0,4))+
  scale_y_continuous(limits = c(-4, 4)) +
  theme(axis.text.y = ggplot2::element_text(size = 12, colour ="black",
                                            hjust = 0,
                                            vjust = 0.5,
                                            angle=0))+
  annotate("text", x=1:3, y=-3, label=paste(estimates_rev[,1],' (',estimates_rev[,2],',',estimates_rev[,3],')', sep=''))+
  annotate("text", x=3.5, y=-3, label="Odds ratio (95% CI)",fontface = 2)

orchy

ggsave('figures/updatedapproach_simplecats_orchardplot.tiff', plot=last_plot(), width=8, height=6)

# UPDATED APPROACH, ORIGINAL SYMPTOM CATEGORIES ----------

#metadata_2018_anyexp$symptom[grep('Ear discharge \\(single symptom case definition\\)|Ear infection \\(sensitive case definition\\)|Earache \\(single symptom case definition\\)',metadata_2018_anyexp$symptom)]<-'Symptoms of ear ailments (sensitive case definitions)'

#main model_H2
model_H2 <- metafor::rma.mv(yi, vi,
                            mod = ~ 0 + symptom, 
                            random = list(~1|studyid, ~1|studyid/observation_id),
                            data = metadata_2018_anyexp,
                            control=list(rel.tol=1e-8))

model_H2
#main model_H2
model_H2_noranf <- metafor::rma.mv(yi, vi,
                                   mod = ~ 0 + symptom, 
                                   data = metadata_2018_anyexp,
                                   control=list(rel.tol=1e-8))

#use Jackson approach to calculate I2 - https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
I2<-c(100 * (vcov(model_H2)[1,1] - vcov(model_H2_noranf)[1,1]) / vcov(model_H2)[1,1])

#use Wolfgang's approach to calculate pseudo-R2 of 3 variance components (2 sigmas and tau) - https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-October/000247.html; https://stackoverflow.com/questions/22356450/getting-r-squared-from-a-mixed-effects-multilevel-model_H2-in-metafor; https://gist.github.com/wviechtb/6fbfca40483cb9744384ab4572639169
R2<-100*(max(0,(sum(model$sigma2, model$tau2) -
                  sum(model_H2$sigma2,model_H2$tau2)) / 
               sum(model$sigma2,model$tau2)))

estimates<-exp(cbind(model_H2$b,model_H2$ci.lb,model_H2$ci.ub))
estimates<-round(estimates, 2)

unique(metadata_2018_anyexp$symptom)
neworder<-c('Symptoms of any illness',
            'Symptoms of ear ailments (sensitive case definitions)',
            'Ear ache',
            'Ear discharge',
            'Symptoms of gastrointestinal illness (sensitive case definitions)',
            'Diarrhoea',
            'Nausea',
            'Stomach ache',
            'Vomiting',
            'Symptoms of gastrointestinal illness (specific case definitions)')

neworder_rev<-rev(neworder)

estimates_rev<-estimates[match(neworder_rev,levels(as.factor(metadata_2018_anyexp$symptom))),]

#neworder_rev<-neworder
#estimates_rev<-estimates

orchy<-orchaRd::orchard_plot(model_H2, 
                             group='studyid', 
                             mod='symptom', 
                             tree.order = neworder_rev,
                             xlab='Odds ratio', 
                             angle=0, k.pos=5.5, legend.pos = 'bottom.left')+
  expand_limits(x=c(-1,12))+
  scale_y_continuous(limits = c(-8, 6)) +
  theme(axis.text.y = ggplot2::element_text(size = 8, colour ="black",
                                            hjust = 0,
                                            vjust = 0.5,
                                            angle=0))+
  annotate("text", x=1:10, y=-6.5, label=paste(estimates_rev[,1],' (',estimates_rev[,2],',',estimates_rev[,3],')', sep=''))+
  annotate("text", x=11, y=-6.5, label="Odds ratio (95% CI)",fontface =2)

orchy

ggsave('figures/updatedapproach_origcats_orchardplot.tiff', plot=last_plot(), width=12, height=6)

# ORIGINAL APPROACH, SIMPLIFIED SYMPTOM CATEGORIES ----------

#metadata_2018_anyexp$symptom[grep('Ear discharge \\(single symptom case definition\\)|Ear infection \\(sensitive case definition\\)|Earache \\(single symptom case definition\\)',metadata_2018_anyexp$symptom)]<-'Symptoms of ear ailments (sensitive case definitions)'

#main model_H2
model_H2 <- metafor::rma.mv(yi, vi,
                            mod = ~ 0 + symptom_simplified, 
                            random = list(~symptom_simplified|studyid),
                            struct = 'DIAG',
                            data = metadata_2018_anyexp,
                            control=list(rel.tol=1e-8))

model_H2
#main model_H2
model_H2_noranf <- metafor::rma.mv(yi, vi,
                                   mod = ~ 0 + symptom_simplified, 
                                   data = metadata_2018_anyexp,
                                   control=list(rel.tol=1e-8))

#use Jackson approach to calculate I2 - https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
I2<-c(100 * (vcov(model_H2)[1,1] - vcov(model_H2_noranf)[1,1]) / vcov(model_H2)[1,1])

#use Wolfgang's approach to calculate pseudo-R2 of 3 variance components (2 sigmas and tau) - https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-October/000247.html; https://stackoverflow.com/questions/22356450/getting-r-squared-from-a-mixed-effects-multilevel-model_H2-in-metafor; https://gist.github.com/wviechtb/6fbfca40483cb9744384ab4572639169
R2<-100*(max(0,(sum(model$sigma2, model$tau2) -
                  sum(model_H2$sigma2,model_H2$tau2)) / 
               sum(model$sigma2,model$tau2)))

estimates<-exp(cbind(model_H2$b,model_H2$ci.lb,model_H2$ci.ub))
estimates<-round(estimates, 2)

unique(metadata_2018_anyexp$symptom_simplified)
neworder<-c('Any',
            'Ear',
            'Gastrointestinal')

neworder_rev<-rev(neworder)

estimates_rev<-estimates[match(neworder_rev,levels(as.factor(metadata_2018_anyexp$symptom_simplified))),]

#neworder_rev<-neworder
#estimates_rev<-estimates

orchy<-orchaRd::orchard_plot(model_H2, 
                             group='studyid', 
                             mod='symptom_simplified', 
                             tree.order = neworder_rev,
                             xlab='Odds ratio', 
                             angle=0, k.pos=3, legend.pos = 'bottom.left')+
  expand_limits(x=c(0,4))+
  scale_y_continuous(limits = c(-4, 4)) +
  theme(axis.text.y = ggplot2::element_text(size = 12, colour ="black",
                                            hjust = 0,
                                            vjust = 0.5,
                                            angle=0))+
  annotate("text", x=1:3, y=-3, label=paste(estimates_rev[,1],' (',estimates_rev[,2],',',estimates_rev[,3],')', sep=''))+
  annotate("text", x=3.5, y=-3, label="Odds ratio (95% CI)",fontface = 2)

orchy

ggsave('figures/origapproach_simplecats_orchardplot.tiff', plot=last_plot(), width=8, height=6)

# ORIGINAL APPROACH, ORIGINAL SYMPTOM CATEGORIES -------------------------------------------

#metadata_2018_anyexp$symptom[grep('Ear discharge \\(single symptom case definition\\)|Ear infection \\(sensitive case definition\\)|Earache \\(single symptom case definition\\)',metadata_2018_anyexp$symptom)]<-'Symptoms of ear ailments (sensitive case definitions)'

#main model_H2
model_H2 <- metafor::rma.mv(yi, vi,
                         mod = ~ 0 + symptom, 
                         random = list(~symptom|studyid),
                         struct='DIAG',
                         data = metadata_2018_anyexp,
                         control=list(rel.tol=1e-8))

#main model_H2
model_H2_noranf <- metafor::rma.mv(yi, vi,
                         mod = ~ 0 + symptom, 
                         data = metadata_2018_anyexp,
                         control=list(rel.tol=1e-8))

#use Jackson approach to calculate I2 - https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
I2<-c(100 * (vcov(model_H2)[1,1] - vcov(model_H2_noranf)[1,1]) / vcov(model_H2)[1,1])

#use Wolfgang's approach to calculate pseudo-R2 of 3 variance components (2 sigmas and tau) - https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2017-October/000247.html; https://stackoverflow.com/questions/22356450/getting-r-squared-from-a-mixed-effects-multilevel-model_H2-in-metafor; https://gist.github.com/wviechtb/6fbfca40483cb9744384ab4572639169
R2<-100*(max(0,(sum(model$sigma2, model$tau2) -
                  sum(model_H2$sigma2,model_H2$tau2)) / 
               sum(model$sigma2,model$tau2)))

estimates<-exp(cbind(model_H2$b,model_H2$ci.lb,model_H2$ci.ub))
estimates<-round(estimates, 2)

unique(metadata_2018_anyexp$symptom)
neworder<-c('Symptoms of any illness',
            'Symptoms of ear ailments (sensitive case definitions)',
            'Ear ache',
            'Ear discharge',
            'Symptoms of gastrointestinal illness (sensitive case definitions)',
            'Diarrhoea',
            'Nausea',
            'Stomach ache',
            'Vomiting',
            'Symptoms of gastrointestinal illness (specific case definitions)')
            
neworder_rev<-rev(neworder)

estimates_rev<-estimates[match(neworder_rev,levels(as.factor(metadata_2018_anyexp$symptom))),]

#neworder_rev<-neworder
#estimates_rev<-estimates

orchy<-orchaRd::orchard_plot(model_H2, 
                             group='studyid', 
                             mod='symptom', 
                             tree.order = neworder_rev,
                             xlab='Odds ratio', 
                             angle=0, k.pos=5.5, legend.pos = 'bottom.left')+
  expand_limits(x=c(-1,12))+
  scale_y_continuous(limits = c(-8, 6)) +
  theme(axis.text.y = ggplot2::element_text(size = 8, colour ="black",
                                             hjust = 0,
                                             vjust = 0.5,
                                             angle=0))+
  annotate("text", x=1:10, y=-6.5, label=paste(estimates_rev[,1],' (',estimates_rev[,2],',',estimates_rev[,3],')', sep=''))+
  annotate("text", x=11, y=-6.5, label="Odds ratio (95% CI)",fontface =2)

orchy
  
ggsave('figures/origapproach_origcats.tiff', plot=last_plot(), width=12, height=6)
