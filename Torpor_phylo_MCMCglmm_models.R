## Code for paper titled:
#"Torpid hummingbirds- energy savings in temperate and tropical sites"
## Paper authors: Anusha Shankar*, Rebecca J Schroeder*, 
# Susan M Wethington, Catherine H Graham, Donald R Powers
## *Equal authors
## Code by: Anusha Shankar, github/nushiamme; 
# contact: anusha<dot>shankar<at>stonybrook<dot>edu or nushiamme<at>gmail<dot>com for questions about code
## Thank you Liliana Davalos and Paige Copenhaver-Parry for help with these models!
## MCMCglmm models, accounting for both the phylogenetic structure and 
# the repeated-measures per species

### Contents
## Setup, read files in, format data
## Table 3: Model for Probability of entry into torpor ~ mass
## Supp Table S1: Summary of multiple NEE models
## Table 3: Plot of best model of nighttime energy expenditure 
# (i.e. Mass-corrected NEE ~ torpor duration)
## Table 3, Supp Table S2: Rewarming analyses

library(MCMCglmm)
library(nlme)
library(ape)
library(geiger) # for treedata() function
library(caper)
library(coda) # only for autocorr function

#### Setup ####
setwd("C:\\Users\\nushi\\Dropbox\\Hummingbird energetics\\July2018\\Data\\")
#GFU wd
#setwd("/Users/anshankar/Dropbox/Hummingbird energetics/Feb2018/Data/")

## Read in torpor data file
torpor <- read.csv("Torpor_individual_summaries_NewRER_Mar2020.csv") #Torpor data file, each row is an individual

## Read in McGuire et al. 2014 hummingbird phylogeny
tree<-read.tree("hum294.tre")

## Generic plot theme
my_theme <- theme_classic(base_size = 30) + 
  theme(panel.border = element_rect(colour = "black", fill=NA))

#### Adding columns ####
#Converting NA's in Hours_torpid column into 0's.
torpor$Hours2 <- torpor$Hours_torpid
torpor$Hours2[is.na(torpor$Hours2==TRUE)] <- 0

#Making a column for energy savings as a proportion of energy expenditure 
#savings = ((normothermic-torpor)/normothermic)*100
torpor$savings <- 100-torpor$Percentage_avg_varRER
torpor$savings2 <- 100-torpor$Percentage_avg_varRER
torpor$savings2[is.na(torpor$savings2)] <- 0

## To make column where savings is a binned, ordinal value.
torpor$savings_quantile <- with(torpor, factor(
  findInterval(savings2, c(-Inf,
                       unique(quantile(savings2, probs=seq(0,1, by=0.25))), Inf)), 
  labels=c("1","2","3","4")
))
torpor$savings_quantile <- as.numeric(torpor$savings_quantile)
torpor$savings_quantile2 <- as.factor(torpor$savings_quantile)

## Rewarming kJ column - remove NA's
torpor$kJ_rewarming <- torpor$kJ_RER0.71_rewarming_BeforeOvershoot
torpor$kJ_rewarming[is.na(torpor$kJ_rewarming==TRUE)] <- 0

#### Phylogenetic components - prune tree ####
#Replace tip names in the tree with those in torpor_subset database
## If using an updated tree check tree names again. Using a hard replace here because there aren't
#too many to be replaced.
#(To show tip names, use: tree$tip.label)
tree$tip.label[1]<-"WNJ"
tree$tip.label[12]<-"TBH"
tree$tip.label[92]<-"EMB"
tree$tip.label[93]<-"FBB"
tree$tip.label[95]<-"GCB"
tree$tip.label[219]<-"BBLH"
tree$tip.label[163]<-"BLUH"
tree$tip.label[156]<-"RIHU"

tips<-data.frame(levels(torpor$Species))
colnames(tips) <- "tips"
rownames(tips)<-tips$tips

#match tree to data, prune tree, species names should be in rownnames of "data" 
tre1<-treedata(tree, tips)$phy
#To check that the reationships between species in the trimmed tree look right
plot(tre1) 

#### Models ####
## Now, to run Bayesian models with repeated measures per species (i.e. multiple individuals per species), 
#we setup an inverse matrix and set up a prior
#Using a Bayesian rather than a maximum likelihood model because with an ML model we could include 
#repeated measures, OR we could include a phylogenetic structure. 
#But to get a hierarchy, with both a phylogeny and then repeated measures 
#within the phylogeny, we need turn the phylogeny into an inverse matrix
inv.phylo<-inverseA(tre1,nodes="TIPS",scale=TRUE)
#set up a prior for a phylogenetic mixed model
#Setting priors to be very uninformative
prior<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002)) 
#run the hierarchical phyogenetic model, the name of the species 
#(repeated across rows of observations) 

plot(torpor$Mass, torpor$Tornor)

#### Model for probability of entry into torpor as a function of mass ####
## Table 3
## Frequency converted into bernoulli individual-level torpor-not
##Try library(phylolm)
mfreq2 <- MCMCglmm(Tornor~Mass, random=~Species, family='ordinal',
                   ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
                   verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mfreq2) ## Table 3
plot(mfreq2) ## Look at the model parameter distribution
plot(mfreq2$VCV)

p1<-predict(mfreq2, interval="confidence")

plot(torpor$Mass,torpor$Tornor)
lines(p1[,1]~torpor$Mass,col="red")
lines(p1[,2]~torpor$Mass,col="red",lty=2)
lines(p1[,3]~torpor$Mass,col="red", lty=2)

x <- torpor$Mass
y <- torpor$Tornor
mod.x <- mean(mfreq2$Sol[,1])
mod.y <- mean(mfreq2$Sol[,2])

## Plot probability of torpor use vs Mass
ggplot(torpor, aes(Mass,Tornor)) + geom_point(size=3) + my_theme +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), 
                          se = T) +
  xlab("Mass (g)") + ylab("Probability of torpor use")

#### Nighttime energy expenditure MCMCglmm models (stepwise) ####
## All these model results are summarized in Supp. Table S1
## All these models used to use mass-corrected nighttime energy expenditure as the dependent variable
## But removed Mass-corrected NEE as a function of Mass on December 8. 2019
## And instead using raw NEE in kJ with mass just as a co-variate
mNEE_mass<-MCMCglmm(NEE_kJ_variableRER~Mass, random=~Species, 
              ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
              verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_mass)

## Nee (kJ) as a function of duration and mass
mNEE_dur_mass<-MCMCglmm(NEE_kJ_variableRER~Hours2+Mass, random=~Species, 
              ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
              verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_dur_mass)

## As a function of duration of torpor - this is the best DIC+most parsimonious model - model used in paper
mNEE_dur <-MCMCglmm(NEE_kJ_variableRER~Hours2, random=~Species, 
                   ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
                   verbose=F, nitt = 5e6, thin = 1000)
summary(mNEE_dur)
plot(mNEEmass_dur)


## Of min chamber temperature
mNEE_tc<-MCMCglmm(NEE_kJ_variableRER~Tc_min_C, random=~Species, 
              ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
              verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_tc)

## This is the best model in terms of lowest DIC value (by <2 points), but not most parsimonious
## Duration + min chamber temperature
mNEE_dur_tc <-MCMCglmm(NEE_kJ_variableRER~Hours2+Tc_min_C, random=~Species, 
              ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
              verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_dur_tc) ## Table 3
par(mar = rep(2, 4))
plot(mNEE_dur_tc) ## Look at model parameter distribution

## Mass + Duration + min chamber temperature
mNEE_mass_dur_tc <-MCMCglmm(NEE_kJ_variableRER~Mass+Hours2+Tc_min_C, random=~Species, 
             ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
             verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_mass_dur_tc)

## As a function of hourly energy savings (as a quantile because otherwise 0's swamp differences
# between non-zeros)
mNEE_sav <- MCMCglmm(NEE_kJ_variableRER~savings_quantile, random=~Species, 
                ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
                verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_sav)
plot(mNEE_sav)

## Duration + min temp + savings
mNEE_dur_tc_sav <- MCMCglmm(NEE_kJ_variableRER~Hours2+Tc_min_C+savings_quantile, random=~Species, 
                ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
                verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_dur_tc_sav)

## Mass + duration + min temp + savings
mNEE_mass_dur_tc_sav <- MCMCglmm(NEE_kJ_variableRER~Mass+Hours2+Tc_min_C+savings_quantile, random=~Species, 
                            ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
                            verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_mass_dur_tc_sav)

## NEE ~ 
### Full model including rewarming, Oct 2017
mNEE_full <-MCMCglmm(NEE_kJ_variableRER~Mass+Hours2+Tc_min_C+savings_quantile+
                       kJ_rewarming, 
             random=~Species, ginverse = list(Species=inv.phylo$Ainv), 
             prior=prior, data=torpor, verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_full)
#plot(mNEE_full) 

## Without any phylogenetic corrections- shows that results have an inflated significance when 
#phylo corrections are not done
mNEE_nophylo <-MCMCglmm(NEE_mass~Mass+Hours2+Tc_min_C+savings,
                        data=torpor[torpor$Hours2!=0,])
summary(mNEE_nophylo)

#### Rewarming ####
## First model for rewarming, only taking mass (g) into account
mrewarm <- MCMCglmm(kJ_rewarming~Mass, 
                       random=~Species, family='gaussian',
                       ginverse=list(Species=inv.phylo$Ainv), prior=prior, 
                       data=torpor[torpor$Torpid_not=="T",],
                       verbose=F,nitt=5e6, thin=1000)
summary(mrewarm) ## In Supp. Table S2
plot(mrewarm)

## Second rewarming model, including mass (g) and chamber temperature (deg C) during rewarming
mrewarm_tc <- MCMCglmm(kJ_rewarming~Mass+Rewarming_Tc, 
                       random=~Species, family='gaussian',
                       ginverse=list(Species=inv.phylo$Ainv), prior=prior, 
                       data=torpor[torpor$Torpid_not=="T",],
                       verbose=F,nitt=5e6, thin=1000)
summary(mrewarm_tc) ## Table 3 and in Supp Table S2
plot(mrewarm_tc) ## Look at model results plotted

anova(lm(torpor$kJ_rewarming_mass~torpor$Mass,na.action = na.omit))


library(phylolm)
m1.fit<-phyloglm(formula = Tornor ~ Mass, torpor, tre1, method = "logistic_IG10",
                 boot = 1000,start.beta = c(1,-10), start.alpha = 0.01, full.matrix = T)
summary(m1.fit)
coef(m1.fit)
vcov(m1.fit)
