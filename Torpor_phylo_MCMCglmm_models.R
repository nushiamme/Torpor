## Code for paper titled:
#"Torpid hummingbirds- energy savings in temperate and tropical sites"
## Paper authors: Anusha Shankar*, Rebecca J Schroeder*, 
# Susan M Wethington, Catherine H Graham, Donald R Powers
## *Equal authors
## Code by: Anusha Shankar, github/nushiamme
## Contact: anusha<dot>shankar<at>stonybrook<dot>edu for questions about code
## Thank you Liliana Davalos for help with this!
## Started Nov 23, 2016
## MCMCglmm models, accounting for both the phylogenetic structure and the repeated-measures per
# species

### Contents
## Setup, read files in, format data
## Table 3 and Supp Figure S4: Model for Probability of entry into torpor ~ mass
## NEE models; summarized in Supp Table S2
## Supp Figure S4: Plot of best model of nighttime energy expenditure 
# (i.e. Mass-corrected NEE ~ torpor duration + min chamber temperature)

library(MCMCglmm)
library(nlme)
library(ape)
library(geiger) # for treedata() function
library(caper)
library(lattice) # for xyplot function

#### Setup ####
setwd("C:\\Users\\ANUSHA\\Dropbox\\Hummingbird energetics\\Submission_Jul2017\\Data\\")

## Read in torpor data file
torpor <- read.csv("Torpor_individual_summaries.csv") #Torpor data file, each row is an individual

## Read in McGuire et al. 2014 hummingbird phylogeny
tree<-read.tree("hum294.tre")

#### Adding columns ####

#Mass-correct nighttime energy expenditure - done here to allow for 
#freedom in changing allometric exponent
torpor$NEE_MassCorrected<- torpor$NEE_kJ/(torpor$Mass^(2/3))

#Converting NA's in Hours_torpid column into 0's.
torpor$Hours2 <- torpor$Hours_torpid
torpor$Hours2[is.na(torpor$Hours2==TRUE)] <- 0

#Making a column for energy savings as a proportion of energy expenditure 
#savings = ((normothermic-torpor)/normothermic)*100
torpor$savings <- 100-torpor$Percentage_avg
torpor$savings2 <- 100-torpor$Percentage_avg
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
torpor$kJ_rewarming2 <- torpor$kJ_rewarming_BeforeOvershoot
torpor$kJ_rewarming2[is.na(torpor$kJ_rewarming2==TRUE)] <- 0

#### Phylogenetic components - prune tree ####
#Replace tip names in the tree with those in torpor database
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
tree$tip.label[156]<-"MAHU"

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
prior<-list(G=list(G1=list(V=1,nu=1)),R=list(V=1,nu=1)) 
#run the hierarchical phyogenetic model, the name of the species 
#(repeated across rows of observations) 

#### Model for probability of entry into torpor as a function of mass ####
## Table 3 and Supp. Figure S4
## Frequency converted into bernoulli individual-level torpor-not
mfreq1 <- MCMCglmm(Tornor~Mass, random=~Species, family='categorical',
                   ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
                   verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mfreq1) ## Table 3
plot(mfreq1) ## Supp. Figure S4

#### Nighttime energy expenditure MCMCglmm models (stepwise) ####
## All these model results are summarized in Supp. Table S2
## All these models use mass-corrected nighttime energy expenditure as the dependent variable
## Mass-corrected NEE as a function of Mass
mNEE_mass<-MCMCglmm(NEE_MassCorrected~Mass, random=~Species, 
              ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
              verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_mass)

## As a function of duration of torpor
mNEE_dur<-MCMCglmm(NEE_MassCorrected~Hours2, random=~Species, 
              ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
              verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_dur)

## Of min chamber temperature
mNEE_tc<-MCMCglmm(NEE_MassCorrected~Tc_min_C, random=~Species, 
              ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
              verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_tc)

## This is the best model; lowest DIC value
## Duration + min chamber temperature
mNEE_dur_tc <-MCMCglmm(NEE_MassCorrected~Hours2+Tc_min_C, random=~Species, 
              ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
              verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_dur_tc) ## Table 4
par(mar = rep(2, 4))
plot(mNEE_dur_tc) ## Supp Figure S8

## Mass + Duration + min chamber temperature
mNEE_mass_dur_tc <-MCMCglmm(NEE_MassCorrected~Mass+Hours2+Tc_min_C, random=~Species, 
             ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
             verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_mass_dur_tc)

## As a function of hourly energy savings (as a quantile because otherwise 0's swamp differences
# between non-zeros)
mNEE_sav <- MCMCglmm(NEE_MassCorrected~savings_quantile, random=~Species, 
                ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
                verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_sav)

## Duration + min temp + savings
mNEE_dur_tc_sav <- MCMCglmm(NEE_MassCorrected~Hours2+Tc_min_C+savings_quantile, random=~Species, 
                ginverse = list(Species=inv.phylo$Ainv), prior=prior, data=torpor, 
                verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_dur_tc_sav)

## NEE ~ 
## This was the full model used until July 2017
mNEE_mass_dur_tc_sav <-MCMCglmm(NEE_MassCorrected~Mass+Hours2+Tc_min_C+savings_quantile,
             random=~Species, ginverse = list(Species=inv.phylo$Ainv), 
             prior=prior, data=torpor, verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_mass_dur_tc_sav)
par(mar = rep(2, 4))
plot(mNEE_mass_dur_tc_sav)
autocorr(mNEE_mass_dur_tc_sav) #To check how autocorrelated the variables are

### Full model including rewarming, Oct 2017
mNEE_full <-MCMCglmm(NEE_MassCorrected~Mass+Hours2+Tc_min_C+savings_quantile+
               kJ_rewarming2, 
             random=~Species, ginverse = list(Species=inv.phylo$Ainv), 
             prior=prior, data=torpor, verbose=FALSE, nitt = 5e6, thin = 1000)
summary(mNEE_full)
plot(mNEE_full) 

## Without any phylogenetic corrections- shows that results have an inflated significance when 
#phylo corrections are not done
mNEE_nophylo <-MCMCglmm(NEE_MassCorrected~Mass+Hours2+Tc_min_C+savings,
                        data=torpor[torpor$Hours2!=0,])
summary(mNEE_nophylo)

#### Rewarming ####
## First model for rewarming, only taking mass (g) into account
mrewarm <- MCMCglmm(kJ_rewarming_BeforeOvershoot~Mass, 
                       random=~Species, family='gaussian',
                       ginverse=list(Species=inv.phylo$Ainv), prior=prior, 
                       data=torpor[torpor$Torpid_not=="T",],
                       verbose=F,nitt=5e6, thin=1000)
summary(mrewarm) ## In Supp. Table S3
plot(mrewarm)

## Second rewarming model, including mass (g) and chamber temperature (deg C) during rewarming
mrewarm_tc <- MCMCglmm(kJ_rewarming_BeforeOvershoot~Mass+Rewarming_Tc, 
                       random=~Species, family='gaussian',
                       ginverse=list(Species=inv.phylo$Ainv), prior=prior, 
                       data=torpor[torpor$Torpid_not=="T",],
                       verbose=F,nitt=5e6, thin=1000)
summary(mrewarm_tc) ## Table 5 and in Supp Table S3
plot(mrewarm_tc) ## Supp. Figure S10
autocorr(mrewarm_tc) ## Checking autocorrelation between variables
