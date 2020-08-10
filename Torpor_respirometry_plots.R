## Code for paper titled:
#"Hummingbird torpor in context: duration, more than temperature, 
# is key to nighttime energy savings"
## Paper authors: Anusha Shankar*, Rebecca J Schroeder*, 
# Susan M Wethington, Catherine H Graham, Donald R Powers
## *Equal authors
## Code by: Anusha Shankar, github/nushiamme
## Contact: anusha<dot>shankar<at>stonybrook<dot>edu or nushiamme<at>gmail<dot>com for questions about code
# This script contains code for torpor/respirometry plots

## Contents
## Setup, add columns, define background plot functions and objects
## Figure 3a: VO2 for broad-billed hummingbirds showing inflection point at 14-15 degC, from laboratory data
## Figure 3b: VO2 for broad-billed hummingbirds from field data, with greater variability than laboratory data
## Figure 4: To test the effect of rewarming costs on total nighttime energy expenditure,
# we plot total nighttime energy expenditure including rewarming (small black dots), and excluding rewarming 
# (coloured dots), as a function of individual mass
## Appendix Figures:
## Appendix Figure A2: Sample graph of energy expenditure (Joules) of a green-crowned brilliant over the course of a night
## Appendix Figure A3: Site-wise measures of a. duration (hours) of torpor and b. nighttime energy expenditure (kJ)
## Appendix Figure A4a: Torpor duration vs. minimum chamber temperature (Tc min) for the night. 
## Appendix Figure A4b: Torpor duration vs. average hourly energy savings in torpor relative to normothermy.
## Appendix Figure A5: Average hourly torpid energy savings relative to normothermy for all individuals that used torpor across all sites
## Appendix Figure A6: The probability of entering torpor is a negative binomial function of the mass of the individual. This is a graphical depiction of model 1 in Table 3 of the main paper


#### Setup #### 
## Libraries and reading in data
library(ggplot2)
library(reshape)
library(polynom)
library(gridExtra)
library(plyr) # only for 'revalue' function


## setwd and read in file
setwd("C:\\Users\\nushi\\Dropbox\\Hummingbird energetics\\July2018\\Data")

## Read in files
torpor <- read.csv("Torpor_individual_summaries.csv") # Torpor summaries per individual
bblh_tnz <- read.csv("BroadBill.csv") ## For Figure 3, BBLH minimum body temperature
bblh_VO2_temp_hourly <- read.csv("BBLH_hourly_VO2_field.csv") ## BBLH hourly temperature and VO2
gcb_0720 <- read.csv("E14_0720_GCB_no_bsln_Rewarm.csv")

# General functions ####
## Saving standard theme  
my_theme <- theme_classic(base_size = 25) + 
  theme(axis.title.y = element_text(color = "black", vjust = 2))

## Theme with slightly smaller font
my_theme2 <- my_theme + theme_classic(base_size = 15)

## Theme for Fig A2 
my_theme_blank <- theme_classic(base_size = 30) + 
  theme(axis.title.y = element_text(vjust = 2),
        panel.border = element_blank())

## Template axis labels
Tc.xlab <- expression(atop(paste("Chamber temperature (", degree,"C)"))) # for chamber temperature
NEE_lab <- 'Nighttime energy expenditure (kJ)' # for mass-corrected nighttime energy expenditure
VO2_lab <- expression(paste(VO[2]~mL~O[2]/min))
Tc_min.xlab <- expression(atop(paste("Minimum chamber temperature (", degree,"C)")))

## Function for adding a regression equation to graphs
## (Where y= table$column for y in the equation and x= table$column for x)
lm_eqn <- function(y, x){
  m <- lm(y ~ x);
  eq <- substitute(italic(y) == 
                     a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[[1]], digits = 2), 
                        b = format(coef(m)[[2]], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

# Function to return sample sizes
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

# Adding/Formatting columns ####
## In BBLH tnz data set
bblh_tnz$N_T <- factor(bblh_tnz$N_T, levels=c('T', 'N')) # Reorder levels BBLH tnz dataset

## In torpor dataset
## Make species a sensible order, for just species that used torpor - used in Supp Fig S7
torpor$Species2 <- factor(torpor$Species,
                          levels = c('BBLH','RIHU', 'BLUH', 'GCB','FBB', 'EMB', 'TBH', 'WNJ'), ordered = T)
## Make a column for Scientific code, as Ecology wants it
torpor$Species_sciname <- torpor$Species2
torpor$Species_sciname <- revalue(torpor$Species_sciname, c("BBLH"="CYLA", "RIHU"="EUFU", "BLUH"="LACL", "GCB"="HEJA", "FBB"="HERU", "EMB"="HEIM",
                                  "TBH"="PHSY", "WNJ"="FLME"))

## Rewarming kJ column - remove NA's
torpor$kJ_rewarming <- torpor$kJ_RER0.71_rewarming_BeforeOvershoot
torpor$kJ_rewarming[is.na(torpor$kJ_rewarming==TRUE)] <- 0


#For Supp Fig A3, order and expand site names
torpor$Site_full <- torpor$Site
torpor$Site_full <- factor(torpor$Site_full, levels=c('HC', 'SC', 'SWRS', 'MQ','SL'))
levels(torpor$Site_full) <- c("Harshaw", "Sonoita", "Southwestern Research Station", "Maquipucuna", "Santa Lucia")

## Savings column to convert percentage energy expended in torpor relative to 
# normothermy, into savings relative to normothermy - used in Supp Fig S7

torpor$savings <- 100-torpor$Percentage_avg_varRER

bblh_controlled_torpor <- bblh_tnz[bblh_tnz$N_T=="T",] # Subsetting just torpor measurements from controlled lab data


## for Fig. A4b, 
# create Hours_torpid2 column and change NA's to 0's for lm analyses; keep original Hours_torpid column with NA's
torpor$Hours_torpid2 <- torpor$Hours_torpid
torpor$Hours_torpid2[is.na(torpor$Hours_torpid2)] <- 0

#### Plots ####
## Figure 3a
# Broadbill torpid energy expenditure with temperature, showing an inflection point
# at 14-15 degC. Only using subset of the dataset- just torpid values.
# These measurements were taken under controlled conditions in 5 degC temperature steps, 
# separately from all the other torpor measurements
polyn.formula_lab <- bblh_controlled_torpor$VO2_all ~ poly(bblh_controlled_torpor$Temp_C, 2, raw = TRUE)
m_pol_lab <- lm(polyn.formula_lab, bblh_controlled_torpor)
polyn.eq_lab <- as.character(signif(as.polynomial(coef(m_pol_lab)), 2))
polyn.text_lab <- paste(gsub("x", "~italic(x)", polyn.eq_lab, fixed = TRUE),
                    paste("italic(R)^2",  
                          format(summary(m_pol_lab)$r.squared, digits = 2), 
                          sep = "~`=`~"),
                    sep = "~~~~")
bblh_VO2_temp <- ggplot(bblh_controlled_torpor, aes(Temp_C, VO2_all)) + 
  geom_point(size=3) + my_theme + 
  geom_smooth(method = "lm", se = T, 
              formula = y~ poly(x,2,raw=T), 
              colour = "black") + 
  annotate(geom = "text", x = 8, y = .2, label = polyn.text_lab, 
           family = "serif", hjust = 0, parse = TRUE, size=10) +
  ylab(VO2_lab) + xlab(Tc.xlab) + xlim(3,25) +
  theme(axis.title.x = element_blank())
plot(bblh_VO2_temp)

## Figure 3b, field measurements of torpor - natural temperature and photoperiod
## Torpid bblh VO2 vs temp
bblh_torpid <- na.omit(bblh_VO2_temp_hourly[bblh_VO2_temp_hourly$Torpid_not=="Torpid",])
polyn.formula <- bblh_torpid$VO2 ~ poly(bblh_torpid$Temperature, 2, raw = TRUE)
m_pol <- lm(polyn.formula, bblh_torpid)
polyn.eq <- as.character(signif(as.polynomial(coef(m_pol)), 2))
polyn.text <- paste(gsub("x", "~italic(x)", polyn.eq, fixed = TRUE),
                    paste("italic(R)^2",  
                          format(summary(m_pol)$r.squared, digits = 2), 
                          sep = "~`=`~"),
                    sep = "~~~~")
ggplot(bblh_torpid, aes(Temperature, VO2)) + 
  geom_point(aes(col=Bird_no), size=3, alpha=0.7) +
  geom_smooth(method = "lm", se = T, 
              formula = y~ poly(x,2,raw=T), 
              colour = "black") + 
  annotate(geom = "text", x = 8, y = .2, label = polyn.text, 
           family = "serif", hjust = 0, parse = TRUE, size=10) +
  ylim(-0.01,0.3) + guides(col=guide_legend(title="Individual ID")) +
  my_theme + theme(legend.key.height = unit(2, 'lines')) + xlim(3,25) +
  xlab(Tc.xlab) + ylab(VO2_lab)


## Old Figure 4 (variable RER)
## Nighttime energy expenditure with and without rewarming costs - not mass-corrected
rewarming_mass_NEE_temp <- ggplot(torpor[torpor$Torpid_not=="T",], aes(Mass, NEE_variableRER_minus_rewarming_kJ)) + 
  my_theme2 + geom_point(aes(col=Rewarming_Tc), size=3, alpha=0.5) +
  scale_colour_gradient(low="blue", high="white", name="Chamber \n temperature") +
  geom_point(aes(Mass, NEE_kJ_variableRER), col='black')  +
  geom_segment(aes(x=Mass,y=NEE_kJ_variableRER,xend=Mass,yend=NEE_variableRER_minus_rewarming_kJ, col=Rewarming_Tc),
               arrow = arrow(length = unit(0.01, "npc"),type = "closed"), alpha=0.5) +
  xlab("Mass") + ylab("Nighttime energy expenditure (kJ)") 
rewarming_mass_NEE_temp

## Final Fig. 4; reviewer's suggestion: leave mass out and plot NEE vs. rewarming
rewarming_NEE <- ggplot(torpor[torpor$Torpid_not=="T",], 
                        aes(kJ_rewarming, NEE_kJ_variableRER)) + 
  my_theme + geom_point(aes(fill=Rewarming_Tc), pch=21, size=4, alpha=0.9, col='black') +
  scale_fill_gradient(low="blue", high="white", name="Chamber \ntemperature") +
  xlab("Rewarming (kJ)") + ylab("Nighttime \n energy expenditure (kJ)") 
rewarming_NEE

rewarming_mass <- ggplot(torpor[torpor$Torpid_not=="T",], aes(Mass, kJ_rewarming)) + 
  my_theme + geom_point(aes(fill=Rewarming_Tc), pch=21, size=4, alpha=0.9, col='black') +
  scale_fill_gradient(low="blue", high="white", name="Chamber \ntemperature") +
  geom_smooth(method='lm') +
  xlab ("Mass (g)") + ylab("\n Rewarming (kJ)")
rewarming_mass

grid.arrange(rewarming_NEE, rewarming_mass, ncol=1, nrow=2)

### Figures for Appendix 
## Fig A2: Sample torpor night; Whole night of energy expenditure for one individual
gcb_0720$Category <- factor(gcb_0720$Category, levels=c("B", "N", "R", "T"), 
                            labels=c("B", "Normothermy", "Rewarming", "Torpor"))
torCol <- c("white", "black", "red", "purple")
names(torCol) <- levels(gcb_0720$Category)
colScale <- scale_colour_manual(name = "Category", values = torCol)
ggplot(NULL, aes(x=SampleNo, y=EE_J, col=Category)) +
  geom_path(data=gcb_0720[gcb_0720$SampleNo<20000 & gcb_0720$Category=="Normothermy",], size=1.25) +
  geom_path(data=gcb_0720[gcb_0720$SampleNo>30000 & gcb_0720$Category=="Normothermy",], size=1.25) +
  geom_path(data=gcb_0720[gcb_0720$Category=="Torpor",], size=1.25) +
  geom_path(data=gcb_0720[gcb_0720$Category=="Rewarming",], size=1.25) +
  my_theme_blank + colScale + 
  theme(axis.text.x = element_text(angle=30, hjust=1, size=20),
        legend.key.height=unit(3,"line"),
        axis.line.x = element_line(colour = "grey50"),
        axis.line.y = element_line(colour = "grey50")) +
  scale_x_continuous(breaks= seq(0,36000,3600)) +
  #geom_hline(yintercept=seq(0,50,2)) +
  #ylim(0,50) + 
  xlab("Hours") + 
  ylab("Energy expenditure (J)")

## Figure A3
## Figure A3a : Duration of torpor by site
Hours_site <- ggplot(torpor[torpor$Torpid_not=="T",], aes(Site_full, Hours_torpid)) + geom_boxplot(fill='lightgrey', outlier.size = 3) +
  my_theme2 + theme(axis.text.x = element_text(angle = 20, size=15, hjust=1, color='black'),
                    axis.title=element_text(size=20)) + 
  ylab("Duration of torpor (hours)") + xlab("Site")
Hours_site

## Figure A3b: NEE plot by site
NEE_site <- ggplot(torpor, aes(Site_full, NEE_kJ_variableRER)) + geom_boxplot(fill='lightgrey', outlier.size = 3) +
  my_theme2 + theme(axis.text.x = element_text(angle = 20, size=15, hjust=1, color='black'),
                    axis.title=element_text(size=20)) + 
  ylab(NEE_lab) + xlab("Site")
NEE_site
## Figure A3:
grid.arrange(Hours_site, NEE_site, ncol=2, nrow=1)

## Figure A4a
## Duration vs. Energy savings
ggplot(torpor, aes(savings, Hours_torpid)) + geom_point(aes(col=Species_sciname), size=3, alpha=0.7) + my_theme +
  theme(legend.key.height = unit(3, 'lines')) + xlab("Hourly energy savings (%)") + ylab("Torpor duration (hours)") +
  scale_color_brewer(palette = "Set1", name="Species") + ylim(0,8.5)


## Figure A4b
## Duration vs. minimum chamber temperature of the night
ggplot(torpor, aes(Tc_min_C, Hours_torpid2)) + geom_point(aes(col=Species_sciname), size=3, alpha=0.7) + my_theme +
  theme(legend.key.height = unit(3, 'lines')) + ylab("Torpor duration (hours)") + xlab(Tc_min.xlab) +
  scale_color_brewer(palette = "Set1", name="Species") + ylim(0,8.5)

## Figure A5
## Savings plot by species
savings_plot <- ggplot(torpor[!is.na(torpor$savings),], aes(Species_sciname, savings)) + 
  geom_boxplot(outlier.shape = 19, fill= "light grey") + xlab("Species") + 
  ylab("Hourly torpid energy savings (%)") + theme(legend.position="none") + my_theme +
  stat_summary(fun.data = give.n, geom = "text", vjust=-1, size=10)
savings_plot

## Duration vs. time of entry
## Figure A6
dur_entrytime <- ggplot(torpor, aes(EntryTime_numeric, Hours_torpid)) + 
  geom_point(size=2) + my_theme2 + xlab("Hour of entry") + ylab("Duration of torpor (hours)") +
  geom_smooth(method=lm, size=1, col="black") +
  geom_text(x = 5.5, y = 6, label = lm_eqn(torpor$Hours_torpid, torpor$EntryTime_numeric), 
            parse=T, size=5)
dur_entrytime

