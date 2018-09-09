## Code for paper titled:
#"Torpid hummingbirds- energy savings in temperate and tropical sites"
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
## Supplementary Figure S2: Site-wise plot comparing total duration and total nighttime energy expenditure
## Supplementary Figure S3: Average hourly torpid energy savings relative to normothermy
## Supplementary Figure S4a: Torpor duration vs. minimum chamber temperature (Tc min) for the night. 
## Supplementary Figure S4b: Torpor duration vs. average hourly energy savings in torpor relative to normothermy.


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

# General functions ####
## Saving standard theme  
my_theme <- theme_classic(base_size = 30) + 
  theme(axis.title.y = element_text(color = "black", vjust = 2),
        panel.border = element_rect(colour = "black", fill=NA))

## Theme with slightly smaller font
my_theme2 <- my_theme + theme_classic(base_size = 20)

## Template axis labels
Tc.xlab <- expression(atop(paste("Chamber temperature (", degree,"C)"))) # for chamber temperature
NEE_masslab <- 'Nighttime energy expenditure (kJ/g)' # for mass-corrected nighttime energy expenditure
VO2_lab <- expression(paste(VO[2]~mL~O[2]/min))
Tc_min.xlab <- expression(atop(paste("Minimum chamber temperature (", degree,"C)")))

## Function for adding a regression equation to graphs
## (Where y= table$column for y in the equation and x= table$column for x)
lm_eqn <- function(y, x){
  m <- lm(y ~ x);
  eq <- substitute(italic(y) == 
                     a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
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

#For Supp Fig S5, order and expand site names
torpor$Site_full <- torpor$Site
torpor$Site_full <- factor(torpor$Site_full, levels=c('HC', 'SC', 'SWRS', 'MQ','SL'))
levels(torpor$Site_full) <- c("Harshaw", "Sonoita", "Southwestern Research Station", "Maquipucuna", "Santa Lucia")

## Making a column for mass-corrected total Nighttime energy expenditure
torpor$NEE_Mass <- torpor$NEE_kJ/torpor$Mass

## Savings column to convert percentage energy expended in torpor relative to 
# normothermy, into savings relative to normothermy - used in Supp Fig S7
torpor$savings <- 100-torpor$Percentage_avg

bblh_controlled_torpor <- bblh_tnz[bblh_tnz$N_T=="T",] # Subsetting just torpor measurements from controlled lab data

## for Supp Fig. S4b, 
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
  annotate(geom = "text", x = 10, y = .2, label = polyn.text_lab, 
           family = "serif", hjust = 0, parse = TRUE, size=10) +
  ylab(VO2_lab) + xlab(Tc.xlab)
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
  annotate(geom = "text", x = 10, y = .2, label = polyn.text, 
           family = "serif", hjust = 0, parse = TRUE, size=10) +
  ylim(-0.01,0.3) + guides(col=guide_legend(title="Individual ID")) +
  my_theme + theme(legend.key.height = unit(3, 'lines')) +
  xlab(Tc.xlab) + ylab(VO2_lab)

## Figure 4
## Nighttime energy expenditure with and without rewarming costs - not mass-corrected
rewarming_NEE <- ggplot(torpor[torpor$Torpid_not=="T",], aes(Mass, NEE_without_rewarming_kJ)) + 
  my_theme2 + geom_point(aes(col=Rewarming_Tc), size=3, alpha=0.5) +
  scale_colour_gradient(low="blue", high="red", name="Chamber \n temperature") +
  geom_point(aes(Mass, NEE_kJ), col='black')  +
  geom_segment(aes(x=Mass,y=NEE_kJ,xend=Mass,yend=NEE_without_rewarming_kJ, col=Rewarming_Tc),
               arrow = arrow(length = unit(0.01, "npc"),type = "closed"), alpha=0.5) +
  xlab("Mass") + ylab("Nighttime energy expenditure (kJ)") 
rewarming_NEE

## Supp Figure S2
## Figure S2a : Duration of torpor by site
Hours_site <- ggplot(torpor[torpor$Torpid_not=="T",], aes(Site_full, Hours_torpid)) + geom_boxplot(fill='lightgrey', outlier.size = 3) +
  my_theme2 + theme(axis.text.x = element_text(angle = 20, size=15, hjust=1, color='black'),
                    axis.title=element_text(size=20)) + 
  ylab("Duration of torpor (hours)") + xlab("Site")
Hours_site

## Figure S2b: NEE plot by site
NEE_site <- ggplot(torpor, aes(Site_full, NEE_Mass)) + geom_boxplot(fill='lightgrey', outlier.size = 3) +
  my_theme2 + theme(axis.text.x = element_text(angle = 20, size=15, hjust=1, color='black'),
                    axis.title=element_text(size=20)) + 
  ylab(NEE_masslab) + xlab("Site")
NEE_site
## Figure S2:
grid.arrange(Hours_site, NEE_site, ncol=2, nrow=1)

## Duration vs. time of entry
## Supp Figure S6
dur_entrytime <- ggplot(torpor, aes(EntryTime_numeric, Hours_torpid)) + 
  geom_point(size=2) + my_theme2 + xlab("Hour of entry") + ylab("Duration of torpor (hours)") +
  geom_smooth(method=lm, size=1, col="black") +
  geom_text(x = 5.5, y = 6, label = lm_eqn(torpor$Hours_torpid, torpor$EntryTime_numeric), 
            parse=T, size=5)
dur_entrytime

## Supp Figure S3
## Savings plot by species
savings_plot <- ggplot(torpor[!is.na(torpor$savings),], aes(Species_sciname, savings)) + 
  geom_boxplot(outlier.shape = 19, fill= "light grey") + xlab("Species") + 
  ylab("Hourly torpid energy savings (%)") + theme(legend.position="none") + my_theme +
  stat_summary(fun.data = give.n, geom = "text", vjust=-1, size=10)
savings_plot

## Supp Figure S4a
## Duration vs. Energy savings
ggplot(torpor, aes(100-Percentage_avg, Hours_torpid)) + geom_point(aes(col=Species_sciname), size=3, alpha=0.7) + my_theme +
  theme(legend.key.height = unit(3, 'lines')) + xlab("Hourly energy savings (%)") + ylab("Torpor duration (hours)") +
  scale_color_brewer(palette = "Set1", name="Species")

## Supp Figure S4b
## Duration vs. minimum chamber temperature of the night
ggplot(torpor, aes(Tc_min_C, Hours_torpid2)) + geom_point(aes(col=Species_sciname), size=3, alpha=0.7) + my_theme +
  theme(legend.key.height = unit(3, 'lines')) + ylab("Torpor duration (hours)") + xlab(Tc_min.xlab) +
  scale_color_brewer(palette = "Set1", name="Species")
