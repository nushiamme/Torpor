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
## Figure 3: VO2 for broad-billed hummingbirds showing inflection point at 14-15 degC
## Supplementary Figure S5: Faceted plot comparing total duration and total nighttime energy expenditure per site
## Supplementary Figure S6: Duration of torpor per individual per night, 
# as a function of the hour of entry
## Supplementary Figure S7: Average hourly torpid energy savings relative to normothermy
## Supplementary Figure S8: Average hourly mass-corrected energy expenditure 
# as a function of minimum hourly chamber temperature
## Supplementary Figure S10: Total nighttime energy expenditure including and 
# excluding rewarming, as a function of individual mass


#### Setup #### 
## Libraries and reading in data
library(ggplot2)
library(reshape)
library(polynom)

## setwd and read in file
setwd("C:\\Users\\ANUSHA\\Dropbox\\Hummingbird energetics\\Submission_Nov2017\\Data")
#GFU wd
setwd("/Users/anshankar/Dropbox/Hummingbird energetics/Submission_Nov2017/Data/")

## Use sep=";" if using a csv format from Europe.
torpor <- read.csv("Torpor_individual_summaries.csv") # Torpor summaries per individual
bblh_tnz <- read.csv("BroadBill.csv") ## For Figure 3, BBLH minimum body temperature


# General functions ####
## Saving standard theme  
my_theme <- theme_classic(base_size = 30) + 
  theme(axis.title.y = element_text(color = "black", vjust = 2),
        panel.border = element_rect(colour = "black", fill=NA))

## Theme with slightly smaller font
my_theme2 <- my_theme + theme_classic(base_size = 15)

## Template axis labels
Tc.xlab <- expression(atop(paste("Chamber Temperature (", degree,"C)"))) # for chamber temperature
NEE_corrlab <- bquote('Nighttime energy expenditure (kJ/' ~M^(0.67)*')') # for mass-corrected nighttime energy expenditure
VO2_lab <- expression(paste(VO[2]~mL~O[2]/min))

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

## Make a function to arrange plots with a shared legend
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}


# Adding/Formatting columns ####

## In BBLH tnz data set
bblh_tnz$N_T <- factor(bblh_tnz$N_T, levels=c('T', 'N')) # Reorder levels BBLH tnz dataset

## In torpor dataset
## Make species a sensible order, for just species that used torpor - used in Supp Fig S6
torpor$Species2 <- factor(torpor$Species,
                          levels = c('BBLH','RIHU','GCB','FBB','TBH', "WNJ"), ordered = T)

#For Supp Fig S5, order and expand site names
torpor$Site_full <- torpor$Site
torpor$Site_full <- factor(torpor$Site_full, levels=c('HC', 'SC', 'SWRS', 'MQ','SL'))
levels(torpor$Site_full) <- c("Harshaw", "Sonoita", "Southwestern Research Station", "Maquipucuna", "Santa Lucia")


## Making a column for mass-corrected total Nighttime energy expenditure - useful for summary tables
torpor$NEE_MassCorrected <- torpor$NEE_kJ/(torpor$Mass^(2/3))

## Savings column to convert percentage energy expended in torpor relative to 
# normothermy, into savings relative to normothermy - used in Supp Fig S6
torpor$savings <- 100-torpor$Percentage_avg

# Adding column to torpor dataset dividing energy expenditure by 2/3*Mass to correct for 
# mass with allometric scaling - can be changed to test other allometries
# for average torpid and normothermic hourly Energy expenditure
torpor$AvgEE_normo_MassCorrected <- torpor$Avg_EE_hourly_normo/(torpor$Mass^(2/3))
torpor$AvgEE_torpid_MassCorrected <- torpor$Avg_EE_hourly_torpid/(torpor$Mass^(2/3))

## Subset just BBLH data (for which we have enough data at multiple sites 
#to do more in-depth analyses)
BBLH_torpor <- subset(torpor, Species=="BBLH")

## Melt BBLH dataframe to put torpid and normo in same column
m_BBLH_tor_nor <- melt(BBLH_torpor, id.vars="Tc_min_C", 
                       measure.vars = c("AvgEE_torpid_MassCorrected", 
                                        "AvgEE_normo_MassCorrected"))
levels(m_BBLH_tor_nor$variable)[levels(m_BBLH_tor_nor$variable)=="AvgEE_normo_MassCorrected"] <- 
  "Avg Normothermic EE"
levels(m_BBLH_tor_nor$variable)[levels(m_BBLH_tor_nor$variable)=="AvgEE_torpid_MassCorrected"] <- 
  "Avg Torpid EE"

m_BBLH_tor_nor$variable <- factor(m_BBLH_tor_nor$variable, 
                                   levels = rev(levels(m_BBLH_tor_nor$variable)))

bblh_torpid <- bblh_tnz[bblh_tnz$N_T=="T",]

#### Plots ####
## Figure 3 
# Broadbill torpid energy expenditure with temperature, showing an inflection point
# at 14-15 degC. Only using subset of the dataset- just torpid values.
# These measurements were taken under controlled conditions in 5 degC temperature steps, 
# separately from all the other torpor measurements
bblh_VO2_temp <- ggplot(bblh_torpid, aes(Temp_C, VO2_all)) + 
  geom_point(size=3) + my_theme2 + geom_smooth(stat='smooth', method='loess', color='black') +
  ylab("Oxygen consumption (ml/min)") + xlab(Tc.xlab)
plot(bblh_VO2_temp)

## Trying out second order polynomial
polyn.formula_lab <- bblh_torpid$VO2_all ~ poly(bblh_torpid$Temp_C, 2, raw = TRUE)
m_pol_lab <- lm(polyn.formula_lab, bblh_torpid)
polyn.eq_lab <- as.character(signif(as.polynomial(coef(m_pol_lab)), 2))
polyn.text_lab <- paste(gsub("x", "~italic(x)", polyn.eq_lab, fixed = TRUE),
                    paste("italic(R)^2",  
                          format(summary(m_pol_lab)$r.squared, digits = 2), 
                          sep = "~`=`~"),
                    sep = "~~~~")

bblh_VO2_temp <- ggplot(bblh_tnz[bblh_tnz$N_T=="T",], aes(Temp_C, VO2_all)) + 
  geom_point(size=3) + my_theme + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), col='black') +
  annotate(geom = "text", x = 10, y = .2, label = polyn.text_lab, 
           family = "serif", hjust = 0, parse = TRUE, size=10) +
  ylab(VO2_lab) + xlab(Tc.xlab) + ylim(-0.01,0.3)
plot(bblh_VO2_temp)


## Duration vs. time of entry
## Supp Figure S5
dur_entrytime <- ggplot(torpor, aes(EntryTime_numeric, Hours_torpid)) + 
  geom_point(size=2) + my_theme2 + xlab("Hour of entry") + ylab("Duration of torpor (hours)") +
  geom_smooth(method=lm, size=1, col="black") +
  geom_text(x = 5.5, y = 6, label = lm_eqn(torpor$Hours_torpid, torpor$EntryTime_numeric), 
            parse=T, size=5)
dur_entrytime

## Supp Figure S5
## NEE plot by site
NEE_site <- ggplot(torpor, aes(Site_full, NEE_MassCorrected)) + geom_boxplot(fill='lightgrey', outlier.size = 3) +
  my_theme2 + theme(axis.text.x = element_text(angle = 20, size=15, hjust=1, color='black'),
                    axis.title=element_text(size=20)) + 
  ylab(NEE_corrlab) + xlab("Site")
NEE_site
## Now hours torpid
Hours_site <- ggplot(torpor[torpor$Torpid_not=="T",], aes(Site_full, Hours_torpid)) + geom_boxplot(fill='lightgrey', outlier.size = 3) +
  my_theme2 + theme(axis.text.x = element_text(angle = 20, size=15, hjust=1, color='black'),
                    axis.title=element_text(size=20)) + 
  ylab("Duration of torpor (hours)") + xlab("Site")
Hours_site
grid.arrange(Hours_site, NEE_site, ncol=2, nrow=1)

## Supp Figure S6
## Savings plot by species
savings_plot <- ggplot(torpor[!is.na(torpor$savings),], aes(Species2, savings)) + 
  geom_boxplot(outlier.shape = 19, fill= "light grey") + xlab("Species") + 
  ylab("Hourly torpid energy savings (%)") + theme(legend.position="none") + my_theme +
  stat_summary(fun.data = give.n, geom = "text", vjust=-1, size=10)
savings_plot

## BBLH and GCB plots, building up to Supp Figure S7
## Both normo and torpid avg mass-corrected EE for BBLH on same graph
BBLH_tor_nor <- ggplot(m_BBLH_tor_nor, aes(as.numeric(Tc_min_C), value, shape=variable)) +
  geom_point(aes(shape=variable), size=6) + my_theme2 +
  geom_smooth(method=lm, size=1, col="black") + 
  scale_shape_manual("Hourly Energy Expenditure\n", 
                     values=c(16,1), labels=c("Normothermic", "Torpid")) +
  theme(legend.key.height=unit(3,"line"), legend.position="bottom", 
        axis.title.y = element_text(vjust = 2)) +
  ylab("Mean BBLH Energy Expenditure (kJ/g)") + xlab(Tc.xlab) + ggtitle("b.")
BBLH_tor_nor 

## Both normo and torpid avg mass-corrected EE for GCB on same graph
GCB_tor_nor <- ggplot(m_GCB_tor_nor, aes(as.numeric(Tc_min_C), value, shape=variable)) +
  my_theme2 +  geom_point(aes(shape=variable), size=6) + 
  geom_smooth(method=lm, size=1, col="black") + 
  scale_shape_manual("Hourly Energy Expenditure:\n", 
                     values=c(16,1), labels=c("Normothermic", "Torpid")) +
  theme(legend.key.height=unit(1,"line"),
        legend.position="bottom", axis.title.y = element_text(vjust = 2)) + ggtitle("a.") +
  ylab("Mean GCB Energy Expenditure (kJ/g)") + xlab(Tc.xlab)
GCB_tor_nor 

## Arrange BBLH and GCB tor_nor plots in a single window
## This is Supp Figure S7
grid_arrange_shared_legend(GCB_tor_nor, BBLH_tor_nor)

## Nighttime energy expenditure with and without rewarming costs - not mass-corrected
## Supp Figure S9
rewarming_NEE <- ggplot(torpor[torpor$Torpid_not=="T",], aes(Mass, NEE_without_rewarming_kJ)) + 
  my_theme2 + geom_point(aes(col=Rewarming_Tc), size=3, alpha=0.5) +
  scale_colour_gradient(low="blue", high="red", name="Chamber \n temperature") +
  geom_point(aes(Mass, NEE_kJ), col='black')  +
  geom_segment(aes(x=Mass,y=NEE_kJ,xend=Mass,yend=NEE_without_rewarming_kJ, col=Rewarming_Tc),
               arrow = arrow(length = unit(0.01, "npc"),type = "closed"), alpha=0.5) +
  xlab("Mass") + ylab("Nighttime energy expenditure (kJ)") 
rewarming_NEE