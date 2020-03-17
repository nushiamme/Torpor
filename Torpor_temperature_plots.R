## Code for paper titled:
#"Torpid hummingbirds- energy savings in temperate and tropical sites"
## Paper authors: Anusha Shankar*, Rebecca J Schroeder*, 
# Susan M Wethington, Catherine H Graham, Donald R Powers
## *Equal authors
## Code by: Anusha Shankar, github/nushiamme
## Contact: anusha<dot>shankar<at>stonybrook<dot>edu or nushiamme<at>gmail<dot>com for questions about code
# This script contains code for the two temperature plots

## Contents:
# Figure 2a: Ambient temp plots by hour, per site
# Figure 2b: Chamber Temp plots by hour, per site

require(ggplot2)
require(reshape)
require(gridExtra)

#### Setup ####
## Set wd and read in files
setwd("C:\\Users\\nushi\\Dropbox\\Hummingbird energetics\\July2018\\Data")
#setwd("Dropbox/Hummingbird energetics/Feb2018/Data") #GFU

## Read in files
tc_summ <- read.csv("Tc_AllSites_summ.csv") # Chamber temperatures
ta_summ <- read.csv("Ta_AllSites_summ.csv") # Ambient temperatures
bblh_tnz <- read.csv("BroadBill.csv") ## lab measurements of MR with temperature
bblh_VO2_temp_hourly <- read.csv("BBLH_hourly_VO2_field.csv") ## field measurements of hourly temp and VO2

## General functions, adding columns, ordering factors
my_theme <- theme_classic(base_size = 30) + 
  theme(panel.border = element_rect(colour = "black", fill=NA))

## X-axis labels
Tc.lab <- expression(atop(paste("Chamber temperature ( ", degree,"C)")))
Ta.lab <- expression(atop(paste("Ambient temperature (", degree,"C)")))

## Ordering hour column
ta_summ$Hour2 <- factor(ta_summ$Hour2, levels= c("19", "20", "21", "22", "23", "24",
                                                 "1", "2", "3", "4", "5", "6", "7"), ordered=T)
tc_summ$Hour2 <- factor(tc_summ$Hour2, levels= c("19", "20", "21", "22", "23", "24", 
                                                 "1", "2", "3", "4", "5", "6", "7"), ordered=T)
## Ordering site column
tc_summ$Site <- factor(tc_summ$Site, levels=c('HC','SC','SWRS','MQ','SL'))

## Labels for the x-axis text (Hour)
Hour_labels <- c("1900", "2000", "2100", "2200","2300", "2400", "100", 
                 "200", "300", "400", "500", "600", "700")

Hour_labels_2 <- c("1900", "", "2100", "","2300", "", "100", 
                 "", "300", "", "500", "", "700")


#### Formatting dataframes ####
## Melting dataframes to make them easier to plot
# First for ambient temperature
m.ta <- melt(ta_summ, id.vars = c("Site", "Hour2"), 
             measure.vars = c("Mean_Ta", "Min_Ta", "Max_Ta"))
names(m.ta) <- c("Site", "Hour", "Variable", "Temperature") #Rename columns
m.ta$Site <- factor(m.ta$Site, levels=c('HC','SC','SWRS','MQ','SL')) # Reorder sites
# Rename sites to full names
levels(m.ta$Site) <- c("Harshaw", "Sonoita", "Southwestern Research Station", "Maquipucuna", "Santa Lucia")

# Then for chamber temperature
m.tc <- melt(tc_summ, id.vars = c("Site", "Hour2"), 
             measure.vars = c("Mean_Tc", "Min_Tc", "Max_Tc"))
names(m.tc) <- c("Site", "Hour", "Variable", "Temperature") #Rename columns
m.tc$Site <- factor(m.tc$Site, levels=c('HC','SC','SWRS','MQ','SL')) # Reorder sites
# Rename sites to full names
levels(m.tc$Site) <- c("Harshaw", "Sonoita", "Southwestern Research Station", "Maquipucuna", "Santa Lucia")

#### Plots ####
## Figure 2a
## Ambient temp plots by hour, per site
AmbTemp <- ggplot(m.ta, aes(Hour,Temperature, alpha=Variable)) + facet_grid(.~Site) +  my_theme +
  facet_grid(~Site, labeller = labeller(Site = label_wrap_gen(10))) +
  scale_x_discrete(labels=Hour_labels_2) +
  geom_line(aes(group=Variable, col=Variable), size=1.5) +
  scale_color_manual(values=c("Black", "Blue", "Red")) +
  scale_alpha_manual(values = c(1, 0.5, 0.5)) +
  theme(axis.text.x = element_text(angle = 60, size=18, hjust=1), 
        legend.position="none", plot.title = element_text(size = 30),
        panel.grid.major.y = element_line(size=.1, color="grey75"), 
        strip.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y=element_text(vjust=-3)) + 
  ylim(0,35) +
  ylab(Ta.lab) #+ ggtitle("a.") 
AmbTemp

# Figure 2b
## Chamber Temp plots by hour, per site
ChambTemp <- ggplot(m.tc, aes(Hour,Temperature, alpha=Variable)) + my_theme + 
  facet_grid(~Site, labeller = labeller(Site = label_wrap_gen(10))) +
  scale_x_discrete(labels=Hour_labels_2) +
  geom_line(aes(group=Variable, col=Variable), size=1.5) +
  scale_color_manual(values=c("Black", "Blue", "Red")) +
  scale_alpha_manual(values = c(1, 0.5, 0.5)) +
  theme(axis.text.x = element_text(angle = 60, size=18, hjust=1), 
        legend.position="none", plot.title = element_text(size = 30),
        panel.grid.major.y = element_line(size=.1, color="grey75"), 
        strip.text.x = element_text(size = 18),
        axis.title.y=element_text(vjust=-3)) + 
  ylim(0,35) +
  xlab("Hour") + ylab(Tc.lab)
ChambTemp

grid.arrange(AmbTemp, ChambTemp, ncol=1, nrow=2)

## Table 4 - NEE in Harshaw and Sonoita with current and with warming temperatures
bblh_normo_eqn <- coef(lm(bblh_tnz$VO2_Normothermic~bblh_tnz$Temp_C))
polyn.formula <- bblh_torpid$VO2 ~ poly(bblh_torpid$Temperature, 2, raw = TRUE)
m_pol <- coef(lm(polyn.formula, bblh_torpid))

m.ta_ordered <- m.ta
levels(m.ta_ordered$Hour) <- c(as.character(seq(1:length(levels(m.ta_ordered$Hour)))))
m.ta_ordered$Hour <- as.numeric(as.character(m.ta_ordered$Hour))
earlyRER <- 21.16
laterRER <- 19.67


## For HC
m.ta_HC <- m.ta_ordered[m.ta_ordered$Site=="Harshaw" & m.ta_ordered$Variable=="Mean_Ta",]
m.ta_HC$O2_max_torpor <- 0
m.ta_HC$O2_min_torpor <- 0
m.ta_HC$kJ_hour_minTorpor <- 0
m.ta_HC$kJ_hour_maxTorpor <- 0
m.ta_HC$TempPlus3C <- m.ta_HC$Temperature+3
m.ta_HC$O2_max_torpor_warming3C <- 0
m.ta_HC$O2_min_torpor_warming3C <- 0

## Min torpor use for harshaw, 2 hours
m.ta_HC$O2_min_torpor[m.ta_HC$Hour %in% c(1:9)] <- bblh_normo_eqn[[1]] + 
    bblh_normo_eqn[[2]]*m.ta_HC$Temperature[m.ta_HC$Hour %in% c(1:9)]
m.ta_HC$O2_min_torpor[m.ta_HC$Hour %in% c(10,11)] <- m_pol[[1]] + 
    m_pol[[2]]*m.ta_HC$Temperature[m.ta_HC$Hour %in% c(4,5)] +
    m_pol[[3]]*(m.ta_HC$Temperature[m.ta_HC$Hour %in% c(4,5)]^2)

## Max torpor use for harshaw (6 hours torpor)
m.ta_HC$O2_max_torpor[m.ta_HC$Hour %in% c(1:5)] <- bblh_normo_eqn[[1]] + 
    bblh_normo_eqn[[2]]*m.ta_HC$Temperature[m.ta_HC$Hour %in% c(1:5)]
m.ta_HC$O2_max_torpor[m.ta_HC$Hour %in% c(6:11)] <- m_pol[[1]] + 
    m_pol[[2]]*m.ta_HC$Temperature[m.ta_HC$Hour %in% c(6:11)] +
    m_pol[[3]]*(m.ta_HC$Temperature[m.ta_HC$Hour %in% c(6:11)]^2)

# Calculating kJ using high RER for first two hours and low for the rest
m.ta_HC$kJ_hour_minTorpor[m.ta_HC$Hour<3] <- m.ta_HC$O2_min_torpor[m.ta_HC$Hour<3]*earlyRER*60/1000
m.ta_HC$kJ_hour_minTorpor[m.ta_HC$Hour>3] <- m.ta_HC$O2_min_torpor[m.ta_HC$Hour>3]*laterRER*60/1000
m.ta_HC$kJ_hour_maxTorpor[m.ta_HC$Hour<3] <- m.ta_HC$O2_max_torpor[m.ta_HC$Hour<3]*earlyRER*60/1000
m.ta_HC$kJ_hour_maxTorpor[m.ta_HC$Hour>3] <- m.ta_HC$O2_max_torpor[m.ta_HC$Hour>3]*laterRER*60/1000



m.ta_SC <- m.ta[m.ta$Site=="Sonoita" & m.ta$Variable=="Mean_Ta",]
m.ta_SC$MR_max_torpor <- "NA"
m.ta_SC$MR_min_torpor <- "NA"
m.ta_SC$TempPlus3C <- "NA"
m.ta_SC$MR_max_torpor_warming3C <- "NA"
m.ta_SC$MR_min_torpor_warming3C <- "NA"


bblh_VO2_temp <- ggplot(bblh_tnz, aes(Temp_C, VO2_Normothermic)) + 
  geom_point(size=3) + my_theme + 
  geom_smooth(method = "lm", se = T, 
              formula = y~ x, 
              colour = "black") + 
  annotate(geom = "text", x = 10, y = .2, label = polyn.text_lab, 
           family = "serif", hjust = 0, parse = TRUE, size=10) +
  ylab(VO2_lab) + xlab(Tc.xlab)
plot(bblh_VO2_temp)


