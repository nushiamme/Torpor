## Code for paper titled:
#"Torpid hummingbirds- energy savings in temperate and tropical sites"
## Paper authors: Anusha Shankar*, Rebecca J Schroeder*, 
# Susan M Wethington, Catherine H Graham, Donald R Powers
# *Equal authors
## Script started on: Sept 24, 2016

require(ggplot2) 
require(reshape)

## Set wd and read in files
setwd("C:\\Users\\ANUSHA\\Dropbox\\Hummingbird energetics\\Submission_Jul2017\\Data")

## Read in files
tc_summ <- read.csv("Tc_AllSites_summ.csv") # Chamber temperatures
ta_summ <- read.csv("Ta_AllSites_summ.csv") # Ambient temperatures

## General functions, adding columns, ordering factors
my_theme <- theme_classic(base_size = 30) + 
  theme(panel.border = element_rect(colour = "black", fill=NA))

## X-axis labels
Tc.lab <- expression(atop(paste("Chamber Temperature ( ", degree,"C)")))
Ta.lab <- expression(atop(paste("Ambient Temperature (", degree,"C)")))

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

## Melting dataframes to make them easier to plot
# First for ambient temperature
m.ta <- melt(ta_summ, id.vars = c("Site", "Hour2"), 
             measure.vars = c("Mean_Ta", "Min_Ta", "Max_Ta"))
names(m.ta) <- c("Site", "Hour", "Variable", "Temperature") #Rename columns
m.ta$Site <- factor(m.ta$Site, levels=c('HC','SC','SWRS','MQ','SL')) # Reorder sites
# Rename sites to full names
levels(m.ta$Site) <- c("Harshaw", "Sonoita", "Southwest Research Station", "Maqui", "Santa Lucia")

# Then for chamber temperature
m.tc <- melt(tc_summ, id.vars = c("Site", "Hour2"), 
             measure.vars = c("Mean_Tc", "Min_Tc", "Max_Tc"))
names(m.tc) <- c("Site", "Hour", "Variable", "Temperature") #Rename columns
m.tc$Site <- factor(m.tc$Site, levels=c('HC','SC','SWRS','MQ','SL')) # Reorder sites
# Rename sites to full names
levels(m.tc$Site) <- c("Harshaw", "Sonoita", "Southwest Research Station", "Maqui", "Santa Lucia")

#### Plots ####
## Figure 2
## Chamber Temp plots by hour, per site
ChambTemp <- ggplot(m.tc, aes(Hour,Temperature, alpha=Variable)) + my_theme + 
  facet_grid(~Site, labeller = labeller(Site = label_wrap_gen(10))) +
  #geom_point(aes(group=Variable, col=Variable), size=1.5) +
  geom_line(aes(group=Variable, col=Variable), size=1.5) +
  scale_color_manual(values=c("Black", "Blue", "Red")) +
  scale_alpha_manual(values = c(1, 0.5, 0.5)) +
  theme(axis.text.x = element_text(angle = 60, size=12, hjust=1), 
        legend.position="none", plot.title = element_text(size = 30),
        panel.grid.major.y = element_line(size=.1, color="grey75"), 
        strip.text.x = element_text(size = 18),
        axis.title.y=element_text(vjust=-3)) + 
  xlab("Hour") + ylab(Tc.lab) + scale_x_discrete(labels=Hour_labels) #+ ggtitle("b.")
ChambTemp

## Supp. figure S2
## Ambient temp plots by hour, per site
AmbTemp <- ggplot(m.ta, aes(Hour,Temperature, alpha=Variable)) + facet_grid(.~Site) +  my_theme +
  facet_grid(~Site, labeller = labeller(Site = label_wrap_gen(10))) +
  scale_x_discrete(labels=Hour_labels) +
  geom_line(aes(group=Variable, col=Variable), size=1.5) +
  scale_color_manual(values=c("Black", "Blue", "Red")) +
  scale_alpha_manual(values = c(1, 0.5, 0.5)) +
  theme(axis.text.x = element_text(angle = 60, size=12, hjust=1), 
        legend.position="none", plot.title = element_text(size = 30),
        panel.grid.major.y = element_line(size=.1, color="grey75"), 
        strip.text.x = element_text(size = 18),
        axis.title.y=element_text(vjust=-3)) + 
  xlab("Hour") + ylab(Ta.lab) #+ ggtitle("a.") 
AmbTemp
