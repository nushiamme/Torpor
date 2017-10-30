## Analyzing and plotting ambient and chamber temperature data for torpor paper
## Script by Anusha Shankar
## Paper authors: Anusha Shankar*, Rebecca J Schroeder*, Susan M Wethington, Catherine H Graham, Donald R Powers
## Script started on: Sept 24, 2016
## Thanks
## If you know how to make the hours into ordered factors that are easy to plot without the messiness below, 
## please email anusha<dot>shankar<at>stonybrook<dot>edu
# * = equal authors

require(ggplot2) 
require(reshape)

## Set wd and read in files
setwd("C://Users//shankar//Dropbox//Hummingbird energetics//Submission_Oct2016")

## Read in files
tatc <- read.csv("TempSummary_AllSites.csv")
## Made these two in R with the aggregating chunk below, so avoid that chunk if reading these two in
tc_summ <- read.csv("Tc_AllSites_summ.csv")
ta_summ <- read.csv("Ta_AllSites_summ.csv")


## General functions, adding columns, ordering factors
my_theme <- theme_classic(base_size = 30) + 
  theme(panel.border = element_rect(colour = "black", fill=NA))

pd <- position_dodge(0.1) # move them .05 to the left and right

tatc$Hour_rounded <- factor(tatc$Hour_rounded, 
                            levels= c("1900", "1930", "2000", "2030", "2100", "2130", "2200", "2230", "2300", "2330", "2400",
                                      "2430", "100", "130", "200", "230", "300", "330", "400", "430", "500", "530",
                                      "600", "630", "700"), ordered=T)

tatc$Hour2 <- factor(tatc$Hour2, levels= c("19", "20", "21", "22", "23", "24", "1", "2", "3", "4", "5", "6", "7"), 
                     ordered=T)
ta_summ$Hour2 <- factor(ta_summ$Hour2, levels= c("19", "20", "21", "22", "23", "24", "1", "2", "3", "4", "5", "6", "7"), 
                        ordered=T)
tc_summ$Hour2 <- factor(tc_summ$Hour2, levels= c("19", "20", "21", "22", "23", "24", "1", "2", "3", "4", "5", "6", "7"), 
                        ordered=T)
tc_summ$Site <- factor(tc_summ$Site, levels=c('HC','SC','SWRS','MQ','SL'))

Hour_labels <- c("1900", "2000", "2100", "2200","2300", "2400", "100", "200", "300", "400", "500", "600", "700")

Tc.lab <- expression(atop(paste("Chamber Temperature ( ", degree,"C)")))
Ta.lab <- expression(atop(paste("Ambient Temperature (", degree,"C)")))

m.ta <- melt(ta_summ, id.vars = c("Site", "Hour2"), measure.vars = c("Mean_Ta", "Min_Ta", "Max_Ta"))
names(m.ta) <- c("Site", "Hour", "Variable", "Temperature")
m.tc <- melt(tc_summ, id.vars = c("Site", "Hour2"), measure.vars = c("Mean_Tc", "Min_Tc", "Max_Tc"))
names(m.tc) <- c("Site", "Hour", "Variable", "Temperature")
levels(m.tc$Site) <- c("Harshaw", "Sonoita", "Southwest Research Station", "Maqui", "Santa Lucia")
m.ta$Site <- factor(m.ta$Site, levels=c('HC','SC','SWRS','MQ','SL'))
levels(m.ta$Site) <- c("Harshaw", "Sonoita", "Southwest Research Station", "Maqui", "Santa Lucia")

tc_summ$Site <- factor(tc_summ$Site, levels=c('HC','SC','SWRS','MQ','SL'))
tatc_summ <- merge(ta_summ, tc_summ, by=c("Site", "Hour2"))
tatc_summ$Site <- factor(tatc_summ$Site, levels=c('HC','SC','SWRS','MQ','SL'))

#### Plots ####
## Figure 1
## Chamber Temp plots by hour, per site
ChambTemp <- ggplot(m.tc, aes(Hour,Temperature, alpha=Variable)) + my_theme + 
  facet_grid(~Site, labeller = labeller(Site = label_wrap_gen(10))) +
  #geom_point(aes(group=Variable, col=Variable), size=1.5) +
  geom_line(aes(group=Variable, col=Variable), size=1.5) +
  scale_color_manual(values=c("Black", "Blue", "Red")) +
  scale_alpha_manual(values = c(1, 0.5, 0.5)) +
  theme(axis.text.x = element_text(angle = 60, size=12, hjust=1), legend.position="none", plot.title = element_text(size = 30),
        panel.grid.major.y = element_line(size=.1, color="grey75"), strip.text.x = element_text(size = 18),
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
  theme(axis.text.x = element_text(angle = 60, size=12, hjust=1), legend.position="none", plot.title = element_text(size = 30),
        panel.grid.major.y = element_line(size=.1, color="grey75"), strip.text.x = element_text(size = 18),
        axis.title.y=element_text(vjust=-3)) + #xlim("1900", "2200", "100", "400", "700") +
  xlab("Hour") + ylab(Ta.lab) #+ ggtitle("a.") 
AmbTemp

## Old ambient temp - Black with "Error-like bars" depicting max and min temps
AmbTemp <- ggplot(ta_summ, aes(Hour2,Mean_Ta)) + my_theme + facet_grid(.~Site) +  
  geom_point(aes(group=Site), size=1.5) +
  geom_line(aes(group=Site)) +
  geom_errorbar(aes(ymin= Min_Ta, ymax= Max_Ta), alpha=0.6, width=.1, position=pd) +
  theme(axis.text.x = element_text(angle = 90, size=15), legend.position="none") +
  xlab("Hour") + ylab(Ta.lab) + ggtitle("Sites") + theme(plot.title = element_text(size = 20)) +
  scale_x_discrete(labels=Hour_labels)
AmbTemp
