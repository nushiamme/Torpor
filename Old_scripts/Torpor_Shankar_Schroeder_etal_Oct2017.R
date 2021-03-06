##Code for paper titled:
#"Torpid hummingbirds- energy savings in temperate and tropical sites"
## Paper authors: Anusha Shankar*, Rebecca J Schroeder*, Susan M Wethington, Catherine H Graham, Donald R Powers
## Code by: Anusha Shankar, github/nushiamme
## Contact: anusha<dot>shankar<at>stonybrook<dot>edu for questions
## Plots and analyses for paper on hummingbird torpor in Arizona and Ecuador

#Contents
# Figure 2: Chamber temperatures
# Figure 3: VO2 for broad-billed hummingbirds showing inflection point at 14-15 degC


#### libraries and reading in data ####
library(ggplot2)
library(reshape)
library(gridExtra)
library(grid)
library(wq)
library(RColorBrewer)

## setwd and read in file
setwd("C:\\Users\\ANUSHA\\Dropbox\\Hummingbird energetics\\Submission_Jul2017\\Data")


## Use sep=";" if using a csv format from Europe.
torpor <- read.csv("Torpor_individual_summaries.csv") # Torpor summaries per individual
freq_sites <- read.csv("Frequency_torpor_sites.csv") # Frequency of torpor use per species, organized by site
freq_sp <- read.csv("Frequency_torpor_sp.csv") # Frequency of torpor use per species, organized by species

## Make the column for 'proportion of night spent torpid' a numeric column
torpor$Prop_hours <- as.numeric(as.character(torpor$Prop_hours))

##### Adding columns #######
## Adding column dividing NEE by 2/3*Mass to correct for mass with allometric scaling - can be changed to test
# other allometries
torpor$NEE_MassCorrected<- torpor$NEE_kJ/(torpor$Mass^(2/3))

# Line to arrange Site facets in sensible order
torpor$Site_new <- factor(torpor$Site, levels=c('HC','SC','SWRS','MQ','SL'))

## New column for time of entry into torpor, ordered in a chronological way
torpor$EntryTime_new <- factor(torpor$Time_of_entry, 
                               levels=c('2000', '2100', '2130', '2200', '2230',
                                        '2300', '2400', '100', '130', '200', '330', '500'))

freq_sites$Site_new <- factor(freq_sites$Site, levels=c('HC','SC','SWRS','MQ','SL'))

## Make species a sensible order, for just species that used torpor
torpor$Species2 <- factor(torpor$Species,
                          levels = c('BBLH','MAHU','GCB','FBB','TBH', "WNJ"), ordered = T)

## Create Hours_torpid2 column and change NA's to 0's for lm analyses; keep original Hours_torpid column with NA's
torpor$Hours_torpid2 <- torpor$Hours_torpid
torpor$Hours_torpid2[is.na(torpor$Hours_torpid2)] <- 0

## Savings column to convert percentage energy expended in torpor relative to normothermy into savings
# relative to normothermy
torpor$savings <- 100-torpor$Percentage_avg

#### Make new data frames ####
## Make dataframe to check site-wise frequency of torpor
rate_site <- data.frame(table(torpor$Site,torpor$Torpid_not))
names(rate_site) <- c("Site", "Torpid_not", "N")

## Subset just BBLH data (for which we have enough data at multiple sites 
#to do more in-depth analyses)
BBLH_torpor <- subset(torpor, Species=="BBLH")

## Subset just GCB data (for which we have enough data at multiple sites
#to do more in-depth analyses)
GCB_torpor <- subset(torpor, Species=="GCB")

## Make table to summarize nighttime energy expenditure
nee.agg <- aggregate(torpor$NEE_kJ, 
                     by=list(torpor$Torpid_not, torpor$Site_new, 
                             torpor$Species), 
                     FUN="mean", na.rm=T)
names(nee.agg) <- c("Torpid_not", "Site", "Species", "NEE_kJ")
#nee.agg #to check

## Summarize hours spent torpid
hours.agg <- aggregate(torpor$Hours_torpid, 
                       by=list(torpor$Site_new, 
                               torpor$Species), 
                       FUN="mean", na.rm=T)
names(hours.agg) <- c("Site", "Species", "Mean torpor duration (hours)")
#hours.agg #to check

## Whole animal O2 consumption summary
o2.agg <- aggregate(torpor$O2_ml_min, 
                    by=list(torpor$Torpid_not, torpor$Site_new, 
                            torpor$Species), 
                    FUN="mean", na.rm=T)
names(o2.agg) <- c("Torpid_not", "Site", "Species", "Oxygen_ml_min")
#o2.agg #to check

## Make table to summarize savings
savings.agg <- aggregate(torpor$savings, 
                         by=list(torpor$Site_new, 
                                 torpor$Species), 
                         FUN="mean", na.rm=T)
names(savings.agg) <- c("Site", "Species", "Hourly Torpid Savings (%)")
savings.agg

## Mass summary just by species
mass.agg <- aggregate(torpor$Mass, 
                      by=list(torpor$Species), 
                      FUN="mean", na.rm=T)
names(mass.agg) <- c("Species", "Mass")
mass.agg

## Frequency table -  add proportion column
freq_sp$prop <- (freq_sp$Torpid/freq_sp$Total)*100
freq_sp$mass <- mass.agg$Mass

## Melt BBLH dataframe to put torpid and normo in same column
m_BBLH_tor_nor <- melt(BBLH_torpor, id.vars="Tc_min_C", 
                       measure.vars = c("AvgEE_torpid_MassCorrected", "AvgEE_normo_MassCorrected"))
levels(m_BBLH_tor_nor$variable)[levels(m_BBLH_tor_nor$variable)=="AvgEE_normo_MassCorrected"] <- 
  "Avg Normothermic EE"
levels(m_BBLH_tor_nor$variable)[levels(m_BBLH_tor_nor$variable)=="AvgEE_torpid_MassCorrected"] <- 
  "Avg Torpid EE"

m_BBLH_tor_nor$variable <- factor(m_BBLH_tor_nor$variable, 
                                  levels = rev(levels(m_BBLH_tor_nor$variable)))

## Melt BBLH dataframe to put torpid and normo in same column
m_GCB_tor_nor <- melt(GCB_torpor, id.vars="Tc_min_C", 
                      measure.vars = c("AvgEE_torpid_MassCorrected", "AvgEE_normo_MassCorrected"))
levels(m_GCB_tor_nor$variable)[levels(m_GCB_tor_nor$variable)=="AvgEE_normo_MassCorrected"] <- 
  "Avg Normothermic EE"
levels(m_GCB_tor_nor$variable)[levels(m_GCB_tor_nor$variable)=="AvgEE_torpid_MassCorrected"] <- 
  "Avg Torpid EE"

m_GCB_tor_nor$variable <- factor(m_GCB_tor_nor$variable,levels = 
                                   rev(levels(m_GCB_tor_nor$variable)))

## New aggregated data frame for rewarming
rewarm.agg <- aggregate(torpor$Rewarming_O2_rate, by=list(torpor$Species), 
                        FUN='mean', na.rm=T)

mass.agg2 <- aggregate(torpor$Mass, by=list(torpor$Species), 
                       FUN='mean', na.rm=T)

rewarm.agg <- merge(rewarm.agg, mass.agg2, by=c('Group.1'))
names(rewarm.agg) <- c("Species", "Rewarming_O2_rate", "Mass")

#### General functions ####
## Saving standard theme  
my_theme <- theme_classic(base_size = 30) + 
  theme(axis.title.y = element_text(color = "black", vjust = 2),
        panel.border = element_rect(colour = "black", fill=NA))

## Theme with slightly smaller font
my_theme2 <- my_theme + theme_classic(base_size = 15)

## To arrange graphs
lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 

# Function to return sample sizes
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

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

## Code for easy axis labels
Tc.xlab <- expression(atop(paste("Chamber Temperature (", degree,"C)")))
Ta.xlab <- expression(atop(paste("Ambient Temperature (", degree,"C)")))
Tc_min.xlab <- expression(atop(paste("Minimum Chamber Temperature (", degree,"C)")))
NEE_corrlab <- bquote('Nighttime energy expenditure (kJ/' ~M^(0.67)*')')
Nec_consump_lab <- bquote('Nectar consumption/' ~M^(0.67)*'')

#### Depth of entry into torpor plots #####
# Full site names for plots included in paper
torpor$Site_full <- torpor$Site_new
levels(torpor$Site_full) <- c("Harshaw", "Sonoita", "Southwest Research Station", "Maqui", "Santa Lucia")

## Supp docs figure S3
## Savings plot by species
savings_plot <- ggplot(torpor[!is.na(torpor$savings),], aes(Species2, savings)) + 
  geom_boxplot(outlier.shape = 19, fill= "light grey") + xlab("Species") + 
  ylab("Hourly torpid energy savings (%)") + theme(legend.position="none") + my_theme +
  stat_summary(fun.data = give.n, geom = "text", vjust=-1, size=10)
savings_plot

#### Exploratory plots ####
## Frequency of torpor use
freqplot <- ggplot(freq_table, aes(Temptrop, prop)) + geom_boxplot(fill= "light grey") + 
  ylab("Frequency of torpor use (%)") +  xlab("Region") + my_theme2 + 
  stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, size=5) +
  theme(axis.title.x = element_blank()) + ggtitle("a.")
freqplot

## Temp-trop Plot for Mass-corrected Nighttime energy expenditure
#heightDetails.gtable <- function(x) sum(x$heights)
energyM_temptrop <- ggplot(torpor, aes(Temptrop, NEE_MassCorrected)) + my_theme2 + 
  geom_boxplot(fill= "light grey") + 
  xlab("Region") + ylab(NEE_corrlab) + 
  theme(axis.title.y = element_text(vjust=-0.1)) +
  #(my_theme %+replace% theme(axis.title.y = element_custom())) +
  stat_summary(fun.data = give.n, geom = "text", vjust=-1, size=5) +
  theme(axis.title.x = element_blank()) + ggtitle("b.")
energyM_temptrop

## Hours torpid temptrop
hours_temptrop <- ggplot(na.omit(torpor[,c("Species","Hours_torpid", "Temptrop")]), 
                         aes(Temptrop, Hours_torpid)) + 
  geom_boxplot(fill= "light grey") + my_theme2 + ylab("Torpor duration (hours)") + xlab("Region") + 
  stat_summary(fun.data = give.n, geom = "text", size=5, vjust=-1) +
  theme(axis.title.x = element_blank()) + ggtitle("c.")
hours_temptrop

## Plot for proportion hours spent torpid
prop_hours_plot <- ggplot(na.omit(torpor[,c("Species","Hours_torpid","Site_new","Temptrop","Prop_hours")]), 
                          aes(Temptrop, as.numeric(as.character((Prop_hours))))) + 
  geom_boxplot(fill= "light grey") + my_theme2 + ylab("Percentage of hours spent torpid") + 
  xlab("Region") + stat_summary(fun.data = give.n, geom = "text", size=5, vjust=-1) +
  theme(axis.title.x = element_blank()) + ggtitle("d.")
prop_hours_plot

## Duration vs. time of entry
dur_entrytime <- ggplot(torpor, aes(EntryTime_numeric, Hours_torpid)) + 
  geom_point(size=2) + my_theme2 + xlab("Hour of entry") + ylab("Duration of torpor (hours)") +
  geom_smooth(method=lm, size=1, col="black") +
  geom_text(x = 5.5, y = 6, label = lm_eqn(torpor$Hours_torpid, torpor$EntryTime_numeric), parse=T, size=5)
dur_entrytime

rewarming <- ggplot(torpor, aes(Mass, Rewarming_O2_rate)) + 
  my_theme2 + geom_point(aes(col=Species), size=2) + xlab("Mass") +
  #ylab(bquote('Rate of rewarming ('~VO[2]~ ml/~min^2*')')) +
  scale_color_brewer(palette='Set1') +
  geom_smooth(method=lm, size=1, col="black") +
  geom_text(x = 5, y = 0.07,
            label = lm_eqn(log(torpor$Rewarming_O2_rate), torpor$Mass), parse=T, size=5)
rewarming

rewarming_sp <- ggplot(rewarm.agg, aes(Mass, Rewarming_O2_rate)) + 
  my_theme2 + geom_point(aes(col=Species), size=2) + xlab("Mass") +
  #ylab(bquote('Rate of rewarming ('~VO[2]~ ml/~min^2*')')) +
  scale_color_brewer(palette='Set1') +
  geom_smooth(method=lm, size=1, col="black")
rewarming_sp

rewarming_kJ <- ggplot(torpor, aes(Mass, kJ_rewarming_BeforeOvershoot)) + 
  my_theme2 + geom_point(aes(col=Species), size=3) + xlab("Mass") + ylab("Rewarming (kJ)") +
  #ylab(bquote('Rate of rewarming ('~VO[2]~ ml/~min^2*')')) +
  scale_color_brewer(palette='Set1') #+
  #geom_smooth(method=lm, size=1, col="black") #+
  #geom_text(x = 5, y = 0.07,
            #label = lm_eqn(log(torpor$Rewarming_VO2_rate), torpor$Mass), parse=T, size=5)
rewarming_kJ

bluecols <- brewer.pal(9, 'Blues')
newcol <- colorRampPalette(bluecols)
ncols <- 42
bluecols2 <- newcol(ncols)

rewarming_NEE <- ggplot(torpor, aes(Mass, NEE_without_rewarming_kJ)) + 
  my_theme2 + geom_point(aes(col=Rewarming_Tc), size=3, col=bluecols2) +
  geom_point(aes(Mass, NEE_kJ), col='black')  +
  geom_segment(aes(x=Mass,y=NEE_kJ,xend=Mass,yend=NEE_without_rewarming_kJ, col=Rewarming_Tc),
               col = bluecols2, 
               arrow = arrow(length = unit(0.01, "npc"),type = "closed"), alpha=0.5) +
  xlab("Mass") + ylab("Nighttime energy expenditure (kJ)") 
rewarming_NEE

rewarming_NEE <- ggplot(torpor[torpor$Torpid_not=="T",], aes(Mass, NEE_without_rewarming_kJ)) + 
  my_theme2 + geom_point(aes(col=Rewarming_Tc), size=3, alpha=0.5) +
  scale_colour_gradient(low="blue", high="red", name="Chamber \n temperature") +
  geom_point(aes(Mass, NEE_kJ), col='black')  +
  geom_segment(aes(x=Mass,y=NEE_kJ,xend=Mass,yend=NEE_without_rewarming_kJ, col=Rewarming_Tc),
               arrow = arrow(length = unit(0.01, "npc"),type = "closed"), alpha=0.5) +
  xlab("Mass") + ylab("Nighttime energy expenditure (kJ)") 
rewarming_NEE

rewarming_NEE_Tc <- ggplot(torpor[torpor$Torpid_not=="T",], 
                           aes(Rewarming_Tc, (NEE_kJ-NEE_without_rewarming_kJ))) + 
  my_theme2 + geom_point() +
  xlab(Tc.xlab) + ylab("Rewarming costs (kJ)")
rewarming_NEE_Tc

Nec_consump_lab <- bquote('Nectar consumption/' ~M^(0.67)*'')
## Savings temptrop
temptrop_savings <- ggplot(m.temptrop[m.temptrop$variable=="Percentage_avg",], 
                           aes(Temptrop, 100-(value))) + 
  geom_boxplot(fill= "light grey") + my_theme2 + ylab("Hourly torpid energy savings (%)") + #xlab("Region") +
  stat_summary(fun.data = give.n, geom = "text", size=5, vjust=-1.7) +
  theme(axis.title.x = element_blank()) + ggtitle("e.")
temptrop_savings

## Compiling Temperate tropical measures, comparison
grid.arrange(freqplot, energyM_temptrop, hours_temptrop, prop_hours_plot, temptrop_savings, 
             nrow=3, ncol=2, bottom = textGrob("Region", gp=gpar(fontsize=25), vjust=-0.5))

## BBLH and GCB plot, supp Figure S4
## Both normo and torpid avg EE for BBLH on same graph
BBLH_tor_nor <- ggplot(m_BBLH_tor_nor, aes(as.numeric(Tc_min_C), value, shape=variable)) +
  geom_point(aes(shape=variable), size=6) + my_theme2 +
  geom_smooth(method=lm, size=1, col="black") + 
  scale_shape_manual("Hourly Energy Expenditure\n", 
                     values=c(16,1), labels=c("Normothermic", "Torpid")) +
  theme(legend.key.height=unit(3,"line"), legend.position="bottom", axis.title.y = element_text(vjust = 2)) +
  ylab("Mean BBLH Energy Expenditure (kJ/g)") + xlab(Tc.xlab) + ggtitle("b.")
BBLH_tor_nor 

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

grid_arrange_shared_legend(GCB_tor_nor, BBLH_tor_nor)