## Torpor paper, A. Shankar, R. Schroeder et al.
## Plots and analyses for paper on hummingbird torpor in Arizona and Ecuador

####### libraries and reading in data ######
library(ggplot2)
library(reshape)
library(gridExtra)
library(grid)
library(wq)
library(gam)
library(foreign)
library(MASS)
library(ggbiplot)
library(dplyr)
library(rgl)

## setwd and read in file
#wdCH
setwd("C:\\Users\\shankar\\Dropbox\\Hummingbird energetics\\Submission_Oct2016")
#wdMS
#setwd("C:\\Users\\ANUSHA\\Dropbox\\Hummingbird energetics\\Tables_for_paper")

## Used a European computer last- remove the sep=";" if using a csv format where commas are separators.
torpor <- read.csv("Torpor_individual_summaries.csv", sep=";")
freq_sites <- read.csv("Frequency_torpor_sites.csv")
freq_sp <- read.csv("Frequency_torpor_sp.csv")

tempsumm <- read.csv("Temp_summary.csv")

torpor$Percentage_avg <- as.numeric(as.character(torpor$Percentage_avg))
torpor$Prop_hours <- as.numeric(as.character(torpor$Prop_hours))

##### Adding columns #######
## Adding column dividing NEE by 2/3*Mass to correct for mass with allometric scaling
torpor$NEE_MassCorrected<- torpor$NEE_kJ/(torpor$Mass^(2/3))

## Adding columns to correct for mass in Avg EE normo, Min EE normo, torpid, etc. 
torpor$AvgEE_normo_MassCorrected <- torpor$Avg_EE_hourly_normo/(torpor$Mass^(2/3))
torpor$MinEE_normo_MassCorrected <- as.numeric(torpor$Min_EE_normo)/(torpor$Mass^(2/3))
torpor$AvgEE_torpid_MassCorrected <- torpor$Avg_EE_hourly_torpid/(torpor$Mass^(2/3))
torpor$MinEE_torpid_MassCorrected <- as.numeric(torpor$Min_EE_torpid)/(torpor$Mass^(2/3))

# Line to arrange Site facets in sensible order
torpor$Site_new <- factor(torpor$Site, levels=c('HC','SC','SWRS','MQ','SL'))

torpor$EntryTime_new <- factor(torpor$Time_of_entry, 
                                   levels=c('2000', '2100', '2130', '2200', '2230',
                                            '2300', '2400', '100', '130', '200', '330', '500'))

freq_sites$Site_new <- factor(freq_sites$Site, levels=c('HC','SC','SWRS','MQ','SL'))

## Create Hours_torpid2 column and change NA's to 0's for lm analyses; keep original Hours_torpid column with NA's
torpor$Hours_torpid2 <- torpor$Hours_torpid
torpor$Hours_torpid2[is.na(torpor$Hours_torpid2)] <- 0

## Savings column to convert percentage energy expended in torpor relative to normothermy into savings relative to normothermy
torpor$savings <- 100-torpor$Percentage_avg

#### Make new data frames ####
## Melt into temperate-tropical format
m.temptrop <- melt(torpor, id.vars = c("Temptrop", "Species", "Site"), 
                   measure.vars =  c("Hours_torpid", "Prop_hours", "NEE_kJ", "Percentage_avg",
                                     "NEE_MassCorrected", "Mass"))

## Make dataframe to check site-wise frequency of torpor
rate_site <- data.frame(table(torpor$Site,torpor$Torpid_not))
names(rate_site) <- c("Site", "Torpid_not", "N")

## Subset just BBLH data
BBLH_torpor <- subset(torpor, Species=="BBLH")

## Subset just GCB data
GCB_torpor <- subset(torpor, Species=="GCB")

## Make table to summarize nighttime energy expenditure
nee.agg <- aggregate(torpor$NEE_kJ, 
                     by=list(torpor$Torpid_not, torpor$Site_new, 
                             torpor$Species), 
                     FUN="mean", na.rm=T)
names(nee.agg) <- c("Torpid_not", "Site", "Species", "NEE_kJ")
nee.agg

## Mass-corrected NEE aggregate
neet.agg <- aggregate(torpor$NEE_MassCorrected, 
                      by=list(torpor$Torpid_not, torpor$Site_new, 
                              torpor$Species), 
                      FUN="mean", na.rm=T)
names(neet.agg) <- c("Torpid_not", "Site", "Species", "NEE Mass corrected (kJ/g^(2/3))")
neet.agg

## Aggregate mass
mass1.agg <- aggregate(torpor$Mass, 
                       by=list(torpor$Torpid_not, torpor$Site_new, 
                               torpor$Species), 
                       FUN="mean", na.rm=T)
names(mass1.agg) <- c("Torpid_not", "Site", "Species", "Mass")
mass1.agg

nee_Mc_agg <- merge(nee.agg, neet.agg,by=c("Torpid_not", "Site", "Species"))
nee_Mc_agg <- merge(nee_Mc_agg, mass1.agg, by=c("Torpid_not", "Site", "Species"))

write.csv(nee_Mc_agg, "NEE and mass-corrected summary.csv")

## Summarize hours spent torpid
hours.agg <- aggregate(torpor$Hours_torpid, 
                       by=list(torpor$Site_new, 
                               torpor$Species), 
                       FUN="mean", na.rm=T)
names(hours.agg) <- c("Site", "Species", "Mean torpor duration (hours)")
hours.agg

## Whole animal O2 consumption summary
o2.agg <- aggregate(torpor$O2_ml_min, 
                    by=list(torpor$Torpid_not, torpor$Site_new, 
                            torpor$Species), 
                    FUN="mean", na.rm=T)
names(o2.agg) <- c("Torpid_not", "Site", "Species", "Oxygen_ml_min")
o2.agg

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

## Codes for easy axis labels
Tc.xlab <- expression(atop(paste("Chamber Temperature (", degree,"C)")))
Ta.xlab <- expression(atop(paste("Ambient Temperature (", degree,"C)")))
Tc_min.xlab <- expression(atop(paste("Minimum Chamber Temperature (", degree,"C)")))
NEE_corrlab <- bquote('Nighttime energy expenditure (kJ/' ~M^(0.67)*')')
Nec_consump_lab <- bquote('Nectar consumption/' ~M^(0.67)*'')

#### Depth plots #####
torpor$Species2 <- factor(torpor$Species,
                          levels = c('BBLH','MAHU','GCB','FBB','TBH', "WNJ"), ordered = T)

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

#### Comparing temperate and tropical species ####
### Each of the 5 facets for Figure 2 
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
  geom_boxplot(fill= "light grey") + my_theme2 + ylab("Percentage of hours spent torpid") + xlab("Region") +
  stat_summary(fun.data = give.n, geom = "text", size=5, vjust=-1) +
  theme(axis.title.x = element_blank()) + ggtitle("d.")
prop_hours_plot

## Savings temptrop
temptrop_savings <- ggplot(m.temptrop[m.temptrop$variable=="Percentage_avg",], 
                           aes(Temptrop, 100-(value))) + 
  geom_boxplot(fill= "light grey") + my_theme2 + ylab("Hourly torpid energy savings (%)") + #xlab("Region") +
  stat_summary(fun.data = give.n, geom = "text", size=5, vjust=-1.7) +
  theme(axis.title.x = element_blank()) + ggtitle("e.")
temptrop_savings

## Compiling Figure 2
grid.arrange(freqplot, energyM_temptrop, hours_temptrop, prop_hours_plot, temptrop_savings, 
             nrow=3, ncol=2, bottom = textGrob("Region", gp=gpar(fontsize=25), vjust=-0.5))



#### Additional plots not in the paper ####
## Testing the effect of mass on hourly energy expenditure in torpor vs. normothermy
summary(lm(torpor$Percentage_avg[!is.na(torpor$Percentage_avg)] ~ torpor$Mass[!is.na(torpor$Percentage_avg)]))

## And plot the regression
savings_mass <- ggplot(torpor[!is.na(torpor$savings),], aes(Mass, savings)) + 
  geom_point(aes(col=Tc_min_C), size=3) + my_theme + geom_smooth(method = lm, col='black') +
  ylab("Hourly torpid energy savings (%)")
savings_mass

#### Comparing temperate and tropical species ####
## Plot for Nighttime energy expenditure, by temperate-tropics
energy_plot <- ggplot(torpor, aes(Temptrop, NEE_kJ)) + my_theme + geom_boxplot() + xlab("Region") +
  ylab("Nighttime energy expenditure (kJ)") + theme(legend.position="none") + 
  stat_summary(fun.data = give.n, geom = "text", vjust=-1, size=10)
energy_plot

#### Basic frequency, NEE and hours plots ####
## Rate of occurrence vs. mass - seeing if mass is more a determinant than just role
freqmass <- ggplot(freq_sp, aes(mass, prop)) + geom_point() + geom_smooth(method = lm, col='black') +
  ylab("Rate of occurrence of torpor (%)") +  xlab("Mass") + my_theme
freqmass

summary(lm(freq_sp$Freq_torpor ~ freq_sp$mass))

## Testing if nectar consumption might have affected duration of torpor
nec_duration <- ggplot(torpor[!is.na(torpor$Hours_torpid2),], 
                       aes(Hours_torpid2, (Nectar_consumption/(Mass^(0.67))))) + my_theme + 
  geom_point(aes(col=Species), size=4) + geom_smooth(method=lm, color="black") + scale_color_brewer(palette = "Set1") +
  geom_text(x = 4, y = 0.32, label = lm_eqn((torpor$Nectar_consumption/(torpor$Mass^0.67)), torpor$Hours_torpid2), 
            parse=T, size=8) + theme(legend.key.height=unit(3,"line")) +
  xlab("Torpor duration") +  ylab(Nec_consump_lab)
nec_duration

## Nectar consumption vs. hour of entry into torpor
nec_time <- ggplot(torpor[!is.na(torpor$Time_numeric),], aes(Time_numeric, (Nectar_consumption/(Mass^(0.67))))) + 
  my_theme + geom_point(aes(col=Species), size=5) + geom_smooth(method=lm, color="black") +
  scale_color_brewer(palette = "Set1") +
  geom_text(x = 4, y = 0.32, label = lm_eqn((torpor$Nectar_consumption/(torpor$Mass^0.67)), torpor$Time_numeric), 
            parse=T, size=8) + theme(legend.key.height=unit(3,"line")) +
  xlab("Hour of entry") +  ylab(Nec_consump_lab)
nec_time

## Frequency of torpor use
freq <- ggplot(freq_table, aes(Species, prop)) + geom_bar(stat="identity") + 
  geom_text(data=freq_table,aes(x=Species,y=prop,label=paste("n = ", Total)), vjust=-0.5, size=10) +
  ylab("Rate of occurrence of torpor (%)") +  xlab("Species") + my_theme + ylim(0, 109)
freq

freq_site_plot <- ggplot(freq_sites, aes(Species, prop)) + 
  geom_bar(stat="identity", width=0.5, position = position_dodge(width=0.5)) + 
  facet_grid(~Site_new, scale="free_x", space="free") +
  geom_text(data=freq_sites, aes(x=Species,y=prop, label=paste("n = ", Total)), vjust=-0.5, size=8) +
  ylab("Rate of occurrence of torpor (%)") +  xlab("Species") + my_theme + ylim(0, 109)
freq_site_plot

## Hours torpid
hours_plot <- ggplot(na.omit(torpor[,c("Species","Hours_torpid","Site_full")]),
                     aes(Species, Hours_torpid)) + 
  geom_boxplot(outlier.size = 3, fill="light grey") + my_theme +
  facet_grid(~Site_full, labeller = labeller(Site_full = label_wrap_gen(10)), scale="free_x", space="free") +
  ylab("Hours Torpid") + theme(legend.position="none") + ggtitle("b.") +
  stat_summary(position = position_nudge(y = 0.98), fun.data = give.n, geom = "text", size=8)
hours_plot

## Plot for Mass-corrected Nighttime energy expenditure, by species
energyM_plot <- ggplot(torpor, aes(Species, NEE_MassCorrected)) +  my_theme +
  geom_boxplot(aes(col=Species), size=1) + facet_grid(.~Site_new, scale="free_x", space="free") + 
  ylab(NEE_corrlab) + theme(legend.position="none") +
  stat_summary(fun.data = give.n, geom = "text", vjust=-5, size=10)
energyM_plot

## Just BBLH, nighttime EE vs. site
energyM_BBLH <- ggplot(torpor[torpor$Species=="BBLH",], aes(Site_new, NEE_MassCorrected)) +  my_theme +
  geom_boxplot(size=1) + 
  ylab(NEE_corrlab) + theme(legend.position="none") +
  stat_summary(fun.data = give.n, geom = "text", vjust=-5, size=10)
energyM_BBLH

## Energy vs. hours torpid, species labeled
energyM_hours <- ggplot(torpor, aes(Hours_torpid2, NEE_MassCorrected)) +  my_theme +
  geom_point(aes(shape = factor(Species)), size=4) + theme_bw(base_size=30) +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  geom_smooth(method=lm, color="black") +
  geom_text(x = 5, y = 4.5, label = lm_eqn(torpor$NEE_MassCorrected, torpor$Hours_torpid2), 
            parse=T, size=10) +
  labs(shape='Species') + scale_color_brewer(palette = "Set1") + theme_bw(base_size=30) +
  ylab(NEE_corrlab) + xlab("Torpor duration")
energyM_hours

## Energy vs. hours torpid, species labeled
BBLH_energyM_hours <- ggplot(torpor[torpor$Species=="BBLH",], aes(Hours_torpid, NEE_MassCorrected)) +  my_theme +
  geom_point(aes(shape = factor(Species)), size=4) + theme_bw(base_size=30) +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  geom_smooth(method=lm, color="black") +
  geom_text(x = 4, y = 2.5, label = lm_eqn(torpor[torpor$Species=="BBLH",], torpor$NEE_MassCorrected, torpor$Hours_torpid), 
            parse=T, size=10) +
  labs(shape='Species') + scale_color_brewer(palette = "Set1") + theme_bw(base_size=30) +
  ylab(NEE_corrlab) + xlab("Torpor duration")
BBLH_energyM_hours

## Energy vs. hours torpid, without species labeled- for retreat
energy_hours_spUnlabeled <- ggplot(torpor, aes(Hours_torpid2, NEE_MassCorrected)) +  
  geom_point(size=4) + my_theme + geom_smooth(method='loess', color="black") +
  geom_text(x = 7, y = 6, label = paste("R^2 :", " 0.51", sep=""), parse=T, size=10) +
  geom_text(x = 7, y = 5, label = paste("N =", " 42", sep=""), parse=F, size=10) +
  labs(shape='Species') + scale_color_brewer(palette = "Set1") +
  ylab(NEE_corrlab) + xlab("Torpor duration")
energy_hours_spUnlabeled

## Comparing NEE and hours plots
grid.arrange(energyM_plot, hours_plot, nrow=1, ncol=2)

## Comparing energy plots with and without mass-correction
grid.arrange(energy_plot, energyM_plot, nrow=1, ncol=2)

## NEE plot by temperature
energy_temp <- ggplot(torpor, aes(as.numeric(Tc_mean_C), NEE_kJ)) + 
  geom_point(aes(shape = factor(Species)), size=4) + 
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  labs(shape='Species') + xlim(4, 30) +
  scale_color_brewer(palette = "Set1") + theme_bw() + 
  geom_text(aes(label=Torpid_not, hjust=2), size=5, show_guide=F,
            fontface="bold") +
  ylab("Nighttime energy expenditure (kJ)") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
energy_temp 

## Mass-corrected NEE plot by temperature
energyM_temp <- ggplot(torpor, aes(as.numeric(Tc_mean_C), NEE_MassCorrected)) + 
  geom_point(aes(shape = factor(Species)), size=4) + 
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  #scale_shape_manual(values=1:nlevels(torpor$Species)) +
  labs(shape='Species') + xlim(0, 30) +
  scale_color_brewer(palette = "Set1") + theme_bw() + 
  geom_text(aes(label=Torpid_not, hjust=2), size=5, show_guide=F,
            fontface="bold") +
  ylab("Nighttime energy expenditure Mass-corrected (kJ/g)") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
energyM_temp 

## Mass-corrected NEE plot by MINIMUM Tc
energyM_Tcmin <- ggplot(torpor, aes(as.numeric(Tc_min_C), NEE_MassCorrected)) + 
  geom_point(aes(shape = factor(Species)), size=4) + 
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  geom_smooth(method=lm, color="black") +
  #geom_text(x = 16, y = 4.5, label = lm_eqn(torpor, torpor$NEE_MassCorrected, 
  #                                         torpor$Tc_min_C), parse=T) +
  labs(shape='Species') + #xlim(5, 25) +
  scale_color_brewer(palette = "Set1") + theme_bw() + 
  #geom_text(aes(label=Torpid_not, hjust=2), size=5, show_guide=F,
  #        fontface="bold") +
  ylab("Nighttime energy expenditure Mass-corrected (kJ/g)") + xlab(Tc_min.xlab) + my_theme 
energyM_Tcmin

## NEE vs. temp, with and without mass-correction
grid.arrange(energy_temp, energyM_temp, nrow=1, ncol=2)

## NEE plot by chamber temperature, facet by site and color by species
energy_temp_site <- ggplot(torpor, aes(as.numeric(Tc_mean_C), NEE_kJ)) + geom_boxplot() +
  geom_point(aes(shape = factor(Species)), size=4) + 
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  #scale_shape_manual(values=1:nlevels(torpor$Species)) +
  labs(shape='Species') + xlim(-7, 45) + theme_bw() + 
  #scale_color_brewer(palette = "Set1") 
  geom_text(aes(label=Torpid_not, hjust=1.75, fontface="bold"),size=5) +
  facet_grid(.~Site_new) +
  ylab("Nighttime energy expenditure (kJ)") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
energy_temp_site

## NEE plot by chamber temperature, facet by site and color by species
NEE_site <- ggplot(torpor, aes(Site_new, NEE_MassCorrected)) + geom_boxplot(outlier.size = 3) +
  theme_classic(base_size=30) +
  ylab("Nighttime energy expenditure (kJ)") + xlab("Site") +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA)) 
NEE_site

Hours_site <- ggplot(torpor, aes(Site_new, Hours_torpid)) + geom_boxplot(outlier.size = 3) +
  theme_classic(base_size=30) +
  ylab("Hours torpid") + xlab("Site") +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA))
Hours_site

## Plot both Site plots together
lay_out(list(NEE_site, 1, 2), 
        list(Hours_site, 1, 1))

## Mass-corrected NEE plot by chamber temperature, facet by site and color by species
energyM_temp_site <- ggplot(torpor, aes(as.numeric(Tc_mean_C), NEE_MassCorrected)) + 
  geom_point(aes(shape = factor(Species)), size=4) + 
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  #scale_shape_manual(values=1:nlevels(torpor$Species)) +
  labs(shape='Species') + xlim(-7, 45) +
  scale_color_brewer(palette = "Set1") + theme_bw() + 
  geom_text(aes(label=Torpid_not, hjust=1.75, fontface="bold"),size=5) +
  facet_grid(.~Site_new) +
  ylab("Nighttime energy expenditure Mass-corrected (kJ/g)") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
energyM_temp_site

## Plot both energy_temp with and without mass-correction together. i.e. 4 plots
lay_out(list(energy_temp, 1, 1), 
        list(energy_temp_site, 1, 2),
        list(energyM_temp, 2, 1),
        list(energyM_temp_site, 2, 2))

#### Non-mass-corrected min and avg graphs ####
## Min normo EE by Tc
min_normo_EE <- ggplot(torpor, aes(as.numeric(Tc_mean_C), as.numeric(Min_EE_normo))) + theme_bw() + 
  geom_point(aes(shape = factor(Species)), size=4) +  labs(shape ='Species') +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + xlim(5, 32) +
  geom_text(aes(label=Torpid_not, hjust=1.75, fontface="bold"),size=5) +
  ylab("Min EE normothermic") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
min_normo_EE

## Minimum hourly energy expenditure while torpid by Tc
min_torpid_EE <- ggplot(torpor, aes(as.numeric(Tc_mean_C), as.numeric(Min_EE_torpid))) + 
  theme_bw() + geom_point(aes(shape = factor(Species)), size=4) + labs(shape ='Species') +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + #xlim(-7, 50) +
  ylab("Min EE torpid") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
min_torpid_EE

## Avg EE normo, with points labeled as Normo or Torpid birds; to see if overall EE is lower or
## higher for birds that tend to go into torpor
avg_normo_EE <- ggplot(torpor, aes(as.numeric(Tc_mean_C), Avg_EE_hourly_normo)) +  theme_bw() + 
  geom_point(aes(shape = factor(Species)), size=4) + labs(shape ='Species') +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + xlim(5, 30) +
  #facet_grid(.~Site,space="free") +
  geom_text(aes(label=Torpid_not, hjust=1.75, fontface="bold"),size=5) +
  ylab("Avg EE normothermic") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
avg_normo_EE

## Average hourly energy expenditure while torpid
avg_torpid_EE <- ggplot(torpor, aes(Tc_mean_C, Avg_EE_hourly_torpid)) +  theme_bw() + 
  geom_point(aes(shape = factor(Species)), size=4) + labs(shape ='Species') +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + #xlim(5, 30) +
  #facet_grid(.~Site,space="free") +
  ylab("Avg EE torpid") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
avg_torpid_EE

#### Mass-corrected Min and avg graphs ####
## Mass-corrected Min normo EE by Tc
m_min_normo_EE <- ggplot(torpor, aes(as.numeric(Tc_mean_C), 
                                     as.numeric(MinEE_normo_MassCorrected))) +  
  theme_bw() + geom_point(aes(shape = factor(Species)), size=4) +  labs(shape ='Species') +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + xlim(5, 32) +
  geom_text(aes(label=Torpid_not, hjust=1.75, fontface="bold"),size=5) +
  ylab("Min EE normothermic (kJ/g)") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
m_min_normo_EE

## Mass-corrected Minimum hourly energy expenditure while torpid by Tc
m_min_torpid_EE <- ggplot(torpor, aes(as.numeric(Tc_mean_C), 
                                      as.numeric(MinEE_torpid_MassCorrected))) + 
  theme_bw() + geom_point(aes(shape = factor(Species)), size=4) + labs(shape ='Species') +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  geom_smooth(method=lm, color="black") +
  geom_text(x = 16, y = 0.6, label = lm_eqn(torpor, torpor$AvgEE_normo_MassCorrected, 
                                            torpor$Tc_min_C), parse=T) +
  scale_color_brewer(palette = "Set1") + #xlim(-7, 50) +
  ylab("Min EE torpid (kJ/g)") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
m_min_torpid_EE

## Mass-corrected Average hourly energy expenditure while normothermic
m_avg_normo_EE <- ggplot(torpor, aes(as.numeric(Tc_mean_C), AvgEE_normo_MassCorrected)) + 
  theme_bw() + geom_point(aes(shape = factor(Species)), size=4) + labs(shape ='Species') +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + xlim(5, 30) +
  #facet_grid(.~Site,space="free") +
  geom_text(aes(label=Torpid_not, hjust=1.75, fontface="bold"),size=5) +
  ylab("Avg EE normothermic (kJ/g)") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
m_avg_normo_EE

## Mass-corrected Average hourly energy expenditure while torpid
m_avg_torpid_EE <- ggplot(torpor, aes(Tc_mean_C, AvgEE_torpid_MassCorrected)) +  theme_bw() + 
  geom_point(aes(shape = factor(Species)), size=4) + labs(shape ='Species') +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + #xlim(5, 30) +
  #facet_grid(.~Site,space="free") +
  ylab("Avg EE torpid (kJ/g)") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
m_avg_torpid_EE

## Arranging mass-corrected graphs
grid.arrange(m_avg_normo_EE, m_avg_torpid_EE, 
             m_min_normo_EE, m_min_torpid_EE, nrow=2, ncol=2) # avg and min EE
lay_out(list(m_avg_normo_EE_labeled, 1, 1), 
        list(m_avg_torpid_EE, 2, 1)) # Only average EEs

lay_out(list(m_avg_normo_EE, 1, 1),
        list(m_avg_torpid_EE, 1, 2),
        list(m_min_normo_EE, 2, 1), 
        list(m_min_torpid_EE, 2, 2))

#### Mass-corrected mins and avgs, plotted against Tc min instead of Tc mean with regressions ####
## Mass-corrected Min normo EE by Tc min
m_min_normo_EE_Tcmin_eq <- ggplot(torpor, aes(as.numeric(Tc_min_C), 
                                              as.numeric(MinEE_normo_MassCorrected))) +  
  theme_bw() + geom_point(aes(shape = factor(Species)), size=4) +  labs(shape ='Species') +
  geom_smooth(method=lm, color="black") +
  geom_text(x = 19, y = 0.38, label = lm_eqn(torpor, torpor$MinEE_normo_MassCorrected, 
                                             torpor$Tc_min_C), parse=T) +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + xlim(0, 32) +
  geom_text(aes(label=Torpid_not, hjust=1.75, fontface="bold"),size=5) +
  ylab("Min EE normothermic (kJ/g)") + xlab(Tc_min.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
m_min_normo_EE_Tcmin_eq

## Mass-corrected Minimum hourly energy expenditure while torpid by Tc
m_min_torpid_EE_Tcmin_eq <- ggplot(torpor, aes(as.numeric(Tc_min_C), 
                                               as.numeric(MinEE_torpid_MassCorrected))) + 
  theme_bw() + geom_point(aes(shape = factor(Species)), size=4) + labs(shape ='Species') +
  geom_smooth(method=lm, color="black") +
  geom_text(x = 15, y = 0.06, label = lm_eqn(torpor, torpor$MinEE_torpid_MassCorrected, 
                                             torpor$Tc_min_C), parse=T) +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + #xlim(-7, 50) +
  ylab("Min EE torpid (kJ/g)") + xlab(Tc_min.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
m_min_torpid_EE_Tcmin_eq

## Mass-corrected Average hourly energy expenditure while normothermic
m_avg_normo_EE_Tcmin_eq <- ggplot(torpor, aes(as.numeric(Tc_min_C), AvgEE_normo_MassCorrected)) + 
  theme_bw(base_size=30) + geom_point(aes(shape = factor(Species)), size=4) + 
  labs(shape ='Species') + geom_smooth(method=lm, color="black") +
  geom_text(x = 16, y = 0.6, label = lm_eqn(torpor, torpor$AvgEE_normo_MassCorrected, 
                                            torpor$Tc_min_C), parse=T, size=10) +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + #xlim(0, 30) +
  #facet_grid(.~Site,space="free") +
  geom_text(aes(label=Torpid_not, hjust=1.75, fontface="bold"),size=5) +
  ylab("Avg EE normothermic (kJ/g)") + xlab(Tc_min.xlab) +
  theme(axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"), legend.key.height=unit(3,"line")) 
m_avg_normo_EE_Tcmin_eq

## Mass-corrected Average hourly energy expenditure while torpid
m_avg_torpid_EE_Tcmin_eq <- ggplot(torpor, aes(Tc_min_C, AvgEE_torpid_MassCorrected)) + theme_bw() + 
  geom_point(aes(shape = factor(Species)), size=4) + labs(shape ='Species') +
  geom_smooth(method=lm, color="black") +
  #geom_text(x = 13, y = 0.095, label = lm_eqn(torpor, torpor$AvgEE_torpid_MassCorrected, 
  #                                         torpor$Tc_min_C), parse=T) +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + 
  #facet_grid(.~Site,space="free") +
  ylab("Avg EE torpid (kJ/g)") + xlab(Tc_min.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
m_avg_torpid_EE_Tcmin_eq

## Arranging mass-corrected vs. Tc min graphs
lay_out(list(m_avg_normo_EE_Tcmin_eq, 1, 1),
        list(m_avg_torpid_EE_Tcmin_eq, 1, 2),
        list(m_min_normo_EE_Tcmin_eq, 2, 1), 
        list(m_min_torpid_EE_Tcmin_eq, 2, 2))

#### Overall regressions ####

## All species, Avg normo vs Tc regression
m_AvgEE_normo_Tcmin_eq <- ggplot(torpor, aes(as.numeric(Tc_mean_C), AvgEE_normo_MassCorrected)) + 
  theme_bw() + geom_point(aes(shape = factor(Species)), size=4) + labs(shape ='Species') +
  geom_smooth(method=lm, color="black") +
  geom_text(x = 16, y = 0.6, label = lm_eqn(torpor, torpor$AvgEE_normo_MassCorrected, 
                                            torpor$Tc_min_C), parse=T) +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + xlim(7, 30) +
  geom_text(aes(label=Torpid_not, hjust=1.75, fontface="bold"),size=5) +
  ylab("Avg EE normothermic (kJ/g)") + xlab(Tc.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
m_AvgEE_normo_Tcmin_eq

#### Species- specific graphs ####
## Comparing HC and SC BBLH energy expenditure in torpor
ggplot(BBLH_torpor, aes(Site, AvgEE_torpid_MassCorrected)) + my_theme +
  geom_boxplot() + geom_point(aes(col=Tc_min_C), size=5) + 
  scale_colour_gradient(low = "blue", high = "red", "Min chamber temperature\n") +
  ylab("Hourly Energy expenditure in torpid birds (J/h*g)")

ggplot(BBLH_torpor, aes(Site, AvgEE_normo_MassCorrected)) + my_theme +
  geom_boxplot() + geom_point(aes(col=Tc_min_C), size=5) + 
  scale_colour_gradient(low = "blue", high = "red", "Min chamber temperature\n") +
  ylab("Hourly Energy expenditure in normothermic birds (J/h*g)")

# BBLH Avg mass-corrected hourly normothermic EE vs. min Tc with regression line
m_BBLH_avgEE_normo_Tcmin_eq <- ggplot(BBLH_torpor, aes(as.numeric(Tc_min_C), 
                                                       AvgEE_normo_MassCorrected)) + 
  theme_bw() + geom_point(size=4) + geom_smooth(method=lm, color="black") +
  geom_text(x = 20, y = 0.38, label = paste("R^2 :", " 0.291",sep=""), parse=T, size=8) +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + #xlim(0, 30) +
  geom_text(aes(label=Torpid_not, hjust=1.75, fontface="bold"),size=5) +
  ylab("Avg BBLH normothermic EE (kJ/g)") + xlab(Tc_min.xlab) +
  theme(axis.title.x = element_text(size=24, face="bold"),
        axis.text.x = element_text(size=22),
        axis.title.y = element_text(size=24, face="bold"), axis.text.y = element_text(size=24)) 
m_BBLH_avgEE_normo_Tcmin_eq

## BBLH Avg mass-corrected hourly torpid EE vs. min Tc with regression line
m_BBLH_avgEE_torpid_Tcmin_eq <- ggplot(BBLH_torpor, aes(as.numeric(Tc_min_C), 
                                                        AvgEE_torpid_MassCorrected)) + 
  my_theme + geom_point(size=4) + geom_smooth(method=lm) + facet_grid(~Site) +
  geom_text(x = 14, y = 0.07, label = lm_eqn(BBLH_torpor, BBLH_torpor$AvgEE_torpid_MassCorrected, 
                                             BBLH_torpor$Tc_min_C), parse=T, size=8) +
  ylab("Avg BBLH EE torpid (kJ/g)") + xlab(Tc_min.xlab) 
m_BBLH_avgEE_torpid_Tcmin_eq

## Melt BBLH dataframe to put torpid and normo in same column
m_BBLH_tor_nor <- melt(BBLH_torpor, id.vars="Tc_min_C", 
                       measure.vars = c("AvgEE_torpid_MassCorrected", "AvgEE_normo_MassCorrected"))
levels(m_BBLH_tor_nor$variable)[levels(m_BBLH_tor_nor$variable)=="AvgEE_normo_MassCorrected"] <- 
  "Avg Normothermic EE"
levels(m_BBLH_tor_nor$variable)[levels(m_BBLH_tor_nor$variable)=="AvgEE_torpid_MassCorrected"] <- 
  "Avg Torpid EE"

m_BBLH_tor_nor$variable <- factor(m_BBLH_tor_nor$variable, 
                                  levels = rev(levels(m_BBLH_tor_nor$variable)))

## Both normo and torpid avg EE for BBLH on same graph
BBLH_tor_nor <- ggplot(m_BBLH_tor_nor, aes(as.numeric(Tc_min_C), value, shape=variable)) +
  geom_point(aes(shape=variable), size=6) + my_theme2 +
  geom_smooth(method=lm, size=1, col="black") + 
  scale_shape_manual("Hourly Energy Expenditure\n", 
                     values=c(16,1), labels=c("Normothermic", "Torpid")) +
  theme(legend.key.height=unit(3,"line"), legend.position="bottom", axis.title.y = element_text(vjust = 2)) +
  ylab("Mean BBLH Energy Expenditure (kJ/g)") + xlab(Tc.xlab) + ggtitle("b.")
BBLH_tor_nor 

grid.arrange(m_BBLH_avgEE_normo_Tcmin_eq, m_BBLH_avgEE_torpid_Tcmin_eq, nrow=1, ncol=2)

## BBLH comparison in color
BBLH_tor_nor_col <- ggplot(m_BBLH_tor_nor, aes(as.numeric(Tc_min_C), value, color=variable)) +
  my_theme + geom_point(aes(shape=variable), size=6) + 
  geom_smooth(method=lm, size=2) + scale_color_manual(values=c("#0099ff", "#ff0000")) +
  geom_text(x = 22, y = 0.32, label = paste("R^2 :", " 0.0302",sep=""), 
            parse=T, size=8, col="black") + 
  geom_text(x = 22, y = 0.12, label = paste("R^2 :", " 0.343",sep=""), 
            parse=T, size=8, col="black") +
  ylab("Avg BBLH Energy Expenditure (kJ/g)") + xlab(Tc.xlab) +
  theme(legend.key.height=unit(3,"line"), legend.title=element_text(size=20))
BBLH_tor_nor_col

## BBLH Min mass-corrected hourly normothermic EE vs. min Tc with regression line
m_BBLH_minEE_normo_Tcmin_eq <- ggplot(BBLH_torpor, aes(as.numeric(Tc_min_C), 
                                                       MinEE_normo_MassCorrected)) + 
  theme_bw() + geom_point() + geom_smooth(method=lm, color="black") +
  geom_text(x = 16, y = 0.32, label = lm_eqn(BBLH_torpor, BBLH_torpor$MinEE_normo_MassCorrected, 
                                             BBLH_torpor$Tc_min_C), parse=T) +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + #xlim(0, 30) +
  geom_text(aes(label=Torpid_not, hjust=1.75, fontface="bold"),size=5) +
  ylab("Minimum BBLH normothermic EE (kJ/g)") + xlab(Tc_min.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
m_BBLH_minEE_normo_Tcmin_eq

## BBLH Min mass-corrected hourly torpid EE vs. min Tc with regression line
m_BBLH_minEE_torpid_Tcmin_eq <- ggplot(BBLH_torpor, aes(as.numeric(Tc_min_C), 
                                                        MinEE_torpid_MassCorrected)) + 
  theme_bw() + geom_point() + geom_smooth(method=lm) +
  geom_text(x = 14, y = 0.024, label = lm_eqn(BBLH_torpor, BBLH_torpor$MinEE_torpid_MassCorrected, 
                                              BBLH_torpor$Tc_min_C), parse=T, size=6) +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + #xlim(0, 30) +
  geom_text(aes(label=Torpid_not, hjust=1.75, fontface="bold"),size=5) +
  ylab("Minimum BBLH EE torpid (kJ/g)") + xlab(Tc_min.xlab) +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16, face="bold"), axis.text.y = element_text(size=14)) 
m_BBLH_minEE_torpid_Tcmin_eq

lay_out(list(m_BBLH_avgEE_normo_Tcmin_eq, 1, 1),
        list(m_BBLH_avgEE_torpid_Tcmin_eq, 1, 2),
        list(m_BBLH_minEE_normo_Tcmin_eq, 2, 1), 
        list(m_BBLH_minEE_torpid_Tcmin_eq, 2, 2))

## GCB
m_GCB_avgEE_normo_Tcmin_eq <- ggplot(GCB_torpor, aes(as.numeric(Tc_min_C), 
                                                     AvgEE_normo_MassCorrected)) +
  theme_bw() + geom_point(size=4) + geom_smooth(method=lm) +
  geom_text(x=21, y=0.32, label=lm_eqn(GCB_torpor, GCB_torpor$AvgEE_normo_MassCorrected,
                                       GCB_torpor$Tc_min_C), parse=T,size=8) +
  ylab("Avg GCB EE normo (kJ/g)") + xlab(Tc_min.xlab) +
  theme(axis.title.x = element_text(size=24, face="bold"),
        axis.text.x = element_text(size=22),
        axis.title.y = element_text(size=24, face="bold"), axis.text.y = element_text(size=22)) 
m_GCB_avgEE_normo_Tcmin_eq

m_GCB_avgEE_torpid_Tcmin_eq <- ggplot(GCB_torpor, aes(as.numeric(Tc_min_C), 
                                                      AvgEE_torpid_MassCorrected)) +
  theme_bw() + geom_point(size=4) + geom_smooth(method=lm) +
  geom_text(x=21, y=0.10, label=lm_eqn(GCB_torpor, GCB_torpor$AvgEE_torpid_MassCorrected,
                                       GCB_torpor$Tc_min_C), parse=T,size=8) +
  ylab("Avg GCB EE torpid (kJ/g)") + xlab(Tc_min.xlab) +
  theme(axis.title.x = element_text(size=24, face="bold"),
        axis.text.x = element_text(size=22),
        axis.title.y = element_text(size=24, face="bold"), axis.text.y = element_text(size=22)) 
m_GCB_avgEE_torpid_Tcmin_eq

grid.arrange(m_GCB_avgEE_normo_Tcmin_eq, m_GCB_avgEE_torpid_Tcmin_eq, nrow=1, ncol=2)

## Melt BBLH dataframe to put torpid and normo in same column
m_GCB_tor_nor <- melt(GCB_torpor, id.vars="Tc_min_C", 
                      measure.vars = c("AvgEE_torpid_MassCorrected", "AvgEE_normo_MassCorrected"))
levels(m_GCB_tor_nor$variable)[levels(m_GCB_tor_nor$variable)=="AvgEE_normo_MassCorrected"] <- 
  "Avg Normothermic EE"
levels(m_GCB_tor_nor$variable)[levels(m_GCB_tor_nor$variable)=="AvgEE_torpid_MassCorrected"] <- 
  "Avg Torpid EE"

m_GCB_tor_nor$variable <- factor(m_GCB_tor_nor$variable,levels = 
                                   rev(levels(m_GCB_tor_nor$variable)))

## Both normo and torpid avg EE for GCB on same graph
GCB_tor_nor_col <- ggplot(m_GCB_tor_nor, aes(as.numeric(Tc_min_C), value, color=variable)) +
  my_theme + geom_point(aes(shape=variable), size=4) + 
  geom_smooth(method=lm, size=2) + scale_color_manual(values=c("#0099ff", "#ff0000")) +
  geom_text(x = 22, y = 0.32, label = paste("R^2 :", " 0.0302",sep=""), 
            parse=T, size=8, col="black") + 
  geom_text(x = 22, y = 0.12, label = paste("R^2 :", " 0.343",sep=""), 
            parse=T, size=8, col="black") +
  ylab("Avg GCB Energy Expenditure (kJ/g)") + xlab(Tc.xlab) +
  theme(legend.key.height=unit(3,"line"), legend.title=element_text(size=20))
GCB_tor_nor_col

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

library(gtable)
g1 <- ggplotGrob(BBLH_tor_nor)
g1 <- gtable_add_cols(g1, unit(0,"mm")) # add a column for missing legend
g2 <- ggplotGrob(GCB_tor_nor)
g <- cbind(g1, g2, size="first") # stack the two plots
g$widths <- unit.pmax(g1$widths, g2$widths) # use the largest widths
# center the legend vertically
g$layout[grepl("guide", g$layout$name),c("t","b")] <- c(1,nrow(g))
grid.newpage()
grid.draw(g)

#### Statistics ####

## Mean and t tests
mean(torpor$Hours_torpid[torpor$Hours_torpid != 0 & torpor$Site %in% c("HC", "SC", "SWRS")], na.rm=T)
mean(torpor$Hours_torpid[torpor$Hours_torpid != 0 & torpor$Site %in% c("MQ", "SL")], na.rm=T)

hc <- torpor$Hours_torpid[torpor$Species=="BBLH"&torpor$Site=="HC"]
sc <- torpor$Hours_torpid[torpor$Species=="BBLH"&torpor$Site=="SC"]

az_h <- torpor$Hours_torpid[torpor$Hours_torpid != 0 & torpor$Site %in% c("HC", "SC", "SWRS")]
ec_h <- torpor$Hours_torpid[torpor$Hours_torpid != 0 & torpor$Site %in% c("MQ", "SL")]
t.azec.hours <- t.test(az,ec)

az_nee <- torpor$NEE_kJ[torpor$Site %in% c("HC", "SC", "SWRS")]
ec_nee <- torpor$NEE_kJ[torpor$Site %in% c("MQ", "SL")]
t.azec.nee <- t.test(az_nee,ec_)


## Linear and non-linear models
## Regressions
multiple_regression_NEE <- lm(NEE_kJ ~ Tc_mean_C + Tc_min_C + Mass, data=torpor)
summary(multiple_regression_NEE) # show results
# Other useful functions 
coefficients(multiple_regression_NEE) # model coefficients
confint(multiple_regression_NEE, level=0.95) # CIs for model parameters 
fitted(multiple_regression_NEE) # predicted values
residuals(multiple_regression_NEE) # residuals
anova(multiple_regression_NEE) # anova table 
vcov(multiple_regression_NEE) # covariance matrix for model parameters 
influence(multiple_regression_NEE) # regression diagnostics

## Avg normo EE anova (non-mass-corrected)
mul_regr_AvgEEnormo <- lm(Avg_EE_hourly_normo ~ Tc_mean_C + Tc_min_C + Mass, data=torpor)
anova(mul_regr_AvgEEnormo)

##### WEIRD- mean Tc and Mass have a significant effect on Average mass-corrected normo hourly EE
mul_regr_m_AvgEEnormo <- lm(AvgEE_normo_MassCorrected ~ Tc_mean_C + Tc_min_C + Mass, data=torpor)
anova(mul_regr_m_AvgEEnormo)

## Avg non-mass-corrected normo EE for BBLH anova
mul_regr_m_AvgEEnormo_BBLH <- lm(AvgEE_normo_MassCorrected ~ Tc_mean_C + Tc_min_C + 
                                   Mass, data=BBLH_torpor)
anova(mul_regr_m_AvgEEnormo_BBLH)

## Avg non-mass-corrected torpid EE for BBLH anova
mul_regr_m_AvgEEtorpid_BBLH <- lm(AvgEE_torpid_MassCorrected ~ Tc_mean_C + Tc_min_C + 
                                    Mass, data=BBLH_torpor)
anova(mul_regr_m_AvgEEtorpid_BBLH)

## GAMs
torpor$Tc_min_C_sq <- (Tc_min_C)^2
##lm
quad_avgEE_torpid <- lm(torpor$AvgEE_torpid_MassCorrected ~ Tc_min_C + I(Tc_min_C_sq))
predictedEE <- predict(quad_avgEE_torpid,list(Temp=Tc_min_C, Temp2=Tc_min_C^2))
summary(quad_avgEE_torpid)

## Using the poly function
fit1 <- lm(torpor$AvgEE_torpid_MassCorrected ~ torpor$Tc_min_C)
summary(fit1)
fit2b <- lm(torpor$AvgEE_torpid_MassCorrected ~ poly(torpor$Tc_min_C, 2, raw=TRUE))
summary(fit2b)
fit3b <- lm(torpor$AvgEE_torpid_MassCorrected ~ poly(torpor$Tc_min_C, 3, raw=TRUE))
summary(fit3b)
## Plot data points
plot(torpor$Tc_min_C, torpor$AvgEE_torpid_MassCorrected, pch=16, xlab = "Min Temp (deg C)", 
     ylab = "Avg torpid EE Mass-corrected", cex.lab = 1.6, cex.axis=1.6, col = "blue")
lines(sort(torpor$Tc_min_C), predictedEE[order(torpor$Tc_min_C)], col='red', type='b', lwd=3) 

## Interesting useless graph that made me think if there were two separate things happening
plot(torpor$Tc_min_C, torpor$AvgEE_torpid_MassCorrected, type="l", lwd=3)

## Separating torpor data by AvgEE_torpid_MassCorrected curve break.
###### NOTE: NOT ACTUALLY LCT AND UCT. Using those terms just for short convenience. 
##Just temps above and
## below 18 dec C
torpor_LCT <- torpor[torpor$Tc_min_C<=18,]
torpor_UCT <- torpor[torpor$Tc_min_C>18,]

## Quadratic Regression equation for temperatures <= 18 deg C. Good fit! R squared is 0.57. 
## See later, without squared ## term is better fit
quad_avgEE_torpidLCT <- lm(AvgEE_torpid_MassCorrected ~ Tc_min_C + I(Tc_min_C_sq), torpor_LCT)
summary(quad_avgEE_torpidLCT)

## Quadratic Regression equation for temperatures above 18 deg C. Not much of a fit, 
## R squared is 0.20
quad_avgEE_torpidUCT <- lm(AvgEE_torpid_MassCorrected ~ Tc_min_C + I(Tc_min_C_sq), torpor_UCT)
summary(quad_avgEE_torpidUCT)

## lm without sqared term of the same as above. Better fit than quadratic for temperatures <= 18
lm_avgEE_torpidLCT <- lm(AvgEE_torpid_MassCorrected ~ Tc_min_C, torpor_LCT)
summary(lm_avgEE_torpidLCT)
## lm without squared term, for temos >18; even worse fit than with squared.
lm_avgEE_torpidUCT <- lm(AvgEE_torpid_MassCorrected ~ Tc_min_C, torpor_UCT)
summary(lm_avgEE_torpidUCT)

## Mean and min Tc fit equally well; both have R2 of about 0.61
lm_avgEE_torpidLCT_Tc_mean <- lm(AvgEE_torpid_MassCorrected ~ Tc_mean_C, torpor_LCT)
summary(lm_avgEE_torpidLCT)

AvgEE_torpidLCT_Tc_min <- ggplot(torpor_LCT, aes(Tc_min_C, AvgEE_torpid_MassCorrected)) + 
  theme_bw() + geom_point(aes(shape = factor(Species)), size=4) + labs(shape ='Species') +
  geom_smooth(method=lm, color="black") +
  geom_text(x = 11, y = 0.092, label = lm_eqn(torpor_LCT, torpor_LCT$AvgEE_torpid_MassCorrected, 
                                              torpor_LCT$Tc_min_C), parse=T, size=8) +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + 
  #facet_grid(.~Site,space="free") +
  ylab("Avg EE torpid (kJ/g)") + xlab(Tc_min.xlab) +
  theme(axis.title.x = element_text(size=24, face="bold"),
        axis.text.x = element_text(size=22),
        axis.title.y = element_text(size=24, face="bold"), axis.text.y = element_text(size=22)) 
AvgEE_torpidLCT_Tc_min

AvgEE_torpidUCT_Tc_min <- ggplot(torpor_UCT, aes(Tc_min_C, AvgEE_torpid_MassCorrected)) + 
  theme_bw() + geom_point(aes(shape = factor(Species)), size=4) + labs(shape ='Species') +
  geom_smooth(method=lm, color="black") +
  geom_text(x = 21, y = 0.094, label = lm_eqn(torpor_UCT, torpor_UCT$AvgEE_torpid_MassCorrected, 
                                              torpor_UCT$Tc_min_C), parse=T, size=8) +
  scale_shape_manual(values=c(3,1,2,0,15,16,17,23)) +
  scale_color_brewer(palette = "Set1") + 
  #facet_grid(.~Site,space="free") +
  ylab("Avg EE torpid (kJ/g)") + xlab(Tc_min.xlab) +
  theme(axis.title.x = element_text(size=24, face="bold"),
        axis.text.x = element_text(size=22),
        axis.title.y = element_text(size=24, face="bold"), axis.text.y = element_text(size=22)) 
AvgEE_torpidUCT_Tc_min

grid.arrange(AvgEE_torpidLCT_Tc_min, AvgEE_torpidUCT_Tc_min, nrow=1, ncol=2)

## lm() With min Tc
quad_avgEE_normo <- lm(torpor$AvgEE_normo_MassCorrected ~ torpor$Tc_min_C + I(torpor$Tc_min_C_sq))
predictedEE_normo <- predict(quad_avgEE_normo,list(Temp=torpor$Tc_min_C, Temp2=torpor$Tc_min_C_sq))
## Plot avg normo EE vs. min temperatures
plot(torpor$Tc_min_C, torpor$AvgEE_normo_MassCorrected, pch=16, xlab = "Min Temp (deg C)", 
     ylab = "Avg normo EE Mass-corrected", cex.lab = 1.3, col = "blue")
lines(sort(torpor$Tc_min_C), predictedEE_normo[order(torpor$Tc_min_C)], col='red', type='b') 

## lm() With mean Tc
torpor$Tc_mean_C_sq <- (torpor$Tc_mean_C)^2
lm_normo_Tc_mean <- lm(torpor$AvgEE_normo_MassCorrected ~ torpor$Tc_mean_C)
plot(lm_normo_Tc_mean)
lines(torpor$Tc_mean_C, predict(lm_normo_Tc_mean), col='red', type='b')
quad_avgEE_normo_Tcmean <- lm(torpor$AvgEE_normo_MassCorrected ~ torpor$Tc_mean_C) + 
  I(torpor$Tc_mean_C_sq)
predictedEE_normo_Tcmean <- predict(quad_avgEE_normo,list(Temp=torpor$Tc_mean_C, 
                                                          Temp2=torpor$Tc_mean_C_sq))
## Plot avg normo EE vs. mean temperatures
plot(torpor$Tc_mean_C, torpor$AvgEE_normo_MassCorrected, pch=16, xlab = "Mean Temp (deg C)", 
     ylab = "Avg normo EE Mass-corrected", cex.lab = 1.3, col = "blue")
lines(sort(torpor$Tc_mean_C), predictedEE_normo[order(torpor$Tc_mean_C)], col='red', type='b') 

#### PCA - Looks Awesome in 3D ######
t.pc <- torpor[,c("NEE_kJ", "Hours_torpid2", "Tc_min_C", "Percentage_avg")]
names(t.pc) <- c("NEE", "Torpor_duration", "Min Tc", "Hourly_energy_savings")
t.pc$Hourly_energy_savings <- as.numeric(t.pc$Hourly_energy_savings)
t.pc$Hourly_energy_savings[is.na(t.pc$Hourly_energy_savings)] <- 0
t.pc$Torpor_duration[is.na(t.pc$Torpor_duration)] <- 0
t.pc$Avg_EE_hourly_torpid[is.na(t.pc$Avg_EE_hourly_torpid)] <- 0
pc.cr <- prcomp(t.pc[,1:4], center=T, scale. = T)
pc.cr2 <- princomp(t.pc[,1:4], cor=T, scores=T)

plot3d(pc.cr2$scores[,1:3], col=as.integer(torpor$Site),size = 10)

text3d(pc.cr2$scores[,1:3],texts=torpor$Species,adj = 1)
text3d(pc.cr2$loadings[,1:3], texts=rownames(pc.cr2$loadings), col="red")
coords <- NULL
for (i in 1:nrow(pc.cr2$loadings)) {
  coords <- rbind(coords, rbind(c(0,0,0),pc.cr2$loadings[i,1:3]))
}
lines3d(coords, col="red", lwd=4)

## 2D
g.nee <- ggbiplot(pc.cr, obs.scale = 1, var.scale = 1, 
                  groups = as.character(torpor$Site), ellipse = TRUE, 
                  circle = TRUE) + scale_color_brewer(palette = "Set1")
g.nee <- g.nee + theme(legend.direction = 'horizontal', 
                       legend.position = 'top') + theme_classic(base_size = 20) +
  theme(axis.title.y = element_text(vjust = 2), 
        panel.border = element_rect(colour = "black", fill=NA))
print(g.nee)

summary(pc.cr)

###Anovas #####
## Anova of torpid energy savings as a function of site and species. Species doesn't matter, site does!
anova(lm(Percentage_avg~Site_new+Species+Tc_min_C+Mass, data = torpor))

## Anova with Avg hourly EE torpid gives slightly different results
anova(lm(AvgEE_torpid_MassCorrected~Site_new+Species+Tc_min_C+Mass, data = torpor))

## Anova with NEE
anova(lm(NEE_MassCorrected~Site_new+Species+Tc_min_C+Mass, data = torpor))


#### T-tests ######
## Subsetting melted dataframe to get just depth values. Then subtracting from 100 to make them hourly savings.
m.savings <- m.temptrop[m.temptrop$variable=="Percentage_avg",]
m.savings$value <- 100-m.savings$value
t.test(m.savings$value[m.savings$Temptrop=="Temperate"], m.savings$value[m.savings$Temptrop=="Tropical"], 
       paired = F)

## Comparing savings, MQ vs. SL, and HC vs. SC
t.test(m.savings$value[m.savings$Site=="MQ"], m.savings$value[m.savings$Site=="SL"], paired=F)

## Subsetting NEE non-masscorrected, and testing tropical vs. temperate
m.nee <- m.temptrop[m.temptrop$variable=="NEE_kJ",]
t.test(m.nee$value[m.nee$Temptrop=="Temperate"], m.nee$value[m.nee$Temptrop=="Tropical"], 
       paired = F)

## NEE mass-corrected
m.nee <- m.temptrop[m.temptrop$variable=="NEE_MassCorrected",]
t.test(m.nee$value[m.nee$Temptrop=="Temperate"], m.nee$value[m.nee$Temptrop=="Tropical"], 
       paired = F)
t.test(m.nee$value[m.nee$Site=="SC"], m.nee$value[m.nee$Site=="HC"], 
       paired = F)
m.tor <- m.temptrop[m.temptrop$variable=="NEE_MassCorrected" & !is.na(m.temptrop$value[m.temptrop$variable=="Hours_torpid"]),]
m.nor <- m.temptrop[m.temptrop$variable=="NEE_MassCorrected" & is.na(m.temptrop$value[m.temptrop$variable=="Hours_torpid"]),]
t.test(m.tor$value[m.nee$Site=="SC"], m.tor$value[m.nee$Site=="HC"], paired = F)

# Difference between normothermic and torpid NEE in Sonoita birds
(mean(m.tor$value[m.tor$Site=="SC"])/mean(m.nor$value[m.nor$Site=="SC"]))*100
# Difference in NEE between SC and HC torpid birds
((mean(m.tor$value[m.nee$Site=="SC"]) - mean(m.tor$value[m.nee$Site=="HC"]))/(mean(m.tor$value[m.nee$Site=="SC"])))*100
## Diference in NEE between  MAHU and BLUH
((mean(5.459978,5.349119) - mean(5.195791, 4.780688, 6.583659, 4.320507))/mean(5.459978,5.349119))*100

## Rate of occurrence
t.test(freq_table$Rate_occurrence[freq_table$Temptrop=="Temperate"], 
       freq_table$Rate_occurrence[freq_table$Temptrop=="Tropical"], paired=F)
t.test(rate_site$N[rate_site$Site=="HC"], rate_site$N[rate_site$Site=="SC"], paired=F)

## Hours torpid
m.hours <- m.temptrop[m.temptrop$variable=="Hours_torpid",]
t.test(m.hours$value[m.hours$Temptrop=="Temperate"],  m.hours$value[m.hours$Temptrop=="Tropical"], paired=F)

## Hours torpid
m.prophours <- m.temptrop[m.temptrop$variable=="Prop_hours",]
t.test(m.prophours$value[m.prophours$Temptrop=="Temperate"],
       m.prophours$value[m.prophours$Temptrop=="Tropical"], paired=F)

## Testing just hourly torpid EE in HC vs. SC birds
t.test(torpor$Avg_EE_hourly_torpid[torpor$Site=="HC"], 
       torpor$Avg_EE_hourly_torpid[torpor$Site=="SC"], paired=F)
