################################################################################
##  NutNet_preliminaryAnalyses.R: Analysis of plant community composition responses in NutNet experiment at Konza Prairie LTER.
##
##  Author: Kimberly Komatsu
##  Date created: May 2, 2023
################################################################################

library(codyn)
library(nlme)
library(emmeans)
library(vegan)
library(readxl)
library(tidyverse)

setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\nutrient network\\nutnet\\La Pierre_Konza_Saline_SGS_data') #desktop


##### functions #####
`%notin%` <- Negate(`%in%`)

#set options
options(contrasts=c('contr.sum','contr.poly'))

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

#homemade functions
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  


##### data #####
trt <- read.csv('KNZ_NutNet_trt.csv')

sp2007 <- read_excel('NutNet_2007_data\\Konza_NutNet_2007.xls', sheet='cover') %>% 
  filter(site=='2C') %>% 
  pivot_longer(cols=c('june_cover', 'august_cover'), names_to='date', values_to='cover') %>% 
  mutate(season=ifelse(date=='june_cover', 'spring', 'fall')) %>% 
  select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum))
sp2008 <- read_excel('NutNet_2008_data\\KNZ_SER_NutNet_2008.xls', sheet='cover') %>% 
  filter(site=='KNZ') %>% select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum))
sp2009 <- read_excel('NutNet_2009_data\\KNZ_SER_SGSNutNet_2009.xls', sheet='cover') %>% 
  filter(site=='KNZ') %>% select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum))
sp2010 <- read_excel('NutNet_2010_data\\KNZ_SER_SGSNutNet_2010.xls', sheet='cover') %>% 
  filter(site=='KNZ') %>% select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum))
sp2011 <- read_excel('NutNet_2011_data\\KNZ_SER_SGSNutNet_2011.xls', sheet='cover') %>% 
  filter(site=='KNZ') %>% select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum))
sp2012 <- read_excel('NutNet_2012_data\\KNZ_SER_SGSNutNet_2012.xls', sheet='cover') %>% 
  filter(site=='KNZ') %>% select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum))
sp2013 <- read_excel('NutNet_2013_data\\KNZ_SER_SGSNutNet_2013.xls', sheet='cover') %>% 
  filter(site=='KNZ') %>% select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum))
sp2014 <- read_excel('NutNet_2014_data\\KNZ_SER_NutNet_2014.xls', sheet='cover') %>% 
  filter(site=='KNZ') %>% select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum))
sp2015 <- read_excel('NutNet_2015_data\\KNZ_SER_NutNet_2015.xls', sheet='cover') %>% 
  filter(site=='KNZ') %>% select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum))
sp2016 <- read_excel('NutNet_2016_data\\KNZ_SER_NutNet_2016.xls', sheet='cover') %>% 
  filter(site=='KNZ') %>% select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum))
sp2017 <- read_excel('NutNet_2017_data\\KNZ_SER_NutNet_2017.xls', sheet='cover') %>% 
  filter(site=='KNZ') %>% select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum)) %>% 
  mutate(season=str_to_lower(season))
sp2018 <- read_excel('NutNet_2018_data\\KNZ_SER_NutNet_2018.xls', sheet='cover') %>% 
  filter(site=='KNZ') %>% select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum))
sp2019 <- read_excel('NutNet_2019_data\\KNZ_NutNet_2019.xls', sheet='cover') %>% 
  filter(site=='KNZ') %>% select(year, season, plot, sppnum, cover) %>% mutate(sppnum=as.numeric(sppnum))
sp2021 <- read_excel('NutNet_2021_data\\KNZ_NutNet_2021.xls', sheet='cover') %>% 
  filter(site=='konz.us') %>% select(year, season, plot, spp_num, cover) %>% 
  rename(sppnum=spp_num)
sp2022 <- read_excel('NutNet_2022_data\\KNZ_NutNet_2022.xlsx', sheet='cover') %>% 
  mutate(sppnum=as.numeric(spp_num),
         season=ifelse(season=='early', 'spring', 'fall')) %>% 
  filter(site=='konz.us') %>% select(year, season, plot, sppnum, cover) 

spAll <- rbind(sp2007, sp2008, sp2009, sp2010, sp2011, sp2012, sp2013, sp2014, sp2015, sp2016, sp2017, sp2018, sp2019, sp2021, sp2022) %>%
  group_by(year, plot, sppnum) %>%
  summarise(max_cover=max(cover)) %>%
  ungroup() %>%
  rename(code=sppnum) %>% 
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\prairie species lists\\PPS011_new KNZ spp list.csv')) %>%
  filter(gen %notin% c('litter', 'rock', 'dung', 'bare_ground', 'bison_trail')) %>%
  mutate(genus_species=paste(genus, species, sep='_'))


##### relative cover #####
totCover <- spAll %>%
  group_by(year, plot) %>%
  summarise(total_cover=sum(max_cover)) %>% #calculate total cover
  ungroup()

relCover <- spAll %>%
  left_join(totCover) %>%
  mutate(rel_cover=100*(max_cover/total_cover)) %>% #calculate relative cover
  select(-total_cover) %>%
  left_join(trt) #%>%
  # filter(!is.na(code)) #remove entries that were unknowns


##### community metrics #####
commMetrics <- community_structure(relCover, time.var='year', abundance.var='rel_cover', replicate.var='plot') %>%
  left_join(trt)


##### richness response #####
summary(richModel <- lme(richness~year*n*p*k,
                         data=subset(commMetrics, exclose!='Fenced'),
                         random=~1|block/plot,
                         correlation=corCompSymm(form=~year|block/plot), 
                         control=lmeControl(returnObject=T)))
anova.lme(richModel, type='sequential') #significant N*year interaction
emmeans(richModel, pairwise~year*n, adjust="tukey")


#figure - richness by year*n
ggplot(data=barGraphStats(data=subset(commMetrics, exclose!='Fenced'), variable="richness", byFactorNames=c('year', 'n')), aes(x=year, y=mean, color=n)) +
  geom_rect(aes(xmin=2006.5, xmax=2007.5, ymin=0, ymax=20), fill='#F0F0F0', alpha=0.1, color='#F0F0F0') +
  geom_rect(aes(xmin=2008.5, xmax=2009.5, ymin=0, ymax=20), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2010.5, xmax=2011.5, ymin=0, ymax=20), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2012.5, xmax=2013.5, ymin=0, ymax=20), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2014.5, xmax=2015.5, ymin=0, ymax=20), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2016.5, xmax=2017.5, ymin=0, ymax=20), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2018.5, xmax=2019.5, ymin=0, ymax=20), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2020.5, xmax=2021.5, ymin=0, ymax=20), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_point() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(0.1)) +
  geom_smooth(method='lm', se=F) +
  ylab('Plant Species Richness') + xlab('Year') +
  scale_color_manual(values=c('#030E4F', '#F49F1C')) +
  coord_cartesian(ylim=c(10,18))
#export at 1400x600

# ggsave("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\grants\\NSF_FY2023\\FY23_NSF_PopComm_lights under the prairie\\figures\\NutNet_temporal richness_nitrogen.png", width=8, height=5, units='in')


##### response ratios ##### 
commMetricsRR <- commMetrics %>%
  filter(treat_other_name=='control') %>% 
  select(year, block, richness, Evar) %>% 
  rename(richness_ctl=richness,
         Evar_ctl=Evar) %>% 
  right_join(commMetrics) %>% 
  filter(treat_other_name!='control') %>% 
  mutate(lnRR_richness=log(richness/richness_ctl),
         lnRR_Evar=log(Evar/Evar_ctl))

summary(richRRModel <- lme(lnRR_richness~year,
                         data=subset(commMetricsRR, n='Nitrogen'),
                         random=~1|block,
                         correlation=corCompSymm(form=~year|block), 
                         control=lmeControl(returnObject=T)))
anova.lme(richRRModel, type='sequential') #significant N*year interaction
emmeans(richRRModel, pairwise~year, adjust="tukey")

#figure - richness RR
ggplot(data=barGraphStats(data=subset(commMetricsRR, n='Nitrogen'), variable="lnRR_richness", byFactorNames=c('year')), aes(x=year, y=mean)) +
  geom_rect(aes(xmin=2006.5, xmax=2007.5, ymin=-1, ymax=1), fill='#F0F0F0', alpha=0.9) +
  geom_rect(aes(xmin=2008.5, xmax=2009.5, ymin=-1, ymax=1), fill='#F0F0F0', alpha=0.9) +
  geom_rect(aes(xmin=2010.5, xmax=2011.5, ymin=-1, ymax=1), fill='#F0F0F0', alpha=0.9) +
  geom_rect(aes(xmin=2012.5, xmax=2013.5, ymin=-1, ymax=1), fill='#F0F0F0', alpha=0.9) +
  geom_rect(aes(xmin=2014.5, xmax=2015.5, ymin=-1, ymax=1), fill='#F0F0F0', alpha=0.9) +
  geom_rect(aes(xmin=2016.5, xmax=2017.5, ymin=-1, ymax=1), fill='#F0F0F0', alpha=0.9) +
  geom_rect(aes(xmin=2018.5, xmax=2019.5, ymin=-1, ymax=1), fill='#F0F0F0', alpha=0.9) +
  geom_rect(aes(xmin=2020.5, xmax=2021.5, ymin=-1, ymax=1), fill='#F0F0F0', alpha=0.9) +
  geom_point(size=2, color='#F49F1C') +
  geom_smooth(method='lm', color='#F49F1C', se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, color='#F49F1C'), width=.1) +
  ylab(expression(paste('lnRR Plant Species Richness'))) + xlab('Year') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(-0.25,0.15))

# ggsave("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\grants\\NSF_FY2023\\FY23_NSF_PopComm_lights under the prairie\\figures\\NutNet_temporal lnRR richness_nitrogen.png", width=8, height=5, units='in')




##### community change #####
racChange <- RAC_change(df=relCover, time.var='year', species.var='genus_species', abundance.var='rel_cover', replicate.var='plot', reference.time=2007) %>% 
  select(-year) %>% 
  rename(year=year2) %>% 
  left_join(trt)

#gains
summary(gainsModel <- lme(gains~year*n*p*k,
                           data=subset(racChange, exclose!='Fenced'),
                           random=~1|block,
                           correlation=corCompSymm(form=~year|block), 
                           control=lmeControl(returnObject=T)))
anova.lme(gainsModel, type='sequential') #significant N*P*K interaction, no year effect
emmeans(gainsModel, pairwise~year, adjust="tukey")


ggplot(data=barGraphStats(data=subset(racChange, exclose!='Fenced'), variable="gains", byFactorNames=c('year','treat_other_name')), aes(x=year, y=mean, color=treat_other_name)) +
  # geom_rect(aes(xmin=2006.5, xmax=2007.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2008.5, xmax=2009.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2010.5, xmax=2011.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2012.5, xmax=2013.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2014.5, xmax=2015.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2016.5, xmax=2017.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2018.5, xmax=2019.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2020.5, xmax=2021.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_point(size=2) +
  geom_smooth(method='lm', se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  ylab(expression(paste('Plant Species Gains'))) + xlab('Year') +
  # theme(legend.position='none') +
  coord_cartesian(ylim=c(0,0.5))

# ggsave("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\grants\\NSF_FY2023\\FY23_NSF_PopComm_lights under the prairie\\figures\\NutNet_temporal gains_all trt.png", width=8, height=5, units='in')


#losses
summary(lossesModel <- lme(losses~year*n*p*k,
                          data=subset(racChange, exclose!='Fenced'),
                          random=~1|block,
                          correlation=corCompSymm(form=~year|block), 
                          control=lmeControl(returnObject=T)))
anova.lme(lossesModel, type='sequential') #significant N*P*K interaction, year*N, year*NP
emmeans(lossesModel, pairwise~year, adjust="tukey")


ggplot(data=barGraphStats(data=subset(racChange, exclose!='Fenced'), variable="losses", byFactorNames=c('year','treat_other_name')), aes(x=year, y=mean, color=treat_other_name)) +
  # geom_rect(aes(xmin=2006.5, xmax=2007.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2008.5, xmax=2009.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2010.5, xmax=2011.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2012.5, xmax=2013.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2014.5, xmax=2015.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2016.5, xmax=2017.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2018.5, xmax=2019.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2020.5, xmax=2021.5, ymin=0, ymax=1), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_point(size=2) +
  geom_smooth(method='lm', se=F) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  ylab(expression(paste('Plant Species Losses'))) + xlab('Year') +
  # theme(legend.position='none') +
  coord_cartesian(ylim=c(0,1))


ggplot(data=subset(racChange, exclose!='Fenced'), aes(x=as.factor(year), y=losses, color=treat_other_name)) +
  geom_boxplot() +
  scale_color_discrete(breaks=c('control', 'N', 'P', 'K', 'NP', 'NK', 'PK', 'NPK')) +
  ylab(expression(paste('Plant Species Losses (%)'))) + xlab('Year') +
  # theme(legend.position='none') +
  coord_cartesian(ylim=c(0,0.75))

# ggsave("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\grants\\NSF_FY2023\\FY23_NSF_PopComm_lights under the prairie\\figures\\NutNet_temporal losses_all trt.png", width=8, height=5, units='in')




##### PERMANOVA #####
relCover2022 <- relCover %>%
  select(year, plot, n, p, k, exclose, treat_other_name, genus_species, rel_cover) %>% 
  pivot_wider(names_from='genus_species', values_from='rel_cover', values_fill=list(rel_cover=0)) %>%
  filter(year==2022)

print(permanova <- adonis2(formula = relCover2022[,8:106]~n*p*k, data=relCover2022, permutations=999, method="bray"))
#treat_other_name F=1.9534, df=9,29, p=0.009
# n         1   1.3474 0.24905 9.8645  0.001 ***
# p         1   0.3093 0.05718 2.2648  0.062 .  
# k         1   0.0891 0.01646 0.6521  0.674    
# n:p       1   0.3706 0.06850 2.7130  0.036 *  
# n:k       1   0.0536 0.00991 0.3923  0.914    
# p:k       1   0.0923 0.01706 0.6759  0.657    
# n:p:k     1   0.1428 0.02639 1.0453  0.360   

#betadisper
veg <- vegdist(relCover2022[,8:106], method = "bray")
dispersion <- betadisper(veg, relCover2022$treat_other_name)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=1.1807, df=9,20, p=0.317

# sppBC <- metaMDS(relCover2022[,8:106])
# 
# plotData <- relCover2022[,1:9]
# 
# #Use the vegan ellipse function to make ellipses
# veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
# {
#   theta <- (0:npoints) * 2 * pi/npoints
#   Circle <- cbind(cos(theta), sin(theta))
#   t(center + scale * t(Circle %*% chol(cov)))
# }
# 
# BC_NMDS = data.frame(MDS1 = sppBC$points[,1], MDS2 = sppBC$points[,2],group=relCover2022$treat_other_name)
# BC_NMDS_Graph <- cbind(plotData,BC_NMDS)
# BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$treat_other_name, display = "sites",
#                              kind = "se", conf = 0.95, label = T)               
# 
# ord3 <- data.frame(plotData,scores(sppBC,display="sites")) %>%
#   group_by(treat_other_name)
# 
# BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$treat_other_name, display = "sites",
#                              kind = "se", conf = 0.95, label = T)
# BC_Ellipses <- data.frame() #Make a new empty data frame called BC_Ellipses  
# for(g in unique(BC_NMDS$group)){
#   BC_Ellipses <- rbind(BC_Ellipses, cbind(as.data.frame(with(BC_NMDS[BC_NMDS$group==g,], 
#                                                              veganCovEllipse(BC_Ord_Ellipses[[g]]$cov,BC_Ord_Ellipses[[g]]$center,BC_Ord_Ellipses[[g]]$scale)))
#                                           ,group=g))
# } #Generate ellipses points
# 
# ggplot(BC_NMDS_Graph, aes(x=MDS1, y=MDS2, color=group,linetype = group, shape = group)) +
#   geom_point(size=6)+ 
#   geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
#   labs(color="", linetype = "", shape = "") +
#   scale_colour_manual(values=c("brown", "brown", "dark green", "dark green", "dark green", "dark green"), name = "") +
#   scale_linetype_manual(values = c("twodash", "solid", "twodash", "solid", "twodash", "solid"), name = "") +
#   xlab("NMDS1")+ 
#   ylab("NMDS2")+ 
#   theme(axis.text.x=element_text(size=24, color = "black"), axis.text.y = element_text(size = 24, color = "black"), legend.text = element_text(size = 24))



##### simper #####
print(sim <- with(relCover2022, simper(relCover2022[,8:106], treat_other_name)))


#RACs
rankAbundance <- relCover%>%
  filter(year==2022) %>%
  filter(!is.na(genus)) %>% 
  mutate(lifeform=ifelse(lifeform=='s', 'g', lifeform)) %>% 
  mutate(spp_name=str_to_sentence(genus_species)) %>%
  separate(spp_name, into=c('genus', 'species'), sep='_') %>% 
  unite(genus_species, genus, species, sep=' ') %>% 
  group_by(treat_other_name, genus_species, growthform, lifeform) %>%
  summarize(avg_cover=mean(rel_cover)) %>%
  ungroup() %>%
  arrange(treat_other_name, -avg_cover) %>%
  group_by(treat_other_name) %>%
  mutate(rank=seq_along(treat_other_name)) %>%
  ungroup()

ggplot(data=subset(rankAbundance, avg_cover>0 & !(treat_other_name %in% c('fence', 'NPK_fence'))), aes(x=rank, y=avg_cover)) +
  geom_line() +
  geom_point(aes(colour=lifeform, shape=growthform), size=3) +
  scale_color_manual(labels=c("Forb", "Graminoid", "Woody"), values=c('#DE8C00', '#009CC0', '#83431E')) +
  scale_shape_discrete(labels=c("Annual", "Perennial")) +
  xlab('') +
  ylab('Relative Percent Cover') +
  # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
  geom_text(aes(y=avg_cover+1.2, x=rank+0.1, label=genus_species), hjust='left', vjust='center', angle=45, size=4) +
  facet_wrap(~treat_other_name, scales='free_x') +
  coord_cartesian(ylim=c(0,100))
#export at 1400x400

# ggsave("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\grants\\NSF_FY2023\\FY23_NSF_PopComm_lights under the prairie\\figures\\NutNet_RAC2022_all trt.png", width=12, height=8, units='in')


#### temporal trend by number of resources manipulated ####
commMetricsFactor <- commMetrics %>% 
  filter(!(treat_other_name %in% c('fence', 'NPK_fence')),
         !is.na(site)) %>% 
  mutate(plot_mani=ifelse(treat_other_name %in% c('N', 'P', 'K'), 1,
                   ifelse(treat_other_name %in% c('NP', 'NK', 'PK'), 2,
                   ifelse(treat_other_name %in% c('NPK'), 3, 0)))) %>% 
  mutate(plot_mani_N=ifelse(treat_other_name %in% c('P', 'K'), '1',
                   ifelse(treat_other_name %in% c('PK'), '2',
                   ifelse(treat_other_name %in% c('NPK'), '3N', 
                   ifelse(treat_other_name %in% c('N'), '1N', 
                   ifelse(treat_other_name %in% c('NP', 'NK'), '2N', 0))))))

#figure - richness over time by plot mani
ggplot(data=commMetricsFactor, aes(x=year, y=richness, color=as.factor(plot_mani))) +
  geom_rect(aes(xmin=2006.5, xmax=2007.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.1, color='#F0F0F0') +
  geom_rect(aes(xmin=2008.5, xmax=2009.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2010.5, xmax=2011.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2012.5, xmax=2013.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2014.5, xmax=2015.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2016.5, xmax=2017.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2018.5, xmax=2019.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2020.5, xmax=2021.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_jitter(width=0.3, height=0) +
  # geom_point() +
  geom_smooth(se=F) +
  ylab('Plant Species Richness') + xlab('Year') +
  scale_color_manual(values=c('#030E4F', '#FFE230', '#F6830C', '#C11414')) +
  coord_cartesian(ylim=c(0,21))

# ggsave("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\grants\\NSF_FY2023\\FY23_NSF_PopComm_lights under the prairie\\figures\\NutNet_richness_plot mani.png", width=9, height=8, units='in')


#figure - richness over time by plot mani with N indicated
ggplot(data=commMetricsFactor, aes(x=year, y=richness, color=as.factor(plot_mani), line_type=as.factor(n))) +
  geom_rect(aes(xmin=2006.5, xmax=2007.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.1, color='#F0F0F0') +
  geom_rect(aes(xmin=2008.5, xmax=2009.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2010.5, xmax=2011.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2012.5, xmax=2013.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2014.5, xmax=2015.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2016.5, xmax=2017.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2018.5, xmax=2019.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_rect(aes(xmin=2020.5, xmax=2021.5, ymin=0, ymax=25), fill='#F0F0F0', alpha=0.9, color='#F0F0F0') +
  geom_jitter(width=0.3, height=0) +
  # geom_point() +
  geom_smooth(se=F, aes(linetype=as.factor(n), color=as.factor(plot_mani))) +
  ylab('Plant Species Richness') + xlab('Year') +
  scale_color_manual(values=c('#030E4F', '#FFE230', '#F6830C', '#C11414')) +
  coord_cartesian(ylim=c(0,21))

# ggsave("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\grants\\NSF_FY2023\\FY23_NSF_PopComm_lights under the prairie\\figures\\NutNet_richness_plot mani.png", width=9, height=8, units='in')



##### light - 2010 #####
light2010 <- read_excel('NutNet_2010_data\\KNZ_SER_SGSNutNet_2010.xls', sheet='par') %>% 
  filter(site=='KNZ', month==8) %>% 
  left_join(trt) %>% 
  filter(exclose!='Fenced') %>% 
  group_by(site, plot, n, p, k, treat_other_name) %>% 
  summarise(across(c(parg:par1.1), mean, na.rm=T)) %>% 
  ungroup() %>% 
  pivot_longer(cols=c(parg:par1.1), names_to='level', values_to='par') %>% 
  separate(col=level, into=c('drop', 'height_m'), sep='r') %>% 
  select(-drop) %>% 
  mutate(height_m=as.numeric(ifelse(height_m=='g', 0, height_m)))

ggplot(data=barGraphStats(data=light2010, variable="par", byFactorNames=c("height_m", "treat_other_name")), aes(x=height_m, y=mean, color=treat_other_name)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.02)) +
  geom_path()

ggplot(data=barGraphStats(data=subset(light2010, treat_other_name %in% c('control', 'N', 'NP', 'NPK')), variable="par", byFactorNames=c("height_m", "treat_other_name")), aes(x=height_m, y=mean, color=treat_other_name)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.02)) +
  geom_path() +
  xlab('Height (m)') + ylab('PAR')
# ggsave("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\grants\\NSF_FY2023\\FY23_NSF_PopComm_lights under the prairie\\figures\\NutNet_light attenuation.png", width=9, height=8, units='in')
