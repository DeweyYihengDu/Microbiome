# ##############################################################################
#
## Checkerboard plots for follow-ups
#
# ##############################################################################

library("here")
library("tidyverse")
library("ggpubr")


AUCs <- read_tsv(here('poster/data',"clean_AUCs_followup.tsv" ))

AUCs_short <- AUCs %>% 
  select('DrugX', 'DrugXConc', 'ABX', 'ABXConc', 
         'Strain', 'Replicate', 'ABXnormAUC') %>% 
  na.omit()

#calculate median ABXnormAUC
median_AUCs_short <- AUCs_short %>% 
  group_by(DrugX, DrugXConc, ABX, ABXConc, Strain) %>% 
  summarise(median=median(ABXnormAUC))

median_AUCs_short <- median_AUCs_short %>% 
  mutate(label=paste(ABX, Strain, DrugX, sep='-'))

#median_AUCs_short <- median_AUCs_short %>% filter(DrugXConc != 0)

median_AUCs_short <- median_AUCs_short %>% 
  mutate(star = case_when(
    median > 0.25 ~ '*',
    median < 0.25 ~ ''))

median_AUCs_short <- median_AUCs_short %>% 
  mutate(median_new = case_when(
    median > 0  ~ median,
    median <= 0  ~ 0))


myplot_ery <- ggplot(filter(median_AUCs_short, ABX=='Erythromycin'), 
                     aes(x= as.factor(DrugXConc), y=label)) + 
  geom_tile(aes(fill=median_new), colour="white",size=0.25) +
  geom_text(aes(label=star), size = 10) +
  theme_bw(base_size = 14) +
  labs(x="Conc. Drug X in µM",y="")+
  scale_fill_gradientn(name='AUC', colours=c('#fff5ec',"#ff8000"),
                       limits = c(0,0.6)) +
  #remove extra space
  scale_y_discrete(expand=c(0,0)) +
  coord_fixed()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


myplot_doxy <- ggplot(filter(median_AUCs_short, ABX=='Doxycycline'), 
                      aes(x= as.factor(DrugXConc), y=label)) + 
  geom_tile(aes(fill=median_new), colour="white",size=0.25) +
  geom_text(aes(label=star), size = 10) +
  theme_bw(base_size = 14) +
  labs(x="Conc. Drug X in µM",y="")+
  scale_fill_gradientn(name='AUC', colours=c("#fff0f0","#f00000"),
                       limits = c(0,0.6)) +
  #remove extra space
  scale_y_discrete(expand=c(0,0)) +
  coord_fixed()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

pdf(here('poster/figures',"EDFig7c.pdf"),width=18,height=12) 
ggarrange(myplot_ery, myplot_doxy, ncol=1, nrow=2)
dev.off()

# export source data
if (dir.exists(here('poster/source_data'))){
  median_AUCs_short %>% 
    rename(median_adj=median_new) %>% 
    write_tsv(here('poster/source_data', 'EDFig7c.tsv'))
}