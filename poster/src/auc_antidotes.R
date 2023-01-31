# ##############################################################################
#
## strain AUC for different concentration of dicumarol/tolfenamic acid
#
# ##############################################################################

library("here")
library("tidyverse")
library("ggpubr")
library("scales")

pathogens <- read_tsv(here('poster/data','CB_pathogens_raw_data.tsv'))

pathogens <- pathogens %>% 
  filter(DrugX %in% c('Tolfenamic Acid','Dicumarol')) %>% 
  filter(ABXConc == 0.65 | ABXConc == 0.14)


pathogens <- pathogens %>% mutate(Strain_name = case_when(
  Strain == 'Enterococcus faecium' ~ 'E. faecium', 
  Strain == 'Enterococcus faecalis' ~ 'E. faecalis', 
  Strain == 'SADSM20231' ~ 'S. aureus'))

pathogens$DrugX[pathogens$DrugX == 'Tolfenamic Acid'] <- 'Tolfenamic acid'

commensals <- read_tsv(here('poster/data','CB_commensals_raw_data.tsv'))

commensals <- commensals %>% 
  filter(DrugX %in% c('Tolfenamic acid', 'Dicumarol')) %>% 
  filter(ABXConc==0.65)


pathogens_commensals <- rbind(commensals,pathogens)


median_commensals <- pathogens_commensals %>%  
  group_by(DrugX, DrugXConc, ABXConc, Strain, ABX, Strain_name) %>% 
  summarise(median_normAUC=median(normAUC), 
            mean_normAUC=mean(normAUC),
            SDAUC=sd(normAUC)) %>%
  ungroup()


pdf(here('poster/figures',"EDFig9c.pdf"),width=6,height=4) 
ggplot(data=filter(median_commensals, DrugX=='Dicumarol'), 
       aes(x=DrugXConc, y=mean_normAUC, color = Strain_name)) +
  geom_errorbar(aes(ymin=mean_normAUC-SDAUC, ymax=mean_normAUC+SDAUC), 
                width=.05, color = 'black', lty = 'solid') +
  geom_line(size=1) + 
  theme_classic(base_size = 12) +
  labs(x=paste0('Dicumarol',' [µM]'), y='norm. AUC') +
  scale_x_continuous(trans="log", breaks=c(0,2.5,5,10, 20, 40,80,160,320,640))+
  scale_color_manual(values=c('slategrey', 'ivory3','azure3','ivory4',
                              'slategray3','#fbb4b9','#f768a1', 'dimgray', 
                              'slategray1','#ae017e')) +
  coord_fixed(ratio = 4)
dev.off()

pdf(here('poster/figures',"EDFig11a.pdf"),width=6,height=4) 
ggplot(data=filter(median_commensals, DrugX=='Tolfenamic acid'), 
       aes(x=DrugXConc, y=mean_normAUC, color = Strain_name)) +
  geom_errorbar(aes(ymin=mean_normAUC-SDAUC, ymax=mean_normAUC+SDAUC), 
                width=.05, color = 'grey', lty = 'solid') +
  geom_line(size=1) + 
  theme_classic(base_size = 12) +
  labs(x=paste0('Tolfenamic Acid',' [µM]'), y='norm. AUC') +
  scale_x_continuous(trans="log", breaks=c(0,2.5,5,10, 20, 40,80,160,320,640)) +
  scale_color_manual(values=c('slategrey', 'ivory3','azure3','ivory4',
                              'slategray3','#fbb4b9','#f768a1', 'dimgray', 
                              'slategray1','#ae017e')) +
  coord_fixed(ratio = 4)
dev.off()

# export source data
if (dir.exists(here('poster/source_data'))){
  
  median_commensals %>% 
    filter(DrugX=='Dicumarol') %>% 
    write_tsv(here('poster/source_data', 'EDFig9c.tsv'))
  
  median_commensals %>% 
    filter(DrugX=='Tolfenamic acid') %>% 
    write_tsv(here('poster/source_data', 'EDFig11a.tsv'))
}