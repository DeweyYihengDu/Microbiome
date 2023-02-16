# ##############################################################################
#
## Rescue of species determined from 16S sequencing
#
# ##############################################################################

library("here")
library("tidyverse")
library('ggpubr')
library('scales')

# Set colors 
color_min <- c('white','white')
color_max <- c("#ff9933","#ff3333")


# load all checkerboards (CB) for commensals
CB <- read_tsv(here('data','CB_community_pathogen_raw_data.tsv'))

# combine replicates
CB_combined <- CB %>% 
  group_by(DrugX, DrugXConc, ABX, ABXConc, Strain) %>% 
  summarise(median_normAUC=median(normAUC), 
            mean_cfu=mean(cfu_ml), 
            sd_cfu=sd(cfu_ml),
            m=mean(log10(cfu_ml)), s=mad(log10(cfu_ml))) %>% 
  mutate(m.top=m+s, m.bottom=m-s)

CB_combined <- CB_combined %>% 
  mutate(cfu = paste(formatC(mean_cfu, format = "e", digits = 2),
                     formatC(sd_cfu, format = "e", digits = 2), sep = '/'))

# Plot all CB for commensals
all_ABX <- unique(CB_combined$ABX)

for (h in seq_along(all_ABX)){
  part_ABX <- all_ABX[h]
  B <- CB_combined %>% filter(ABX == part_ABX)
  all_DrugX <- unique(B$DrugX)
  for (i in seq_along(all_DrugX)){
    part_DrugX = all_DrugX[i]
    A <- B %>% filter(DrugX==part_DrugX)
    all_strains <- unique(A$Strain)
    for (j in seq_along(all_strains)){
      C <- A %>% filter(Strain == all_strains[j])  
      part_Strain <- C$Strain[1]
      myplot <- ggplot(data = C, aes(x = as.factor(DrugXConc), 
                                     y = as.factor(ABXConc))) +
              geom_tile(aes(fill = median_normAUC), colour = 'white') +
              scale_fill_gradientn(name='AUC', 
                                   colours=c(color_min[h], color_max[h]), 
                                   limits = c(0,1.1))+
              geom_text(aes(label=cfu), size=2)+
              theme_bw(base_size = 12) +
              labs(fill='norm.\n AUC') +
              coord_fixed() +
              xlab(paste(part_DrugX, '[µM]', sep=" ")) + 
              ylab(paste(part_ABX, '[µM]', sep=" ")) +
              theme(axis.text.x= element_text(size = 12, angle = 90, hjust=1), 
                    axis.text.y= element_text(size = 12)) +
              ggtitle(part_Strain)
      myplot_alt <- ggplot(data = C, aes(x = as.factor(DrugXConc), 
                                         y = as.factor(ABXConc))) +
        geom_tile(aes(fill = median_normAUC), colour = 'white') +
        scale_fill_gradientn(name='AUC', 
                             colours=c(color_min[h], color_max[h]), 
                             limits = c(0,1.1))+
        geom_point(aes(size=m), pch=21) + 
        geom_point(aes(size=m.top), pch=21, fill=NA, colour='grey') +
        geom_point(aes(size=m.bottom), pch=21, fill=NA, colour='grey') + 
        scale_size_continuous(range = c(1, 20), trans = 'exp') + 
        theme_bw(base_size = 12) +
        labs(fill='norm.\n AUC') +
        coord_fixed() +
        xlab(paste(part_DrugX, '[µM]', sep=" ")) + 
        ylab(paste(part_ABX, '[µM]', sep=" ")) +
        theme(axis.text.x= element_text(size = 12, angle = 90, hjust=1), 
              axis.text.y= element_text(size = 12)) +
        ggtitle(part_Strain)
      assign(paste(part_ABX , part_DrugX, part_Strain, sep='_'), myplot)
      assign(paste(part_ABX , part_DrugX, part_Strain, 'alt', sep='_'), 
             myplot_alt)
    }
  }
}


pdf(here('figures', 'Fig3c.pdf'),width=6,height=4) 
ggarrange(`Erythromycin_Dicumarol_community + Enterococcus faecalis_alt`,
          ncol=1, nrow=1)
dev.off()


pdf(here('figures', 'EDFig10c.pdf'),width=6,height=4) 
ggarrange(`Erythromycin_Benzbromarone_community + Enterococcus faecalis_alt`,
          ncol=1, nrow=1)
dev.off()


# export source data
if (dir.exists(here('source_data'))){
  # Fig3c
  CB_combined %>% 
    filter(DrugX=='Dicumarol') %>% 
    select(DrugX, DrugXConc, ABX, ABXConc, median_normAUC, mean_cfu, sd_cfu) %>% 
    write_tsv(here('source_data', 'Fig3c.tsv'))
  
  CB_combined %>% 
    filter(DrugX=='Benzbromarone') %>% 
    select(DrugX, DrugXConc, ABX, ABXConc, median_normAUC, mean_cfu, sd_cfu) %>%
    write_tsv(here('source_data', 'EDFig10c.tsv'))
}