library("here")
library("tidyverse")
library('ggpubr')
library('scales')

# Set colors 
color_min <- c('white','white')
color_max <- c("#156077","#f46f20")

# load all checkerboards (CB) for commensals
CB <- read_tsv(here('poster/data','CB_community_raw_data.tsv'))

# combine replicates
CB_combined <- CB %>% 
  group_by(DrugX, DrugXConc, ABX, ABXConc, Strain) %>% 
  summarise(median_normAUC=median(normAUC))

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
              theme_bw(base_size = 12) +
              labs(fill='norm.\n AUC') +
              coord_fixed() +
              xlab(paste(part_DrugX, '[µM]', sep=" ")) + 
              ylab(paste(part_ABX, '[µM]', sep=" ")) +
              theme(axis.text.x= element_text(size = 12, angle = 90, hjust=1), 
                    axis.text.y= element_text(size = 12)) +
              ggtitle(part_Strain)
      assign(paste(part_ABX , part_DrugX, part_Strain, sep='_'), myplot)
    }
  }
}


pdf(here('poster/figures', 'Fig3a.pdf'),width=3,height=2) 
ggarrange(Erythromycin_Dicumarol_community,
          ncol=1, nrow=1)
dev.off()

pdf(here("poster/figures", 'EDFig10a.pdf'),width=3,height=2) 
ggarrange(Erythromycin_Benzbromarone_community,
          ncol=1, nrow=1)
dev.off()

pdf(here('poster/figures', 'EDFig11b.pdf'),width=3,height=2) 
ggarrange(`Erythromycin_Tolfenamic Acid_community`,
          ncol=1, nrow=1)
dev.off()

# export source data
if (dir.exists(here('poster/source_data'))){
  # Fig3a
  CB_combined %>% 
    filter(DrugX=='Dicumarol') %>% 
    select(-Strain) %>% 
    write_tsv(here('poster/source_data', 'Fig3a.tsv'))
  
  CB_combined %>% 
    filter(DrugX=='Benzbromarone') %>% 
    select(-Strain) %>% 
    write_tsv(here('poster/source_data', 'EDFig10a.tsv'))
  
  CB_combined %>% 
    filter(DrugX=='Tolfenamic Acid') %>% 
    select(-Strain) %>% 
    write_tsv(here('poster/source_data', 'EDFig11b.tsv'))
}
