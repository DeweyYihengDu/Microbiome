# ##############################################################################
# 
## Checkerboard plots
#
# ##############################################################################

library("here")
library("tidyverse")
library("ggpubr")
library("scales")
library("cowplot")


# Set colors 
color_min <- c('white','white')
color_max <- c("#ff9933","#ff3333")


# load all checkerboards (CB) for commensals
CB <- read_tsv(here('data','CB_commensals_raw_data.tsv'))

# combine replicates
CB_combined <- CB %>% 
  group_by(DrugX, DrugXConc, ABX, ABXConc, Strain_name, Strain) %>% 
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
    all_strains <- unique(A$Strain_name)
    for (j in seq_along(all_strains)){
      C <- A %>% filter(Strain_name == all_strains[j])  
      part_Strain <- C$Strain_name[1]
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
              theme(axis.text.x= element_text(size = 12, angle = 90, 
                                              hjust=1, vjust=0.5),
                    axis.text.y= element_text(size = 12)) +
              ggtitle(part_Strain)
      assign(paste(part_ABX , part_DrugX, part_Strain, sep='_'), myplot)
    }
  }
}



pdf(here('figures', 'EDFig7d.pdf'),width=18,height=12) 
ggarrange(`Erythromycin_Dicumarol_B. vulgatus`,
          `Erythromycin_Benzbromarone_B. vulgatus`, 
          `Erythromycin_Tolfenamic acid_B. vulgatus`, 
          `Erythromycin_Diethylstilbestrol_B. uniformis`, 
          `Erythromycin_Hexestrol_B. uniformis`, 
          `Erythromycin_Lorglumide sodium salt_B. uniformis`, 
          `Doxycycline_Diflunisal_B. vulgatus`, 
          `Doxycycline_Tolfenamic acid_B. vulgatus`, 
          `Doxycycline_Iopanoic acid_B. vulgatus`, 
          `Doxycycline_Tiratricol_B. uniformis`, 
          ncol=6, nrow=2, align='v')
dev.off()


pdf(here('figures', 'EDFig9a.pdf'),width=18,height=12) 
ggarrange(`Erythromycin_Dicumarol_B. fragilis NT`,
          `Erythromycin_Benzbromarone_B. fragilis NT`,
          `Erythromycin_Tolfenamic acid_B. fragilis NT`, 
          `Doxycycline_Diflunisal_B. fragilis NT`, 
          `Doxycycline_Tolfenamic acid_B. fragilis NT`,
          
          `Erythromycin_Dicumarol_B. thetaiotaomicron`, 
          `Erythromycin_Benzbromarone_B. thetaiotaomicron`, 
          `Erythromycin_Tolfenamic acid_B. thetaiotaomicron`,
          `Doxycycline_Diflunisal_B. thetaiotaomicron`, 
          `Doxycycline_Tolfenamic acid_B. thetaiotaomicron`,
          
          `Erythromycin_Dicumarol_P. copri`, 
          `Erythromycin_Benzbromarone_P. copri`, 
          `Erythromycin_Tolfenamic acid_P. copri`,
          `Doxycycline_Diflunisal_P. copri`,
          `Doxycycline_Tolfenamic acid_P. copri`,
          
          `Erythromycin_Dicumarol_B. ovatus`, 
          `Erythromycin_Benzbromarone_B. ovatus`, 
          `Erythromycin_Tolfenamic acid_B. ovatus`,
          `Doxycycline_Diflunisal_B. ovatus`, 
          `Doxycycline_Tolfenamic acid_B. ovatus`,
          
          `Erythromycin_Dicumarol_B. caccae`, 
          `Erythromycin_Benzbromarone_B. caccae`, 
          `Erythromycin_Tolfenamic acid_B. caccae`,
          `Doxycycline_Diflunisal_B. caccae`, 
          `Doxycycline_Tolfenamic acid_B. caccae`,
          
          `Erythromycin_Dicumarol_P. distasonis`, 
          `Erythromycin_Benzbromarone_P. distasonis`, 
          `Erythromycin_Tolfenamic acid_P. distasonis`, 
          `Doxycycline_Diflunisal_P. distasonis`, 
          `Doxycycline_Tolfenamic acid_P. distasonis`,
          
          ncol=5, nrow=6, align='v')
dev.off()


# export source data
if (dir.exists(here('source_data'))){
  
  CB_combined %>% 
    filter(Strain_name %in% c('B. vulgatus','B. uniformis')) %>% 
    write_tsv(here('source_data', 'EDFig7d.tsv'))
  
  CB_combined %>% 
    filter(Strain_name %in% c('B. fragilis NT', 'B. thetaiotaomicron', 
                              'P. copri', 'B. ovatus', 'B. caccae', 
                              'P. distasonis')) %>% 
    write_tsv(here('source_data', 'EDFig9a.tsv'))
}


# ##############################################################################
# also, for pathogens


# load all checkerboards (CB) for commensals
CB <- read_tsv(here('data','CB_pathogens_raw_data.tsv'))

# combine replicates
CB_combined <- CB %>% 
  group_by(DrugX, DrugXConc, ABX, ABXConc, Strain) %>% 
  summarise(median_normAUC=median(normAUC))

# Plot all CB for pathogens
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
        theme(axis.text.x= element_text(size = 12, angle = 90, 
                                        hjust=1, vjust=0.5), 
              axis.text.y= element_text(size = 12)) +
        ggtitle(part_Strain)
      assign(paste(part_ABX , part_DrugX, part_Strain, sep='_'), myplot)
    }
  }
}


pdf(here('figures', 'EDFig9b.pdf'),width=18,height=12) 
ggarrange(`Erythromycin_Dicumarol_Enterococcus faecalis`, 
          `Erythromycin_Benzbromarone_Enterococcus faecalis`, 
          `Erythromycin_Tolfenamic Acid_Enterococcus faecalis`, 
          `Erythromycin_Dicumarol_Enterococcus faecium`, 
          `Erythromycin_Benzbromarone_Enterococcus faecium`, 
          `Erythromycin_Tolfenamic Acid_Enterococcus faecium`, 
          Erythromycin_Dicumarol_SADSM20231, 
          Erythromycin_Benzbromarone_SADSM20231, 
          `Erythromycin_Tolfenamic Acid_SADSM20231`, 
          Doxycycline_Diflunisal_SADSM20231, 
          `Doxycycline_Tolfenamic Acid_SADSM20231`,
          ncol=6, nrow=2, align = 'v')
dev.off()

# export source data
if (dir.exists(here('source_data'))){
  CB_combined %>% 
    write_tsv(here('source_data', 'EDFig9b.tsv'))
}