# ##############################################################################
#
## Z-scores for Dox and Ery
#
# ##############################################################################

library("here")
library("tidyverse")

finalAUCs <- read_tsv(here('poster/data', "finalAUCs.tsv"))

# calculate z-scores per Experiment
finalAUCs$z <- ave(finalAUCs$normAUC, finalAUCs$Exp, FUN=scale)

#  Remove replicates that are very different
finalAUCs_wide <- finalAUCs %>% 
  select(c(Drug, ABX, Strain, Replicate, z)) %>%  
  spread(Replicate, z) %>% 
  mutate(ratio=R1/R2)

# Remove data if data point are x fold different 
#   (between x=5 and x=8 no effect on final hitlist)
x = 8

finalAUCs_wide <- finalAUCs_wide %>% 
  mutate(outlier = case_when(ratio > x ~ TRUE,
                             ratio < 1/x ~ TRUE))

finalAUCs_clean <-  finalAUCs_wide  %>% 
  filter(is.na(outlier)) %>% 
  mutate(z_mean=((R1+R2)/2))

# For Erythromycin, NT5001, z-mean is equal to z of R2

final_mean_z_scores <- finalAUCs_clean %>% 
  mutate(mean_z_score = case_when(is.na(z_mean) ~ R2, 
                                  !is.na(z_mean) ~ z_mean)) %>% 
  select(c(Drug, ABX, Strain, mean_z_score))

final_mean_z_scores <- final_mean_z_scores %>% 
  mutate(Exp = case_when(
    ABX=='Erythromycin' & Strain=='NT5001' ~ 'ERY_B.vulgatus',
    ABX=='Erythromycin' & Strain=='NT5002' ~ 'ERY_B.uniformis',
    ABX=='Doxycycline' & Strain=='NT5001' ~ 'DOX_B.vulgatus',
    ABX=='Doxycycline' & Strain=='NT5002' ~ 'DOX_B.uniformis'))

# hitlist if final_mean_z_scores > 3

hitlist <- final_mean_z_scores %>% 
  filter(mean_z_score >3) 

hitlist <- hitlist %>% 
  mutate(followup = case_when(
    Drug == 'Diethylstilbestrol' &	ABX== 'Doxycycline'	& 
      Strain == 'NT5001' ~ 'No',
    Drug == 'Efavirenz'	&	ABX== 'Doxycycline'	& Strain == 'NT5001'	~ 'No', 
    Drug == 'Diclazuril'	&	ABX== 'Doxycycline'	& Strain == 'NT5001'	~ 'No', 
    Drug == 'Clonixin Lysinate'	&	ABX== 'Doxycycline' & 
      Strain == 'NT5001'~ 'No', 
    Drug == 'Meclofenamic acid sodium salt monohydrate'	&	
      ABX == 'Doxycycline' & Strain == 'NT5001' ~ 'No', 
    TRUE ~ 'Yes'))

hitlist_flag <- hitlist %>% 
  mutate(hit= TRUE)

final_mean_z_scores <- left_join(final_mean_z_scores,hitlist_flag )

o = ggplot(data= hitlist, aes(x=Exp, y=mean_z_score)) +
  geom_point(color=alpha('#999999',0.5), size=3, aes(shape=followup)) +
  geom_boxplot(data=final_mean_z_scores,
               aes(x=Exp, y=mean_z_score, color=ABX, fill=ABX), 
               outlier.shape=NA, lwd=0.75) +
  coord_flip() +
  theme_linedraw(base_size = 20) +
  theme(title = element_text(face = "bold"), 
        axis.text.y = element_text(face = "bold", size = 12), 
        axis.title.y = element_blank()) +
  ylab('z-scores') +
  # change axis label
  theme(aspect.ratio = 0.5, legend.title = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values=c("#ff3333","#ff9933")) +
  scale_fill_manual(values=alpha(c("#ff3333","#ff9933"),0.5)) +
  scale_shape_manual(values=c(4,16))

ggsave(here('poster/figures', 'EDFig7b.pdf'),width=8, height=6, units='in')

# export source data
if (dir.exists(here("poster/source_data"))) {
  write_tsv(final_mean_z_scores, file = here('poster/source_data', 'EDFig7b.tsv'))
}

hitlist_export <- hitlist %>% 
  mutate(Strain=case_when(Strain=='NT5001' ~ 'B.vulgatus',
                          Strain=='NT5002' ~ 'B.uniformis')) %>% 
  arrange(desc(mean_z_score))
hitlist_export$Exp <- NULL
hitlist_export$mean_z_score <- round(hitlist_export$mean_z_score, digits = 3)

write_tsv(hitlist_export, file = here("poster/files","hitlist.tsv"))
