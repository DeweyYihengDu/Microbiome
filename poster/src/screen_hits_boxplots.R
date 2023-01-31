# ##############################################################################
#
## Boxplots with screen hits
#
# ##############################################################################

library("tidyverse")
library("ggpubr")
library("here")

abx.colors.df <- read_delim(here('poster/data','ABX_color_code.csv'),delim=';')
abx.colors <- abx.colors.df$Color
abx.colors.light <- alpha(abx.colors, 0.5)
names(abx.colors) <- abx.colors.df$ABX_class
names(abx.colors.light) <- abx.colors.df$ABX_class

combined_hits <- read_tsv(here('poster/data', "combined_hits.tsv"), 
                          col_types = cols(
                            target_species = col_character(),
                            prestwick_ID = col_character(),
                            human_use = col_logical(),
                            n_hit = col_double(),
                            hit = col_logical()))

antibiotic_annotation <- read_csv(here('poster/files','ABX_annotation.csv'))


species_annotation <-  read_tsv(here('poster/files', "species_annotation.tsv"))

combined_pv <- read_tsv(here('poster/data', "combined_pv.tsv"), 
                        col_types = cols(
                          NT_code = col_character(),
                          prestwick_ID = col_character(),
                          pv_comb_fdr_BH = col_double(),
                          AUC = col_double(),
                          hit = col_logical(),
                          target_species = col_character()))

  
annotated_antibiotics <- left_join(antibiotic_annotation,combined_hits, 
                                   by=c('prestwick_ID','target_species'))

antibiotics_per_class <- annotated_antibiotics %>% 
  group_by(class) %>% 
  count()
annotated_antibiotics <- left_join(annotated_antibiotics, 
                                   antibiotics_per_class, 
                                   by=('class'))
annotated_antibiotics <- annotated_antibiotics %>% 
  mutate(label= paste0(class,' (', as.character(n), ')'))



# ##############################################################################
# focus on tetracyclines, macrolides, quinolones and beta-lactams

annotated_main_ABX_classes <- annotated_antibiotics %>% 
  filter(EUCAST_Comparison %in% 
           c("beta-lactams", "(Fluoro-)quinolones", 'Sulfonamide', 
             'Tetracyclines', 'Aminoglycosides',
             'Macrolides, lincosamides and streptogramins'))

pdf(here('poster/figures',"EDFig4a.pdf"),width=18,height=12)
  ggplot(data=annotated_main_ABX_classes, 
         aes(x=reorder(EUCAST_Comparison, n_hit, median), y=n_hit)) +
  geom_boxplot(aes(color = EUCAST_Comparison, fill=EUCAST_Comparison),
               show.legend = FALSE, outlier.shape=NA) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, fill='black', 
               color='black', stackratio = 0.75, position='identity', 
               binwidth=1) +
  theme_linedraw(base_size = 8) +
  # add labels
  labs(title = "Main classes - Screen", 
       y = "number of inhibited strains", x = NULL) +
  # change font size
  theme(title = element_text(face = "bold"), 
        axis.title = element_text(face = "bold", size = 8)) +
  # change axis label
  theme(aspect.ratio = 1,legend.title = element_blank(), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=abx.colors) + 
  scale_fill_manual(values=abx.colors.light)
dev.off()
  

# ##############################################################################
# focus on quinolones

annotated_quinolones <- annotated_antibiotics %>% 
  filter(EUCAST_Comparison =='(Fluoro-)quinolones')

quinolones_per_class <- annotated_quinolones %>% 
  group_by(more_info) %>% 
  count()
annotated_quinolones <- left_join(annotated_quinolones, 
                                  quinolones_per_class, by=('more_info'))
annotated_quinolones <- annotated_quinolones %>% 
  mutate(label= paste0(more_info,' (', as.character(n.y), ')')) %>% 
  filter(!is.na(more_info))

pdf(here('poster/figures',"EDFig4b.pdf"),width=18,height=12)
  ggplot(data=annotated_quinolones, 
         aes(x=reorder(label, n_hit, median), y=n_hit)) +
  geom_boxplot(outlier.shape=NA, lwd=0.75, 
               color=abx.colors.df %>% 
                 filter(str_detect(ABX_class, 'quinolones')) %>% 
                 pull(Color), 
               fill=abx.colors.df %>% 
                 filter(str_detect(ABX_class, 'quinolones')) %>% 
                 pull(Color) %>% alpha(., 0.5)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, fill='black', 
               color='black', stackratio = 0.75, 
               position='identity', binwidth=1)+
  theme_linedraw(base_size = 8) +
  # add labels
  labs(title = "Quinolones - Screen", y = "number of inhibited strains", 
       x = NULL) +
  # change font size
  theme(aspect.ratio = 1,title = element_text(face = "bold"), 
        axis.title = element_text(face = "bold", size = 8)) +
  # change axis label
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()

# ##############################################################################
# focus - Bacteroides and beta-lactams

beta_lactams <- antibiotic_annotation %>% 
  filter(EUCAST_Comparison =="beta-lactams")
beta_lactams <- left_join(beta_lactams, combined_pv, by='prestwick_ID')

Bacteroidetes <- species_annotation %>% 
  filter(Genus =='Bacteroides' & Screen == TRUE)
Bacteroidetes_betalactams <- inner_join(Bacteroidetes, beta_lactams, 
                                        by='NT_code')
Bacteroidetes_betalactams_hits <- Bacteroidetes_betalactams %>% 
  filter(hit==TRUE)
Bacteroidetes_betalactams_hits <- Bacteroidetes_betalactams_hits %>% 
  group_by(prestwick_ID) %>% 
  count()
beta_lactams <- antibiotic_annotation %>% 
  filter(class == 'Beta-lactam')
Bacteroidetes_betalactams_hits <- left_join(beta_lactams, 
                                            Bacteroidetes_betalactams_hits, 
                                            by='prestwick_ID')
Bacteroidetes_betalactams_hits$n[is.na(Bacteroidetes_betalactams_hits$n)] <- 0
Bacteroidetes_betalactams_counts <- beta_lactams %>% 
  group_by(subgroup) %>% 
  count()
Bacteroidetes_betalactams_hits <- left_join(Bacteroidetes_betalactams_hits, 
                                            Bacteroidetes_betalactams_counts, 
                                            by='subgroup')
Bacteroidetes_betalactams_hits <- Bacteroidetes_betalactams_hits %>% 
  mutate(label = paste0(subgroup,' (', as.character(n.y), ')'))


pdf(here('poster/figures',"EDFig4c.pdf"),width=18,height=12)
  ggplot(data=Bacteroidetes_betalactams_hits, 
         aes(x=reorder(label, n.x, median), y=n.x)) +
  geom_boxplot(outlier.shape=NA, lwd=1, 
               color=abx.colors.df %>% 
                 filter(ABX_class=='beta-lactams') %>% 
                 pull(Color), 
               fill = abx.colors.df %>% 
                 filter(ABX_class=='beta-lactams') %>% 
                 pull(Color) %>% 
                 alpha(., 0.5)) + 
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize=0.5, fill='#bdbdbd', color='#bdbdbd') +
  theme_linedraw(base_size = 12) +
  ylim(0,8) +
  # add labels
  labs(y = "number of inhibited Bacteroides (8)", x = NULL) +
  # change font size
  theme(axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
  # change axis label
  theme(aspect.ratio = 1, legend.title = element_blank(), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()

# export source data
if (dir.exists(here('source_data'))){
  
  annotated_main_ABX_classes %>% 
    select(prestwick_ID, chemical_name, EUCAST_Comparison, n_hit, n) %>% 
    write_tsv(here('source_data', 'EDFig4a.tsv'))
  
  annotated_quinolones %>% 
    select(prestwick_ID, chemical_name, EUCAST_subgroup, n_hit, n.y) %>% 
    write_tsv(here('source_data', 'EDFig4b.tsv'))
  
  Bacteroidetes_betalactams_hits %>% 
    select(prestwick_ID, chemical_name, EUCAST_subgroup, n.x, n.y) %>% 
    write_tsv(here('source_data', 'EDFig4c.tsv'))
}


