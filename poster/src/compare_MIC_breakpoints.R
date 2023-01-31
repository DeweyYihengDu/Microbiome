# ##############################################################################
# 
## Compare MICs to clinical breakpoints
#
# ##############################################################################

library("here")
library("tidyverse")
library("readxl")
library("ggpubr")

abx.colors.df <- read_delim(here('poster/files','ABX_color_code.csv'),delim=';')
abx.colors <- abx.colors.df$Color
names(abx.colors) <- abx.colors.df$ABX_class

# ##############################################################################
# Import breakpoints and make some formatting adjustments

breakpoints <- read_tsv(here('poster/data','eucast_2019_v9_breakpoints.tsv'))

breakpoints <- breakpoints %>% 
  mutate(Anaerobier = case_when(
    species == 'Anaerobes, Grampositive' ~ 'Anaerobes, Grampositive', 
    species == 'Anaerobes, Gramnegative' ~ 'Anaerobes, Gramnegative', 
    species == 'PK PD breakpoints' ~ 'PK PD breakpoints',
    species != 'Anaerobes, Grampositive' & 
      species !='Anaerobes, Gramnegative' ~ ''))

breakpoints <- breakpoints %>% 
  mutate(ABX_class = case_when(
    drug_group == 'Penicillins' ~ 'beta-lactams', 
    drug_group == 'Cephalosporins' ~ 'beta-lactams', 
    drug_group == 'Carbapenems' ~ 'beta-lactams', 
    drug_group == 'Monobactams' ~ 'beta-lactams', 
    drug_group =='Glycopeptides' ~ 'Glycopeptides and lipoglycopeptides',
    drug_group =='Macrolides' ~ 'Macrolides, lincosamides and streptogramins',
    drug_group =='Macrolides and lincosamides' ~ 
      'Macrolides, lincosamides and streptogramins',
    is.na(drug_group) ~ 'Miscellaneous agents',
    TRUE ~ drug_group))

breakpoints$ABX_class[breakpoints$ABX_class == "Fluoroquinolones"] <- 
  "(Fluoro-)quinolones"

breakpoints <- breakpoints %>% 
  mutate(ABX_subclass = case_when(
    drug_group == 'Penicillins' ~ 'Penicillins', 
    drug_group == 'Cephalosporins' ~ 'Cephalosporins', 
    drug_group == 'Carbapenems' ~ 'Carbapenems', 
    drug_group == 'Monobactams' ~ 'Monobactams', 
    is.na(drug_group) ~ 'Miscellaneous agents',
    TRUE ~ ABX_class))


# ##############################################################################
# combine with our MICs

breakpoints_main <- breakpoints %>% 
  filter(Anaerobier != 'PK PD breakpoints') %>% 
  select(c(ABX_class, sensitive)) %>% 
  mutate(study = 'clinical breakpoints') %>% 
  mutate(mean_qualifier = '=') %>% 
  filter(ABX_class %in% 
           c("beta-lactams", "(Fluoro-)quinolones", 'Sulfonamide', 
             'Tetracyclines', 'Aminoglycosides',
             'Macrolides, lincosamides and streptogramins')) %>% 
  mutate(shp ='grey')

MICs <- read_xlsx(here('poster/data','MICs.xlsx'), sheet = 3)

ABX <- read_xlsx(here('poster/data','MICs.xlsx'), sheet = 1)

species_annotation <-  read_tsv(here('poster/files','species_annotation.tsv'))

MICs <- full_join(MICs, ABX, by='Drug_Abbreviation') %>% 
  mutate(mean_MICs = (R1_value + R2_value)/2) %>% 
  mutate(mean_qualifier = case_when(
    R1_qualifier == '=' & R2_qualifier ==  '=' ~ "=",
    R1_qualifier == '>' & R2_qualifier ==  '>' ~ ">",
    R1_qualifier == '<' & R2_qualifier ==  '<' ~ "<",
    R1_qualifier == '>' & R2_qualifier ==  '=' ~ "=",
    R1_qualifier == '=' & R2_qualifier ==  '>' ~ "=",
    R1_qualifier == '<' & R2_qualifier ==  '=' ~ "=",
    R1_qualifier == '=' & R2_qualifier ==  '<' ~ "=")) %>% 
  mutate(shp = case_when(
    mean_qualifier == '=' ~ 'grey', 
    mean_qualifier == '>' ~ 'black',
    mean_qualifier == '<' ~ 'black'))

MICs <- left_join(MICs, species_annotation, by='NT_code')


antibiotics_per_class <- ABX %>% 
  group_by(EUCAST_comparison) %>% 
  count()
MICs <- left_join(MICs, antibiotics_per_class, by='EUCAST_comparison')

######################################################################################################################################
# MICs without expansion for Bacteroides!)

MICs_Screen <- MICs %>% 
  filter(Screen == TRUE)

MICs_Screen_main <- MICs_Screen %>% 
  filter(EUCAST_comparison %in% 
           c("beta-lactams", "(Fluoro-)quinolones", 'Sulfonamide', 
             'Tetracyclines', 'Aminoglycosides',
             'Macrolides, lincosamides and streptogramins'))

names(MICs_Screen_main)[names(MICs_Screen_main)=="mean_MICs"] <- "sensitive"
names(MICs_Screen_main)[names(MICs_Screen_main)=="EUCAST_comparison"] <- 
  "ABX_class"
MICs_Screen_main <- MICs_Screen_main %>% 
  mutate(study = 'MICs')

MICs_short <- MICs_Screen_main %>% 
  select(c(sensitive, ABX_class, study, mean_qualifier, shp))

# Combine breakpoints and MICs
breakpoints_MICs <- rbind(MICs_short, breakpoints_main)

# Set levels manually to determine order
lev <- MICs_short %>% 
  group_by(ABX_class) %>% 
  summarise(MICs_median = median(sensitive)) %>% 
  mutate(rank = dense_rank(MICs_median)) %>% 
  select(c(ABX_class, rank)) %>% 
  arrange(desc(rank))

breakpoints_MICs$ABX_class <- factor(breakpoints_MICs$ABX_class, 
                                     levels=lev$ABX_class)

breakpoints_MICs$study <- as.factor(breakpoints_MICs$study)
breakpoints_MICs$study <- factor(breakpoints_MICs$study,
                                 levels(breakpoints_MICs$study)[c(2,1)])


pdf(here('poster/figures',"Fig1c.pdf"),width=8.27,height=11.69)

  ggplot(data=breakpoints_MICs, aes(x=ABX_class, y=sensitive, fill=study)) +
  geom_boxplot(show.legend = FALSE, outlier.shape=NA, lwd=0.75, 
               aes(color=ABX_class), position=position_dodge(1)) +
  theme_linedraw(base_size = 12) +
  scale_y_continuous(trans='log2')+
  # add labels
  labs(title = "Main classes - MICs", 
       y = "MIC/breakpoint in µg/ml", x = NULL) +
  # change font size
  theme(title = element_text(face = "bold"), 
        axis.title = element_text(face = "bold", size = 12)) +
  # change axis label
  theme(legend.title = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_manual(values=abx.colors) + 
  scale_fill_manual(values=alpha(c('white', '#bdbdbd'),0.5)) +
  coord_fixed(ratio=0.5)

dev.off()

# export source data
if (dir.exists(here('poster/source_data'))){
  breakpoints_MICs %>% 
    transmute(value=sensitive, ABX_class, type=study) %>% 
    write_tsv(here('poster/source_data', 'Fig1c.tsv'))
}