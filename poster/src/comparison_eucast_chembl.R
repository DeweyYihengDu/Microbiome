# ##############################################################################
#
## Comparisons with EUCAST and ChEMBL
#
# ##############################################################################

library("here")
library("tidyverse")
library("readxl")
library("glue")
library("ggpattern")
library("stringi")


chembl <- read_tsv(here("poster/data/chembl_assays.tsv"), col_types=cols(
  name = col_character(),
  ATC = col_character(),
  tax_id = col_character(),
  species_name = col_character(),
  assay_type = col_character(),
  value = col_double(),
  unit = col_character()
))

chembl <- chembl %>% filter(!stringi::stri_detect_regex(species_name, " sp\\."))

compounds <- read_tsv(here("poster/data/compounds.txt"))

eucast <- read_tsv(here("poster/data/median_eucast.tsv"))

abundance <- read_tsv(here("poster/data/most_abundant_bugs_full.tsv")) %>% 
  select(-NCBI_id) %>% 
  mutate(specI = stri_sub(mOTU, -5, -2))

taxon_map <- read_tsv(here("poster/data/prok-refdb-v11.0.0_specI-v2_clusters-v1.map"), 
                      col_names = c("specI", "NCBI_id", "representative")) %>% 
  separate(NCBI_id, "tax_id", sep = "\\.", extra = "drop") %>% 
  mutate(specI = stri_sub(specI, -4, -1)) %>% 
  unique()

commensals <- taxon_map %>% select(-representative) %>% unique() %>%
  left_join(abundance, by="specI") %>% 
  group_by(tax_id) %>% 
  summarise(in_healthy = any(!is.na(prevalence)), 
            prevalence = median(prevalence, na.rm = TRUE), 
            average_rel_ab = median(average_rel_ab, na.rm = TRUE), 
            median_rel_ab = median(median_rel_ab, na.rm = TRUE))

eucast_names <- eucast %>% select(Species) %>% unique()

eucast_names <- eucast_names %>% 
  mutate(name = stri_replace_all_regex(
    Species, 
    " *(DIN|etest|EUCAST|spp.*|group.*|coagulase .*|ATCC .*| MRS.| MSS.|H37RV)$", 
    "", opts_regex=stri_opts_regex(case_insensitive=TRUE)))
eucast_names <- eucast_names %>% 
  mutate(name = stri_replace_all_fixed(name, "Staphylococcus xylosis", 
                                       "Staphylococcus xylosus"))
eucast_names <- eucast_names %>% 
  mutate(name = stri_replace_all_fixed(name, "Streptococcus, viridans", 
                                       "Streptococcus viridans"))
eucast_names <- eucast_names %>% 
  filter(!name %in% c("Anaerobic cocci Gram positive", 
                      "Candida paratropicalis", "Pichia anomala"))

species <- bind_rows(
  chembl %>% select(tax_id),
  commensals %>% select(tax_id)
) %>% unique()

species_class <- read_tsv(here("poster/data/species_class.tsv"), col_types = "cccc")

eucast_mapping <- bind_rows(
  species_class %>% filter(id == query) %>% right_join(eucast_names) %>% 
    filter(!is.na(id)),
  eucast_names %>% 
    inner_join(species_class %>% filter(rank == "species"), 
               by = c("name" = "query")) %>% 
    select(-name, name = name.y)
)

eucast <- eucast %>% 
  filter(!is.na(chemical_name)) %>% 
  right_join(eucast_mapping %>% 
               select(-query)) 

commensals <- commensals %>% 
  left_join(species_class %>% 
              filter(rank == "species"), by = c("tax_id" = "query")) %>% 
  select(name, id, in_healthy, prevalence, average_rel_ab, median_rel_ab) %>% 
  unique()

commensals <- commensals %>% 
  group_by(name, id) %>% 
  summarise(in_healthy = any(in_healthy), 
            prevalence = median(prevalence, na.rm = TRUE), 
            average_rel_ab = median(average_rel_ab, na.rm = TRUE), 
            median_rel_ab = median(median_rel_ab, na.rm = TRUE) ) 

bacteria <- species_class %>% 
  group_by(query) %>% 
  filter(any(name == "Bacteria"))

# map subspecies to species
d_chembl <- chembl %>% semi_join(bacteria, by=c("tax_id" = "query")) %>% 
  inner_join(species_class %>% filter(rank == "species"), 
             by=c("tax_id" = "query")) %>% 
  select(-species_name, -rank, -tax_id, species_name = name.y, 
         compound_name = name.x, tax_id = id)

d_eucast <- eucast %>% filter(rank == "species") %>%
  inner_join(compounds, by=c("chemical_name"="Drug_Name")) %>% 
  semi_join(bacteria, by=c("id" = "query")) 

abx_species <- bind_rows(
  d_chembl %>% select(compound_name, ATC, species_name, tax_id),
  d_eucast %>% select(compound_name = chemical_name, ATC, 
                      species_name = name, tax_id = id)
) %>% unique()

abx_species <- abx_species %>% left_join(commensals, by = c("tax_id" = "id"))

# abx_species %>% select(species_name, in_healthy) %>% unique() %>% 
# group_by(in_healthy) %>% count() 

lm_ann_ <- read_tsv(here("poster/data/species_with_counts-LM_Disease association.tsv"))
lm_ann <- lm_ann_ %>% filter(!stri_detect_regex(
  species_name, "(?<!complex), sp\\."))
lm_ann <- lm_ann %>% 
  left_join(abx_species %>% select(species_name, tax_id) %>% unique()) 

lm_ann <- bind_rows(
  lm_ann %>% filter(!is.na(NT_code), is.na(tax_id)) %>% select(-tax_id) %>% 
    left_join(species_class, by = c("species_name"= "name")) %>% 
    select(-query) %>% rename(tax_id = id) %>% unique(),
  lm_ann %>% filter(is.na(NT_code) | !is.na(tax_id))
)

lm_ann <- lm_ann %>% 
  left_join(commensals %>% ungroup() %>% select(-name), by=c("tax_id"="id"))


our_abx <- read_xlsx(here("poster/data/MICs.xlsx")) %>% 
  select(Drug_Abbreviation, compound_name = Drug_Name, Class, MG_gmol) %>% 
  left_join(compounds, by = c("compound_name" = "Drug_Name"))

our_abx <- our_abx %>% 
  mutate(Class = replace_na(Class, "")) %>% 
  mutate(Class_ = fct_lump(Class, n = 5)) %>% arrange(Class_) 

our_mic <-  read_xlsx(here("poster/data/MICs.xlsx"), sheet = "MICs")

our_mic_species <- our_mic %>% 
  select(Drug_Abbreviation, NT_code) %>% 
  unique() %>% 
  left_join(our_abx, by = "Drug_Abbreviation")

our_mic_species <- lm_ann %>% filter(!is.na(NT_code)) %>% 
  separate_rows(NT_code) %>% 
  right_join(our_mic_species, by = "NT_code") %>% 
  select(-NT_code) %>% unique() %>% 
  select(compound_name, ATC, species_name, tax_id, in_healthy, 
         prevalence, average_rel_ab, median_rel_ab) 

combined_species <- bind_rows(
  our_mic_species %>%   mutate(our_mic = TRUE),
  abx_species %>% mutate(existing_mic = TRUE) %>% select(-name)
) %>% group_by(compound_name, ATC, species_name, tax_id, in_healthy, 
               prevalence, average_rel_ab, median_rel_ab) %>% 
  summarise( our_mic = any(!is.na(our_mic)), 
             existing_mic = any(!is.na(existing_mic)) ) %>% 
  ungroup()

combined_species <- combined_species %>% 
  mutate(new_mic = our_mic & !existing_mic) %>% 
  filter(!stri_detect_fixed(compound_name, "Sulfamethoxazole"))


lm_ann__ <- lm_ann %>% mutate(prevalent = prevalence > 0.01) %>% 
  select(species_name, prevalent, `Pathobiont/Pathogen/Disease-associated`)

lm_ann__ <- lm_ann__ %>% 
  mutate( prevalent = ifelse( is.na(prevalent), FALSE, prevalent) )

labels <- tribble(
  ~prevalent, ~`Pathobiont/Pathogen/Disease-associated`, 
  ~new_mic, ~category, ~species_category,
  FALSE, FALSE, FALSE, "previously known: rare non-pathogenic species", 
  "rare non-pathogenic species", 
  TRUE, FALSE, FALSE, "previously known: common non-pathogenic species", 
  "common non-pathogenic species",
  TRUE, FALSE, TRUE, "new data: common non-pathogenic species", 
  "common non-pathogenic species",
  FALSE, TRUE, FALSE, "previously known: pathogenic species", 
  "pathogenic species",
  TRUE, TRUE, FALSE, "previously known: potentially pathogenic species (pathobionts occurring in >1% of people)", 
  "potentially pathogenic species",
  TRUE, TRUE, TRUE, "new data: for potentially pathogenic species", 
  "potentially pathogenic species", 
) %>% mutate(category = fct_rev(fct_inorder(category)), 
             species_category = fct_rev(fct_inorder(species_category)))

combined_species__ <- combined_species %>% 
  left_join(lm_ann__, by="species_name") %>% 
  left_join(labels, by = c("new_mic", "prevalent", 
                           "Pathobiont/Pathogen/Disease-associated")) %>% 
  ungroup() %>% 
  select(compound_name, species_name, category, 
         `Pathobiont/Pathogen/Disease-associated`) %>% 
  left_join(our_abx %>% select(compound_name, Class = Class_)) %>% 
  ungroup() %>% 
  arrange(Class, compound_name) %>% 
  mutate(compound_name = fct_inorder(compound_name)) %>% 
  unique()


combined_species___ <- combined_species__ %>% arrange(category,) %>% 
  group_by(Class, compound_name, 
           `Pathobiont/Pathogen/Disease-associated`, category) %>% 
  count() %>% 
  group_by(Class, compound_name, 
           `Pathobiont/Pathogen/Disease-associated`) %>%  
  arrange(desc(category)) %>% 
  mutate( x2 = cumsum(n), x1 = x2-n, 
          y1 = -as.integer(Class)-as.integer(compound_name)+
            as.integer(!`Pathobiont/Pathogen/Disease-associated`)/2.5, 
          y2 = y1 - 1/2.5 ) %>% 
  ungroup()

breaks_compound_name <- combined_species___ %>% ungroup() %>% 
  filter(!`Pathobiont/Pathogen/Disease-associated`) %>% 
  select(compound_name, y2) %>% unique()

palette_ <- c("#c51b7d", "#f1b6da", "#fde0ef", "#000000", "#aaaaaa", "#cccccc")

combined_species___ <- combined_species__ %>% arrange(category,) %>% 
  group_by(Class, compound_name, 
           `Pathobiont/Pathogen/Disease-associated`, category) %>% 
  count() %>% 
  group_by(Class, compound_name, `Pathobiont/Pathogen/Disease-associated`) %>%  
  arrange(desc(category)) %>% 
  mutate( x2 = cumsum(n), x1 = x2-n, 
          y1 = -as.integer(Class)-as.integer(compound_name)+
            as.integer(!`Pathobiont/Pathogen/Disease-associated`)/2.5, 
          y2 = y1 - 1/2.5 ) %>% 
  ungroup() %>% 
  mutate(COLUMN = as.ordered(Class) >= "Peptide")

combined_species___ %>% ggplot(aes(fill = category)) +
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2)) +
  geom_text(
    data = combined_species___ %>% select(COLUMN, Class, y2) %>% 
      group_by(Class) %>% filter(y2 == max(y2)) %>% unique(),
    aes(x = 0, y = y2+1, label = Class), color = "black", size = 4, hjust = 0,
    inherit.aes = F
  ) +
  scale_fill_manual(values = palette_, guide = guide_legend(reverse=TRUE)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(breaks = breaks_compound_name$y2, 
                     labels = breaks_compound_name$compound_name, 
                     expand = c(0,0.5)) +
  theme_minimal() +
  xlab("Number of species with known MIC") +
  ylab("") +
  theme(legend.position = "bottom", legend.direction = "vertical", 
        # axis.text.y = element_blank(),
        # panel.spacing.y = unit(0.1, "lines"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text.y = element_text(angle = 0, hjust=1), 
        legend.title = element_blank(),
        strip.text = element_blank()
  ) +
  facet_wrap(~COLUMN, scale = "free_y")

ggsave(here("poster/figures", "EDFig3c.pdf"), width = 10, height = 5)

# export source data
if (dir.exists(here('poster/source_data'))){
  combined_species___ %>% 
    select(1:5) %>% 
    write_tsv(here('poster/source_data', 'EDFig3c.tsv'))
}


combined_species__2 <- combined_species %>%
  left_join(lm_ann__, by="species_name") %>% 
  left_join(labels, by = c("new_mic", "prevalent", 
                           "Pathobiont/Pathogen/Disease-associated")) %>% 
  ungroup() %>% 
  select(compound_name, species_name, category, species_category, 
         `Pathobiont/Pathogen/Disease-associated`, existing_mic, our_mic, 
         new_mic, prevalent) %>% 
  left_join(our_abx %>% select(compound_name, Class = Class_)) %>% 
  ungroup() %>% 
  arrange(Class, compound_name) %>% 
  mutate(compound_name = fct_inorder(compound_name)) %>% 
  unique()


d <- combined_species__2 %>% 
  group_by(species_category, 
           `Pathobiont/Pathogen/Disease-associated`, existing_mic, 
           our_mic, new_mic, prevalent) %>% count()

dd <- d %>% group_by(species_category) %>% 
  arrange(!existing_mic, our_mic) %>% 
  mutate( x2 = cumsum(n), x1 = x2 - n, 
          y1 = as.integer(species_category)/2.5, y2 = y1 - 1/2.5 + 1/8)

breaks_species_category <- dd %>% mutate( y = (y1+y2)/2 ) %>% 
  select(species_category, y) %>% unique()

ddt <- d %>% group_by(`Pathobiont/Pathogen/Disease-associated`, prevalent) %>% 
  summarise( increase = sum(n[new_mic]) / sum(n[!new_mic]), 
             label = paste0("+", round(100*increase), "% new MICs" )) %>% 
  filter(increase > 0)

ddt <- ddt %>% left_join(dd %>% filter(new_mic))

dd %>% 
  ggplot(aes(fill = new_mic)) +
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2)) +
  geom_rect_pattern(data = dd %>% filter(existing_mic, our_mic),
                    aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
                    pattern_size = 0, pattern_fill = "#444444", 
                    pattern_spacing = 0.01, pattern_angle = 45, 
                    pattern_density = 0.25,
                    show.legend = F) +
  geom_text(data = ddt, aes(x = x2+10, y = (y1+y2)/2, label = label), 
            hjust = 0, size = 3) +
  scale_fill_manual(values = c("#bbbbbb", "#444444"), 
                    labels = c("MICs exist in EUCAST or ChEMBL    ", 
                               "MICs from this screen                                                                  ")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(breaks = breaks_species_category$y, 
                     labels = breaks_species_category$species_category, 
                     expand = c(0,0.1)) +
  theme_minimal() +
  xlab("Number of MICs (drug-species pairs)") +
  ylab("") +
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        axis.text.y = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.title = element_blank(),
  )

ggsave(here("poster/figures/Fig1b.pdf"), width = 5, height = 2)

# export source data
if (dir.exists(here('poster/source_data'))){
  dd %>% 
    arrange(species_category) %>% 
    select(1:7) %>% 
    write_tsv(here('poster/source_data', 'Fig1b.tsv'))
}



chembl_mics <- d_chembl %>% filter(assay_type == "MIC", unit == "ug.mL-1") %>% 
  semi_join(our_mic_species) %>% 
  group_by(compound_name, species_name, tax_id) %>% 
  summarise(MIC_ug_ml = median(value)) %>% 
  ungroup() %>% mutate(source = "ChEMBL") %>% 
  filter(!stri_detect_fixed(compound_name, "Sulfamethoxazole"))

eucast_mics <- d_eucast %>% rename(compound_name = chemical_name, 
                                   species_name = name, tax_id = id) %>% 
  semi_join(our_mic_species) %>% 
  group_by(compound_name, species_name, tax_id) %>% 
  summarise(MIC_ug_ml = median(Conc_ug_ml)) %>% ungroup() %>% 
  mutate(source = "EUCAST") %>% 
  filter(!stri_detect_fixed(compound_name, "Sulfamethoxazole"))

our_mics <- our_mic %>% left_join(our_abx, by = "Drug_Abbreviation") %>% 
  left_join(lm_ann %>% filter(!is.na(NT_code)) %>% separate_rows(NT_code)) %>% 
  unique() %>% 
  filter(R1_qualifier == "=", R2_qualifier == "=") %>% rowwise() %>% 
  mutate(MIC_ug_ml = mean(R1_value, R2_value)) %>% ungroup() %>% 
  group_by(compound_name, species_name, tax_id, ATC) %>% 
  summarise(MIC_ug_ml = median(MIC_ug_ml)) %>% ungroup() %>% 
  select(compound_name, ATC, species_name, MIC_ug_ml, tax_id) %>% 
  mutate(source = "new data") %>% 
  filter(!stri_detect_fixed(compound_name, "Sulfamethoxazole"))

mics <- bind_rows( 
  inner_join(chembl_mics, eucast_mics, 
             by = c("compound_name", "species_name", "tax_id")),
  inner_join(our_mics, chembl_mics, 
             by = c("compound_name", "species_name", "tax_id")),
  inner_join(our_mics, eucast_mics, 
             by = c("compound_name", "species_name", "tax_id"))
)

mics <- mics %>% mutate(combo = as.character(glue("{source.x} vs. {source.y}")))

mics <- bind_rows(
  mics %>% mutate(dataset = "All MICs"),
  mics %>% group_by(compound_name, tax_id) %>% filter(n() == 3) %>% 
    mutate(dataset = "Shared across all datasets")
)

mic_corr <- mics %>% group_by(dataset, combo) %>% 
  summarise( corr = cor(MIC_ug_ml.x, MIC_ug_ml.y, method = "s"), n = n(), 
             cor.test(MIC_ug_ml.x, MIC_ug_ml.y, method = "s")$p.value )

bb <- 10**(-2:2)
ll <- unlist(lapply(bb, as.character))

mics %>% ggplot(aes(MIC_ug_ml.x, MIC_ug_ml.y)) + geom_point(size = 0.5) +
  geom_text( data = mic_corr, 
             aes(x = min(mics$MIC_ug_ml.x), y = max(mics$MIC_ug_ml.y),
                 label = paste0("r[s] == ",round(corr, digits = 2))), 
             hjust = 0, vjust = 1, parse=T ) +
  geom_text( data = mic_corr, 
             aes(x = min(mics$MIC_ug_ml.x), y = max(mics$MIC_ug_ml.y),
                 label = paste0("N == ", n)), hjust = 0, vjust = 3, 
             parse=TRUE ) +
  facet_grid(dataset ~ combo) + 
  scale_x_log10(name = "MIC (µg/ml)", breaks = bb, labels = ll) + 
  scale_y_log10(name = "MIC (µg/ml)", breaks = bb, labels = ll) +
  coord_equal() + theme_minimal() +
  theme( 
    strip.background = element_rect(fill = "lightgrey", color = "lightgrey"),
    panel.grid.minor = element_blank()
  )

ggsave(here("poster/figures", "EDFig3b.pdf"), height = 4, width = 6)

# export source data
if (dir.exists(here('poster/source_data'))){
  mics %>% 
    write_tsv(here('poster/source_data', 'EDFig3b.tsv'))
}
