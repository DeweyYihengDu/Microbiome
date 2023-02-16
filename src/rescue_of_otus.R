# ##############################################################################
#
## Rescue of species determined from 16S sequencing
#
# ##############################################################################

library("here")
library("tidyverse")
library("reshape2")
library("vegan")
library("gridExtra")
library("egg")
library("DESeq2")

samples <- read_tsv(here("data", "antidotes_samples.tsv"), col_types = cols(
  .default = col_character(),
  `16s_primer_index` = col_double(),
  dilution = col_double(),
  time_h = col_double(),
  erythromycin = col_double(),
  col = col_double(),
  count = col_double()
))

samples <- samples %>% 
  mutate(condition = fct_relevel(condition, c("no treatment", "no antidote")))

counts <- read_tsv(here("data", "antidotes_16S_counts.tsv.gz"), 
                   col_types = cols(
                     otu = col_character(),
                     sample = col_character(),
                     count = col_double(),
                     rel_abundance = col_double(),
                     size = col_double(),
                     domain = col_character(),
                     phylum = col_character(),
                     class = col_character(),
                     order = col_character(),
                     family = col_character(),
                     genus = col_character(),
                     species = col_character(),
                     consensus_taxonomy = col_character(),
                     otu_name = col_character()))

otus <- counts %>% 
  select(-sample, -count, -rel_abundance, -size) %>% 
  unique()

count_matrix <- acast(counts, sample ~ otu, value.var = "count")


processSubject <- function(samples) {
  dds_tmt <- DESeq2::DESeqDataSetFromMatrix(
    countData = t(count_matrix[ samples %>% pull(sample), ]),
    colData = samples %>% select(condition),
    design = ~condition)
  
  dds_tmt <- dds_tmt[rowSums(DESeq2::counts(dds_tmt)) >= 10,]
  
  dds_tmt$condition <- relevel(dds_tmt$condition, ref = "no treatment")
  
  dds_tmt <- DESeq2::DESeq(dds_tmt)
  res_tmt <- DESeq2::results(dds_tmt, cooksCutoff=FALSE)
  
  dfRes <- function(s, vs = "no treatment") {
    data.frame(DESeq2::results(dds_tmt, contrast=c("condition", vs, s), 
                               cooksCutoff=FALSE)) %>% 
      rownames_to_column(var = "otu") %>% 
      mutate(antidote = s) %>% tibble()
  }
  
  resNoAntidote <- dfRes("no antidote") %>% 
    select(-antidote) %>% 
    rename_at(vars(-otu), function(x) paste0(x,"_no_antidote"))
  
  resAntidotes <- bind_rows(
    dfRes("20 benzbromarone"), 
    dfRes("20 dicoumarol"),
    dfRes("40 tolfenamic acid")
  ) %>% 
    rename_at(vars(-otu, -antidote), function(x) paste0(x,"_with_antidote"))

  resAntidotes2 <- bind_rows(
    dfRes("20 benzbromarone", "no antidote"), 
    dfRes("20 dicoumarol", "no antidote"),
    dfRes("40 tolfenamic acid", "no antidote")
  ) %>% 
    rename_at(vars(-otu, -antidote), function(x) paste0(x,"_vs_no_antidote"))
  
  resCombined <- 
    resNoAntidote %>% 
    left_join(resAntidotes, by = "otu") %>% 
    left_join(resAntidotes2, by = c("otu", "antidote")) 
  
  resCombined %>% 
    mutate(sample_name = as.character(samples[1,"sample_name"]))
}

resCombined <- samples %>% 
  group_by(sample_name) %>% 
  do(processSubject(.))

padj_cutoff <- 0.1
min_log2_fc_difference <- 1

resCombined <- resCombined %>% 
  filter(!is.na(padj_no_antidote)) %>% 
  mutate(reduced = log2FoldChange_no_antidote > min_log2_fc_difference & 
           padj_no_antidote < padj_cutoff)

resCombined <- resCombined %>% 
  ungroup() %>% 
  group_by(sample_name, reduced) %>% 
  mutate(padj_vs_no_antidote_only_reduced = ifelse(
    reduced, p.adjust(pvalue_vs_no_antidote[reduced], "BH"), NA)) %>% 
  ungroup() %>% 
  mutate(
    interesting = reduced & 
      log2FoldChange_no_antidote - 
      log2FoldChange_with_antidote > min_log2_fc_difference & 
      padj_vs_no_antidote_only_reduced < padj_cutoff,
    reduced_by_antidote = log2FoldChange_with_antidote > 
      min_log2_fc_difference & 
      padj_vs_no_antidote < padj_cutoff,
    outcome = ifelse(
      reduced, ifelse(interesting, 
                      "reduced by erythromycin, rescued by antidote", 
                      "reduced by erythromycin, but not rescued"), 
      ifelse(reduced_by_antidote, 
             "reduced by antidote but not by erythromycin", 
             "not reduced" )))

interesting_OTU <- resCombined %>% 
  group_by(otu) %>% 
  filter(any(interesting)) %>% 
  left_join(otus, by = "otu")

interesting_counts <- counts %>% 
  semi_join(interesting_OTU, by = "otu") %>% 
  left_join(samples %>% select(-count), by = "sample")


d <- resCombined %>% 
  group_by(antidote, otu) %>% 
  filter(any(reduced)) %>% 
  ungroup() %>% 
  left_join(otus) %>%
  mutate(order = fct_rev(
    fct_lump_n(fct_infreq(replace_na(order, "unknown")), 2, 
               other_level = "other"))) %>%
  arrange(desc(as.character(order))) %>% 
  mutate(order = fct_inorder(order)) %>% 
  group_by(antidote, order, sample_name) %>% 
  summarise(n_reduced = sum(reduced), n_rescued = sum(interesting))

p <- d %>% filter(n_reduced > 0) %>% 
  ggplot(aes(sample_name, order)) + 
  geom_point(aes(size = n_reduced, color = n_rescued / n_reduced) ) +
  scale_color_gradient(low = "#dddddd", high = "black", 
                       labels = scales::percent, 
                       name = "Fraction of rescued OTUs") +
  scale_size_area(name = "Number of affected OTUs") +
  facet_wrap(~antidote, ncol = 1) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "darkgrey"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6), 
        legend.text = element_text(size = 5),  
        legend.title = element_text(size = 6), 
        legend.position = "bottom", legend.direction = "horizontal", 
        legend.box = "vertical", legend.box.just = "left")

pp <- grid.arrange(set_panel_size(p, width = unit(6, "cm"), 
                                  height = unit(2, "cm")))
ggsave(here("figures", "Fig3d.pdf"), pp, height = 16, width = 8, units = "cm")

# export source data
if (dir.exists(here('source_data'))){
  d %>% 
    filter(n_reduced > 0) %>% 
    write_tsv(here('source_data', 'Fig3d.tsv'))
}

dd <- d %>% group_by(antidote, order) %>% 
  summarise(Rescued = sum(n_rescued), 
            `Not rescued` = sum(n_reduced) - sum(n_rescued), 
    f = sum(n_rescued)/sum(n_reduced)) %>% 
  mutate(lbl = paste0(" ", round(100*(f),0), "%", 
                      ifelse(Rescued > 10, " rescued", "")))

p <- dd %>% 
  pivot_longer(3:4) %>% 
  ggplot(aes(order, value, fill = name)) + 
  geom_col(width = 0.75) + 
  geom_text(data = dd %>% filter(Rescued >= 5), 
            aes(x = order, y = 0, label = lbl), 
            color = "white", size = 2, hjust = 0, inherit.aes = F) +
  geom_text(data = dd %>% filter(Rescued < 5), 
            aes(x = order, y = Rescued, label = lbl), 
            color = "black", size = 2, hjust = 0, inherit.aes = F) +
  coord_flip() +
  ylab("Number of OTUs") +
  scale_fill_manual(values = c("#dddddd", "black"), name = "") +
  facet_wrap(~antidote, ncol = 1) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 5), 
        axis.title.x = element_text(size = 6), 
        legend.position = "none")

pp <- grid.arrange(set_panel_size(p, 
                                  width = unit(6, "cm"), 
                                  height = unit(2, "cm")))
ggsave(here("figures", "Fig3e.pdf"), pp, height = 16, width = 8, units = "cm")

# export source data
if (dir.exists(here('source_data'))){
  dd %>% 
    select(-lbl) %>% 
    write_tsv(here('source_data', 'Fig3e.tsv'))
}


d <- resCombined %>% 
  group_by(antidote, otu) %>% 
  filter(any(reduced)) %>% 
  ungroup() %>% 
  left_join(otus) %>% 
  filter(order == "Bacteroidales")  %>% 
  mutate(genus = fct_rev(
    fct_lump_n(fct_infreq(replace_na(genus, "unknown")), 2))) %>%
  group_by(antidote, genus) %>% 
  summarise(n_reduced = sum(reduced), n_rescued = sum(interesting)) %>% 
  ungroup() %>% 
  arrange(desc(as.character(genus))) %>% 
  mutate(genus = fct_inorder(genus)) 

dd <- d %>% group_by(antidote, genus) %>% 
  summarise(Rescued = sum(n_rescued), 
            `Not rescued` = sum(n_reduced) - sum(n_rescued), 
    f = sum(n_rescued)/sum(n_reduced)) %>% 
  mutate(lbl = paste0(" ", round(100*(f),0), "%", 
                      ifelse(Rescued > 10, " rescued", "")))

p <- dd %>% 
  pivot_longer(3:4) %>% 
  ggplot(aes(genus, value, fill = name)) + 
  geom_col(width = 0.75) + 
  geom_text(data = dd, aes(x = genus, y = Rescued, label = lbl), 
            color = "black", size = 2, hjust = 0, inherit.aes = F) +
  coord_flip() +
  ylab("Number of OTUs") +
  scale_fill_manual(values = c("#dddddd", "black"), name = "") +
  facet_wrap(~antidote, ncol = 1) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5), 
        axis.text.y = element_text(size = 6), 
        axis.title.x = element_text(size = 6), 
        legend.position = "none")

pp <- grid.arrange(set_panel_size(p, width = unit(3, "cm"), 
                                  height = unit(2, "cm")))

ggsave(here("figures", "Fig3f.pdf"), pp, height = 16, width = 8, units = "cm")

# export source data
if (dir.exists(here('source_data'))){
  dd %>% 
    select(-lbl) %>% 
    write_tsv(here('source_data', 'Fig3f.tsv'))
}