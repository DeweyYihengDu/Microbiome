# ##############################################################################
#
## Screen Heatmap for the supplement
#
# ##############################################################################

# packages
library("here")
library("tidyverse")
library("readxl")
library("ape")
library("ComplexHeatmap")
library("circlize")
library("ggthemes")


# ##############################################################################
# data

# abx info
abx.info <- read_csv(here('poster/files', 'ABX_annotation.csv')) 

# genome info
genome.info <- read_tsv(here('poster/files','species_annotation.tsv'))

# screen data
screen.data <- read_tsv(here('poster/data', 'combined_pv.tsv')) %>% 
  filter(prestwick_ID %in% abx.info$prestwick_ID) %>% 
  select(NT_code, prestwick_ID, hit) %>% 
  mutate(hit=as.numeric(hit)) %>% 
  spread(key=prestwick_ID, value=hit) %>% 
  as.data.frame()
rownames(screen.data) <- screen.data$NT_code
screen.data$NT_code <- NULL

# read in tree
tree <- read.tree(here('poster/files', 'phylogenetic_tree.tre'))

# colors
#abx.colors.df <- read_delim(here('poster/files','ABX_color_code.csv'), delim=';')
#abx.colors <- abx.colors.df$Color
#names(abx.colors) <- abx.colors.df$ABX_class

# ##############################################################################
# hetamap for all classes

# order by hits per class
abx.class.sorted <- abx.info %>% 
  mutate(no_hits=colSums(screen.data[, prestwick_ID], na.rm = TRUE)) %>% 
  group_by(EUCAST_Comparison) %>% 
  summarise(m=median(no_hits)) %>% 
  arrange(desc(m))

abx.info <- abx.info %>% 
  mutate(EUCAST_Comparison = 
           factor(EUCAST_Comparison, 
                  levels = abx.class.sorted$EUCAST_Comparison)) %>% 
  arrange(EUCAST_Comparison)


abx.info %>% 
  mutate(no_hits=colSums(screen.data[, prestwick_ID], na.rm = TRUE)) %>% 
  ggplot(aes(x=EUCAST_Comparison, y=no_hits, fill=EUCAST_Comparison)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1) +
    #scale_fill_manual(values=abx.colors, guide='none') + 
    coord_flip() + 
    xlab('') + 
    ylab('Number of hits') + 
    theme_few() + 
    theme(panel.grid.major.x = element_line(colour='lightgrey'))


df.plot <- t(screen.data[,abx.info$prestwick_ID])
rownames(df.plot) <- abx.info$chemical_name[match(rownames(df.plot), 
                                                  abx.info$prestwick_ID)]

# annotations
df.row <- data.frame(Group=abx.info$EUCAST_Comparison)
side.annot <- rowAnnotation(df=df.row, 
                            show_legend=FALSE)
#side.annot <- rowAnnotation(df=df.row, col=list(Group=abx.colors),show_legend=FALSE)
# order by phylogeny
species.clustering <- as.hclust(chronos(tree))
df.plot <- df.plot[,species.clustering$labels]
colnames(df.plot) <- genome.info$Species_short[match(colnames(df.plot), 
                                                     genome.info$NT_code)]


pdf(here('poster/figures','EDFig1.pdf'), 
    width = 8.27, height = 11.7, useDingbats = FALSE)
h <- Heatmap(df.plot,
             cluster_columns = species.clustering,
             clustering_method_rows = 'ward.D2',
             clustering_distance_rows = 'binary',
             col=c('grey90','grey10'), na_col='white',
             split=abx.info$EUCAST_Comparison,
             heatmap_legend_param = list(title='Abx Sensitivity', 
                                         color_bar='discrete', 
                                         labels=c('resistant', 'sensitive'), 
                                         at=c(0,1)),
             show_row_dend = FALSE,
             row_names_gp = gpar(fontsize = 8))

side.annot + h
dev.off()

# source data
if (dir.exists(here('poster/source_data'))){
  df.plot %>% 
    as_tibble(rownames = 'Drug') %>% 
    pivot_longer(-Drug, names_to = 'strain', values_to = 'screen.hit') %>% 
    left_join(abx.info %>% transmute(Drug=chemical_name, EUCAST_Comparison)) %>% 
    write_tsv(here('poster/source_data', 'EDFig1.tsv'))
}
# ##############################################################################
# same heatmap for MIC


# Revisions
# show all strains in the heatmap
# strains
sp.hm <- c("NT5001", "NT5002", "NT5003", "NT5004", "NT5009", "NT5011",
           "NT5019", "NT5025", "NT5026", "NT5032", "NT5033", "NT5050",
           "NT5054", "NT5078", "NT5083", "NT5065", "NT5057", "NT5049", 
           "NT5051", "NT5052", "NT5053", "NT5055", "NT5056", "NT5058", 
           "NT5059", "NT5066", "NT5064")

# get MIC info
mic.info <- read_excel(here('poster/data','MICs.xlsx'), sheet = 1)

# get MIC data
mic.data <- read_excel(here('poster/data','MICs.xlsx'), sheet = 3) %>% 
  mutate(MIC=exp(rowMeans(as.matrix((select(., R1_value, R2_value) %>% 
                                       mutate_all(log)))))) %>% 
  select(Drug_Abbreviation, NT_code, MIC) %>%
  mutate(MIC=log2(MIC)) %>% 
  spread(key=Drug_Abbreviation, value=MIC) %>% 
  as.data.frame()
rownames(mic.data) <- mic.data$NT_code
mic.data$NT_code <- NULL

# normalize MIC data
mic.info <- mic.info %>% 
  separate(Range, sep='\\s*- ', into = c('low', 'high'), 
           remove=FALSE, convert=TRUE) %>% 
  mutate(log.low=log2(low)) %>% 
  mutate(log.high=log2(high))

# re-order by EUCast Class
mic.info <- mic.info %>% 
  mutate(EUCAST_comparison = 
           factor(EUCAST_comparison, 
                  levels = abx.class.sorted$EUCAST_Comparison)) %>% 
  arrange(EUCAST_comparison)

df.plot <- mic.data
for (d in colnames(df.plot)){
  low <- mic.info %>% filter(Drug_Abbreviation==d) %>% pull(log.low)
  high <- mic.info %>% filter(Drug_Abbreviation==d) %>% pull(log.high)
  df.plot[,d] <- (df.plot[,d] - low)/(high-low)
}
# df.plot.mic <- df.plot[rowSums(is.na(df.plot)) == 0,]
df.plot.mic <- 1 - df.plot
colnames(df.plot.mic) <- mic.info$Drug_Name[match(colnames(df.plot.mic), 
                                                  mic.info$Drug_Abbreviation)]

# take strains
df.plot.mic <- df.plot.mic[sp.hm,]

# re-order by phylogeny
tree_2 <- read.tree(here('poster/files','iq_tree.nw'))
tree_2 <- keep.tip(tree_2, rownames(df.plot.mic))
species.clustering.2 <- as.hclust(chronos(tree_2))
df.plot.mic <- t(df.plot.mic[species.clustering.2$labels,])
colnames(df.plot.mic) <- genome.info$Species_short[
  match(colnames(df.plot.mic), genome.info$NT_code)]

# re-order by EUCAST class
df.plot.mic <- df.plot.mic[mic.info$Drug_Name,]

# prepare heatmap
df.row <- data.frame(Group=mic.info$EUCAST_comparison)
side.annot <- rowAnnotation(df=df.row,
                            show_legend=FALSE)
#side.annot <- rowAnnotation(df=df.row, col=list(Group=abx.colors),show_legend=FALSE)
# prepare matrix for cell annotation
small_mat <- 2**mic.data
rownames(small_mat) <- genome.info$Species_short[match(rownames(small_mat), 
                                                       genome.info$NT_code)]
colnames(small_mat) <- mic.info$Drug_Name[match(colnames(small_mat),
                                                mic.info$Drug_Abbreviation)]
small_mat <- t(small_mat[colnames(df.plot.mic), rownames(df.plot.mic)])


# remove negative values
df.plot.mic[df.plot.mic < 0] <- 0


# plot heatmap
h.1 <- Heatmap(df.plot.mic,
               cluster_columns = species.clustering.2,
               cluster_rows = FALSE,
               col=c('grey90','grey10'), na_col='white',
               split=mic.info$EUCAST_comparison,
               show_row_dend = FALSE,
               row_names_side = 'right',
               row_names_gp = gpar(fontsize = 8),
               name = 'MIC',
               cell_fun = function(j, i, x, y, width, height, fill) {
                 temp <- small_mat[i, j]
                 if (is.na(temp)) {
                   string <- ''
                  } else if (temp > 20){
                   string <- sprintf('%.0f', temp)
                  } else if (temp > 2){
                    string <- sprintf("%.1f", temp)
                  } else {
                    string <- sprintf("%.2f", temp)
                  }
                 grid.text(string, 
                           x, y, gp = gpar(fontsize = 6, col='white'))
                 })

# save heatmap
pdf(here('poster/figures','EDFig2.pdf'), 
    width = 8.27, height = 11.7, useDingbats = FALSE)
side.annot + h.1
dev.off()

# export source data
if (dir.exists(here('poster/source_data'))){
df.plot.mic %>% 
  as_tibble(rownames = 'Drug') %>% 
  pivot_longer(-Drug, names_to = 'strain', values_to = 'relative.MIC') %>% 
  left_join(mic.info %>% transmute(Drug=Drug_Name, EUCAST_comparison)) %>% 
  left_join(genome.info %>% transmute(strain=Species_short, NT_code)) %>% 
  full_join(
    mic.data %>% as_tibble(rownames = 'NT_code') %>% 
      pivot_longer(-NT_code, names_to='Drug_abbreviation', 
                   values_to = 'real.MIC') %>% 
      left_join(mic.info %>% transmute(Drug_abbreviation=Drug_Abbreviation, 
                                       Drug=Drug_Name)) %>% 
      select(real.MIC, Drug, NT_code)) %>% 
  transmute(Drug, strain, relative.MIC, real.MIC, Group=EUCAST_comparison) %>% 
  write_tsv(here('poster/source_data', 'EDFig2.tsv'))
}
