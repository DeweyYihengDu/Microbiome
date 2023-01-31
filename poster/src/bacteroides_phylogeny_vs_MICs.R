# ##############################################################################
#
## Compare phylogenetic and MIC-based similarity in Bacteroides
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
library("vegan")
library("dendextend")


# ##############################################################################
# genome info
genome.info <- read_tsv(here('poster/files','species_annotation.tsv'))

# get MIC info
mic.info <- read_excel(here('poster/data','MICs.xlsx'), sheet = 1)

tree <- read.tree(here('poster/files','iq_tree.nw'))

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

# restrict to BL and Bacteroides
df.plot <- mic.data[genome.info %>% 
                      filter(Genus=='Bacteroides') %>% 
                      pull(NT_code),
                    mic.info %>% 
                      filter(EUCAST_comparison == 'beta-lactams') %>% 
                      pull(Drug_Abbreviation)]
df.plot <- df.plot[,colSums(is.na(df.plot)) == 0]
org.mat <- df.plot
for (d in colnames(df.plot)){
  low <- mic.info %>% filter(Drug_Abbreviation==d) %>% pull(log.low)
  high <- mic.info %>% filter(Drug_Abbreviation==d) %>% pull(log.high)
  df.plot[,d] <- (df.plot[,d] - low)/(high-low)
}
df.plot <- 1 - df.plot

tree <- keep.tip(tree, rownames(df.plot))

rownames(df.plot) <- genome.info$Species_short[match(rownames(df.plot), 
                                                     genome.info$NT_code)]

# ##############################################################################
# plot heatmap
org.mat <- 2**org.mat
h <- Heatmap(df.plot,
             clustering_method_rows = 'ward.D2',
             clustering_method_columns = 'mcquitty',
             col=c('grey90','grey10'), na_col='white',
             row_names_gp = gpar(fontsize = 8),
             name = ' ',
             cell_fun = function(j, i, x, y, width, height, fill) {
               temp <- org.mat[i, j]
               if (temp > 20){
                 string <- sprintf('%.0f', temp)
               } else if (temp > 2){
                 string <- sprintf("%.1f", temp)
               } else {
                 string <- sprintf("%.2f", temp)
               }
               grid.text(string, 
                         x, y, gp = gpar(fontsize = 6, col='white'))
             })

print(h)


# ##############################################################################
# print tree
tree$tip.label <- genome.info$Species_short[match(tree$tip.label, 
                                                  genome.info$NT_code)]

plot(tree)


# ##############################################################################
# pairwise distances

phylo.dist <- as.matrix(cophenetic.phylo(tree))

mic.dist <- as.matrix(vegan::vegdist(log2(org.mat), method='euclidean'))
rownames(mic.dist) <- colnames(mic.dist) <- genome.info$Species_short[
  match(rownames(mic.dist), genome.info$NT_code)]

mic.dist <- mic.dist[rownames(phylo.dist), colnames(phylo.dist)]
diag(mic.dist) <- NA
mic.dist[upper.tri(mic.dist)] <- NA


df.plot <- expand.grid(rownames(mic.dist), colnames(mic.dist)) %>% 
  as_tibble() %>% transmute(Strain_1=as.character(Var1),
                            Strain_2=as.character(Var2))
df.plot <- df.plot %>% 
  mutate(MIC.distance=c(mic.dist), phylo.dist=c(phylo.dist)) %>% 
  filter(!is.na(MIC.distance))

df.plot %>% 
  ggplot(aes(x=phylo.dist, y=1-MIC.distance)) + 
    geom_point(col=alpha("black", .33)) + 
    xlab('Phylogenetic distance') + 
    ylab('MIC similarity') + 
    theme_base(base_size = 10)


# ##############################################################################
# compare trees

mic.dist <- vegdist(log2(org.mat), method='euclidean')
mic.dend <- hclust(mic.dist, method = 'ward.D2')
mic.dend$labels <- genome.info$Species_short[
  match(mic.dend$labels, genome.info$NT_code)]
mic.dend <- as.dendrogram(mic.dend)

tree.dend <- as.dendrogram(as.hclust(chronos(tree)))

dendlist(mic.dend, tree.dend) %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram() 
dendlist(mic.dend, tree.dend) %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  entanglement() 

# ##############################################################################
# plot pairwise distances as double heatmap, see Maier et al. ED 9 b

phylo.clustering <- as.hclust.phylo(chronos(tree))

mat.1 <- phylo.dist
mat.1 <- mat.1[phylo.clustering$labels, rev(phylo.clustering$labels)]


mat.2 <- as.matrix(mic.dist)
rownames(mat.2) <- genome.info$Species_short[
  match(rownames(mat.2), genome.info$NT_code)]
colnames(mat.2) <- genome.info$Species_short[
  match(colnames(mat.2), genome.info$NT_code)]
mat.2 <- mat.2[phylo.clustering$labels, rev(phylo.clustering$labels)]


pdf(here('poster/figures','EDFig4e.pdf'), width = 8, height = 7)
plot(tree)
Heatmap(mat.1, cluster_rows = FALSE, cluster_columns = FALSE, 
        row_names_side = 'left', column_names_side = 'top',
        col = rev(RColorBrewer::brewer.pal(9, 'Purples')),
        name='Phylogenetic\ndistance')
Heatmap(mat.2, cluster_rows = FALSE, cluster_columns = FALSE, 
        row_names_side = 'left', column_names_side = 'top',
        col=rev(RColorBrewer::brewer.pal(9, 'Oranges')),
        name='Abx\ndistance')  
dev.off()

# export source data
if (dir.exists(here('poster/source_data'))){
  mat.1 <- mat.1[,rownames(mat.1)]
  mat.1[lower.tri(mat.1)] <- NA
  diag(mat.1) <- NA
  a <- mat.1 %>% 
    as_tibble(rownames = 'species_1') %>% 
    pivot_longer(-species_1, names_to = 'species_2', 
                 values_to = 'phylo.dist') %>% 
    filter(!is.na(phylo.dist))
  mat.2 <- mat.2[rownames(mat.1), rownames(mat.1)]
  mat.2[lower.tri(mat.2)] <- NA
  diag(mat.2) <- NA
  b <- mat.2 %>% 
    as_tibble(rownames = 'species_1') %>% 
    pivot_longer(-species_1, names_to = 'species_2', 
                 values_to = 'mic.dist') %>% 
    filter(!is.na(mic.dist))
  full_join(a, b) %>% 
    write_tsv(here('poster/source_data', 'EDFig4e.tsv'))
}
