# ##############################################################################
#
## PCA of screen
#
# ##############################################################################

# packages
library("labdsv")
library("tidyverse")
library("readxl")
library("ggthemes")
library("here")

# ##############################################################################
# data

# abx info
abx.info <- read_csv(here('poster/files','ABX_annotation.csv'))

# genome info
genome.info <- read_tsv(here('poster/files', 'species_annotation.tsv'))

# screen data
screen.data <- read_tsv(here('poster/data', 'combined_pv.tsv')) %>% 
  filter(prestwick_ID %in% abx.info$prestwick_ID)
  
hit.mat <- screen.data %>% 
  select(NT_code, prestwick_ID, hit) %>% 
  mutate(hit=as.numeric(hit)) %>% 
  spread(key=prestwick_ID, value=hit) %>% 
  as.data.frame()
rownames(hit.mat) <- hit.mat$NT_code
hit.mat$NT_code <- NULL

auc.mat <- screen.data %>% 
  select(NT_code, prestwick_ID, AUC) %>% 
  spread(key=prestwick_ID, value=AUC) %>% 
  as.data.frame()
rownames(auc.mat) <- auc.mat$NT_code
auc.mat$NT_code <- NULL

# colors
abx.colors.df <- read_delim(here('poster/files', 'ABX_color_code.csv'), delim=';')
abx.colors <- abx.colors.df$Color
names(abx.colors) <- abx.colors.df$ABX_class

### REVISIONS!
### Remove Sulfonamide antibiotics
auc.mat <- auc.mat[,!(colnames(auc.mat) %in% (abx.info %>% 
                        filter(EUCAST_Comparison == 'Sulfonamide') %>% 
                        pull(prestwick_ID)))]

# ##############################################################################
# compute PCA on AUCs

auc.mat[is.na(auc.mat)] <- 1
pca.results <- pca(auc.mat, cor=TRUE)

df.plot <- tibble(NT_code=rownames(pca.results$scores), 
                  PC1=pca.results$scores[,'PC1'],
                  PC2=pca.results$scores[,'PC2'])
df.plot <- df.plot %>% 
  mutate(Phylum=genome.info$Phylum[match(NT_code, genome.info$NT_code)])

df.plot$PC1 <-  (df.plot$PC1 - mean(df.plot$PC1))/sd(df.plot$PC1)
df.plot$PC2 <-  (df.plot$PC2 - mean(df.plot$PC2))/sd(df.plot$PC2)

var.explained <- (pca.results$sdev**2/pca.results$totdev)[1:2]


df.loadings <- tibble(prestwick_ID=rownames(pca.results$loadings),
                      PC1=pca.results$loadings[,'PC1'],
                      PC2=pca.results$loadings[,'PC2'])
df.loadings <- df.loadings %>% 
  mutate(Class=abx.info$EUCAST_Comparison[match(prestwick_ID, 
                                                abx.info$prestwick_ID)])
df.loadings$PC1 <-  (df.loadings$PC1 - 
                       mean(df.loadings$PC1))/sd(df.loadings$PC1)
df.loadings$PC2 <-  (df.loadings$PC2 - 
                       mean(df.loadings$PC2))/sd(df.loadings$PC2)

df.loadings <- df.loadings %>% 
  group_by(Class) %>% 
  summarise(PC1=mean(PC1), PC2=mean(PC2), n=n()) %>% 
  filter(n > 4)

g <- df.plot %>% 
  ggplot(aes(x=PC1, y=PC2)) + 
    geom_point(aes(shape=Phylum)) + 
    scale_shape_manual(values=c(5, 0, 6, 4, 21, 8)) +
    theme_base(base_size = 10) + 
    geom_segment(data=df.loadings, aes(x=0, y=0, 
                                       xend=PC1, yend=PC2, col=Class), 
                 inherit.aes = FALSE) + 
    scale_color_manual(values=abx.colors) + 
    xlab(paste0('PC1 [', sprintf(fmt='%.2f', var.explained[1]*100),'%]')) +
    ylab(paste0('PC2 [', sprintf(fmt='%.2f', var.explained[2]*100),'%]')) + 
    theme(aspect.ratio = var.explained[2]/var.explained[1],
          plot.background = element_blank())
ggsave(g, filename = here('poster/figures', 'pca.pdf'),
       width = 10, height = 4, useDingbats=FALSE)

g2 <- df.plot %>% 
  filter(!Phylum%in%c('Fusobacteria', 'Verrucomicrobia')) %>% 
  ggplot(aes(x=Phylum, y=PC1)) + 
    geom_boxplot(col='darkgrey', fill=NA, outlier.shape = NA) +
    geom_jitter(aes(shape=Phylum), width = 0.075) + 
    coord_flip() + 
    theme_base(base_size = 10) + 
    scale_shape_manual(values=c(5, 0, 6, 21), guide='none')
g2
