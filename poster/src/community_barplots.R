# ##############################################################################
#
## Antidote and Abx in Community 
#
# ##############################################################################

library('here')
library("tidyverse")

# data had been pre-processed with DADA2
# ASVs were then assigned to species-level by blasting against reference 16S
# sequences
# Sequences were assigned to a given species if the best hit matched 
# with more than 97% identity
df.plot <- read_tsv(here('poster/data', 'iv_community.tsv'))

df.plot.sum <- df.plot %>% 
  filter(Time=='24h') %>% 
  group_by(plot_group, species) %>% 
  summarise(m=median(abundance), s=mad(abundance),
            .groups = 'drop') %>% 
  group_by(plot_group) %>%  
  mutate(m.rel=m/sum(m)) %>% 
  filter(species!='Other')

sp.order <- df.plot %>% 
  filter(plot_group=='no_abx-no_antidote') %>% 
  group_by(species) %>% 
  summarise(m=median(abundance), .groups = 'drop') %>% 
  arrange(desc(m))

atds <- setdiff(unique(df.plot$Antidote), 'none')
for (atd in atds){
  g <- df.plot %>% 
    filter(Time!='11h') %>% 
    # filter(Antidote==atd) %>% 
    mutate(Replicate=paste0(Row, Column)) %>% 
    filter(species!='Other') %>% 
    mutate(species=factor(species, levels = sp.order$species)) %>% 
    mutate(plot_group=factor(plot_group, 
                             levels = c('no_abx-no_antidote',
                                        paste0('no_abx-', atd),
                                        'abx-no_antidote',
                                        paste0('abx-', atd)))) %>% 
    filter(!is.na(plot_group)) %>% 
    ggplot(aes(x=Replicate, y=abundance, fill=species)) + 
    geom_bar(stat='identity') + 
    facet_grid(~plot_group, scales = 'free', space = 'free') + 
    theme_bw() + 
    ylab("Abundance") + 
    xlab("") + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank()) + 
    scale_fill_manual(values=c( 'ivory3', 'ivory4', 'slategray3',
                                'slategrey',  'azure3', 'dimgray',
                                'slategray1')) +
    ggtitle(atd) +
    NULL
  assign(atd, g)
}

# save
pdf(here('poster/figures',"Fig3b.pdf"),width=5,height=4) 
dicoumarol
dev.off()


pdf(here('poster/figures',"EDFig10b.pdf"),width=5,height=4) 
benzbromarone
dev.off()

pdf(here('poster/figures',"EDFig11c.pdf"),width=5,height=4) 
`tolfenamic acid`
dev.off()

# export source data
if (dir.exists(here('poster/source_data'))){
  # Fig3b
  df.plot %>% 
    filter(Time!='11h') %>% 
    filter(group_2 %in% c('no_antidote', 'dicoumarol')) %>% 
    mutate(Replicate=paste0(Row, Column)) %>% 
    select(species, abundance, group, group_2, Replicate) %>% 
    write_tsv(here('poster/source_data', 'Fig3b.tsv'))
  
  df.plot %>% 
    filter(Time!='11h') %>% 
    filter(group_2 %in% c('no_antidote', 'benzbromarone')) %>% 
    mutate(Replicate=paste0(Row, Column)) %>% 
    select(species, abundance, group, group_2, Replicate) %>% 
    write_tsv(here('poster/source_data', 'EDFig10b.tsv'))
  
  df.plot %>% 
    filter(Time!='11h') %>% 
    filter(group_2 %in% c('no_antidote', 'tolfenamic acid')) %>% 
    mutate(Replicate=paste0(Row, Column)) %>% 
    select(species, abundance, group, group_2, Replicate) %>% 
    write_tsv(here('poster/source_data', 'EDFig11c.tsv'))
}
