# ##############################################################################
#
## B.vulgatus survival with antidotes
#
# ##############################################################################

library('here')
library("tidyverse")
library("ggpubr")
library("scales")


Raw_data_survival <- read_tsv(here('data','survival_antidote_raw_data.tsv'))


Ery <- ggplot(data=filter(Raw_data_survival, ABX=='Ery'), 
              aes(x=Label, y=X..survival)) +
  geom_boxplot(aes(color = ABX, fill=ABX),outlier.shape=NA, lwd=1) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, 
               fill='#ff9933', color='#ff9933') +
  # add labels
  labs(title = "B. vulgatus [5 x MIC]", y = "% survival", x = NULL) +
  theme_linedraw(base_size = 12) +
  geom_hline(yintercept=100, linetype="dashed", color = "#bdbdbd")+
  geom_hline(yintercept=0.01, linetype="dashed", color = "#bdbdbd")+
  scale_y_continuous(trans=log10_trans(), 
                     labels=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000), 
                     limits = c(0.00001,1000), 
                     breaks =c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000))+
  # change font size
  theme(title = element_text(face = "bold"),legend.position = "none", 
        axis.title = element_text(face = "bold", size = 12)) +
  # change axis label
  theme(axis.text.x = element_text(face = "bold",angle=45, vjust=0.5, size=10), 
        axis.text.y = element_text(size=10), aspect.ratio = 1, 
        legend.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values=c('#ff9933')) +
  scale_fill_manual(values=alpha(c('#ff9933'),0.5))


Doxy <- ggplot(data=filter(Raw_data_survival, ABX=='Doxy'), 
               aes(x=reorder(Label, X..survival, median), y=X..survival)) +
  geom_boxplot(aes(color = ABX, fill=ABX),outlier.shape=NA, lwd=1) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill="#ff3333", 
               color="#ff3333") +
  # add labels
  labs(title = "B. vulgatus [5 x MIC]", y = "% survival", x = NULL) +
  theme_linedraw(base_size = 12) +
  geom_hline(yintercept=100, linetype="dashed", color = "#bdbdbd")+
  geom_hline(yintercept=0.01, linetype="dashed", color = "#bdbdbd")+
  scale_y_continuous(trans=log10_trans(), 
                     labels=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000), 
                     limits = c(0.00001,1000), 
                     breaks =c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000))+
  # change font size
  theme(title = element_text(face = "bold"),legend.position = "none", 
        axis.title = element_text(face = "bold", size = 12)) +
  # change axis label
  theme(axis.text.x = element_text(face = "bold",angle=45, vjust=0.5, size=10), 
        axis.text.y = element_text(size=10), aspect.ratio = 1, 
        legend.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values=c("#ff3333")) +
  scale_fill_manual(values=alpha(c("#ff3333"),0.5))

pdf(here('figures',"EDFig7e.pdf"),width=8.27,height=11.69)
ggarrange(Ery, Doxy, ncol = 2)
dev.off()

# export source data
if (dir.exists(here('source_data'))){
  Raw_data_survival %>% 
    transmute(ABX, Drug, survival=X..survival, Replicate) %>% 
    write_tsv(here('source_data', 'EDFig7e.tsv'))
}

