# ##############################################################################
#
## Animal experiments
#
# ##############################################################################

library("here")
library("tidyverse")
library("ggpubr")
library("scales")

Raw_data_animal_experiment <- read_tsv(
  here('poster/data','Raw_data_animal_experiments.tsv'))


# analysis for dicumarol
annotated_data_dicum <- Raw_data_animal_experiment %>% 
  filter(Treatment %in% c('erythromycin', 'erythromycin+dicoumarol'))

annotated_data_dicum_long <- annotated_data_dicum %>% 
  mutate(time = as.factor(Time)) %>% 
  filter(time != 0.375) %>%
  filter(time !=0.25)

pdf(here('poster/figures',"Fig3g.pdf"),width=8.27,height=11.69)

ggplot(data=annotated_data_dicum_long, 
       aes(x=time, y=copies_g_feces,fill=Treatment, color=Treatment)) +
  geom_boxplot(outlier.shape=NA, lwd=0.75, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), 
               dotsize=0.5, aes(color = Treatment)) +
  theme_linedraw(base_size = 12) +
  scale_y_continuous(trans='log10')+
  # add labels
  labs(title='Dicoumarol',y = "B. vulgatus copies per g feces", 
       x = 'day post treatment') +
  # change font size
  theme(title = element_text(face = "bold"), 
        axis.title = element_text(size = 12)) +
  # change axis label
  theme(legend.title = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y= element_text(size = 12),
        axis.text.x= element_text(size = 12)) +
  scale_fill_manual(values=alpha(c('#FF9933', 'grey'),0.5)) +
  scale_color_manual(values=c('#FF9933', 'grey')) +
  coord_fixed(ratio=1)

dev.off()


annotated_data_dicum_short <- annotated_data_dicum %>% 
  mutate(time = as.factor(Time)) %>% 
  filter(Time < 3)  %>%  
  filter(!is.na(Ery_Conc_nmol_g))

df.summary2 <- annotated_data_dicum_short %>%
  group_by(Treatment, Time) %>%
  summarise(
    sd = sd(Ery_Conc_nmol_g),
    Ery_Conc_nmol_g = mean(Ery_Conc_nmol_g))

pdf(here('poster/figures',"Fig3h.pdf"),width=8.27,height=11.69)

ggplot(annotated_data_dicum_short, aes(Time, Ery_Conc_nmol_g)) +
  geom_jitter(aes(color = Treatment),position = position_jitter(0)) + 
  geom_line(aes(color = Treatment), data = df.summary2, size=1.5) +
  theme_linedraw(base_size = 12) +
  # add labels
  labs(title='Fecal erythromycin concentration',
       y = "erythromycin [nmol/g]", x = 'day post treatment') +
  # change font size
  theme(title = element_text(face = "bold"), 
        axis.title = element_text(size = 12)) +
  # change axis label
  theme(legend.title = element_blank(), aspect.ratio = 1, 
        axis.text.x= element_text(size = 12),  
        axis.text.y= element_text(size = 12)) +
  scale_color_manual(values=c("#FF9933", "grey"))

dev.off()

# analysis for benzbromarone
annotated_data_benz <- Raw_data_animal_experiment %>% 
  filter(Treatment %in% c('erythromycin', 'erythromycin+benzbromaron')) %>% 
  filter(Exp_ID != 'TV228')

annotated_data_benz_long <- annotated_data_benz %>% 
  mutate(time = as.factor(Time)) %>% 
  filter(time != 0.375) %>% 
  filter(time !=0.25)

pdf(here('poster/figures',"EDFig10g.pdf"),width=8.27,height=11.69)

ggplot(data=annotated_data_benz_long, 
       aes(x=time, y=copies_g_feces,fill=Treatment, color=Treatment)) +
  geom_boxplot(outlier.shape=NA, lwd=0.75, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), 
               dotsize=0.5, aes(color = Treatment)) +
  theme_linedraw(base_size = 12) +
  scale_y_continuous(trans='log10')+
  # add labels
  labs(title='Benzbromarone',y = "B. vulgatus copies per g feces", 
       x = 'day post treatment') +
  # change font size
  theme(title = element_text(face = "bold"), 
        axis.title = element_text(size = 12)) +
  # change axis label
  theme(legend.title = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y= element_text(size = 12),
        axis.text.x= element_text(size = 12)) +
  scale_fill_manual(values=alpha(c('#FF9933', 'grey'),0.5)) +
  scale_color_manual(values=c('#FF9933', 'grey')) +
  coord_fixed(ratio=1)

dev.off()


annotated_data_benz_short <- annotated_data_benz %>% 
  mutate(time = as.factor(Time)) %>% 
  filter(Time < 3)  %>%  
  filter(!is.na(Ery_Conc_nmol_g))

df.summary2_benz <- annotated_data_benz_short %>%
  group_by(Treatment, Time) %>%
  summarise(
    sd = sd(Ery_Conc_nmol_g),
    Ery_Conc_nmol_g = mean(Ery_Conc_nmol_g))

pdf(here('poster/figures',"EDFig10h.pdf"),width=8.27,height=11.69)

ggplot(annotated_data_benz_short, aes(Time, Ery_Conc_nmol_g)) +
  geom_jitter(aes(color = Treatment),position = position_jitter(0)) + 
  geom_line(aes(color = Treatment), data = df.summary2_benz, size=1.5) +
  theme_linedraw(base_size = 12) +
  # add labels
  labs(title='Fecal erythromycin concentration',
       y = "erythromycin [nmol/g]", x = 'day post treatment') +
  # change font size
  theme(title = element_text(face = "bold"), 
        axis.title = element_text(size = 12)) +
  # change axis label
  theme(legend.title = element_blank(), aspect.ratio = 1, 
        axis.text.x= element_text(size = 12),  
        axis.text.y= element_text(size = 12)) +
  scale_color_manual(values=c("#FF9933", "grey"))

dev.off()



# Mann-Whitney U test

# Ery+Dicum: all time points

TPs <- unique(annotated_data_dicum$Time)

statistics <- data.frame(matrix(ncol = 3, nrow =length(TPs))) 
x <- c("Time", "p_value_Bvulgatus", "p_value_Ery_Conc") 
colnames(statistics) <- x

for (i in seq_along(TPs)){
  a <- annotated_data_dicum %>% filter(Time == TPs[i]) 
  b <- wilcox.test(a$copies_g_feces ~ a$Treatment, exact=FALSE,
                   correct=FALSE, alternative = "two.sided", paired = FALSE)
  c <- try(wilcox.test(a$Ery_Conc_nmol_g ~ a$Treatment, exact=FALSE,
                       correct=FALSE, alternative = "two.sided", 
                       paired = FALSE))
  statistics$Time[i] <- i
  statistics$p_value_Bvulgatus[i] <- b$p.value
  statistics$p_value_Ery_Conc[i] <- c$p.value
}


# Ery+Benz: all time points

TPs <- unique(annotated_data_benz$Time)

statistics2 <- data.frame(matrix(ncol = 3, nrow =length(TPs))) 
x <- c("Time", "p_value_Bvulgatus", "p_value_Ery_Conc") 
colnames(statistics2) <- x

for (i in seq_along(TPs)){
  a <- annotated_data_benz %>% filter(Time == TPs[i]) 
  b <- wilcox.test(a$copies_g_feces ~ a$Treatment, exact=FALSE,
                   correct=FALSE, alternative = "two.sided", paired = FALSE)
  c <- try(wilcox.test(a$Ery_Conc_nmol_g ~ a$Treatment, exact=FALSE, 
                       correct=FALSE, alternative = "two.sided", 
                       paired = FALSE))
  statistics2$Time[i] <- i
  statistics2$p_value_Bvulgatus[i] <- b$p.value
  statistics2$p_value_Ery_Conc[i] <- c$p.value
}

# export source data
if (dir.exists(here('poster/source_data'))){
  # Fig3g
  annotated_data_dicum_long %>% 
    select(Mouse_ID, Time, copies_g_feces, Treatment, Exp_ID) %>% 
    write_tsv(here('poster/source_data', 'Fig3g.tsv'))
  # Fig3h
  annotated_data_dicum_short %>% 
    select(Mouse_ID, Time, Ery_Conc_nmol_g, Treatment) %>% 
    write_tsv(here('poster/source_data', 'Fig3h.tsv'))
  
  
  annotated_data_benz_long %>% 
    select(Mouse_ID, Time, copies_g_feces, Treatment, Exp_ID) %>% 
    write_tsv(here('poster/source_data', 'EDFig10g.tsv'))
  
  annotated_data_benz_short %>% 
    select(Mouse_ID, Time, Ery_Conc_nmol_g, Treatment) %>% 
    write_tsv(here('poster/source_data', 'EDFig10h.tsv'))
}

