# ##############################################################################
#
## Show congruence between screen and MIC via ROCs
#
# ##############################################################################

# packages
library("here")
library("tidyverse")
library("readxl")
library("pROC")
library("ggthemes")

# ##############################################################################
# load data
abx.colors.df <- read_delim(here('files','ABX_color_code.csv'),delim=';')
abx.colors <- abx.colors.df$Color
names(abx.colors) <- abx.colors.df$ABX_class
abx.colors['all'] <- 'black'

# abx info
abx.info <- read_csv(here('files','ABX_annotation.csv'))

# mic info
mic.info <- read_excel(here('data','MICs.xlsx'), sheet = 1)

overlap.abx <- intersect(abx.info$prestwick_ID, mic.info$prestwick_ID) 

# screen AUC
screen.data <- read_tsv(here('data','combined_pv.tsv')) %>% 
  filter(prestwick_ID %in% overlap.abx) %>% 
  select(NT_code, prestwick_ID, AUC) %>% 
  spread(key=prestwick_ID, value=AUC) %>% 
  as.data.frame()
rownames(screen.data) <- screen.data$NT_code
screen.data$NT_code <- NULL

# mic data
mic.data <- read_excel(here('data','MICs.xlsx'), sheet = 3) %>% 
  mutate(MIC=exp(rowMeans(as.matrix((select(., R1_value, R2_value) %>% 
                                     mutate_all(log)))))) %>% 
  select(Drug_Abbreviation, NT_code, MIC)

mic.data <- right_join(mic.data, mic.info %>% 
                         select(Drug_Abbreviation, prestwick_ID))
mic.data <- mic.data %>% 
  filter(prestwick_ID %in% overlap.abx) %>% 
  filter(NT_code %in% rownames(screen.data)) %>% 
  select(NT_code, prestwick_ID, MIC) %>%
  spread(key=prestwick_ID, value=MIC) %>% 
  as.data.frame()
rownames(mic.data) <- mic.data$NT_code
mic.data$NT_code <- NULL

screen.data <- screen.data[rownames(mic.data),]

# ##############################################################################
# load info about concentrations

conc.info <- read_excel(here('files','Suppl_Table_1.xlsx')) %>% 
  filter(prestwick_ID %in% overlap.abx)

df.plot <- tibble(prestwick_ID=character(0), auroc=double(0))
df.roc.all <- tibble(screen=double(0), mic=double(0), mic.bin=double(0),
                     prestwick_ID=character(0))

modif=2

# binarize the mic.data
mic.data.bin <- mic.data
for (i in colnames(mic.data.bin)){
  conc <- conc.info %>% filter(prestwick_ID==i) %>% 
    pull(`screen conc. (20 µM as µg/ml)`) %>% 
    round
  mic <- mic.data.bin[,i]
  auc <- screen.data[,i]
  
  df.temp <- tibble(screen=auc, mic=mic) %>% 
    na.omit()
  df.temp <- df.temp %>% 
    mutate(mic.bin=case_when(mic > conc*modif ~ 1, mic < conc/modif ~ 0)) %>% 
    na.omit() %>% 
    mutate(prestwick_ID=i)
  df.roc.all <- df.roc.all %>% 
    bind_rows(df.temp)
  
  if (length(unique(df.temp$mic.bin)) == 1){
    next()
  } else {
    auroc <- roc(predictor=df.temp$screen, response=df.temp$mic.bin,
                 levels = c(0,1), direction='<')$auc
    df.plot <- df.plot %>% 
      add_row(prestwick_ID=i, auroc=as.numeric(auroc))
  }
}

# ##############################################################################
# calculate AUROCs
roc.overall <- roc(predictor=df.roc.all$screen, response=df.roc.all$mic.bin,
                   levels = c(0,1), direction = '<')

df.plot <- tibble(Sensitivity=roc.overall$sensitivities,
                  Specificity=roc.overall$specificities,
                  type='all')

auroc.values <- c('all'=roc.overall$auc)

# add for all classes
df.roc.all <- df.roc.all %>% 
  mutate(Class=abx.info$EUCAST_Comparison[match(prestwick_ID, 
                                                abx.info$prestwick_ID)])

for (cl in unique(df.roc.all$Class)){
  df.temp <- df.roc.all %>% 
    filter(Class==cl)
  if (length(unique(df.temp$mic.bin))==1){
    next()
  } else {
    roc.class <- roc(predictor=df.temp$screen, response=df.temp$mic.bin,
                     levels = c(0,1), direction = '<')
    auroc.values[cl] <- roc.class$auc
    df.plot <- df.plot %>% 
      bind_rows(tibble(Sensitivity=roc.class$sensitivities,
                       Specificity=roc.class$specificities,
                       type=cl))
  }
}


# ##############################################################################
# plot everything

types <- unique(df.plot$type)

g <- df.plot %>% 
  arrange(Sensitivity) %>% 
  mutate(type=factor(type, levels=types)) %>% 
  ggplot(aes(x=Specificity, y=Sensitivity)) + 
  geom_abline(slope = 1, intercept = 1, colour='grey', lty=3) + 
  geom_line(aes(group=type, col=type)) + 
  theme_base() + 
  scale_x_reverse() +
  scale_color_manual(values=abx.colors[types], name='',
                     labels=paste0(types, ' (AUROC: ',
                                   sprintf(fmt='%.2f', 
                                           auroc.values[types]), ')')) +
  theme(legend.position = c(0.65, 0.23),
        legend.background = element_blank())
ggsave(g, filename = here('figures','EDFig3a.pdf'),
       width = 8, height = 8, useDingbats=FALSE)

# ##############################################################################
# get table with FP/FN/TP/TN/N
J <- roc.overall$sensitivities+roc.overall$specificities-1
threshold <- roc.overall$thresholds[which(J==max(J))]



df.table <- tibble(TP=df.roc.all %>% filter(screen > threshold) %>% 
                     pull(mic.bin) %>% sum, 
                   FP=df.roc.all %>% filter(screen > threshold) %>% 
                     pull(mic.bin) %>% -1 %>% abs %>% sum, 
                   FN=df.roc.all %>% filter(screen < threshold) %>% 
                     pull(mic.bin) %>% sum, 
                   TN=df.roc.all %>% filter(screen < threshold) %>% 
                     pull(mic.bin) %>% -1 %>% abs %>% sum, 
                   n=length(unique(df.roc.all$prestwick_ID)),
                   cond.pos=sum(df.roc.all$mic.bin),
                   cond.neg=sum(df.roc.all$mic.bin == 0),
                   type='all')

for (cl in unique(df.roc.all$Class)){
  df.temp <- df.roc.all %>% 
    filter(Class==cl)
  df.table <- df.table %>% 
    add_row(TP=df.temp %>% filter(screen > threshold) %>% 
              pull(mic.bin) %>% sum, 
            FP=df.temp %>% filter(screen > threshold) %>% 
              pull(mic.bin) %>% -1 %>% abs %>% sum, 
            FN=df.temp %>% filter(screen < threshold) %>% 
              pull(mic.bin) %>% sum, 
            TN=df.temp %>% filter(screen < threshold) %>% 
              pull(mic.bin) %>% -1 %>% abs %>% sum, 
            n=length(unique(df.temp$prestwick_ID)),
            cond.pos=sum(df.temp$mic.bin),
            cond.neg=sum(df.temp$mic.bin == 0),
            type=cl)
}


df.table$Sensitivity <- df.table$TP/df.table$cond.pos
df.table$Specificity <- df.table$TN/df.table$cond.neg

df.table <- full_join(df.table, 
                      enframe(auroc.values, name='type', value='auroc'))
df.table

# export source data
if (dir.exists(here('source_data'))){
  df.plot %>% 
    write_tsv(here('source_data', 'EDFig3a.tsv'))
  df.table %>% 
    write_tsv(here('source_data', 'EDFig3a-1.tsv'))
}
