# ##############################################################################
#
## This script finds the abundance and prevalence of microbes in healthy human 
## gut.
#
# ##############################################################################

library("here")

# load data --------------------------------------------------------------------
# THIS NEED TO BE UPDATED
base_p = "https://www.embl.de/download/zeller/milanese/lisa_paper/"
# load files
labels = read.csv(paste0(base_p,"labels.tsv"), sep ="\t",
                  stringsAsFactors = F, head = T,skip = 1)
motus_raw = read.csv(paste0(base_p,"all_files.motus"), sep ="\t",
                  stringsAsFactors = F, head = T,
                  skip = 2, row.names = 1)

# select only healthy ("-1" = "healthy") ---------------------------------------
healthy_patients = colnames(labels)[which(labels==-1)]
motus = motus_raw[,healthy_patients ]

cat(paste0("number of samples: ",length(healthy_patients),"\n") )


# filtering --------------------------------------------------------------------
# remove -1 from mOTUs table
motus = motus[1:7726,]

# calcualte relative abundance
motus = t(t(motus)/colSums(motus))

# put to zero everything that is less than 10^-4
motus[which(motus < 10^-4)] = 0 

# prevalence filter - 5 samples
pos = which (rowSums( motus > 0 ) < 5 & rowSums( motus > 0 ) > 0)
motus[pos,] = 0

# log transform
motus_log = log10(motus+0.00001)


# calculate relative abundance and prevalence ----------------------------------
average_rel_ab = rowMeans(motus)

median_rel_ab = apply(motus, 1, median) 

average_rel_ab_log = rowMeans(motus_log)
average_rel_ab_log = rev(sort(average_rel_ab_log))

median_rel_ab_log = apply(motus_log, 1, median) 

prevalence = rowSums(motus>0) / dim(motus)[2]

# find mOTUS that are different from zero
motus_diff_from_zero = average_rel_ab[which(average_rel_ab > 0)]

# merge into one matrix --------------------------------------------------------
# we use order of avg. log rel ab
to_print = matrix(0,length(motus_diff_from_zero),5)
rownames(to_print) = names(average_rel_ab_log)[1:length(motus_diff_from_zero)]
colnames(to_print) = c("average_log_rel_ab","prevalence",
                       "average_rel_ab","median_rel_ab",
                       "median_log_rel_ab")

to_print[,1] = average_rel_ab_log[rownames(to_print)]
to_print[,2] = prevalence[rownames(to_print)]
to_print[,3] = average_rel_ab[rownames(to_print)]
to_print[,4] = median_rel_ab [rownames(to_print)]
to_print[,5] = median_rel_ab_log[rownames(to_print)]

write.table(to_print, here('files',"most_abundant_bugs.tsv"),
            sep = "\t",quote = F)
