
## Preparation of neighbourhood data for linear mixed effects model ###

## Megan Cole 

## Scripts required to run previous to this:
# -Normalisation, concatenation and scaling of single cell data
# -Rphenograph clustering of normalised and scaled data
# -Febe script of retrieving domain/treatment/metacluster info
# -Preparation for neighbouRhood permutation analysis
# -neighbouRhood permutation analysis 
# neighbouRhood permutation downstream analysis 

# Clear working environment
rm(list = ls())

# Load packages
library(dplyr)

experiment_path = "Neighbourhood_analysis/output/"

output_dir = "Stats_input/"

normal_log2fc = read.csv(
  paste(
    experiment_path,
    "output_stats/neighbouRhood_dat_perm_mean_log2FC_normal.csv",
    sep = ""
  )
)
normal_log2fc$domain = 'Normal'
interface_log2fc = read.csv(
  paste(
    experiment_path,
    "output_stats/neighbouRhood_dat_perm_mean_log2FC_interface.csv",
    sep = ""
  )
)
interface_log2fc$domain = 'Inteface'
tumour_log2fc = read.csv(
  paste(
    experiment_path,
    "output_stats/neighbouRhood_dat_perm_mean_log2FC_tumour.csv",
    sep = ""
  )
)
tumour_log2fc$domain = 'Tumour'

all_log2fc = rbind(normal_log2fc, interface_log2fc, tumour_log2fc)

subset_log2fc = all_log2fc[which(
  all_log2fc$FirstLabel == 'Macrophages type 1' |
    all_log2fc$FirstLabel == 'Macrophages type 2'
), ]

subset_log2fc$MouseID = ""

subset_log2fc[which(subset_log2fc$group == 1), 'MouseID'] = 'BRAC3529.2d_MRTX'
subset_log2fc[which(subset_log2fc$group %in% c(2, 3)), 'MouseID'] = 'BRAC3495.3f_Vehicle'
subset_log2fc[which(subset_log2fc$group == 4), 'MouseID'] = 'BRAC3529.2b_MRTX'
subset_log2fc[which(subset_log2fc$group == 5), 'MouseID'] = 'BRAC3326.4e_Vehicle'
subset_log2fc[which(subset_log2fc$group %in% c(6:8)), 'MouseID'] = 'BRAC3438.6f_Vehicle'
subset_log2fc[which(subset_log2fc$group %in% c(9:12)), 'MouseID'] = 'BRAC4002.3c_MRTX'

subset_log2fc = subset_log2fc %>% select(log2, FirstLabel, SecondLabel, domain, group, MouseID)

write.csv(
  subset_log2fc,
  paste(output_dir, "neighbouRhood_figure4d_lme_input.csv", sep = ""),
  row.names = F
)

#############################################################################################################################################