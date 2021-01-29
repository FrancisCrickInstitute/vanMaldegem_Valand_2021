
#### Bodenmiller neighBouRhood analysis method to look at spatial relationships between cells ####
#### compared to a random permutation of the data ####
## Run on neighbour data for each treatment group separately ##

## Megan Cole

## Scripts required to run previous to this:
# -Normalisation, concatenation and scaling of single cell data
# -Rphenograph clustering of normalised and scaled data
# -Febe script of retrieving domain/treatment/metacluster info
# -Preparation for neighbourhood analysis

# Clear working environment
rm(list=ls())

# Install Bodenmiller neighbouRhood package
devtools::install_github("BodenmillerGroup/neighbouRhood")

# Open required packages
library(data.table)
library(dplyr)
library(magrittr)
library(dtplyr)
library(ggplot2)
library(parallel)
library(neighbouRhood)
library(gplots)
library(RColorBrewer)
library(colorspace)

##############################
#### SET GLOBAL VARIABLES ####
##############################

experiment_path = "~path/Neighbourhood_input/"

#treatment = "Vehicle"
treatment = "MRTX"

# Load objectTable and Object Relationships datasets
clusterdata = read.csv(paste(experiment_path, "objectTable/objectTable_", treatment, ".csv", sep = ""))
neighbourdata = read.csv(paste(experiment_path, "object_relationships/Object_relationships_", treatment, ".csv", sep = ""))

n_perm = 5000
ncores = 10

##############################################

# Convert clusterdata to a data table 
clusterdata = data.table(clusterdata)
clusterdata = clusterdata[,-c("X")]
names(clusterdata)

# Rename column names to ensure they match what the function looks for
names(neighbourdata) = gsub(".Object.Name", " Object Name", names(neighbourdata))
names(neighbourdata) = gsub(".Image.Number", " Image Number", names(neighbourdata))
names(neighbourdata) = gsub(".Object.Number", " Object Number", names(neighbourdata))
names(neighbourdata) = gsub("Module.Number", "Module Number", names(neighbourdata))
neighbourdata = neighbourdata[,-1]
names(neighbourdata)
neighbourdata = data.table(neighbourdata)

# Prepare the data - combine clusterdata and neighbourdata into one variable 
data = prepare_tables(clusterdata, neighbourdata, objname = "AllCellsMask")

# Calculate baseline statistics
data_baseline = apply_labels(data[[1]], data[[2]]) %>%
  aggregate_histo()

# Save baseline statistics
write.csv(data_baseline, paste(experiment_path, "output_stats/neighbouRhood_data_baseline_", treatment, ".csv", sep = ""))

# Calculate permutation statistics
set.seed(12312)
dat_perm = rbindlist(mclapply(1:n_perm, function(x){
  dat_labels = shuffle_labels(data[[1]])
  apply_labels(dat_labels, data[[2]]) %>%
    aggregate_histo()
},mc.cores = ncores
), idcol = 'run')

# Save permutation statistics  
write.csv(dat_perm, paste(experiment_path, "output_stats/neighbouRhood_dat_perm_", treatment, ".csv", sep = ""))

# Calculate pvalues
dat_p <- calc_p_vals(data_baseline, dat_perm, n_perm = 5000, p_tresh = 0.01)

# Save p-values
write.csv(dat_p, paste(experiment_path, "output_stats/neighbouRhood_dat_pvalues_", treatment, ".csv", sep = ""))

# Enrichment/depletion score for each neighbour pairing when compared to permutation statistics
pmat = dcast(dat_p, 'FirstLabel ~ SecondLabel', value.var = 'sigval', fun.aggregate = sum,
             fill=0, drop=F) 

# Make FirstLabel names = row names & remove FirstLabel column
rname = pmat$FirstLabel

pmat = pmat %>%
  select(-c('FirstLabel')) %>%
  as.matrix()

row.names(pmat) <- rname

# Select clusters for plotting
select = c("Epithelium", "Fibroblasts", "Lymphocytes", "Myeloid: dendritic cells", "Myeloid: macrophages type 1",
           "Myeloid: macrophages type 2", "Myeloid: neutrophils", "Normal", "Tumour", "Vessels")
pmat = pmat[select, select]

# Save pmat
write.csv(pmat, paste(experiment_path, "output_stats/neighbouRhood_pmat_", treatment, ".csv", sep = ""))


#### Supplemental Figure 6a ####

## Plot heatmap of permutation results and save as pdf ##
pdf(paste(experiment_path, "output_plots/heatmap/neighbouRhood_permutation_histogram_", treatment, ".pdf", sep = ""), width = 9, height = 9)

hr <- hclust(dist(pmat), method="ward.D")
heatmap.2(pmat,
          dendrogram = "none",
          Colv = FALSE,
          Rowv = FALSE,
          trace = "none",
          margins = c(15,15),
          col = diverge_hsv(11),
          density.info ='none',
          keysize = 1,
          breaks = seq(-6,5, 1) #Range will change based on number of images in the analysis
)
dev.off()

#############################################################################################################################################


