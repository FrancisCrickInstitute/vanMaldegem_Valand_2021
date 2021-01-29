
#### Downstream analysis of neighbourhood permutation test ####

## Megan Cole 

## Scripts required to run previous to this:
# -Normalisation, concatenation and scaling of single cell data
# -Rphenograph clustering of normalised and scaled data
# -Preparation for neighbourhood analysis
# -NeighbouRhood permutation analysis

# Clear working environment
rm(list=ls())

# Load required packages 
library(dplyr)
library(ggplot2)

############################
##### GLOBAL VARIABLES #####
############################
experiment_path = "~path/Neighbourhood_input/"
output_path = experiment_path

pmat_Cont = read.csv(paste(experiment_path, "output_stats/neighbouRhood_pmat_Vehicle.csv", sep = ""))
pmat_MRTX = read.csv(paste(experiment_path, "output_stats/neighbouRhood_pmat_MRTX.csv", sep = ""))

###############################################################################################################

#### Supplemental Figure 6b ####

#############################################################
##### NEIGHBOUR ANALYSIS BASED ON PERMUTATION STATISTIS ##### 
#############################################################

#### Preparing pmat ####
rname = pmat_Cont$X

pmat_Cont = pmat_Cont %>%
  select(-c('X')) %>%
  as.matrix()
pmat_MRTX = pmat_MRTX %>%
  select(-c('X')) %>%
  as.matrix()

row.names(pmat_Cont) <- rname
colnames(pmat_Cont) <- rname
row.names(pmat_MRTX) <- rname
colnames(pmat_MRTX) <- rname

# Transpose the data frame so that cells of interest are cols and neighbours are rows
pmat_t_Cont = as.data.frame(t(pmat_Cont))
pmat_t_MRTX = as.data.frame(t(pmat_MRTX))

# Add cell type names column
pmat_t_Cont$cellType = colnames(pmat_t_Cont)
pmat_t_MRTX$cellType = colnames(pmat_t_MRTX)

# Add a treatment group column
pmat_t_Cont$treatment = "Vehicle"
pmat_t_MRTX$treatment = "MRTX"

# Combine Cont and MRTX pmat info 
pmat = rbind(pmat_t_Cont, pmat_t_MRTX)

# Create a list of bar charts for neighbour interaction comparison across treatment groups, one plot for each cell type 
p <- list()

for(i in 1:length(unique(pmat$cellType))){
  
  print(unique(pmat$cellType)[i])
  
  # Remove the spaces in metacluster names
  cellType_fix = gsub(" ", "_", unique(pmat$cellType)[i]) 
  
  p[[i]] = local({
    i=i
    # Create bar plot 
    ggplot(pmat,
           aes(x = cellType, y = pmat[,i], fill = factor(treatment, levels = c("Vehicle", "MRTX")), width = 0.75))+
      geom_bar(stat = "identity", position = "dodge", colour = "black")+ 
      labs(x = "Metacluster", y = "Enrichment score",
           title = paste(unique(pmat$cellType)[i]),
           fill = "treatment")+
      scale_y_continuous(breaks = seq(-6,6, by = 2))+
      theme_minimal()+
      theme(plot.margin = margin(4,4,4,4,"cm"))+
      theme(axis.text.x = element_text(angle = 45,size =10, hjust = 1), plot.title = element_text(size = 14), axis.title=element_text(size=12)) 
  })
}

## Plotting from the list of bar charts ## 
pdf(paste(experiment_path,"output_plots/bar/Epithelium_neighbours_barchart", ".pdf", sep = ""), width = 9, height = 7)
p[1]
dev.off()

pdf(paste(experiment_path,"output_plots/bar/Fibroblasts_neighbours_barchart", ".pdf", sep = ""), width = 9, height = 7)
p[2]
dev.off()

pdf(paste(experiment_path,"output_plots/bar/Lymphocytes_metacluster_neighbours_barchart", ".pdf", sep = ""), width = 9, height = 7)
p[3]
dev.off()

pdf(paste(experiment_path,"output_plots/bar/Myeloid_dendritic_cells_metacluster_neighbours_barchart", ".pdf", sep = ""), width = 9, height = 7)
p[4]
dev.off()

pdf(paste(experiment_path,"output_plots/bar/Myeloid_macrophages_type1_metacluster_neighbours_barchart", ".pdf", sep = ""), width = 9, height = 7)
p[5]
dev.off()

pdf(paste(experiment_path,"output_plots/bar/Myeloid_macrophages_type2_metacluster_neighbours_barchart", ".pdf", sep = ""), width = 9, height = 7)
p[6]
dev.off()

pdf(paste(experiment_path,"output_plots/bar/Myeloid_neutrophils_metacluster_neighbours_barchart", ".pdf", sep = ""), width = 9, height = 7)
p[7]
dev.off()

pdf(paste(experiment_path,"output_plots/bar/Normal_metacluster_neighbours_barchart", ".pdf", sep = ""), width = 9, height = 7)
p[8]
dev.off()

pdf(paste(experiment_path,"output_plots/bar/Tumour_metacluster_neighbours_barchart", ".pdf", sep = ""), width = 9, height = 7)
p[9]
dev.off()

pdf(paste(experiment_path,"output_plots/bar/Vessels_metacluster_neighbours_barchart", ".pdf", sep = ""), width = 9, height = 7)
p[10]
dev.off()


#### Figure 4e and Supplemental Figure 6c ####

###########################################################
#### Comparing baseline and permutation neighbour data ####

# Load baseline and permutation data for each treatment group
data_baseline_Cont = read.csv(paste(experiment_path, "output_stats/neighbouRhood_data_baseline_Vehicle.csv", sep = ""))
data_baseline_MRTX = read.csv(paste(experiment_path, "output_stats/neighbouRhood_data_baseline_MRTX.csv", sep = ""))

dat_perm_Cont = read.csv(paste(experiment_path, "output_stats/neighbouRhood_dat_perm_Vehicle.csv", sep = ""))
dat_perm_MRTX = read.csv(paste(experiment_path, "output_stats/neighbouRhood_dat_perm_MRTX.csv", sep = ""))

# Add column of treatment name to each dataset 
data_baseline_Cont$treatment = "Vehicle"
dat_perm_Cont$treatment = "Vehicle"

data_baseline_MRTX$treatment = "MRTX"
dat_perm_MRTX$treatment = "MRTX"

# Merge data from treatment groups
data_baseline = rbind(data_baseline_Cont, data_baseline_MRTX)
dat_perm = rbind(dat_perm_Cont, dat_perm_MRTX)

# Create unique ID based on pairing neighbour combination, group and treatment
data_baseline$ID = paste(data_baseline$FirstLabel, "-", data_baseline$SecondLabel, "_", data_baseline$group, "_", data_baseline$treatment, sep = "")
dat_perm$ID = paste(dat_perm$FirstLabel, "-", dat_perm$SecondLabel, "_", dat_perm$group, "_", dat_perm$treatment, sep = "")

# Get an average of the ct values from each permutation run 
dat_perm_mean = dat_perm[,-c(1,2)] %>% 
  group_by(ID) %>%
  summarise_each(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))

# Check if dat_perm_mean and dat_baseline have the same number of objects 
difference = setdiff(dat_perm_mean[,c("ID","group", "FirstLabel", "SecondLabel", "treatment")], data_baseline[,c("ID","group", "FirstLabel", "SecondLabel", "treatment")])

# If there are missing baseline values, add ct = 0 for those missing values
if(nrow(difference) > 0){
  # List the values that are missing 
  print(difference$ID)
  # Create a vector of zeros for the length of values that are missing
  difference$ct = 0
  difference = difference[,c("group", "FirstLabel", "SecondLabel", "ct", "treatment", "ID")]
  data_baseline = rbind(data_baseline[,-c(1)], difference)
  if(nrow(data_baseline) == nrow(dat_perm_mean)){
    print("Dimensions of baseline and mean permutation statistics are now equal")
  } else {
    print(paste("Dimensions of baseline and mean permutation statistics are not equal - nrow baseline = ", nrow(data_baseline),
                ", nrow mean permutation = ", nrow(dat_perm_mean), sep = ""))
  }
} else{
  print("Dimensions of baseline and mean permutation statistics are equal")
}

# Order baseline data by ID, so it matches mean permutation data
data_baseline = data_baseline[order(data_baseline$ID),]

# Check that IDs from baseline and mean permutation statistics match 
match = dat_perm_mean$ID == data_baseline$ID
table(match)["FALSE"]
table(match)["TRUE"]

# Add baseline ct values to mean permutation data frame 
dat_perm_mean$ct_baseline = data_baseline$ct
# Determine the difference between baseline and permutation ct values
dat_perm_mean$ct_difference = dat_perm_mean$ct_baseline-dat_perm_mean$ct

#### Add p-value data to dat_perm for plotting ####
# Load dat_p from bodenmiller permutation analysis
dat_p_Cont = read.csv(paste(experiment_path, "output_stats/neighbouRhood_dat_pvalues_Vehicle.csv", sep = ""))
dat_p_MRTX = read.csv(paste(experiment_path, "output_stats/neighbouRhood_dat_pvalues_MRTX.csv", sep = ""))

# Add column of treatment names 
dat_p_Cont$treatment = "Vehicle"
dat_p_MRTX$treatment = "MRTX"

# Create a unique ID based on first and second labels, image number and treatment 
dat_p_Cont$ID = paste(dat_p_Cont$FirstLabel, "-", dat_p_Cont$SecondLabel, "_", dat_p_Cont$group, "_", dat_p_Cont$treatment, sep = "")
dat_p_MRTX$ID = paste(dat_p_MRTX$FirstLabel, "-", dat_p_MRTX$SecondLabel, "_", dat_p_MRTX$group, "_", dat_p_MRTX$treatment, sep = "")

# Combine p-value data for treatment groups 
dat_p = rbind(dat_p_Cont, dat_p_MRTX)
dim(dat_p)

# Reorder p-value data to match dat_perm_mean
dat_p = dat_p[order(dat_p$ID),]

# Check that the ordering of ID values for dat_p and dat_perm_mean are the same
match = dat_p$ID == dat_perm_mean$ID
table(match)["FALSE"]
table(match)["TRUE"]

# Add significance values from dat_p to dat_perm_mean
dat_perm_mean$sig = dat_p$sig

# Select cell types of interest for plotting
select = c("Epithelium", "Fibroblasts", "Lymphocytes", "Myeloid: dendritic cells", "Myeloid: macrophages type 1",
           "Myeloid: macrophages type 2", "Myeloid: neutrophils", "Normal", "Tumour", "Vessels")

dat_perm_mean = subset(dat_perm_mean, FirstLabel %in% select & SecondLabel %in% select)
dim(dat_perm_mean)


##################################
## Calculating log2 fold change ##

# Calculate log2 fold change between baseline and permutation ct values
dat_perm_mean$log2 = log2(dat_perm_mean$ct_baseline/dat_perm_mean$ct)

# Find the second smallest value in the dataset, as smallest = -Inf
small = unique(sort(dat_perm_mean$log2))[2]

# Round to the nearest decimal point 
small = round(small, 1)

# Replace the -Inf values with this value to plot data that would otherwise be excluded 
dat_perm_mean = dat_perm_mean %>%
  mutate(log2 = ifelse(FirstLabel == "Epithelium" & log2 == "-Inf", small, log2))

#Create an empty list
p = list()

for(i in 1:length(unique(dat_perm_mean$FirstLabel))){
  
  print(unique(dat_perm_mean$FirstLabel)[i])
  
  #Subset the data associated with the cellType of interest
  select_metacluster = dat_perm_mean[which(dat_perm_mean$FirstLabel == unique(dat_perm_mean$FirstLabel)[i]),]
  
  p[[i]] = local({
    i=i
    ggplot(select_metacluster, aes(x = SecondLabel, y = log2, 
                                   colour = factor(treatment, levels = c("Vehicle", "MRTX")), 
                                   fill = factor(treatment, levels = c("Vehicle", "MRTX"))))+
      geom_hline(yintercept=0)+
      geom_dotplot(data = select_metacluster[which(select_metacluster$sig ==FALSE),], binaxis = 'y', stackdir='center', dotsize=1,  fill = "white") +
      geom_dotplot(data = select_metacluster[which(select_metacluster$sig ==TRUE),], binaxis = 'y', stackdir='center', dotsize=1) +
      labs(y = "Log2FC enrichment",
           title = paste(unique(select_metacluster$FirstLabel)))+
      theme_minimal()+
      theme(plot.margin = margin(4,4,4,4,"cm"))+
      theme(axis.text.x = element_text(angle = 45,size =10, hjust = 1), plot.title = element_text(size = 14), axis.title=element_text(size=12), 
            legend.title = element_blank(), axis.title.x=element_blank())
  })
}


## Saving plots ## 
pdf(paste(experiment_path,"output_plots/point/Epithelium_neighbours_pointPlot_log2significance", ".pdf", sep = ""), width = 9, height = 7)
p[1]
dev.off()

pdf(paste(experiment_path,"output_plots/point/Fibroblasts_neighbours_pointPlot_log2significance", ".pdf", sep = ""), width = 9, height = 7)
p[2]
dev.off()

pdf(paste(experiment_path,"output_plots/point/Lymphocytes_neighbours_pointPlot_log2significance", ".pdf", sep = ""), width = 9, height = 7)
p[3]
dev.off()

pdf(paste(experiment_path,"output_plots/point/Myeloid_dendritic_cells_neighbours_pointPlot_log2significance", ".pdf", sep = ""), width = 9, height = 7)
p[4]
dev.off()

pdf(paste(experiment_path,"output_plots/point/Myeloid_macrophages_type1_neighbours_pointPlot_log2significance", ".pdf", sep = ""), width = 9, height = 7)
p[5]
dev.off()

pdf(paste(experiment_path,"output_plots/point/Myeloid_macrophages_type2_neighbours_pointPlot_log2significance", ".pdf", sep = ""), width = 9, height = 7)
p[6]
dev.off()

pdf(paste(experiment_path,"output_plots/point/Myeloid_neutrophils_neighbours_pointPlot_log2significance", ".pdf", sep = ""), width = 9, height = 7)
p[7]
dev.off()

pdf(paste(experiment_path,"output_plots/point/Normal_neighbours_pointPlot_log2significance", ".pdf", sep = ""), width = 9, height = 7)
p[8]
dev.off()

pdf(paste(experiment_path,"output_plots/point/Tumour_neighbours_pointPlot_log2significance", ".pdf", sep = ""), width = 9, height = 7)
p[9]
dev.off()

pdf(paste(experiment_path,"output_plots/point/Vessels_neighbours_pointPlot_log2significance", ".pdf", sep = ""), width = 9, height = 7)
p[10]
dev.off()


################################################################



