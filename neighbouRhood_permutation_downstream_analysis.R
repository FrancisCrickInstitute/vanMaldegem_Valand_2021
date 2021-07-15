

#### Downstream analysis of neighbourhood permutation test ####

## Megan Cole
## 2021

## Scripts required to run previous to this:
# -Normalisation, concatenation and scaling of single cell data
# -Rphenograph clustering of normalised and scaled data
# -Preparation for neighbourhood permuation
# -NeighbouRhood permutation

# Clear working environment
rm(list = ls())

# Load required packages
library(dplyr)
library(ggplot2)
library(textshape)
library(tibble)

############################
##### GLOBAL VARIABLES #####
############################

experiment_path = "Neighbourhood_analysis/output/"

###############################################################################################################

###########################################################
#### Comparing baseline and permutation neighbour data ####

data_prep = function(domain, path, select) {
  # Load baseline and permutation data for each treatment group
  data_baseline_Cont = read.csv(
    paste(
      path,
      "output_stats/neighbouRhood_data_baseline_",
      domain,
      "_Vehicle.csv",
      sep = ""
    )
  )
  data_baseline_MRTX = read.csv(
    paste(
      path,
      "output_stats/neighbouRhood_data_baseline_",
      domain,
      "_MRTX.csv",
      sep = ""
    )
  )
  
  dat_perm_Cont = read.csv(
    paste(
      path,
      "output_stats/neighbouRhood_Vehicle/neighbouRhood_dat_perm_wholeTissue_Vehicle.csv",
      sep = ""
    )
  )
  dat_perm_MRTX = read.csv(
    paste(
      path,
      "output_stats/neighbouRhood_MRTX/neighbouRhood_dat_perm_wholeTissue_MRTX.csv",
      sep = ""
    )
  )
  
  # Add column of treatment name to each dataset
  data_baseline_Cont$treatment = "Vehicle"
  dat_perm_Cont$treatment = "Vehicle"
  
  data_baseline_MRTX$treatment = "MRTX"
  dat_perm_MRTX$treatment = "MRTX"
  
  # Merge data from treatment groups
  data_baseline = rbind(data_baseline_Cont, data_baseline_MRTX)
  dat_perm = rbind(dat_perm_Cont, dat_perm_MRTX)
  
  # Create unique ID based on pairing neighbour combination, group and treatment
  data_baseline$ID = paste(
    data_baseline$FirstLabel,
    data_baseline$SecondLabel,
    data_baseline$group,
    data_baseline$treatment,
    sep = "_"
  )
  dat_perm$ID = paste(
    dat_perm$FirstLabel,
    dat_perm$SecondLabel,
    dat_perm$group,
    dat_perm$treatment,
    sep = "_"
  )
  
  # Get an average of the ct values from each permutation run
  dat_perm_mean = dat_perm[,-c(1, 2)] %>%
    group_by(ID) %>%
    summarise_each(funs(if (is.numeric(.))
      mean(., na.rm = TRUE)
      else
        first(.)))
  
  # Check if dat_perm_mean and dat_baseline have the same number of objects
  difference = setdiff(dat_perm_mean[, c("ID", "group", "FirstLabel", "SecondLabel", "treatment")], data_baseline[, c("ID", "group", "FirstLabel", "SecondLabel", "treatment")])
  
  # If there are missing baseline values, add ct = 0 for those missing values
  if (nrow(difference) > 0) {
    # List the values that are missing
    print(difference$ID)
    # Create a vector of zeros for the length of values that are missing
    difference$ct = 0
    # Order difference columns to match order of data_baseline, ready for binding
    difference = difference[, c("group",
                                "FirstLabel",
                                "SecondLabel",
                                "ct",
                                "treatment",
                                "ID")]
    data_baseline = rbind(data_baseline[,-c(1)], difference)
    #data_baseline = rbind(data_baseline, difference)
    if (nrow(data_baseline) == nrow(dat_perm_mean)) {
      print("Dimensions of baseline and mean permutation statistics are now equal")
    } else {
      print(
        paste(
          "Dimensions of baseline and mean permutation statistics are not equal - nrow baseline = ",
          nrow(data_baseline),
          ", nrow mean permutation = ",
          nrow(dat_perm_mean),
          sep = ""
        )
      )
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
  
  #### Add p-value data to dat_perm for plotting ####
  # Load dat_p from bodenmiller permutation analysis
  dat_p_Cont = read.csv(
    paste(
      path,
      "output_stats/neighbouRhood_Vehicle/neighbouRhood_dat_pvalues_wholeTissue_Vehicle.csv",
      sep = ""
    )
  )
  dat_p_MRTX = read.csv(
    paste(
      path,
      "output_stats/neighbouRhood_MRTX/neighbouRhood_dat_pvalues_wholeTissue_MRTX.csv",
      sep = ""
    )
  )
  
  # Add column of treatment names
  dat_p_Cont$treatment = "Vehicle"
  dat_p_MRTX$treatment = "MRTX"
  
  # Create a unique ID based on first and second labels, image number and treatment
  dat_p_Cont$ID = paste(
    dat_p_Cont$FirstLabel,
    dat_p_Cont$SecondLabel,
    dat_p_Cont$group,
    dat_p_Cont$treatment,
    sep = "_"
  )
  dat_p_MRTX$ID = paste(
    dat_p_MRTX$FirstLabel,
    dat_p_MRTX$SecondLabel,
    dat_p_MRTX$group,
    dat_p_MRTX$treatment,
    sep = "_"
  )
  
  # Combine p-value data for treatment groups
  dat_p = rbind(dat_p_Cont, dat_p_MRTX)
  dim(dat_p)
  
  # Save dat_p for combined treatment groups
  if (domain == 'wholeTissue') {
    write.csv(dat_p,
              paste(path, "output_stats/dat_p_Vehicle_MRTX.csv", sep = ""))
  }
  
  # Reorder p-value data to match dat_perm_mean
  dat_p = dat_p[order(dat_p$ID),]
  
  # Check that the ordering of ID values for dat_p and dat_perm_mean are the same
  match = dat_p$ID == dat_perm_mean$ID
  table(match)["FALSE"]
  table(match)["TRUE"]
  
  # Add significance values from dat_p to dat_perm_mean
  dat_perm_mean$sig = dat_p$sig
  
  dat_perm_mean = subset(dat_perm_mean, FirstLabel %in% select &
                           SecondLabel %in% select)
  dim(dat_perm_mean)
  
  print("Data prep function complete")
  
  # Call log2FC function
  log2FC(domain, path, dat_perm_mean)
}


##################################
## Calculating log2 fold change ##

log2FC = function(domain, path, data) {
  # Calculate log2 fold change between baseline and permutation ct values
  data$log2 = log2(data$ct_baseline / data$ct)
  
  # Save data with calculations of difference between baseline and permutation ct values
  write.csv(
    data,
    paste(
      path,
      "output_stats/neighbouRhood_dat_perm_mean_log2FC_",
      domain,
      ".csv",
      sep = ""
    ),
    row.names = F
  )
  
  for (i in 1:length(unique(data$FirstLabel))) {
    select_metacluster = data[which(data$FirstLabel == unique(data$FirstLabel)[i]),]
    print(unique(data$FirstLabel)[i])
    
    p = ggplot(select_metacluster,
               aes(
                 x = SecondLabel,
                 y = log2,
                 colour = factor(treatment, levels = c("Vehicle", "MRTX")),
                 fill = factor(treatment, levels = c("Vehicle", "MRTX"))
               )) +
      geom_hline(yintercept = 0) +
      geom_dotplot(
        data = select_metacluster[which(select_metacluster$sig == FALSE),],
        binaxis = 'y',
        stackdir = 'center',
        dotsize = 1,
        fill = "white"
      ) +
      geom_dotplot(
        data = select_metacluster[which(select_metacluster$sig == TRUE),],
        binaxis = 'y',
        stackdir = 'center',
        dotsize = 1
      ) +
      theme_minimal() +
      labs(title = unique(data$FirstLabel)[i],
           x = "",
           y = "Log2FC enrichment") +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = 14
        ),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)
      ) +
      theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
      scale_x_discrete(
        limits = c(
          "CD4 T cells",
          "Regulatory T cells",
          "CD8 T cells",
          "B cells",
          "Dendritic cells cDC1",
          "Dendritic cells other",
          "Macrophages type 1",
          "Macrophages type 2",
          "Neutrophils",
          "Endothelium",
          "Epithelium",
          "Fibroblasts",
          "Tumour"
        )
      )
    
    if (unique(data$FirstLabel)[i] == "Macrophages type 1" |
        unique(data$FirstLabel)[i] == "Macrophages type 2") {
      print("Setting y limits")
      p = p + ylim(-4, 2.5)
    }
    
    filename = paste(
      unique(data$FirstLabel)[i],
      "_neighbours_pointPlot_log2significance_perTreatment.png",
      sep = ""
    )
    ggsave(
      plot = p,
      device = "png",
      width = 7,
      height = 5,
      dpi = 300,
      bg = 'White',
      path = paste(
        path,
        "output_plots/point/new/",
        domain,
        "_v_wholeTissue/",
        sep = ""
      ),
      filename = filename
    )
  }
  print("log2FC plotting complete")
}

#######################
#### Run functions ####

# Cell types of interest
cellTypes = c(
  "B cells",
  "CD4 T cells",
  "CD8 T cells",
  "Dendritic cells cDC1",
  "Dendritic cells other",
  "Endothelium",
  "Epithelium",
  "Fibroblasts",
  "Macrophages type 1",
  "Macrophages type 2",
  "Neutrophils",
  "Regulatory T cells",
  "Tumour"
)

# Call the function to prep the data & then run log2FC calculation and plotting
data_prep("tumour", experiment_path, cellTypes)
# Domain options: 'normal', 'interface', 'tumour', 'wholeTissue'


## END 
################################################################





