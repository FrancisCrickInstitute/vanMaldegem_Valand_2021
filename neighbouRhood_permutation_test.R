
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
rm(list = ls())

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
library(textshape)
library(tibble)

##############################
#### SET GLOBAL VARIABLES ####
##############################

experiment_path = "Neighbourhood_analysis/output/"

##############################################

neighbourhood_permutation = function(path, n_perm, ncores, treatment, domain) {
  # Load objectTable and Object Relationships datasets
  clusterdata = read.csv(paste(
    path,
    "objectTable/objectTable_wholeTissue_",
    treatment,
    ".csv",
    sep = ""
  ))
  neighbourdata = read.csv(
    paste(
      path,
      "object_relationships/Object_relationships_",
      domain,
      "_",
      treatment,
      ".csv",
      sep = ""
    )
  )
  
  # Drop columns not required
  drop = c("X", "filename", 'treatment', 'domain', 'cellID')
  clusterdata = clusterdata[, !names(clusterdata) %in% drop]
  # Convert clusterdata to a data table
  clusterdata = data.table(clusterdata)
  names(clusterdata)
  
  # Rename column names to ensure they match what the function looks for
  names(neighbourdata) = gsub(".Object.Name", " Object Name", names(neighbourdata))
  names(neighbourdata) = gsub(".Image.Number", " Image Number", names(neighbourdata))
  names(neighbourdata) = gsub(".Object.Number", " Object Number", names(neighbourdata))
  names(neighbourdata) = gsub("Module.Number", "Module Number", names(neighbourdata))
  neighbourdata = neighbourdata[, -1]
  names(neighbourdata)
  neighbourdata = data.table(neighbourdata)
  
  # Prepare the data - combine clusterdata and neighbourdata into one variable
  data = prepare_tables(clusterdata, neighbourdata, objname = "AllCellsMask")
  
  # Calculate baseline statistics
  data_baseline = apply_labels(data[[1]], data[[2]]) %>%
    aggregate_histo()
  
  # Save baseline statistics
  write.csv(
    data_baseline,
    paste(
      path,
      "output_stats/neighbouRhood_",
      treatment,
      "/neighbouRhood_data_baseline_",
      domain,
      "_",
      treatment,
      ".csv",
      sep = ""
    )
  )
  
  print('Baseline statistics complete')
  
  # Only run permutation statistics on data for whole tissue
  if (domain == 'wholeTissue') {
    print('Running permutation statistics')
    # Calculate permutation statistics
    set.seed(12312)
    dat_perm = rbindlist(mclapply(1:n_perm, function(x) {
      dat_labels = shuffle_labels(data[[1]])
      apply_labels(dat_labels, data[[2]]) %>%
        aggregate_histo()
    }, mc.cores = ncores), idcol = 'run')
    
    # Save permutation statistics
    write.csv(
      dat_perm,
      paste(
        path,
        "output_stats/neighbouRhood_",
        treatment,
        "/neighbouRhood_dat_perm_",
        domain,
        "_",
        treatment,
        ".csv",
        sep = ""
      )
    )
    
    # Calculate pvalues
    dat_p <-
      calc_p_vals(data_baseline,
                  dat_perm,
                  n_perm = 5000,
                  p_tresh = 0.01)
    
    # Save p-values
    write.csv(
      dat_p,
      paste(
        path,
        "output_stats/neighbouRhood_",
        treatment,
        "/neighbouRhood_dat_pvalues_",
        domain,
        "_",
        treatment,
        ".csv",
        sep = ""
      )
    )
    
    # Enrichment/depletion score for each neighbour pairing when compared to permutation statistics
    pmat = dcast(
      dat_p,
      'FirstLabel ~ SecondLabel',
      value.var = 'sigval',
      fun.aggregate = sum,
      fill = 0,
      drop = F
    )
    
    # Make FirstLabel names = row names & remove FirstLabel column
    pmat = pmat %>% remove_rownames %>% column_to_rownames(var = "FirstLabel")
    
    
    
    # Select clusters for plotting
    select = c(
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
    
    pmat = pmat[names(pmat) %in% select, names(pmat) %in% select]
    
    #pmat = pmat[select, select]
    pmat = as.matrix(pmat)
    
    
    # Save pmat
    write.csv(
      pmat,
      paste(
        path,
        "output_stats/neighbouRhood_",
        treatment,
        "/neighbouRhood_pmat_",
        domain,
        "_",
        treatment,
        ".csv",
        sep = ""
      )
    )
    
    print("neighbourhood permutation analysis complete")
  }
}


## Run neighbourhood_permutation function ##
neighbourhood_permutation(experiment_path, 5000, 10, "Vehicle", "wholeTissue")

# treatment = "Vehicle", "MRTX"
# domain = "normal", "interface", "tumour", "wholeTissue"

#############################################################################################################################################

