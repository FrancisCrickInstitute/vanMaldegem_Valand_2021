
#### Preparation of objectTable and Object relationships for Bodenmiller neigbouRhood permutation script ####

## Megan Cole

## Scripts required to run previous to this:
# -Normalisation, concatenation and scaling of single cell data
# -Rphenograph clustering of normalised and scaled data 
# -Febe script of retrieving treatment/metacluster info

# Clear working environment
rm(list = ls())

library(dplyr)

##############################
#### SET GLOBAL VARIABLES ####
##############################

# Set path for loading objectTable prep data
data_path = "Figures_input/"

filename = "input_celldata_complete.csv"

# Set path for loading object relationships datasets
experiment_path = "Neighbourhood_analysis/input/"

# Set path for saving outputs to
output_path = "Neighbourhood_analysis/output/"

# Read in marker expression data
celldata = read.csv(paste(data_path, filename, sep = ""))


################################################################


##############################
##### PREP - objectTable #####

create_objectTable = function(cellNames) {
  # cellNames = column name with cell type labels
  
  print(names(celldata))
  
  # Select columns required to generate objectTable data
  # (ObjectNumber, ImageNumber, cluster, filename, treatment, domain)
  celldata = celldata[, c("ObjectNumber",
                          "ImageNumber",
                          cellNames,
                          "filename",
                          "treatment",
                          "domain")]
  names(celldata)
  
  print(unique(celldata$filename))
  print(unique(celldata$ImageNumber))
  
  # Re-label ImageNumber based on order of each unique image in the dataset
  celldata = celldata %>% mutate(
    ImageNumber = ifelse(
      filename == "AllCellsMask_BRAC3495.3f_ROI1_t1_Vehicle.csv",
      "2",
      ImageNumber
    ),
    ImageNumber = ifelse(
      filename == "AllCellsMask_BRAC3495.3f_ROI1_t2_Vehicle.csv",
      "3",
      ImageNumber
    ),
    ImageNumber = ifelse(
      filename == "AllCellsMask_BRAC3529.2b_ROI1_t2_MRTX.csv",
      "4",
      ImageNumber
    ),
    ImageNumber = ifelse(
      filename == "AllCellsMask_BRAC3326.4e_ROI1_t1_Vehicle.csv",
      "5",
      ImageNumber
    ),
    ImageNumber = ifelse(
      filename == "AllCellsMask_BRAC3438.6f_ROI1_t1_Vehicle.csv",
      "6",
      ImageNumber
    ),
    ImageNumber = ifelse(
      filename == "AllCellsMask_BRAC3438.6f_ROI2_t1_Vehicle.csv",
      "7",
      ImageNumber
    ),
    ImageNumber = ifelse(
      filename == "AllCellsMask_BRAC3438.6f_ROI3_t1_Vehicle.csv",
      "8",
      ImageNumber
    ),
    ImageNumber = ifelse(
      filename == "AllCellsMask_BRAC4002.3c_ROI1_t1_MRTX.csv",
      "9",
      ImageNumber
    ),
    ImageNumber = ifelse(
      filename == "AllCellsMask_BRAC4002.3c_ROI2_t1_MRTX.csv",
      "10",
      ImageNumber
    ),
    ImageNumber = ifelse(
      filename == "AllCellsMask_BRAC4002.3c_ROI2_t2_MRTX.csv",
      "11",
      ImageNumber
    ),
    ImageNumber = ifelse(
      filename == "AllCellsMask_BRAC4002.3c_ROI3_t1_MRTX.csv",
      "12",
      ImageNumber
    )
  )
  print(unique(celldata$ImageNumber))
  
  # Add a unique cellID to each object
  celldata$cellID = paste(celldata$ImageNumber, "_", celldata$ObjectNumber, sep = "")
  
  # Create variables to seperate the vehicle and treatment groups
  print(unique(celldata$treatment))
  celldata_control = celldata %>% filter(treatment == "Vehicle")
  celldata_MRTX = celldata %>% filter(treatment == "MRTX")
  
  treatment_sep = list("Control" = celldata_control, "MRTX" = celldata_MRTX)
  
  # Use just the columns needed for saving in the correct order - ObjectNumber, ImageNumber, label
  celldata_control_order = celldata_control[, c("ObjectNumber", "ImageNumber", cellNames)]
  celldata_MRTX_order = celldata_MRTX[, c("ObjectNumber", "ImageNumber", cellNames)]
  
  # Change the cluster column to be named 'label' in prep for Bodenmiller neighbouRhood package
  names(celldata_control_order)[names(celldata_control_order) == cellNames] <-
    'label'
  names(celldata_MRTX_order)[names(celldata_MRTX_order) == cellNames] <-
    'label'
  
  names(celldata_control_order)
  names(celldata_MRTX_order)
  
  # Save clusterdata
  write.csv(
    celldata_control_order,
    paste(
      output_path,
      "objectTable/objectTable_wholeTissue_Vehicle.csv",
      sep = ""
    ),
    row.names = FALSE
  )
  write.csv(
    celldata_MRTX_order,
    paste(
      output_path,
      "objectTable/objectTable_wholeTissue_MRTX.csv",
      sep = ""
    ),
    row.names = FALSE
  )
  
  # Return separated treatment data for separating object tables by domain
  return(treatment_sep)
  
  print('Object tables created')
}


#######################################
##### PREP - Object relationships #####

create_object_relationships = function() {
  # Vector of filename order
  order = c(
    "Object relationships_BRAC3529.2d_ROI1_t1_MRTX.csv",
    "Object relationships_BRAC3495.3f_ROI1_t1_Vehicle.csv",
    "Object relationships_BRAC3495.3f_ROI1_t2_Vehicle.csv",
    "Object relationships_BRAC3529.2b_ROI1_t2_MRTX.csv",
    "Object relationships_BRAC3326.4e_ROI1_t1_Vehicle.csv",
    "Object relationships_BRAC3438.6f_ROI1_t1_Vehicle.csv",
    "Object relationships_BRAC3438.6f_ROI2_t1_Vehicle.csv",
    "Object relationships_BRAC3438.6f_ROI3_t1_Vehicle.csv",
    "Object relationships_BRAC4002.3c_ROI1_t1_MRTX.csv",
    "Object relationships_BRAC4002.3c_ROI2_t1_MRTX.csv",
    "Object relationships_BRAC4002.3c_ROI2_t2_MRTX.csv",
    "Object relationships_BRAC4002.3c_ROI3_t1_MRTX.csv"
  )
  csvfileNames = c()
  for (o in order) {
    print(o)
    c = paste(experiment_path, o, sep = "/")
    csvfileNames = append(csvfileNames, c)
  }
  csvfileNames
  
  # Merge the object relationship files for images of interest together
  for (i in 1:length(csvfileNames)) {
    # Read the csv files in the folder
    datatemp = read.csv(csvfileNames[i])
    
    # Select only data where there is a neighbour relationship
    datatemp = datatemp[which(datatemp$Relationship == "Neighbors"), ]
    unique(datatemp$Relationship)
    
    filename = gsub(paste(experiment_path, "/", sep = ""), "", csvfileNames[i])
    print(filename)
    datatemp$filename = filename
    
    # Change image number relative to its position in the list of ROIs
    datatemp$First.Image.Number <- i
    datatemp$Second.Image.Number <- i
    
    if (i == 1) {
      data <- datatemp
    }
    else {
      # Combine data files
      tryCatch({
        data <- rbind(data, datatemp)
      },
      error = function(error_message) {
        message("There was an error concatenating the datasets together due to:")
        message(error_message)
      })
      print("DONE")
    }
  }
  
  # Split the data based on treatment group
  Control = data[grep("Vehicle", data$filename), ]
  MRTX = data[grep("MRTX", data$filename), ]
  
  object_relationships = list('Control' = Control, "MRTX" = MRTX)
  
  drop = c("filename")
  
  # Save object relationship files
  write.csv(
    Control[, !(names(Control) %in% drop)],
    paste(
      output_path,
      "object_relationships/Object_relationships_wholeTissue_Vehicle.csv",
      sep = ""
    )
  )
  write.csv(
    MRTX[, !(names(MRTX) %in% drop)],
    paste(
      output_path,
      "object_relationships/Object_relationships_wholeTissue_MRTX.csv",
      sep = ""
    )
  )
  
  return(object_relationships)
  
  print("object relationships created")
}


######################################################
###### Seperating object relationships by domain #####

object_relationship_domains = function(Control,
                                       MRTX,
                                       objectTable_cont,
                                       objectTable_mrtx) {
  ## Control ##
  # Create a cellID column based on First Object Number and First Image Number
  Control$cellID = paste(Control$First.Image.Number,
                         Control$First.Object.Number,
                         sep = "_")
  
  ## Subset the object relationshp data by domain using domain info from objectTables
  cont_split = objectTable_cont %>% group_split(domain)
  Control_interface = Control %>%
    filter(cellID %in% cont_split[[1]]$cellID)
  
  Control_normal = Control %>%
    filter(cellID %in% cont_split[[3]]$cellID)
  
  Control_tumour = Control %>%
    filter(cellID %in% cont_split[[4]]$cellID)
  
  drop = c("filename", "cellID")
  
  # Save object tables separated by domain
  write.csv(
    Control_interface[, !(names(Control_interface) %in% drop)],
    paste(
      output_path,
      "object_relationships/Object_relationships_interface_Vehicle.csv",
      sep = ""
    ),
    row.names = FALSE
  )
  write.csv(
    Control_normal[, !(names(Control_normal) %in% drop)],
    paste(
      output_path,
      "object_relationships/Object_relationships_normal_Vehicle.csv",
      sep = ""
    ),
    row.names = FALSE
  )
  write.csv(
    Control_tumour[, !(names(Control_tumour) %in% drop)],
    paste(
      output_path,
      "object_relationships/Object_relationships_tumour_Vehicle.csv",
      sep = ""
    ),
    row.names = FALSE
  )
  
  ## MRTX ##
  MRTX$cellID = paste(MRTX$First.Image.Number,
                      MRTX$First.Object.Number,
                      sep = "_")
  
  mrtx_split = objectTable_mrtx %>% group_split(domain)
  MRTX_interface = MRTX %>%
    filter(cellID %in% mrtx_split[[1]]$cellID)
  
  MRTX_normal = MRTX %>%
    filter(cellID %in% mrtx_split[[3]]$cellID)
  
  MRTX_tumour = MRTX %>%
    filter(cellID %in% mrtx_split[[4]]$cellID)
  
  write.csv(
    MRTX_interface[, !(names(MRTX_interface) %in% drop)],
    paste(
      output_path,
      "object_relationships/Object_relationships_interface_MRTX.csv",
      sep = ""
    ),
    row.names = FALSE
  )
  write.csv(
    MRTX_normal[, !(names(MRTX_normal) %in% drop)],
    paste(
      output_path,
      "object_relationships/Object_relationships_normal_MRTX.csv",
      sep = ""
    ),
    row.names = FALSE
  )
  write.csv(
    MRTX_tumour[, !(names(MRTX_tumour) %in% drop)],
    paste(
      output_path,
      "object_relationships/Object_relationships_tumour_MRTX.csv",
      sep = ""
    ),
    row.names = FALSE
  )
  
  print("Domain object relationships created")
}


#############################################################################


## Prep object tables for whole tissue - separated by treatment
clusterdata = create_objectTable("clustername3")

# Prep object relationships data for whole tissue - separated by treatment
object_relationships = create_object_relationships()

# Prep object relationships data - separated by domains - separated by treatment
object_relationship_domains(object_relationships[['Control']],
                            object_relationships[['MRTX']],
                            clusterdata[['Control']],
                            clusterdata[['MRTX']])

## END
#############################################################################################################################################


