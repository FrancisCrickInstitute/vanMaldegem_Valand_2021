
#### Preparation of objectTable and Object relationships for Bodenmiller neigbouRhood permutation script ####

## Megan Cole

## Scripts required to run previous to this:
# -Normalisation, concatenation and scaling of single cell data
# -Rphenograph clustering of normalised and scaled data 
# -Cell_data_postprocessing.R
# -Code_for_figures.R
# Or load the complete dataset

# Clear working environment
rm(list=ls())

##############################
#### SET GLOBAL VARIABLES ####
##############################

# Set path for loading objectTable prep data
data_path = "~path/Figure_3_input/"

filename = "input_celldata_complete.csv"

# Set path for loading object relationships datasets
experiment_path = "~path/Neighbourhood_input/"

# Set path for saving outputs to 
output_path = "~path/Neighbourhood_output/"

# Read in marker expression data 
celldata = read.csv(paste(data_path, filename, sep = ""))

################################################################

############################
#### PREP - objectTable ####
############################

names(celldata)

# Select columns required to generate objectTable data
# (ObjectNumber, ImageNumber, cluster, filename, treatment)
celldata = celldata[,c("ObjectNumber", "ImageNumber", "metacluster2", "filename", "treatment")]
names(celldata)

unique(celldata$filename)

# Re-label ImageNumber based on order of each unique image in the dataset
for(i in 1:length(unique(celldata$filename))){
  
  print(unique(celldata$filename)[i])
  
  # Seperate the data associated with unique(celldata$filename)[i]
  celldatai = celldata[which(celldata$filename == unique(celldata$filename)[i]),]
  
  # Change the image number based on the order of the tissue files in the dataset
  celldatai$ImageNumber <- i
  print(unique(celldatai$ImageNumber))
  
  # Merge the edited celldata 
  if(i == 1) {
    celldata_numberOrder <- celldatai
  }
  else {
    # Combine data files into one
    tryCatch(
      {
        celldata_numberOrder <- rbind(celldata_numberOrder,celldatai) 
      },
      error = function(error_message){
        message("There was an error concatenating the datasets together due to:")
        message(error_message)
      }
    )
    print("DONE")
  }
}

# Add a unique cellID to each object
celldata_numberOrder$cellID = paste(celldata_numberOrder$ImageNumber, "_", celldata_numberOrder$ObjectNumber, sep = "")

# Create variables to seperate the vehicle and treatment groups
unique(celldata_numberOrder$treatment)
clusterdata_control = celldata_numberOrder[which(celldata_numberOrder$treatment == "Vehicle"),]
rownames(clusterdata_control)<-1:nrow(clusterdata_control)
clusterdata_MRTX = celldata_numberOrder[which(celldata_numberOrder$treatment == "MRTX"),]
rownames(clusterdata_MRTX)<-1:nrow(clusterdata_MRTX)

# Use just the columns needed for saving in the correct order - ObjectNumber, ImageNumber, label
clusterdata_control_order = clusterdata_control[,c("ObjectNumber","ImageNumber","metacluster2")]
clusterdata_MRTX_order = clusterdata_MRTX[,c("ObjectNumber","ImageNumber","metacluster2")]

# Change the cluster column to be named 'label' in prep for Bodenmiller neighbouRhood package 
names(clusterdata_control_order)[names(clusterdata_control_order) == 'metacluster2'] <- 'label'
names(clusterdata_MRTX_order)[names(clusterdata_MRTX_order) == 'metacluster2'] <- 'label'

names(clusterdata_control_order)
names(clusterdata_MRTX_order)

# Save clusterdata
write.csv(clusterdata_control_order, paste(output_path, "objectTable/objectTable_Vehicle.csv", sep = ""))
write.csv(clusterdata_MRTX_order, paste(output_path, "objectTable/objectTable_MRTX.csv", sep = ""))

#############################################################################

#####################################
#### PREP - Object relationships ####
#####################################

# Vector of filename order
order = c("Object relationships_BRAC3529.2d_ROI1_t1_MRTX.csv", "Object relationships_BRAC3495.3f_ROI1_t1_Vehicle.csv", "Object relationships_BRAC3495.3f_ROI1_t2_Vehicle.csv",
          "Object relationships_BRAC3529.2b_ROI1_t2_MRTX.csv", "Object relationships_BRAC3326.4e_ROI1_t1_Vehicle.csv", "Object relationships_BRAC3438.6f_ROI1_t1_Vehicle.csv",
          "Object relationships_BRAC3438.6f_ROI2_t1_Vehicle.csv", "Object relationships_BRAC3438.6f_ROI3_t1_Vehicle.csv", "Object relationships_BRAC4002.3c_ROI1_t1_MRTX.csv",
          "Object relationships_BRAC4002.3c_ROI2_t1_MRTX.csv", "Object relationships_BRAC4002.3c_ROI2_t2_MRTX.csv", "Object relationships_BRAC4002.3c_ROI3_t1_MRTX.csv")
csvfileNames = c()
for (o in order){ 
  print(o)
  c = paste(experiment_path, o, sep = "/")
  csvfileNames = append(csvfileNames, c)
}
csvfileNames

# Merge the object relationship files for images of interest together 
for(i in 1:length(csvfileNames)){
  
  print(csvfileNames[i])
  # Read the csv files in the folder
  datatemp = read.csv(csvfileNames[i])
  print(dim(datatemp))
  
  # Select only data where there is a neighbour relationship
  datatemp = datatemp[which(datatemp$Relationship == "Neighbors"),]
  print(dim(datatemp))
  unique(datatemp$Relationship)
  
  filename = gsub(paste(experiment_path, "/", sep = ""), "", csvfileNames[i])
  print(filename)
  datatemp$filename = filename
  names(datatemp)
  
  # Change image number relative to its position in the list of ROIs
  datatemp$First.Image.Number <- i
  datatemp$Second.Image.Number <- i
  
  print(unique(datatemp$First.Image.Number))
  print(unique(datatemp$Second.Image.Number))
  
  if(i == 1) {
    data <- datatemp
  }
  else {
    # Combine data files
    tryCatch(
      {
        data <- rbind(data,datatemp) 
      },
      error = function(error_message){
        message("There was an error concatenating the datasets together due to:")
        message(error_message)
      }
    )
    print("DONE")
  }
}

# Split the data based on treatment group
Control = data[grep("Vehicle", data$filename),]
MRTX = data[grep("MRTX", data$filename),]

drop = c("filename") 

# Save object relationship files 
write.csv(Control[,!(names(Control) %in% drop)], paste(output_path, "object_relationships/Object_relationships_Vehicle.csv", sep = ""))
write.csv(MRTX[,!(names(MRTX) %in% drop)], paste(output_path, "object_relationships/Object_relationships_MRTX.csv", sep = ""))

#############################################################################################################################################

