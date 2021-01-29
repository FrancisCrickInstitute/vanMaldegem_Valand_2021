
## Script for loading raw datasets, normlising each individually, concatenating them together and then scaling the concatenated dataset

## Segmentation must be complete before running this script

## Megan Cole

# Clear working environment
rm(list=ls())

# Load necessary packages
library("dplyr")
library("readr")

##############################
#### SET GLOBAL VARIABLES ####
##############################

# Path for data to load for scaling, concatenation and normalisation 
experiment_path = "~path/Normalisation_input"

# Path for saving output 
output_path = "~path/Normalisation_output/"


# Vector of filename order
order = c("AllCellsMask_BRAC3529.2d_ROI1_t1_MRTX.csv", "AllCellsMask_BRAC3495.3f_ROI1_t1_Vehicle.csv", "AllCellsMask_BRAC3495.3f_ROI1_t2_Vehicle.csv",
          "AllCellsMask_BRAC3529.2b_ROI1_t2_MRTX.csv", "AllCellsMask_BRAC3326.4e_ROI1_t1_Vehicle.csv", "AllCellsMask_BRAC3438.6f_ROI1_t1_Vehicle.csv",
          "AllCellsMask_BRAC3438.6f_ROI2_t1_Vehicle.csv", "AllCellsMask_BRAC3438.6f_ROI3_t1_Vehicle.csv", "AllCellsMask_BRAC4002.3c_ROI1_t1_MRTX.csv",
          "AllCellsMask_BRAC4002.3c_ROI2_t1_MRTX.csv", "AllCellsMask_BRAC4002.3c_ROI2_t2_MRTX.csv", "AllCellsMask_BRAC4002.3c_ROI3_t1_MRTX.csv")
csvfileNames = c()
for (o in order){ 
  print(o)
  c = paste(experiment_path, o, sep = "/")
  csvfileNames = append(csvfileNames, c)
  }
# Produce a list of the data files under the correct order
csvfileNames


####################################################################################################

#################################
## FUNCTIONS FOR NORMALISATION ##
#################################

# Extract the normalisation factor based on a percentile
extract.normfactor <- function(csvt, normCh, probs) {
  factor1 <- as.numeric(1/quantile(csvt[[normCh]], probs=probs))
  
  if (is.na(factor1)) {
    print("factor1 set to 1 because the percentile was found to be zero. Try increasing the threshold (thres) or the percentile (perctl).")
    factor1 = 1
  }
  factor1
}

#### NORMALISATION OF MARKER EXPRESSION WITHIN ONE IMAGE STACK ####
# Normalise a channel by multiplying with the normalisation factor
normalise.csv = function(csvt, normCh, probs){
  factor1 = extract.normfactor(csvt, normCh, probs)
  print(paste("normalisation factor based on", normCh, "is factor: ", factor1))
  for (Ch in c(1:ncol(csvt))){
    csvt[,Ch] = csvt[,Ch]*factor1
  }
  csvt
}

#### SCALING BETWEEN CONCATENATED IMAGES STACKS ####
normalise.stack = function(csvt, thres, perctl){
  for (Ch in c(1:ncol(csvt))){
    print(names(csvt[Ch]))
    m = subset(csvt, csvt[,Ch]>thres)
    factor1 = extract.normfactor(m, normCh=Ch, probs=perctl)
    print(factor1)
    csvt[,Ch] = csvt[,Ch]*factor1
  }
  csvt
}

###########################################
#############  NORMALISATION  #############
###########################################

for(i in 1:length(csvfileNames)){
  
  # Print csv image name 
  print(csvfileNames[i])
  
  # Read csv file
  celldatatemp = read.csv(csvfileNames[i])
  
  # Print dataset dimensions
  print(dim(celldatatemp))
  print(names(celldatatemp))
  
  # Shorten columns names
  colnames(celldatatemp) <- gsub("_c1", "", colnames(celldatatemp))
  colnames(celldatatemp) <- gsub("Intensity_MeanIntensity", "MI", colnames(celldatatemp))
  print(colnames(celldatatemp))
  
  # Add filename column 
  celldatatemp$filename = csvfileNames[i]
  celldatatemp$filename = gsub(experiment_path, "", celldatatemp$filename)
  celldatatemp$filename = gsub("/", "", celldatatemp$filename)
  unique(celldatatemp$filename)
  
  # Assign treatment groups bassed on keywords in the filename 
  if(length(grep("mrtx", celldatatemp$filename, ignore.case = TRUE))) {  ##Changed on 20200921, made case insensitive
    celldatatemp$treatment = "MRTX"
  } else {
    celldatatemp$treatment = "Vehicle"
  }
  
  # Set apart the columns that do not need to be normalised:
  subset_data = c(grep("ObjectNumber$", names(celldatatemp), value = TRUE),
                  grep("Location_Center_X", names(celldatatemp), value = TRUE),
                  grep("Location_Center_Y", names(celldatatemp), value = TRUE),
                  grep("TumourImage", names(celldatatemp), value = TRUE),
                  grep("TumourminusInterface", names(celldatatemp), value = TRUE),
                  grep("NormalminusInterface", names(celldatatemp), value = TRUE),
                  grep("NormalImage", names(celldatatemp), value = TRUE),
                  grep("InterfaceImage", names(celldatatemp), value = TRUE),
                  grep("StructuralImage", names(celldatatemp), value = TRUE),
                  grep("Number_Object_Number", names(celldatatemp), value = TRUE),
                  grep("AreaShape_Area", names(celldatatemp), value = TRUE),
                  grep("AreaShape_MajorAxisLength", names(celldatatemp), value = TRUE),
                  grep("AreaShape_Perimeter", names(celldatatemp), value = TRUE),
                  grep("ImageNumber", names(celldatatemp), value = TRUE),
                  grep("filename", names(celldatatemp), value = TRUE),
                  grep("treatment", names(celldatatemp), value = TRUE))
  
  celldata_subset <- celldatatemp[,subset_data]
  
  # Remove subset_data columns from celldatatemp
  celldatatemp = celldatatemp[,!names(celldatatemp) %in% subset_data]
  print(colnames(celldatatemp))
  
  # Select mean intensity marker values for normalisation
  celldatatemp <- cbind(celldatatemp[,grep("MI_", names(celldatatemp))])
  
  # NORMALISATION - Call normalise.csv function for within stack normalisation
  celldatatemp <- normalise.csv(celldatatemp, "MI_Xenon134", c(0.5))
  
  # Add in subset data 
  celldatatemp = cbind(celldatatemp, celldata_subset)
  
  ## CONCATENATE THE DATA FILES ##
  if(i == 1) {
    celldata <- celldatatemp
  }  else {
    # Combine multiple data files into one variable using rbind function
    tryCatch(
      {
        celldata <- rbind(celldata,celldatatemp)
      },
      error = function(error_message){
        message("There was an error concatenating the datasets together due to:")
        message(error_message)
      }
    )
    print("DONE")
  }
}

names(celldata)

# Save normalised and concatenated file 
write.csv(celldata, paste(output_path, "normalised_concatenated_data.csv", sep = ""))
print("Normalised and concatenated file saved")


###############################
##########  SCALING  ##########
###############################

## SCALING - Call normalise.stack function for between stack scaling
celldata_subset <- celldata[,subset_data]
celldatatemp = celldata[,!names(celldata) %in% subset_data]
names(celldatatemp)

# Scale the data to 99th percentile
celldatatemp <- normalise.stack(celldatatemp, 0.001, c(0.99)) 
celldata = cbind(celldatatemp, celldata_subset)

# Save normalised, concatenated and scaled file 
write.csv(celldata, paste(output_path, "normalised_concatenated_scaled_data_001threshold.csv", sep = ""))
print("Normalised, concatenated and scaled file saved")


