
# Rphenograph clustering of normalised, concatenated and scaled data 

## Scripts required to run previous to this:
# -Normalisation, concatenation and scaling of single cell data

## Megan Cole 

#Load packages
library(ggplot2)
library(plyr)
library(cytofkit)
library(Rphenograph)
print("Phenograph package loaded")

##############################
#### SET GLOBAL VARIABLES ####
##############################

# Name of data file to read in 
filename = "20210115_normalised_concatenated_scaled_data_001threshold.csv"

experiment_path = "~path/Clustering_input/"

output_path = "~path/Clustering_output/"

# Load data 
celldata = read.csv(paste(experiment_path, filename, sep = ""))

####################################################################

# Obtain information about the loaded dataset
dim(celldata)
names(celldata)
unique(celldata$filename)


# Remove image areas with a drop in counts
celldata = subset(celldata,
                    !(celldata[,"filename"] == 
                          grep("BRAC3326.4e_ROI1", unique(celldata$filename), ignore.case = TRUE, value = TRUE) &
                      celldata[,"Location_Center_Y"] >= 1887 &
                      celldata[,"Location_Center_Y"] <= 1924)
                  &
                    !(celldata[,"filename"] == 
                          grep("BRAC3326.4e_ROI1", unique(celldata$filename), ignore.case = TRUE, value = TRUE) &
                      celldata[,"Location_Center_Y"] >= 2380 &
                      celldata[,"Location_Center_Y"] <= 2445)

                  &
                      !(celldata[,"filename"] == 
                        grep("BRAC3438.6f_ROI2", unique(celldata$filename), ignore.case = TRUE, value = TRUE) &
                        celldata[,"Location_Center_Y"] >= 0 &
                        celldata[,"Location_Center_Y"] <= 317))

dim(celldata)


# Indicate expected length of selectcolumns and select columns for clustering
columns_of_interest = 17
selectcolumns = c(grep("MI_B220", names(celldata), ignore.case = TRUE), 
                  grep("MI_CD103", names(celldata), ignore.case = TRUE),
                  grep("MI_CD11c", names(celldata), ignore.case = TRUE),
                  grep("MI_CD3$", names(celldata), ignore.case = TRUE),
                  grep("MI_CD44", names(celldata), ignore.case = TRUE),
                  grep("MI_CD45", names(celldata), ignore.case = TRUE),
                  grep("MI_CD4$", names(celldata), ignore.case = TRUE),
                  grep("MI_CD68", names(celldata), ignore.case = TRUE),
                  grep("MI_CD8$", names(celldata), ignore.case = TRUE), 
                  grep("MI_EPCAM", names(celldata), ignore.case = TRUE),
                  grep("MI_F480", names(celldata), ignore.case = TRUE),
                  grep("MI_LY6G", names(celldata), ignore.case = TRUE),
                  grep("MI_MHCcII", names(celldata), ignore.case = TRUE),
                  grep("MI_NKp46", names(celldata), ignore.case = TRUE),
                  grep("MI_PECAM", names(celldata), ignore.case = TRUE),
                  grep("MI_PVR", names(celldata), ignore.case = TRUE),
                  grep("MI_aSMA", names(celldata),ignore.case = TRUE))

selectcolumns
length(selectcolumns)
# Check if correct columns have been selected
if(length(selectcolumns) < columns_of_interest){
  print("WARNING: there is a problem with the data selection - missing column name(s)")
} else if(length(selectcolumns) > columns_of_interest){
  # Remove repeat selection of columns 
  selectcolumns = unique(selectcolumns)
  length(unique(selectcolumns))
  #CD86 might get picked up with the CD8 selection - need to remove it if it does 
  MI_CD86 = which(names(celldata[,selectcolumns]) == "MI_CD86")
  if(isEmpty(MI_CD86)){
    length(selectcolumns)
  } else {
    selectcolumns = selectcolumns[-c(MI_CD86)]
    length(selectcolumns) 
  }
  if(length(selectcolumns) > columns_of_interest){
    # If the length is still above 17, print a warning message that data selection is wrong
    stop(paste("There is a problem wtith data selection - after an attempt to fix, there are still too many columnns selected, expected:", 
                columns_of_interest, "found:", length(selectcolumns)))
  }
  else if(length(selectcolumns) < columns_of_interest){
    stop(paste("There is a problmen with the data selection - missing column name(s), expected:", columns_of_interest,
                "found:", length(selectcolumns)))
  } else{
    print(paste("The list of selected columns has been corrected - there are now", columns_of_interest, "selected columns"))
  }
} else {
  print(paste("The correct number of columns were found - ", columns_of_interest))
}

# Retrieve names of select columns 
names(celldata[,selectcolumns])

# Run Rphenograph on data of select columns
Rphenograph_out <- Rphenograph(celldata[,selectcolumns], k=20)
print(class(Rphenograph_out))

celldata$cluster <- as.numeric(membership(Rphenograph_out[[2]]))

# Save data with clustering information added 
write.csv(celldata, paste(output_path, "Rphenograph_output_", max(celldata$cluster), "clusters_k20_", columns_of_interest, "markers.csv", sep = ""))
print("file saved")

############################################################################################################
