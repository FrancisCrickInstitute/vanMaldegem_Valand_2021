# Febe van Maldegem, The Francis Crick Institute, 2021

# Script to read in the datafile containing all ROIs after normalisation, concatenation, scaling, and clustering, 
# and add columns that describe treatment, tissue domains, cluster names and metaclusters.


# Load packages
library("readr")
library("ggplot2")
library("dplyr")

# Set directory where files will be saved
# Set output_dir as path/outputs/
output_dir = "~path/outputs/"

# Set filenm as path/Figure_3_input/input_celldata_postprocessing.csv"
filenm = "~path/Figure_3_input/input_celldata_postprocessing.csv"
# Read in dataset
celldata = read.csv(filenm)


# Set column for treatment
celldata$treatment = 1
celldata$treatment[grep("mrtx", celldata$filename, ignore.case = TRUE)] = "MRTX"
celldata$treatment[grep("veh", celldata$filename, ignore.case = TRUE)] = "Vehicle"
unique(celldata$treatment)

# Make columns with ROI_name
celldata$ROI_name = substr(celldata$filename, start = regexpr("BRAC",celldata$filename[]), stop = regexpr("\\.csv",celldata$filename[]))

# ROIname2 just as number
unique_ROInames = unique(celldata[order(factor(celldata$treatment, levels = c("Vehicle", "MRTX")), celldata$ROI_name, celldata$ObjectNumber),"ROI_name"])
for (u in unique_ROInames){
  celldata[which(celldata$ROI_name == u), "ROI_name2"] = stringr::str_pad(grep(u, unique_ROInames),2, side="left", pad="0")
}
# ROI_name3 as number_treatment
celldata$ROI_name3 = paste(celldata$ROI_name2,celldata$treatment, sep = "_" )

# Generate columns to annotate domain ("Normal", "Tumour", "Interface")
celldata$domain = as.character(celldata$MI_NormalminusInterface)
celldata$domain[which(celldata$MI_NormalminusInterface>0)] = "Normal"
celldata$domain[which(celldata$MI_TumourminusInterface>0)] = "Tumour"
celldata$domain[which(celldata$MI_InterfaceImage>0)] = "Interface"
celldata$domain[which(celldata$domain == "0")] = "n/a" #There is a small proportion
unique(celldata$domain)

# Column domain2 also includes the "Structural" domain
celldata$domain2 = celldata$domain
celldata$domain2[which(celldata$MI_StructuralImage>0)] = "Structural"



####################################################################################################
# Manual cluster splitting ("expert gating") and cluster annotation
####################################################################################################

# Generate a columnm "cluster2" that makes a string out of the cluster numbers
celldata$cluster2 = stringr::str_pad(celldata$cluster, 2, side="left", pad="0") 

# Split clusters based on thresholds for marker intensities
celldata[which(celldata$cluster == 21
               & celldata$MI_CD206 >= 0.5), "cluster2"] = "21A"
celldata[which(celldata$cluster == 21
               & celldata$MI_CD68 >= 0.5), "cluster2"] = "21A"
celldata[which(celldata$cluster2 == "21"), "cluster2"] = "21B"

celldata[which(celldata$cluster == 05
               & celldata$MI_F480 >= 0.5), "cluster2"] = "05A"
celldata[which(celldata$cluster == 05
               & celldata$MI_CD11c >= 0.6), "cluster2"] = "05B"
celldata[which(celldata$cluster2 == "05"), "cluster2"] = "05C"

celldata[which(celldata$cluster == 29
               & celldata$MI_EPCAM >= 0.7), "cluster2"] = "29A"
celldata[which(celldata$cluster2 == "29"), "cluster2"] = "29B"

celldata[which(celldata$cluster == 15
               & celldata$MI_NKp46 >= 1.5 
               & celldata$MI_pS6 <= 0.1
               & celldata$MI_Foxp3 <= 0.2), "cluster2"] = "15A"
celldata[which(celldata$cluster2 == "15"), "cluster2"] = "15B"


# Manually assign cluster labels 
celldata[which(celldata$cluster2 == "02"), "clustername"] = "02_Tumour"
celldata[which(celldata$cluster2 == "10"), "clustername"] = "10_Tumour"
celldata[which(celldata$cluster2 == "17"), "clustername"] = "17_Tumour"
celldata[which(celldata$cluster2 == "22"), "clustername"] = "22_Tumour"
celldata[which(celldata$cluster2 == "21A"), "clustername"] = "21A_Macrophages"
celldata[which(celldata$cluster2 == "21B"), "clustername"] = "21B_Tumour"
celldata[which(celldata$cluster2 == "26"), "clustername"] = "26_Macrophages"
celldata[which(celldata$cluster2 == "16"), "clustername"] = "16_Neutrophils"
celldata[which(celldata$cluster2 == "01"), "clustername"] = "01_Tumour"
celldata[which(celldata$cluster2 == "15A"), "clustername"] = "15A_NK cells"
celldata[which(celldata$cluster2 == "15B"), "clustername"] = "15B_Tumour"
celldata[which(celldata$cluster2 == "24"), "clustername"] = "24_Fibroblasts"
celldata[which(celldata$cluster2 == "05A"), "clustername"] = "05A_Macrophages"
celldata[which(celldata$cluster2 == "05B"), "clustername"] = "05B_Dendritic cells"
celldata[which(celldata$cluster2 == "05C"), "clustername"] = "05C_Tumour"
celldata[which(celldata$cluster2 == "27"), "clustername"] = "27_Dendritic cells"
celldata[which(celldata$cluster2 == "11"), "clustername"] = "11_Macrophages"
celldata[which(celldata$cluster2 == "07"), "clustername"] = "07_Endothelium"
celldata[which(celldata$cluster2 == "08"), "clustername"] = "08_Macrophages"
celldata[which(celldata$cluster2 == "06"), "clustername"] = "06_T cells"
celldata[which(celldata$cluster2 == "25"), "clustername"] = "25_T cells"
celldata[which(celldata$cluster2 == "13"), "clustername"] = "13_Endothelium"
celldata[which(celldata$cluster2 == "19"), "clustername"] = "19_Endothelium"
celldata[which(celldata$cluster2 == "12"), "clustername"] = "12_Normal lung"
celldata[which(celldata$cluster2 == "20"), "clustername"] = "20_Neutrophils"
celldata[which(celldata$cluster2 == "29A"), "clustername"] = "29A_Epithelium"
celldata[which(celldata$cluster2 == "29B"), "clustername"] = "29B_Endothelium"
celldata[which(celldata$cluster2 == "09"), "clustername"] = "09_B cells"
celldata[which(celldata$cluster2 == "14"), "clustername"] = "14_Endothelium"
celldata[which(celldata$cluster2 == "04"), "clustername"] = "04_Normal lung"
celldata[which(celldata$cluster2 == "30"), "clustername"] = "30_Neutrophils"
celldata[which(celldata$cluster2 == "03"), "clustername"] = "03_Dendritic cells"
celldata[which(celldata$cluster2 == "23"), "clustername"] = "23_Macrophages"
celldata[which(celldata$cluster2 == "18"), "clustername"] = "18_Epithelium"
celldata[which(celldata$cluster2 == "28"), "clustername"] = "28_Fibroblasts"


unique(celldata$clustername)

celldata$metacluster = as.character(celldata$clustername)
celldata[grep("Tumour", celldata$clustername), "metacluster"] = "Tumour"
celldata[grep("Normal", celldata$clustername), "metacluster"] = "Normal"
celldata[grep("Endothelium", celldata$clustername), "metacluster"] = "Vessels"
celldata[grep("Macrophage", celldata$clustername), "metacluster"] = "Myeloid: macrophages"
celldata[grep("T cells", celldata$clustername), "metacluster"] = "Lymphocytes"
celldata[grep("B cells", celldata$clustername), "metacluster"] = "Lymphocytes"
celldata[grep("NK cells", celldata$clustername), "metacluster"] = "Lymphocytes"
celldata[grep("Neutrophils", celldata$clustername), "metacluster"] = "Myeloid: neutrophils"
celldata[grep("Fibroblasts", celldata$clustername), "metacluster"] = "Fibroblasts"
celldata[grep("Dendritic", celldata$clustername), "metacluster"] = "Myeloid: dendritic cells"
celldata[grep("Epithelium", celldata$clustername), "metacluster"] = "Epithelium"
unique(celldata$metacluster)



# Split macrophage metacluster into three subtypes, two major subtypes

celldata$metacluster2 = as.character(celldata$metacluster)
celldata[which(celldata$cluster2 %in% c("11")), "metacluster2"] = "Myeloid: macrophages type 1"
celldata[which(celldata$cluster2 %in% c("05A", "08", "26" )), "metacluster2"] = "Myeloid: macrophages type 2"
unique(celldata$metacluster2)

# Write file to be used as input for the figures, or use "input_celldata_complete.csv" to continue
filenm = paste(output_dir, "celldata_complete.csv", sep = "")
write.csv(celldata, file = filenm)
