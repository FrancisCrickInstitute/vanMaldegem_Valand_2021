# Febe van Maldegem, The Francis Crick Institute, 2021

# These scripts takes a raw data files (no normalisation or scaling), from segmenting with different strategies


library("readr")
library("dplyr")
library(gplots) # For heatmap.2()
library(ggplot2)
library("RColorBrewer")

# Set global variables
experiment_path = "~path/Figure_2_input"
output_dir = "~path/Figure_2_output/"


# Load all files from experiment_path into list of dataframes: "my_data", from a folder that ends with "input/"
extract.filename = function(file){
  unlist(strsplit(file, "input/"))[2]
}
open.files = function(experiment_path){
  ROI_path = dir(experiment_path, recursive = FALSE, include.dirs = TRUE, full.names = TRUE, pattern = ".csv")
  my_data <- lapply(ROI_path, read.csv)
  names(my_data) = lapply(ROI_path, extract.filename)
  print(names(my_data))
  my_data
}

# Open all the files into a list "my_data":
my_data = open.files(experiment_path)
# Change order of appearance of files
my_data = my_data[c(4,1,3,2)]

names(my_data)

# Rename columns for convenience
for (i in 1:length(my_data)){
  colnames(my_data[[i]]) = gsub("_c1", "", colnames(my_data[[i]]))
  colnames(my_data[[i]]) = gsub("Intensity_MeanIntensity", "MI", colnames(my_data[[i]]))
}


##############################################################################################
# Generate heatmap that plots signal enrichment.
#
# Figure 2d

#  Function to calculate signal enrichment
heatmap_signal_noise_top500 = function(my_data, markerlist = c(), n = 500){
  heatmap_df = data.frame()
  for (i in 1:length(my_data)){
    print(i)
    for (marker in markerlist){
      md = my_data[[i]]
      md = md[order(md[,paste("MI_", marker, sep = "")], decreasing = TRUE),]
      md = md[1:n,]
      md_rest = my_data[[i]][!(rownames(my_data[[i]]) %in% rownames(md)),]
      ratio = mean(md[,paste("MI_", marker, sep = "")])/mean(md_rest[,paste("MI_", marker, sep = "")])
      heatmap_df[i, marker] = ratio
    }
  }
  heatmap_df = as.matrix(heatmap_df)
  row.names(heatmap_df) = names(my_data)
  return(heatmap_df)
}

markerlist = c("aSMA","PVR", "CD45","B220","CD3","CD4", "CD8", "CD44", "CD103", "MHCcII",  "CD11c", "NKp46", "F480", "CD68","LY6G","EPCAM", "PECAM")
heatmap_df = heatmap_signal_noise_top500(my_data, markerlist = markerlist, n = 500)

pdf(file=paste(output_dir, "_heatmap_signalenrichment.pdf", sep = ""), width=20, height=10)

heatmap(heatmap_df, scale = "column", col= colorRampPalette(brewer.pal(8, "Blues"))(25), margins=c(10,0), cexRow = 1, cexCol = 1, Rowv = NA, Colv = NA)
# Or with the scale included: 
# heatmap.2(heatmap_df, col= colorRampPalette(brewer.pal(8, "Blues"))(25), scale = "column", cexRow = 1.9, cexCol = 1.9,
#           density.info="none",  # turns off density plot inside color legend
#           trace="none", keysize=0.6, dendrogram ="none",        
#           Colv = NULL, Rowv = NULL,
#           margins =c(15,15))
dev.off()


##############################################################################################
# Generate heatmap that plots pollution of mutually exclusive markers within selected population.
#
# Fig 2e

heatmap_signal_noise_within_top500 = function(my_data, inclusion_markerlist = c(), exclusion_markerlist = c(), n=500){
  heatmap_df = data.frame()
  for (i in 1:length(my_data)){
    print(i)
    for (marker in inclusion_markerlist){
      md = my_data[[i]]
      md = md[order(md[,paste("MI_", marker, sep = "")], decreasing = TRUE),]
      md = md[1:n,]
      for (marker2 in exclusion_markerlist){
        heatmap_df[i, paste(marker, "/", marker2, sep = "")] = mean(md[,paste("MI_", marker, sep = "")])/mean(md[,paste("MI_", marker2, sep = "")])
      }
    }
  }
  heatmap_df = as.matrix(heatmap_df)
  row.names(heatmap_df) = names(my_data)
  return(heatmap_df)
}


heatmap_df = heatmap_signal_noise_within_top500(my_data, 
                                                inclusion_markerlist = c( "CD4", "aSMA" , "PECAM", "F480"), 
                                                exclusion_markerlist = c( "F480",  "CD103", "aSMA", "PVR", "PECAM"), 
                                                n=500)

colnames(heatmap_df)
minlist = c("aSMA/CD8", "F480/F480", "aSMA/aSMA", "PECAM/PECAM", "CD45/CD45", "PVR/PVR", "CD103/CD103") #, #, "CD4/aSMA") #, "CD4/F480", "F480/B220", "PECAM/F480", "aSMA/CD103" )
heatmap_df = heatmap_df[,-which(colnames(heatmap_df) %in% minlist)]

pdf(file=paste(output_dir, "_heatmap_signalnoise.pdf", sep = ""), width=20, height=10)
heatmap(heatmap_df, 
        scale = "column", col= colorRampPalette(brewer.pal(8, "Blues"))(25), margins=c(17,0), cexRow = 1, cexCol = 0.6, Colv = NA, Rowv = NA)
# Or with the scale included: 
# heatmap.2(heatmap_df, col= colorRampPalette(brewer.pal(8, "Blues"))(25), scale = "column", cexRow = 1.9, cexCol = 1.9,
#           density.info="none",  # turns off density plot inside color legend
#           trace="none", keysize=0.6, dendrogram ="none",         
#           Colv = NA, Rowv = NA,
#           margins =c(15,15))
dev.off()


##############################################################################################
# Comparison to manual annotations

# Fig 2f

# This function takes in the reference file, as well as the single cell data. 
# It ranks the single cell data to the top 500 cells, and determines how many of those cells "overlap" with the reference set.
# Overlap is defined by the centre of cell A being within a 5px distance of the centre of cell B.
find_overlap_top500 = function(input_df1 = reference, input_df2 = my_data[[i]], marker = marker, radius = radius, n=n) {
  df1 = input_df1
  df2 = data.frame()
  df2 = round(input_df2[,"Location_Center_X"])
  df2 = as.data.frame(df2)
  df2$Y = round(input_df2[,"Location_Center_Y"])
  df2[,marker] = input_df2[,paste("MI_", marker, sep = "")]
  df2 = df2[order(df2[,marker], decreasing = TRUE),]
  df2 = df2[1:n,]
  names(df2) = c("X", "Y", marker)
  rownames(df2) = c(1:nrow(df2))
  for (x in 1:nrow(df1)){
    c = ((df2$X - df1[x,"X"])^2 + (df2$Y - df1[x,"Y"])^2)
    df1$distance[x] = sqrt(min(((df2$X - df1[x,"X"])^2 + (df2$Y - df1[x,"Y"])^2), na.rm = TRUE))
    f = order(c)[1]
    df1[x,"X_df2"] = df2[f,"X"]
    df1[x,"Y_df2"] = df2[f,"Y"]
    rm(c, f)  
  }
  overlap = df1[which(df1$distance <= radius),]
  outputs = list(overlap = overlap)
  return(outputs)
}

# This function loops the "find_overlap_top500" function over the different segmentation files. 
find_overlap_with_reference_top500 = function(input_df1 = reference, my_data = my_data, marker = marker, radius = radius, n=n){
  outputs_list = list()
  print(dim(reference[which(reference[,marker] == 1),]))
  for (i in names(my_data)){
    local({
      i=i
      print(i)
      outputs = find_overlap_top500(input_df1 = reference, input_df2 = my_data[[i]], marker = marker, radius = radius, n=n)
      outputs_list[[i]] <<- outputs[[1]]
      print(dim(outputs[[1]]))
      
    })
  }
  return(outputs_list)
}
# Choose the relevant reference file
reference_file = "~path/Figure_2_input/manual_annotation/manual_annotation_CD4Tcells.csv"
reference_file = "~path/Figure_2_input/manual_annotation/manual_annotation_CD8Tcells.csv"
reference_file = "~path/Figure_2_input/manual_annotation/manual_annotation_DCs.csv"
reference_file = "~path/Figure_2_input/manual_annotation/manual_annotation_F480.csv"
reference_file = "~path/Figure_2_input/manual_annotation/manual_annotation_aSMA.csv"


reference = read.csv(reference_file, sep = ";")

names(reference)[1] = "X"
marker = "CD4"  #match marker to reference file
marker = "CD8"
marker = "CD103"
marker = "F480"
marker = "aSMA"


radius = 5 #pixels
n=500 #top500 cells
outputs_list = find_overlap_with_reference_top500(input_df1 = reference, my_data = my_data, marker = marker, radius = radius, n=n)

# Make new table first time, then disable so that new measurements will be added to the table
results_table = as.data.frame(names(outputs_list)) #disable after first use
cn = paste("overlap_", marker, "_", n, "n_", radius, "px", sep = "")
results_table[[cn]]  = lapply(outputs_list, dim)
results_table[[cn]] = lapply(results_table[[cn]], "[[", 1)
results_table


##############################################################################################
# Example domain segmentation

# Fig 2g

# Plotting the cells as dots using their X-Y coordinates, coloured by domain
cd = my_data[[1]]
cd$domain = "n/a"
cd$domain[which(cd$MI_NormalminusInterface>0)] = "Normal"
cd$domain[which(cd$MI_TumourminusInterface>0)] = "Tumour"
cd$domain[which(cd$MI_InterfaceImage>0)] = "Interface"
cd$domain[which(cd$MI_StructuralImage>0)] = "Structural"

cd$MI_PECAM = cd$MI_PECAM*100000
cd$MI_CD44 = cd$MI_CD44*100000
cd$MI_EPCAM = cd$MI_EPCAM*100000

domain_cols <- c("Normal" = "red", "Interface" = "green", "Tumour" = "purple", "Structural" = "cyan")
p = ggplot(cd, aes(x= Location_Center_X, y = -Location_Center_Y, color = factor(domain, levels=c("Normal", "Interface", "Tumour", "Structural")))) +
  geom_point(size = 0.2) +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_rect(size = 0.5),
        panel.spacing = unit(0, "mm"),
        plot.background = element_blank()) +
  scale_color_manual(values = domain_cols)
p

filename = "domain_segmentation.png"
ggsave(plot = p, device = "png", width=4, height=5, dpi=300, path = output_dir, filename = filename)

# Violin plots for PECAM expression split by domain
p = ggplot(cd[which(cd$domain != "n/a"),], 
       aes(x= factor(domain, levels=c("Normal", "Interface", "Tumour", "Structural")), 
           y=MI_PECAM, fill = factor(domain, levels=c("Normal", "Interface", "Tumour", "Structural")))) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "none",    
    axis.title.x = element_blank(),
    axis.title.y = element_text(),
    axis.ticks = element_blank(),
    legend.title = element_blank()) +
  ylab("PECAM")+
  scale_fill_manual(values = domain_cols)
p
filename = "PECAM_domain_distribution.png"
ggsave(plot = p, device = "png", width=4, height=5, dpi=300, path = output_dir, filename = filename)

# Violin plots for CD44 expression split by domain
p = ggplot(cd[which(cd$domain != "n/a"),], 
       aes(x= factor(domain, levels=c("Normal", "Interface", "Tumour", "Structural")), 
           y=MI_CD44, fill = factor(domain, levels=c("Normal", "Interface", "Tumour", "Structural")))) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "none",    
    axis.title.x = element_blank(),
    axis.title.y = element_text(),
    axis.ticks = element_blank(),
    legend.title = element_blank()) +
  ylab("CD44")+
  scale_fill_manual(values = domain_cols)
p
filename = "CD44_domain_distribution.png"
ggsave(plot = p, device = "png", width=4, height=5, dpi=300, path = output_dir, filename = filename)
