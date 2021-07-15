# Febe van Maldegem, The Francis Crick Institute, 2021

# Script to generate the figures from Van Maldegem, Valand, et al. 2021


# Clear loaded packages
# lapply(paste("package:",names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

# Load packages
# install::packages("devtools")
library(tidyverse)
library("readr")
library("ggplot2")
library("dplyr")
library(gplots) #library for plotting data
library(tools)
# install.packages('HSAUR') ### statistical analysis package
library(HSAUR)
# install.packages('pheatmap') ### pheatmap function for k-means clustering visualisation
library('pheatmap')
# install.packages('Rtsne') ### installing packages and libraries for tSNE plots in R
library(Rtsne)
library(scales)
# install.packages('colorRamps') ### installing a package to produce plots with different colours
library(colorRamps)
# install.packages("corrplot")
library(corrplot)
# install.packages("Hmisc")
library(Hmisc)
library("RColorBrewer")
# install.packages("plotly")
library("plotly")
# install.packages("umap")
library(umap)
# install.packages("ggpubr")
library(ggpubr)
library(factoextra)


#
# Set working directory to path where files are stored

# Set output_dir as ~working_dir/outputs/
output_dir = "output/"

filenm = "Figures_input/input_celldata_complete.csv"
# This file contains all the additional columns and calculations that are made in this script

# Read in dataset
celldata = read.csv(filenm)


##############################################################################################
# Visualising the clustering results with a heatmap
#
# Figure 3b

# Function to generate a heatmap with markers versus clusters, labeled with clustername2
# lineage markers only (17):
markerlist = c("CD45", "CD3", "CD4", "Foxp3","CD8", "B220", "NKp46", "MHCcII","CD11c","CD103", "F480", "CD68",  "LY6G", "PECAM", "EPCAM","aSMA","PVR", "CD44")
heatmap_markers_clusters = function(celldata, markerlist){
  heatmap_df = data.frame()
  unique_clusters =  c("CD4 T cells","Regulatory T cells","CD8 T cells","B cells","NK cells",
                       "Dendritic cells cDC1", "Dendritic cells other", "Macrophages","Neutrophils",
                       "Endothelium","Epithelium","Fibroblasts","Tumour")
  
    for (cl in unique_clusters){
    md = celldata[which(celldata$clustername2 == cl),]
    cl2 = unique(md$clustername2)
    for (marker in markerlist){
      heatmap_df[cl, marker] =  mean(md[,paste("MI_", marker, sep = "")])
      heatmap_df[cl, "cl2"] = cl2
    }
  }
  cl2name = c()
  cl2name = heatmap_df$cl2
  heatmap_df = heatmap_df[,-2]
  heatmap_df = as.matrix(heatmap_df)
  row.names(heatmap_df) = cl2name
  return(heatmap_df)
}
heatmap_df = heatmap_markers_clusters(celldata[which(celldata$clustername2 != "Unclassified"),], markerlist)
dev.off()
pdf(file=paste(output_dir, "heatmap_clustername2.pdf", sep = ""), width=12, height=10)
heatmap.2(heatmap_df,
          col= colorRampPalette(brewer.pal(8, "Blues"))(25), scale = "row",Rowv = NA,Colv = NA, cexRow = 2.5, cexCol = 2.5,
          density.info="none",  # turns off density plot inside color legend
          trace="none", keysize=0.6,
          lmat = rbind(c(3,4),c(2,1)), #lwid = c(1.5,4), lhei = c(1.5,4,1), #This moves the key
          margins =c(10,25)
          )
dev.off()



##############################################################################################

# Figure 3c
# Stacked bargraph of the cell types coloured by domain
# This makes two bar plots, one for Vehicle, and one for MRTX treatment.

# Statistics per ROI
# Calculate cell type distribution across domains for each ROI "cd_prop_d":
cd_prop_d = data.frame()
n = 0

for (r in sort(unique(celldata$ROI_name3))){
  print(r)
  treatment = as.character(unique(celldata[which(celldata$ROI_name3 == r), "treatment"]))
  MouseID = as.character(unique(celldata[which(celldata$ROI_name3 == r), "MouseID"]))
  for (cl in sort(unique(celldata$clustername2))){
    print(cl)
    cd = celldata[which(celldata$ROI_name3 == r & celldata$domain2 != "n/a" & celldata$clustername2 == cl),]
    a = nrow(cd)
    for (d in c("Normal", "Interface", "Tumour", "Structural")){
      print(d)
      n = n+1
      cd = celldata[which(celldata$ROI_name3 == r & celldata$domain != "n/a" & 
                            celldata$clustername2 == cl & celldata$domain2 == d),]
      b = nrow(cd)
      cd_prop_d[n,"treatment"] = treatment
      cd_prop_d[n,"MouseID"] = MouseID
      cd_prop_d[n,"ROI_name3"] = r
      cd_prop_d[n,"clustername2"] = cl
      cd_prop_d[n,"domain2"] = d
      cd_prop_d[n,"prop"] = b/a
      
    }
  }  
}
cd_prop_d$domain2 = factor(cd_prop_d$domain2, levels = c("Normal", "Interface", "Tumour", "Structural"))

cluster_order = c("CD4 T cells","Regulatory T cells","CD8 T cells","B cells", "NK cells",
             "Dendritic cells cDC1", "Dendritic cells other", "Macrophages","Neutrophils",
             "Endothelium","Epithelium","Fibroblasts","Tumour")
domain_cols <- c("Normal" = "red", "Interface" = "green", "Tumour" = "purple", "Structural" = "cyan")
p = ggbarplot(cd_prop_d[which(cd_prop_d$treatment == "Vehicle"),], x = "clustername2", y="prop", add = "mean_se", fill = "domain2") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, len = 5)) +
  scale_x_discrete(position = "top", limits = cluster_order) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,size =16, hjust = 0), 
        axis.text.y = element_text(size =16),
        plot.title = element_text(size = 6), 
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "top",
        axis.title=element_blank(),
        strip.text = element_text(size = 12),
        panel.spacing = unit(0, "mm"),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA, size = 1),
        strip.background = element_rect(size = 0.5),
        plot.background = element_blank()) +
  scale_fill_manual(values = domain_cols)
  
filename = "Stacked_metaclusters_by_domain_Vehicle.pdf"
ggsave(plot = p, device = "pdf", width=6, height=4.5, dpi=300, path = output_dir, filename = filename)

p = ggbarplot(cd_prop_d[which(cd_prop_d$treatment == "MRTX"),], x = "clustername2", y="prop", add = "mean_se", fill = "domain2") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, len = 5)) +
  scale_x_discrete(position = "top", limits = cluster_order) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,size =16, hjust = 0),
        axis.text.y = element_text(size =16),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title= element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0, "mm"),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA, size = 1),
        strip.background = element_rect(size = 0.5),
        plot.background = element_blank())+
  scale_fill_manual(values = domain_cols)
p
filename = "Stacked_metaclusters_by_domain_MRTX.pdf"
ggsave(plot = p, device = "pdf", width=6, height=4.5, dpi=300, path = output_dir, filename = filename)

#####################################

##########################
# Figure 3g
# Proportion of cell types in tumour domain
p = ggplot(celldata[which(celldata$domain == "Tumour" & celldata$clustername2 != "Unclassified"),], 
       aes(x = factor(treatment, levels = c("Vehicle", "MRTX")), fill = factor(clustername2, levels = cluster_order))) +
  geom_bar(position = "fill", colour = "black") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,size =14, hjust = 1), 
        axis.text.y = element_text(size =14),
        plot.title = element_text(size = 6), 
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title=element_blank()) 
p
filename = "Stacked_metaclusters_by_treatment.pdf"
ggsave(plot = p, device = "pdf", width=4.2, height=5, dpi=300, path = output_dir, filename = filename)



##############################################################################################
# Figure 4 - Macrophages 


# Figure 4a
# Run uMAP on macrophage metacluster, or skip this step and load the dataset with uMAP coordinates

selectcolumns = c("MI_MHCcII", "MI_F480",
                  "MI_PDL1", "MI_CD86", "MI_CD68",
                  "MI_CD206", "MI_pS6", "MI_CD11c") #, "MI_LY6G", "MI_CD103", "MI_CD44")
cd = celldata[which(celldata$clustername2 == "Macrophages"),]
umap_results = umap(cd[,selectcolumns], verbose = TRUE, n_neighbors = 10)
cd$umap1 = umap_results$layout[,1]
cd$umap2 = umap_results$layout[,2]
macrophage_uMAP = cd

umap_total = cd
# Load  macrophage dataset including uMAP coordinates
filenm = "Figures_input/macrophage_umap.csv"
cd = read.csv(filenm)
macrophage_uMAP = cd


# Load macrophage uMAP including uMAP coordinates
filenm = 
  filenm = "Figures_input/macrophage_umap.csv"
cd = read.csv(filenm)


# Figure 4a
# uMAP coloured by cluster
cluster_cols = c("yellow", "#00BA38", "#619CFF", "cyan", "#E76BF3", "orange") # "#F8766D")
p =ggplot(cd[sample(nrow(cd), nrow(cd)),], aes(x=umap1, y=umap2, colour = as.factor(cluster2))) + 
  geom_point(size = 0.3) +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line =  element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank()) +
  scale_color_manual(values = cluster_cols)
p
filename = "uMAP_macrophage_clusters.png"
ggsave(plot = p, device = "png", width=4, height=4, dpi=300, path = output_dir, filename = filename)

# uMAP coloured by domain
domain_cols <- c("Normal" = "red", "Interface" = "green", "Tumour" = "purple", "Structural" = "cyan")
p = ggplot(cd[which(cd$domain %in% c( "Normal", "Tumour", "Interface")),], 
           aes(x=umap1, y=umap2, color = factor(domain, levels=c("Normal", "Interface", "Tumour")))) +
  geom_point(size = 0.3) +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line =  element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank()) +
  scale_color_manual(values = domain_cols)
p
filename = "uMAP_macrophage_domains.png"
ggsave(plot = p, device = "png", width=4, height=4, dpi=300, path = output_dir, filename = filename)

# uMAP coloured by treatment
p = ggplot(cd[sample(nrow(cd), nrow(cd)),], 
           aes(x=umap1, y=umap2, color = factor(treatment, levels=c("Vehicle", "MRTX")))) +
  geom_point(size = 0.3) +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line =  element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank())
p
filename = "uMAP_macrophage_treatment.png"
ggsave(plot = p, device = "png", width=4, height=4, dpi=300, path = output_dir, filename = filename)



# uMAPs coloured by marker intesities
selectcolumns = c("MI_MHCcII", "MI_F480",
                  "MI_PDL1", "MI_CD86", "MI_CD68",
                  "MI_CD206", "MI_pS6", "MI_CD11c")

for (i in selectcolumns){
  print(i)
  p = ggplot(cd[which(cd$treatment == "Vehicle"),], aes(x=umap1, y=umap2, color = get(i))) +
    geom_point(size = 0.3) +
    theme_classic() +
    theme(plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.line =  element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position =  "none") +
    scale_color_distiller(palette = "Spectral", limits = c(0,3))
  
  filename = paste("uMAP_MF_Vehicle_", i, ".png", sep = "")
  ggsave(plot = p, device = "png", width=4, height=4, dpi=300, path = output_dir, filename = filename)
}


# Figure 4a
# uMAP coloured by cluster
cd = cd[which(cd$cluster2 %in% c("05A", "08", "11", "26")),]
cluster_cols = c("yellow", "#00BA38", "#619CFF", "orange") # "#F8766D")
p =ggplot(cd[sample(nrow(cd), nrow(cd)),], aes(x=umap1, y=umap2, colour = as.factor(cluster2))) + 
  geom_point(size = 0.3) +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line =  element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank()) +
  scale_color_manual(values = cluster_cols)
p
filename = "uMAP_macrophage_clusters.png"
ggsave(plot = p, device = "png", width=4, height=4, dpi=300, path = output_dir, filename = filename)



# Assign macrophage subsets based on uMAP

cd$umap1 = macrophage_uMAP$umap1
cd$umap2 = macrophage_uMAP$umap2

cd[which(cd$umap1 >-1 & cd$umap2 >2.8 & cd$umap1 <5), "clustername3"] = "Macrophages type 1"
cd[which(cd$clustername2 != "Macrophages type 1" & cd$cluster != 23), "clustername3"] = "Macrophages type 2"

celldata$clustername3 = celldata$clustername2
celldata[which(celldata$X %in% cd[which(cd$clustername2 == "Macrophages type 1"),"X"]), "clustername3"] = "Macrophages type 1"
celldata[which(celldata$X %in% cd[which(cd$clustername2 == "Macrophages type 2"),"X"]), "clustername3"] = "Macrophages type 2"

unique(celldata$clustername3)


 

# Figure 4c
# Stacked bar graph of domain distribution in macrophage types 
domain_cols <- c("Normal" = "red", "Interface" = "green", "Tumour" = "purple", "Structural" = "cyan")

cd_prop_d = data.frame()
n = 0

for (r in sort(unique(celldata$ROI_name3))){
  print(r)
  treatment = as.character(unique(celldata[which(celldata$ROI_name3 == r), "treatment"]))
  MouseID = as.character(unique(celldata[which(celldata$ROI_name3 == r), "MouseID"]))
  for (cl in sort(unique(celldata$clustername3))){
    print(cl)
    cd = celldata[which(celldata$ROI_name3 == r & celldata$domain2 != "n/a" & celldata$clustername3 == cl),]
    a = nrow(cd)
    for (d in c("Normal", "Interface", "Tumour", "Structural")){
      print(d)
      n = n+1
      cd = celldata[which(celldata$ROI_name3 == r & celldata$domain2 != "n/a" & 
                            celldata$clustername3 == cl & celldata$domain2 == d),]
      b = nrow(cd)
      cd_prop_d[n,"treatment"] = treatment
      cd_prop_d[n,"MouseID"] = MouseID
      cd_prop_d[n,"ROI_name3"] = r
      cd_prop_d[n,"clustername3"] = cl
      cd_prop_d[n,"domain2"] = d
      cd_prop_d[n,"prop"] = b/a
      
    }
  }  
}
cd_prop_d$domain2 = factor(cd_prop_d$domain2, levels = c("Normal", "Interface", "Tumour", "Structural"))

#  With error bars per ROI
p = ggbarplot(cd_prop_d[which(cd_prop_d$clustername3 %in% c("Macrophages type 1", "Macrophages type 2")),], 
              x = "treatment", y = "prop", facet.by = "clustername3", panel.labs = list(clustername3 = c("Type 1", "Type2")), 
                  add = "mean_se", fill = "domain2") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, len = 5)) +
  scale_x_discrete(position = "bottom", limits = c("Vehicle", "MRTX")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,size =14, hjust = 0),
        axis.text.y = element_text(size =10),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title= element_blank(),
        strip.text = element_text(size = 14),
        panel.spacing = unit(0, "mm"),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA, size = 1),
        strip.background = element_rect(size = 0.5),
        plot.background = element_blank())+
  scale_fill_manual(values = domain_cols)

p
filename = "Domain_distribution_macrophages.pdf"
ggsave(plot = p, device = "pdf", width=4.4, height=3, dpi=300, path = output_dir, filename = filename)




# Figure 4e and f
# Visualising the cell outlines

#  2 sets of images

# Set1 - Vehicle
path = "/Users/vanmalf/Documents/Publication/Methods_IMC_paper/Submission/Data_availability/"
TIFFnm = paste(path, "Supplemental_figure_5_input/all_cells_mask_BRAC3438.6f_ROI1_t1_Vehicle.tiff", sep = "")
TIFFol = paste(path, "Supplemental_figure_5_input/Cells_outline_BRAC3438.6f_ROI1_t1_Vehicle.tiff", sep = "")
cd = celldata[which(celldata$ROI_name == "BRAC3438.6f_ROI1_t1_Vehicle."),]

# Set2 - MRTX
TIFFnm = paste(path, "Supplemental_figure_5_input/all_cells_mask_BRAC4002.3c_ROI1_t1_MRTX.tiff", sep = "")
TIFFol = paste(path, "Supplemental_figure_5_input/Cells_outline_BRAC4002.3c_ROI1_t1_MRTX.tiff", sep = "")
cd = celldata[which(celldata$ROI_name == "BRAC4002.3c_ROI1_t1_MRTX."),]

# For each plot, reload the tiff files and start script from here
# Read TIFFs of the relevant set
TIFF = readTIFF(TIFFnm)
TIFF2 = readTIFF(TIFFol)
to_long_df <- function(TIFF) {
  names(TIFF)<- c(1:length(TIFF))
  TIFF <- melt(TIFF, id = c(row(TIFF), names(TIFF)))
  names(TIFF)[1] <- "y"
  names(TIFF)[2] <- "x"
  TIFF
}  
Outline = TIFF2
ROI1 = TIFF
ROI1 = to_long_df(ROI1)
Outline = to_long_df(Outline)
ROI1$unique_px_ID = c(1:nrow(ROI1))
ROI2 = ROI1
n=2

# Now choose whether to plot images for Figure 4e or 4f

# For Figure 4e
clusters = c("Macrophages type 1", "CD4 T cells", "Dendritic cells other","Dendritic cells cDC1" )
labels = c("","", "Macrophages type 1", "CD4 T cells", "Dendritic cells other","Dendritic cells cDC1")
cols = c("Black", "White",  "Green", "Cyan", "Magenta", "Purple")

# For Figure 4f
clusters = c("Fibroblasts", "Macrophages type 2")
labels = c("","", "Fibroblasts", "Macrophages type 2")
cols = c("Black", "White",  "Yellow", "Red")

# Continue here for both 4e and 4f
for (cl in clusters){
  cluster_xy = cd[which(cd$clustername3 == cl),c("ObjectNumber","Location_Center_X","Location_Center_Y")]
  cluster_xy$Location_Center_X = round(cluster_xy$Location_Center_X)
  cluster_xy$Location_Center_Y = round(cluster_xy$Location_Center_Y)
  names(cluster_xy) = c("ObjectNumber","x","y")
  colours_in_mask <- inner_join(ROI1, cluster_xy[,c("x","y")]) 
  min = min(unique(colours_in_mask$value))
  colours_in_mask = colours_in_mask[-(which(colours_in_mask$value == min)),]
  ROI2[which(ROI1$value %in% colours_in_mask$value), "value"] = n
  n = n+1
}

background = ROI2[which(ROI2$value < 1),"unique_px_ID"]
ROI2[which(Outline$value == 1 ),"value"] = 1
ROI2[which(ROI2$unique_px_ID %in% background),"value"] = 0


p = ggplot(ROI2,aes(x=x,y=-y, fill=as.factor(value))) +
  geom_raster() +
  theme_void() +
  theme(legend.title=element_blank(),
        legend.text = element_text(colour = "black", size = 14),
        legend.justification = c(0, 1)) +
  scale_fill_manual(values = alpha(cols,  1), labels = labels)
p

filename = "Cell outlines clusters.pdf"
ggsave(plot = p, device = "pdf", width=6, height=5, dpi=300, path = output_dir, filename = filename)



###########################################################################################
# Figure 5 - T cells
###########################################################################################

# For Figure 5a
# Run uMAP on T cell metacluster, or skip this step and load the dataset with uMAP coordinates

selectcolumns = c("MI_CD3", "MI_CD4",
                  "MI_CD8", "MI_TCRgd",
                  "MI_Foxp3","MI_PD1", "MI_CD103")

cd = celldata[which(celldata$cluster2 %in% c("25","06A", "06B")),]
umap_results = umap(cd[,selectcolumns], verbose = TRUE, n_neighbors = 10)
cd$umap1 = umap_results$layout[,1]
cd$umap2 = umap_results$layout[,2]
Tcell_uMAP = cd

# Load T cell dataset including uMAP coordinates
filenm = "Figures_input/Tcell_umap.csv"
Tcell_uMAP = read.csv(filenm)
cd = Tcell_uMAP

# uMAP coloured by cluster
cluster_cols = c("#00BA38", "orange", "blue")
p =ggplot(cd[sample(nrow(cd), nrow(cd)),], aes(x=umap1, y=umap2, colour = as.factor(cluster2))) + 
  geom_point(size = 0.3) +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line =  element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_color_manual(values = cluster_cols) +
  guides(colour = guide_legend(override.aes = list(size=5)))
p
filename = "uMAP_Tcell_clusters.png"
ggsave(plot = p, device = "png", width=4, height=4, dpi=300, path = output_dir, filename = filename)


# uMAP coloured by domain
domain_cols <- c("Normal" = "red", "Interface" = "green", "Tumour" = "purple", "Structural" = "cyan")
p = ggplot(cd[which(cd$domain %in% c( "Normal", "Tumour", "Interface")),], 
           aes(x=umap1, y=umap2, color = factor(domain, levels=c("Normal", "Interface", "Tumour")))) +
  geom_point(size = 0.3) +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line =  element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_color_manual(values = domain_cols) +
  guides(colour = guide_legend(override.aes = list(size=5)))
p
filename = "uMAP_Tcells_domain.png"
ggsave(plot = p, device = "png", width=4, height=4, dpi=300, path = output_dir, filename = filename)

# uMAP coloured by treatment
p = ggplot(cd[sample(nrow(cd), nrow(cd)),], 
           aes(x=umap1, y=umap2, color = factor(treatment, levels=c("Vehicle", "MRTX")))) +
  geom_point(size = 0.3) +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line =  element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=5)))
p
filename = "uMAP_Tcells_treatment.png"
ggsave(plot = p, device = "png", width=4, height=4, dpi=300, path = output_dir, filename = filename)


selectcolumns = c("MI_CD3", "MI_CD4",
                  "MI_CD8", "MI_TCRgd",
                  "MI_Foxp3","MI_PD1", "MI_CD103")

# uMAPs coloured by markers, Vehicle
for (i in selectcolumns){
  print(i)
  p = ggplot(cd[which(cd$treatment == "Vehicle"),], aes(x=umap1, y=umap2, color = get(i))) +
    geom_point(size = 0.3) +
    theme_classic() +
    theme(plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.line =  element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          legend.title = element_blank()) +
    scale_color_distiller(palette = "Spectral", limits = c(0,3))
  
  filename = paste("uMAP_Tcells_Vehicle_", i, ".png", sep = "")
  ggsave(plot = p, device = "png", width=4, height=4, dpi=300, path = output_dir, filename = filename)
}

# uMAPs coloured by markers, MRTX
for (i in selectcolumns){
  print(i)
  p = ggplot(cd[which(cd$treatment == "MRTX"),], aes(x=umap1, y=umap2, color = get(i))) +
    geom_point(size = 0.3) +
    theme_classic() +
    theme(plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.line =  element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          legend.title = element_blank()) +
    scale_color_distiller(palette = "Spectral", limits = c(0,3))
  
  filename = paste("uMAP_Tcells_MRTX_", i, ".png", sep = "")
  ggsave(plot = p, device = "png", width=4, height=4, dpi=300, path = output_dir, filename = filename)
}








# Figure 5b
#  Histograms to compare expression of PD1 on CD8 T cells across domains or treatments
for (t in c("CD8 T cells", "CD4 T cells", "Regulatory T cells")){
  p = ggplot(celldata[which(celldata$clustername3 == t),], aes(x = MI_PD1, colour = factor(treatment, levels = c("Vehicle", "MRTX")))) + # 
    geom_density(aes(y = stat(count / sum(count)))) +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 14)) +
    scale_x_continuous(limits = c(0, 2)) +
    ylab("Normalised counts") +
    xlab(paste("Mean intensity PD1 on ", t, sep = ""))+
    guides(colour = guide_legend(override.aes = list(size=1)))
  filename = paste(t, "_PD1_treatment_histogram.pdf", sep = "")
  ggsave(plot = p, device = "pdf", width=5.2, height=3, dpi=300, path = output_dir, filename = filename)
}

for (t in c("CD8 T cells", "CD4 T cells", "Regulatory T cells")){
  p = ggplot(celldata[which(celldata$clustername3 == t & 
                              celldata$treatment == "MRTX" &
                              celldata$domain2 != "n/a"),], aes(x = MI_PD1, colour = factor(domain2, levels=c("Normal", "Interface", "Tumour", "Structural")))) + # 
    geom_density(aes(y = stat(count / sum(count)))) +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 14)) +
    scale_x_continuous(limits = c(0, 2)) +
    ylab("Normalised counts") +
    xlab(paste("Mean intensity PD1 on", t, "(MRTX)", sep = " ")) +
    scale_color_manual(values = domain_cols) +
    guides(colour = guide_legend(override.aes = list(size=1)))
  filename =  paste(t, "_PD1_domains_histogram.pdf", sep = "")
  ggsave(plot = p, device = "pdf", width=5.2, height=3, dpi=300, path = output_dir, filename = filename)
}
  


# Figure 5c
# Calculate the distance to the nearest cell of a cell type, 
# Or skip this step and load the dataset with distances
cd = celldata[,c("Location_Center_X", "Location_Center_Y", "clustername3", "ROI_name", "treatment")]
names(cd)
# For cell n from cd, calculate distance to nearest cell of clusterID
cl_v = sort(unique(cd$clustername3))
cd_result = data.frame()
for (ROI in unique(cd$ROI_name)){
  cd_ROI = cd[which(cd$ROI_name == ROI),]
  for (clusterID in cl_v) {
    cd_cl = cd_ROI[which(cd_ROI$clustername3 == clusterID, arr.ind = FALSE),]
    for (n in 1:nrow(cd_ROI)) {
      cd_ROI$temp[n] = sqrt(min((cd_cl[,1]- as.numeric(cd_ROI[n,1]))^2 + (cd_cl[,2]- as.numeric(cd_ROI[n,2]))^2))
    }
    column_name = paste("dist_cluster_", clusterID, sep = "")
    names(cd_ROI)[names(cd_ROI) == "temp"] = column_name
  }
  cd_result = rbind(cd_result, cd_ROI) 
}

celldata = cbind(celldata, cd_result[,grep("dist_cluster", names(cd_result))])

# Load celldata dataset including cluster distances
filenm = "~path/Figures_input/input_celldata_complete.csv"
celldata = read.csv(filenm)

# For Figure 5c
cluster_order = c("CD4 T cells","Regulatory T cells","B cells",
                  "Dendritic cells cDC1", "Dendritic cells other", "Macrophages type 1", "Macrophages type 2","Neutrophils",
                  "Endothelium","Epithelium","Fibroblasts","Tumour")
p = ggplot(celldata[which(!celldata$clustername3 %in% c("Macrophages", "Unclassified", "CD8 T cells", "NK cells")),], 
           aes(x = clustername3, y=`dist_clustername3_CD8.T.cells`, fill = factor(treatment, levels = c("Vehicle", "MRTX")))) + 
  
   
    geom_boxplot() +
  scale_color_gradientn(colours = rainbow(5)) +
  scale_x_discrete(limits = cluster_order) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right") +
  ylab("Distance to nearest CD8 T cell (px)") +
  scale_y_continuous(trans='log2') 
p
filename = "distance_metacluster2_to_CD8.pdf"
ggsave(plot = p, device = "pdf", width=9, height=5, dpi=300, path = output_dir, filename = filename)


# Figure 5d
# Heatmap of T cell neighbours
unique_clusters =  c("CD4 T cells","Regulatory T cells","CD8 T cells","B cells",#"NK cells",
                     "Dendritic cells cDC1", "Dendritic cells other", "Macrophages type 1", "Macrophages type 2","Neutrophils",
                     "Endothelium","Epithelium","Fibroblasts","Tumour")

bins = list(c(0,20), c(20,40), c(40,60), c(60,80), c(80,100))
my_palette = colorRampPalette(c("yellow", "purple"))(n = 199)
for (celltype in c("CD4 T cells", "CD8 T cells", "Regulatory T cells")){
  dist = paste("dist_clustername3_", celltype, sep = "")
  for (domain in c("Normal", "Interface", "Tumour")){
    Neighbours_matrix = data.frame()
    for (tr in c("Vehicle", "MRTX")){
      for (bin in bins){
        print(bin)
        for (cl in unique_clusters){
         Neighbours_matrix[cl,paste(tr, " ", bin[1], "-", bin[2], "px", sep = "")] =  tally(celldata[which(celldata[,dist]  > bin[1] & celldata[,dist] < bin[2] &
                                                                                                              celldata$treatment == tr & celldata$clustername3 == cl 
                                                                                                            & celldata$domain == domain
          ),])
        }
        
      }
    }
    Neighbours_matrix = as.matrix(Neighbours_matrix)
    pdf(file=paste(output_dir, celltype, "_", domain, "_Neighbours_matrix.pdf", sep = ""), width=10, height=15)
    heatmap.2(Neighbours_matrix, 
              col= my_palette, scale = "col",Rowv = NA,Colv = NA, cexRow = 2.5, cexCol = 2.5,
              density.info="none",  # turns off density plot inside color legend
              trace="none", adjCol = c(1,0.5), 
              colsep = c(5), sepwidth = c(0.05,0),
              key = TRUE,  keysize=0.6,
              # lmat = rbind(c(3,4),c(2,1)), #lwid = c(1.5,4), lhei = c(1.5,4,1), #This moves the key
              margins =c(25,25)
    )
    dev.off()
  }
}



##################

# Figure 6


# Figure 6a, PCA variables

# PCA Principal component analysis, per mouse
selectcolumns = c("ROI_name3", "treatment", "MouseID", "MI_CD45", "MI_aSMA", "MI_MHCcII", "MI_Vimentin", "MI_PECAM", "MI_F480", "MI_CD68", "MI_EPCAM", "MI_CD44",
                  "MI_LY6G", "MI_CD3", "MI_PDL1", "MI_CD103", "MI_Foxp3", "MI_TCRgd", "MI_PVR", "MI_CD86",
                  "MI_CD8", "MI_CD206", "MI_pS6", "MI_CD4", "MI_casp3", "MI_Ki67", "MI_B220", "MI_CD11c") #"CD11b",CD24
cd = celldata[,selectcolumns]
PCA_input = data.frame()
MI_Ch = grep("MI_", selectcolumns, value = TRUE)
for (n in MI_Ch){
  print(n)
  for (rn in unique(cd$MouseID)){
    md = cd[which(cd$MouseID == rn),]
    PCA_input[rn,substr(n, start = 4, stop = nchar(n))] = mean(md[[n]])
    PCA_input[rn,"MouseID"] = unique(md$MouseID)
    PCA_input[rn,"treatment"] = unique(md$treatment)
  }
}
MI_Ch = sapply(X = MI_Ch, FUN = substr, start = 4, stop = 20)
PCA_results = prcomp(PCA_input[,MI_Ch])

p = fviz_pca_var(PCA_results,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
)

filename = "PCA_variables.pdf"
ggsave(plot = p, device = "pdf", width=5, height=4, dpi=300, path = output_dir, filename = filename)



# Figure 6b
stat_summary_ROI = celldata %>%
  group_by(ROI_name3, MouseID, treatment) %>%
  summarise(MI_Vimentin = mean(MI_Vimentin))

p = ggplot(stat_summary_ROI, aes(x = factor(treatment, levels = c("Vehicle", "MRTX")), y = MI_Vimentin)) +  
  geom_dotplot(binaxis='y', stackdir='center', stackgroups = TRUE, binpositions="all", aes(fill = factor(MouseID))) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="crossbar", color="grey", width=0.2) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        legend.text = element_text(size = 14))  +
  ylab("Vimentin") 
p
filename = "Vimentin in ROIs.pdf"
ggsave(plot = p, device = "pdf", width=5, height=4, dpi=300, path = output_dir, filename = filename)



# Figure 6c
# Expression of vimentin across metaclusters
cluster_order = c("CD4 T cells","Regulatory T cells","CD8 T cells", "B cells",
                  "Dendritic cells cDC1", "Dendritic cells other", "Macrophages type 1", "Macrophages type 2","Neutrophils",
                  "Endothelium","Epithelium","Fibroblasts","Tumour")

p = ggplot(celldata[which(!celldata$clustername3 %in% c("Macrophages", "Unclassified", "NK cells")),], 
           aes(x = clustername3, y=MI_Vimentin, fill=factor(treatment, levels = c("Vehicle", "MRTX")))) + 
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(),
        legend.title = element_blank())  +
  scale_x_discrete(limits = cluster_order) +
  ylab("Vimentin") +
  scale_y_continuous(trans='log2') +
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="black", position=position_dodge(0.9)) 
p
filename = "Vimentin_metacluster.pdf"
ggsave(plot = p, device = "pdf", width=7, height=3.5, dpi=300, path = output_dir, filename = filename)


# Figure 6e
# Correlating vimentin with proportions of cell types 
cd_prop = data.frame()
for (r in unique(celldata$ROI_name3)){
  print(r)
  cd = celldata[which(celldata$ROI_name3 == r),] 
  cd_prop[r, "treatment"] = unique(cd$treatment)
  cd_prop[r, "ROI_name3"] = unique(cd$ROI_name3)
  cd_prop[r, "MouseID"] = unique(cd$MouseID)
  cd_prop[r, "Vimentin"] = mean(cd$MI_Vimentin)
  total_cc = nrow(cd)
  for (cl in unique(celldata$clustername3)){
    #     print(cl)
    cd_prop[r, cl] = nrow(cd[which(cd$clustername3 == cl),])/total_cc
  }
}


for (i in unique(celldata$clustername3)){
  p = ggplot(cd_prop, aes(x = !!sym(i), y = `Vimentin`, colour = factor(treatment, levels = c("Vehicle", "MRTX")))) +  
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = FALSE, size = 0.5) +
    theme_classic() +
    theme(legend.title = element_blank(),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          #axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14)) +
    stat_cor()
  
  filename = paste(i, "vimentin_corrs.pdf", sep = " ")
  ggsave(plot = p, device = "pdf", width=5, height=4, dpi=300, path = output_dir, filename = filename)
}  


