# Script to generate te supplemental figures from Van Maldegem, Valand, et al. 2021


# Load packages
library("readr")
library("ggplot2")
library("dplyr")
library(gplots)
library("RColorBrewer")
library("EBImage")
library("tiff")
library("reshape")



# Set directory where files will be saved
# Set output_dir as path/outputs/
output_dir = "~path/outputs/"

# Set filenm as path/input_celldata_complete.csv"
filenm = "~path/Supplemental_figure_5_input/input_celldata_complete.csv"
# Read in dataset
celldata = read.csv(filenm)



##############################################################################################
# Supplemental Figure 4 

# Supplemental Figure 4a
# Counts in clusters coloured by treatment
p = ggplot(celldata, aes(x = clustername, fill = factor(treatment, levels = c("Vehicle", "MRTX")))) +
  geom_bar(stat = "count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,size =8, hjust = 1),
        axis.text.y = element_text(size =6),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        plot.title = element_text(size = 8), 
        axis.title=element_text(size=8))  +
  labs(y = "count")
p
filename = "cluster counts treatment.png"
ggsave(plot = p, device = "png", width=9, height=4, dpi=300, path = output_dir, filename = filename)


# Supplemental Figure 4b
# Counts in metaclusters coloured by treatment
p = ggplot(celldata, aes(x = metacluster, fill = factor(treatment, levels = c("Vehicle", "MRTX")))) +
  geom_bar(stat = "count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,size =8, hjust = 1),
        axis.text.y = element_text(size =6),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        plot.title = element_text(size = 8), 
        axis.title=element_text(size=8))  +
  labs(y = "count")
p
filename = "metacluster counts treatment.png"
ggsave(plot = p, device = "png", width=4.5, height=3, dpi=300, path = output_dir, filename = filename)

# Supplemental Figure 4c
# Principal component analysis
selectcolumns = c("ROI_name", "treatment", "MI_CD45", "MI_aSMA", "MI_MHCcII", "MI_Vimentin", "MI_PECAM", "MI_F480", "MI_CD68", "MI_EPCAM", "MI_CD44",
                  "MI_LY6G", "MI_CD3", "MI_PDL1", "MI_CD103", "MI_Foxp3", "MI_TCRgd", "MI_PVR", "MI_CD86",
                  "MI_CD8", "MI_CD206", "MI_pS6", "MI_CD4", "MI_casp3", "MI_Ki67", "MI_B220", "MI_CD11c") #"CD11b",CD24
cd = celldata[,selectcolumns]
PCA_input = data.frame()
MI_Ch = grep("MI_", selectcolumns, value = TRUE)
for (n in MI_Ch){
  print(n)
  for (rn in unique(cd$ROI_name)){
    md = cd[which(cd$ROI_name == rn),]
    PCA_input[rn,n] = mean(md[[n]])
    PCA_input[rn,"ROI_name"] = unique(md$ROI_name)
    PCA_input[rn,"treatment"] = unique(md$treatment)
  }
}

PCA_results = prcomp(t(PCA_input[,MI_Ch]))
PCA_results2 = as.data.frame(PCA_results$rotation)
PCA_results2$treatment = PCA_input$treatment
PCA_results2$ROI_name = PCA_input$ROI_name

PCA_results2 = cbind(PCA_results2, PCA_input[,MI_Ch])
PCA_results2

p = ggplot(PCA_results2, aes(x = PC1, y = PC2, colour = factor(treatment, levels = c("Vehicle", "MRTX")))) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  geom_point()
p
filename = "PCA plot.png"
ggsave(plot = p, device = "png", width=4, height=3, dpi=300, path = output_dir, filename = filename)

# Supplemental Figure 4d
# Metaclusters plotted using X-Y coordinates
p = ggplot(celldata, aes(x= Location_Center_X, y = -Location_Center_Y, color = as.factor(metacluster))) +
  geom_point(size = 0.2) +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_rect(size = 0.5),
        panel.spacing = unit(0, "mm"),
        plot.background = element_blank()) +
  facet_wrap(~ ROI_name3, nrow = 2, ncol = 6)
p

##############################################################################################
# Supplemental Figure 5 - Macrophages and T cells


# Load macrophage uMAP including uMAP coordinates
filenm = "~path/Figure_3_input/macrophage_umap.csv"
cd = read.csv(filenm)


# Supplemental Figure 5a
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

# Supplemental Figure 5b and c
# Visualising the cell outlines

#  2 sets of images

# Set1 - Vehicle
TIFFnm = "~path/Supplemental_figure_5_input/all_cells_mask_BRAC3438.6f_ROI1_t1_Vehicle.tiff"
TIFFol = "~path/Supplemental_figure_5_input/Cells_outline_BRAC3438.6f_ROI1_t1_Vehicle.tiff"
cd = celldata[which(celldata$ROI_name == "BRAC3438.6f_ROI1_t1_Vehicle."),]

# Set2 - MRTX
TIFFnm = "~path/Supplemental_figure_5_input/all_cells_mask_BRAC4002.3c_ROI1_t1_MRTX.tiff"
TIFFol = "~path/Supplemental_figure_5_input/Cells_outline_BRAC4002.3c_ROI1_t1_MRTX.tiff"
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

# Now choose whether to plot images for Sup Figure 5b or 5c

# For Supplemental Figure 5b
clusters = c(06,03,25,11)
labels = c("","", "CD4 T cells", "CD103 DCs", "CD8 T cells", "Macrophages, CD68+ CD11c+")
for (cl in clusters){
  cluster_xy = cd[which(cd$cluster == cl),c("ObjectNumber","Location_Center_X","Location_Center_Y")]
  cluster_xy$Location_Center_X = round(cluster_xy$Location_Center_X)
  cluster_xy$Location_Center_Y = round(cluster_xy$Location_Center_Y)
  names(cluster_xy) = c("ObjectNumber","x","y")
  colours_in_mask <- inner_join(ROI1, cluster_xy[,c("x","y")]) 
  min = min(unique(colours_in_mask$value))
  colours_in_mask = colours_in_mask[-(which(colours_in_mask$value == min)),]
  ROI2[which(ROI1$value %in% colours_in_mask$value), "value"] = n
  n = n+1
}

# For Supplemental Figure 5c
metaclusters = c("Fibroblasts", "Myeloid: macrophages type 2")
labels = c("","", "fibroblasts", "Macrophages type 2")
for (cl in metaclusters){
  cluster_xy = cd[which(cd$metacluster2 == cl),c("ObjectNumber","Location_Center_X","Location_Center_Y")]
  cluster_xy$Location_Center_X = round(cluster_xy$Location_Center_X)
  cluster_xy$Location_Center_Y = round(cluster_xy$Location_Center_Y)
  names(cluster_xy) = c("ObjectNumber","x","y")
  colours_in_mask <- inner_join(ROI1, cluster_xy[,c("x","y")]) 
  min = min(unique(colours_in_mask$value))
  colours_in_mask = colours_in_mask[-(which(colours_in_mask$value == min)),]
  ROI2[which(ROI1$value %in% colours_in_mask$value), "value"] = n
  n = n+1
}

# Continue here for both 5b and 5c
background = ROI2[which(ROI2$value < 1),"unique_px_ID"]
ROI2[which(Outline$value == 1 ),"value"] = 1
ROI2[which(ROI2$unique_px_ID %in% background),"value"] = 0

cols = c("Black", "White",  "Cyan", "Magenta", "Yellow", "Green")
p = ggplot(ROI2,aes(x=x,y=-y, fill=as.factor(value))) +
  geom_raster() +
  theme_void() +
  theme(legend.title=element_blank(),
        legend.text = element_text(colour = "black")) +
  scale_fill_manual(values = alpha(cols,  1), labels = labels)
p

filename = "Cell outlines clusters.pdf"
ggsave(plot = p, device = "pdf", width=5.6, height=5, dpi=300, path = output_dir, filename = filename)


# Supplemental Figure 5d
# Visualise T cells onto x-y coordinates
ggplot(celldata, aes(x= Location_Center_X, y = -Location_Center_Y)) +
  geom_point(size = 0.2, color = "black") +
  geom_point(data = celldata[which(celldata$cluster %in% c(6,25)),], size = 0.2, color = "yellow") +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_rect(size = 0.5),
        panel.spacing = unit(0, "mm"),
        plot.background = element_blank()) +
  facet_wrap(~ ROI_name3, nrow = 2, ncol = 6)





# For Supplemental Figure 5e
# Run uMAP on T cell metacluster, or skip this step and load the dataset with uMAP coordinates

selectcolumns = c("MI_CD3", "MI_CD4",
                  "MI_CD8", "MI_TCRgd",
                  "MI_Foxp3","MI_PD1")

cd = celldata[which(celldata$cluster %in% c(25,6)),]
umap_results = umap(cd[,selectcolumns], verbose = TRUE, n_neighbors = 10)
cd$umap1 = umap_results$layout[,1]
cd$umap2 = umap_results$layout[,2]
Tcell_uMAP = cd

# Load T cell dataset including uMAP coordinates
filenm = "~path/Figure_3_input/Tcell_umap.csv"
cd = read.csv(filenm)


# Supplemental Figure 5e
# uMAP coloured by cluster
cluster_cols = c("#00BA38", "orange")
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
        legend.title = element_blank()) +
  scale_color_manual(values = domain_cols)
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
        legend.title = element_blank())
p
filename = "uMAP_Tcells_treatment.png"
ggsave(plot = p, device = "png", width=4, height=4, dpi=300, path = output_dir, filename = filename)


selectcolumns = c("MI_CD3", "MI_CD4",
                  "MI_CD8", "MI_TCRgd",
                  "MI_Foxp3","MI_PD1")

# uMAPs coloured by markers
for (i in selectcolumns){
  print(i)
  p = ggplot(cd, aes(x=umap1, y=umap2, color = get(i))) +
    geom_point(size = 0.3) +
    theme(plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.line =  element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.title = element_blank()) +
    scale_color_distiller(palette = "Spectral", limits = c(0,3))
  
  filename = paste("uMAP_Tcells_", i, ".png", sep = "")
  ggsave(plot = p, device = "png", width=4, height=4, dpi=300, path = output_dir, filename = filename)
}

