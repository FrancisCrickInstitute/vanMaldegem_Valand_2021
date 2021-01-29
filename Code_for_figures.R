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

# Download the slingshot package (requires R v3.6.2 or higher)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# The following initializes usage of Bioc devel
# BiocManager::install(version='devel')
# BiocManager::install("slingshot")
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("destiny")
library(SingleCellExperiment)
library(slingshot, quietly = TRUE)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(gam)



# Set directory where files will be saved
# Set output_dir as path/outputs/
output_dir = "~path/outputs/"

# Set filenm as path/input_celldata_complete.csv"
filenm = "~path/Figure_3_input/input_celldata_complete.csv"
# Read in dataset
celldata = read.csv(filenm)


##############################################################################################
# Visualising the clustering results with a heatmap
#
# Figure 3a

# Function to generate a heatmap with markers versus clusters, labeled with clustername
# lineage markers only (17):
markerlist = c("CD45", "aSMA", "MHCcII", "F480", "CD68", "EPCAM", "CD44", "LY6G", "CD3", "CD103", "CD8", "CD4", "B220", "CD11c", "NKp46", "PECAM", "PVR")
heatmap_markers_clusters = function(celldata, markerlist){
  heatmap_df = data.frame()
  unique_clusters = sort(unique(celldata$clustername))
  for (cl in unique_clusters){
    md = celldata[which(celldata$clustername == cl),]
    cl2 = unique(md$clustername)
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
heatmap_df = heatmap_markers_clusters(celldata, markerlist)
pdf(file=paste(output_dir, "_heatmap_clustername.pdf", sep = ""), width=14, height=12)
heatmap.2(heatmap_df, col= colorRampPalette(brewer.pal(8, "Blues"))(25), scale = "column", cexRow = 1.9, cexCol = 1.9,
          density.info="none",  # turns off density plot inside color legend
          trace="none", keysize=0.6,
          margins =c(8,25))
dev.off()




##############################################################################################
# Figure 3 - Exploring clusters with heatmap and tSNE


# Figure 3a
# Run tSNE, or skip this step and load the dataset with tSNE coordinates
selectcolumns = c("MI_CD45", "MI_aSMA", "MI_MHCcII", "MI_PECAM", "MI_F480", "MI_CD68", "MI_EPCAM", 
                  "MI_CD44", "MI_LY6G", "MI_CD3", "MI_CD103", "MI_PVR",
                  "MI_CD8", "MI_CD4", "MI_B220", "MI_CD11c", "MI_NKp46") 
tSNE_results = Rtsne(celldata[,selectcolumns], num_threads = 4, perplexity = 100, verbose = TRUE, check_duplicates = FALSE, max_iter = 1500)
celldata$tSNE1 = tSNE_results$Y[,1]
celldata$tSNE2 = tSNE_results$Y[,2]

# Load celldata including tSNE coordinates
filenm = "~path/Figure_3_input/input_celldata_complete.csv"
celldata = read.csv(filenm)


# For figure 3b
# Plot tSNE coloured by Phenograph cluster2 & clustername
points_in_cluster_centers = data.frame()
for (cl in unique(celldata$cluster2)){
  x = median(celldata$tSNE1[which(celldata$cluster2 == cl)])
  y = median(celldata$tSNE2[which(celldata$cluster2 == cl)])
  points_in_cluster_centers[cl, "x"] = x
  points_in_cluster_centers[cl, "y"] = y
  points_in_cluster_centers[cl, "cl"] = cl
}
p = ggplot(celldata[sample(nrow(celldata), 50000), ], aes(x=tSNE1, y=tSNE2, color = as.factor(clustername))) +
  geom_point(size = 0.1) +
  geom_text(data = points_in_cluster_centers, stat = "identity", mapping = aes(x = x, y = y, label = cl, size = 12), colour = "black") +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(size=5)))
p
filename = "celldata_tSNE_cluster.pdf"
ggsave(plot = p, device = "pdf", width=8, height=5, dpi=300, path = output_dir, filename = filename)


# Figure 3c
# tSNE coloured by metacluster
p = ggplot(celldata[sample(nrow(celldata), 50000), ], aes(x=tSNE1, y=tSNE2, color = as.factor(metacluster))) +
  geom_point(size = 0.1) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(size=5)))
p
filename = "celldata_tSNE_metacluster.pdf"
ggsave(plot = p, device = "pdf", width=6.8, height=5, dpi=300, path = output_dir, filename = filename)


# Figure 3d
# tSNE coloured by treatment
p = ggplot(celldata[sample(nrow(celldata), 50000), ], aes(x=tSNE1, y=tSNE2, color = factor(treatment, levels = c("Vehicle", "MRTX")))) +
  geom_point(size = 0.1) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(size=5)))
p
filename = "celldata_tSNE_treatment.pdf"
ggsave(plot = p, device = "pdf", width=6.8, height=5, dpi=300, path = output_dir, filename = filename)


# Figure 3e
# tSNE coloured by domain2
domain_cols <- c("Normal" = "red", "Interface" = "green", "Tumour" = "purple", "Structural" = "cyan")
cd = celldata[which(celldata$domain2 != "n/a"),] # Some cells were not assigned a domain, they are excluded from this analysis
p = ggplot(cd[sample(nrow(cd), 50000), ], aes(x=tSNE1, y=tSNE2, color = factor(domain2, levels=c("Normal", "Interface", "Tumour", "Structural")))) +
  geom_point(size = 0.1) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(values = domain_cols)
p
filename = "celldata_tSNE_domains.pdf"
ggsave(plot = p, device = "pdf", width=6.8, height=5, dpi=300, path = output_dir, filename = filename)


# Figure 3f
# Stacked bargraph of the metaclusters coloured by domain
# This makes two bar plots, one for Vehicle, and one for MRTX treatment.
domain_cols <- c("Normal" = "red", "Interface" = "green", "Tumour" = "purple", "Structural" = "cyan")
p = ggplot(celldata[which(celldata$treatment == "Vehicle" & celldata$domain2 != "n/a"),], aes(x = metacluster, fill = factor(domain2, levels=c("n/a", "Normal", "Interface", "Tumour", "Structural")))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,size =12, hjust = 1), 
        plot.title = element_text(size = 6), 
        legend.title = element_blank(),
        axis.title=element_blank())  +
  scale_fill_manual(values = domain_cols)
p
filename = "Stacked_metaclusters_by_domain_Vehicle.png"
ggsave(plot = p, device = "png", width=6, height=3.3, dpi=300, path = output_dir, filename = filename)

p = ggplot(celldata[which(celldata$treatment == "MRTX" & celldata$domain2 != "n/a"),], aes(x = metacluster, fill = factor(domain2, levels=c("n/a", "Normal", "Interface", "Tumour", "Structural")))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,size =12, hjust = 1), 
        plot.title = element_text(size = 6), 
        legend.title = element_blank(),
        axis.title=element_blank())  +
  scale_fill_manual(values = domain_cols)
p
filename = "Stacked_metaclusters_by_domain_MRTX.png"
ggsave(plot = p, device = "png", width=6, height=3.3, dpi=300, path = output_dir, filename = filename)


# Figure 3g
# Proportion of metaclusters in tumour domain
p = ggplot(celldata[which(celldata$domain == "Tumour"),], 
       aes(x = factor(treatment, levels = c("Vehicle", "MRTX")), fill = as.factor(metacluster))) +
  geom_bar(position = "fill", colour = "black") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,size =14, hjust = 1), 
        plot.title = element_text(size = 6), 
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title=element_blank()) 
p
filename = "Stacked_metaclusters_by_metacluster_treatment.png"
ggsave(plot = p, device = "png", width=3.7, height=3.5, dpi=300, path = output_dir, filename = filename)



##############################################################################################
# Figure 4 - Macrophages and T cells


# Figure 4a
# Run uMAP on macrophage metacluster, or skip this step and load the dataset with uMAP coordinates

selectcolumns = c("MI_MHCcII", "MI_F480",
                  "MI_PDL1", "MI_CD86", "MI_CD68",
                  "MI_CD206", "MI_pS6", "MI_CD11c") #, "MI_LY6G", "MI_CD103", "MI_CD44")
cd = celldata[which(celldata$metacluster == "Myeloid: macrophages"),]
umap_results = umap(cd[,selectcolumns], verbose = TRUE, n_neighbors = 10)
cd$umap1 = umap_results$layout[,1]
cd$umap2 = umap_results$layout[,2]
macrophage_uMAP = cd


# Load  macrophage dataset including uMAP coordinates
filenm = "~path/Figure_3_input/macrophage_umap.csv"
cd = read.csv(filenm)
macrophage_uMAP = cd

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


# Figure 4b
# macrophage types differ in CD68 and F480 expression
metacluster_cols = c("magenta", "orange") 
p = ggplot(cd[which(cd$metacluster2 %in% c("Myeloid: macrophages type 1", "Myeloid: macrophages type 2")),],
           aes(x = MI_CD68, y = MI_F480, colour = metacluster2)) + 
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  xlab("CD68")  +
  ylab("F480") +
  scale_x_continuous(limits = c(0, 3.5)) +
  scale_y_continuous(limits = c(0, 3.5)) +
  scale_color_manual(values = metacluster_cols)
p
filename = "F480_CD68_macrophage_types.png"
ggsave(plot = p, device = "png", width=5, height=4, dpi=300, path = output_dir, filename = filename)


# Figure 4c
# Stacked bar graph of domain distribution in macrophage types 
domain_cols <- c("Normal" = "red", "Interface" = "green", "Tumour" = "purple", "Structural" = "cyan")
p = ggplot(cd[which(cd$metacluster2 %in% c("Myeloid: macrophages type 1", "Myeloid: macrophages type 2") & cd$domain2 != "n/a"),], aes(x = factor(treatment, levels = c("Vehicle", "MRTX")), fill = factor(domain2, levels=c("n/a", "Normal", "Interface", "Tumour", "Structural")))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,size =12, hjust = 1), 
        plot.title = element_text(size = 6), 
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title=element_blank(),
        strip.text = element_text(size = 12),
        panel.spacing = unit(0, "mm"),
        strip.background = element_rect(size = 0.5),
        plot.background = element_blank())  +
  scale_fill_manual(values = domain_cols) +
  facet_wrap(~metacluster2)
filename = "Domain_distribution_macrophages.png"
ggsave(plot = p, device = "png", width=4.4, height=3, dpi=300, path = output_dir, filename = filename)



#' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Figure 4d
# uMAP of Type 1 macrophages coloured by treatment, with pseudotime trajectory

#' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


umap = macrophage_uMAP[which(macrophage_uMAP$metacluster2 == "Myeloid: macrophages type 2"),]
# Determine lineages 
startClus = "26"
lin1 <-  getLineages(umap[,c("umap1","umap2")], umap$cluster2, start.clus = startClus)
lin1

# Convert the lineages into a smooth curve 
crv1 <- getCurves(lin1)
crv1

colours = c("#00BFC4", "#F8766D")

# Randomise order of objects for plotting evenness per group 
umap_sampled = umap[sample(nrow(umap), nrow(umap)),]

# UMAP with slingshot pseudotime curve, coloured by treatment
pdf(paste(output_path, "MacType2_pseudotime_startClus", startClus, "_colourByTreatment_crv1.pdf", sep = ""))
plot(umap_sampled[,c("umap1", "umap2")], col = colours[as.factor(umap_sampled$treatment)], asp = 1, pch = 16, cex = 0.4, axes = FALSE, xlab = "", ylab = "")
lines(crv1, lwd = 2, col = 'black')
dev.off()


# For Figure 4e please see separate script on neighbouRhood analysis


# Figure 4f
# Domain distribution CD4 and CD8 T cells

domain_cols <- c("Normal" = "red", "Interface" = "green", "Tumour" = "purple", "Structural" = "cyan")
p = ggplot(celldata[which(celldata$cluster %in% c(6,25)  
          &   celldata$domain %in% c( "Normal", "Tumour", "Interface")),], 
          aes(x = factor(treatment, levels = c("Vehicle", "MRTX")), 
          fill = factor(domain2, levels=c("n/a", "Normal", "Interface", "Tumour", "Structural")))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,size =12, hjust = 1), 
        plot.title = element_text(size = 6), 
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title=element_blank(),
        strip.text = element_text(size = 12),
        panel.spacing = unit(0, "mm"),
        strip.background = element_rect(size = 0.5),
        plot.background = element_blank())  +
  scale_fill_manual(values = domain_cols) +
  facet_wrap(~clustername)
p
filename = "Domain_distribution_Tcells.png"
ggsave(plot = p, device = "png", width=4.4, height=3, dpi=300, path = output_dir, filename = filename)



# Figure 4g
#  Histograms to compare expression of PD1 on CD8 T cells across domains or treatments
p = ggplot(celldata[which(celldata$cluster == 25),], aes(x = MI_PD1, colour = factor(treatment, levels = c("Vehicle", "MRTX")))) + # 
  geom_density(aes(y = stat(count / sum(count)))) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_x_continuous(limits = c(0, 2)) +
  ylab("Normalised counts") +
  xlab("Mean intensity PD1 on CD8 T cells")
filename = "PD1_CD8_treatment_histogram.pdf"
ggsave(plot = p, device = "pdf", width=4.6, height=3, dpi=300, path = output_dir, filename = filename)


p = ggplot(celldata[which(celldata$cluster == 25 & 
                            celldata$treatment == "MRTX" &
                            celldata$domain2 != "n/a"),], aes(x = MI_PD1, colour = factor(domain2, levels=c("Normal", "Interface", "Tumour", "Structural")))) + # 
  geom_density(aes(y = stat(count / sum(count)))) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  ylab("Normalised counts") +
  xlab("Mean intensity PD1 on CD8 T cells (MRTX)") +
  scale_color_manual(values = domain_cols)
filename = "PD1_CD8_domains_histogram.pdf"
ggsave(plot = p, device = "pdf", width=4.6, height=3, dpi=300, path = output_dir, filename = filename)

# Figure 4h
# Calculate the distance to the nearest cell of a cluster, 
# Or skip this step and load the dataset with distances
cd = celldata[,c("Location_Center_X", "Location_Center_Y", "cluster2", "clustername", "ROI_name", "treatment")]
names(cd)
# For cell n from cd, calculate distance to nearest cell of clusterID
cl_v = sort(unique(cd$cluster2))
cd_result = data.frame()
for (ROI in unique(cd$ROI_name)){
  cd_ROI = cd[which(cd$ROI_name == ROI),]
  for (clusterID in cl_v) {
    cd_cl = cd_ROI[which(cd_ROI$cluster2 == clusterID, arr.ind = FALSE),]
    for (n in 1:nrow(cd_ROI)) {
      cd_ROI$temp[n] = sqrt(min((cd_cl[,1]- as.numeric(cd_ROI[n,1]))^2 + (cd_cl[,2]- as.numeric(cd_ROI[n,2]))^2))
    }
    column_name = paste("dist_cluster_", clusterID, sep = "")
    names(cd_ROI)[names(cd_ROI) == "temp"] = column_name
  }
  cd_result = rbind(cd_result, cd_ROI) 
}
cd = cd_result
celldata = cbind(celldata, cd_result[,grep("dist_cluster", names(cd_result))])

# Load celldata dataset including cluster distances
filenm = "~path/Figure_3_input/input_celldata_complete.csv"
celldata = read.csv(filenm)
# write.csv(celldata, filenm)

# Distance from CD8 T cell cluster to metaclusters, separated by treatments 
# For figure 4G
p = ggplot(celldata[which(celldata$metacluster2 != "Myeloid: macrophages"),], aes(x = reorder(metacluster2, `dist_cluster_25`), y=`dist_cluster_25`, fill = factor(treatment, levels = c("Vehicle", "MRTX")))) + 
  geom_boxplot() +
  scale_color_gradientn(colours = rainbow(5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right") +
  ylab("Distance to nearest CD8 T cell") +
  scale_y_continuous(trans='log2') 
p
filename = "distance_metacluster2_to_CD8.pdf"
ggsave(plot = p, device = "pdf", width=9, height=5, dpi=300, path = output_dir, filename = filename)


# Figure 4i
# Expression of vimentin across metaclusters
p = ggplot(celldata[which(celldata$metacluster2 != "Myeloid: macrophages"),], aes(x = metacluster2, y=MI_Vimentin, fill=factor(treatment, levels = c("Vehicle", "MRTX")))) + 
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(),
        legend.title = element_blank())  +
  ylab("Vimentin") +
  scale_y_continuous(trans='log2') +
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="black", position=position_dodge(0.9)) 
p
filename = "Vimentin_metacluster.pdf"
ggsave(plot = p, device = "pdf", width=7, height=4, dpi=300, path = output_dir, filename = filename)

