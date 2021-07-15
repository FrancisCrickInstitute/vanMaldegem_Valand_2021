# Script to generate te Supplementary figures from Van Maldegem, Valand, et al. 2021


# Load packages
library(tidyverse)
library("readr")
library("ggplot2")
library("dplyr")
library(gplots)
library("RColorBrewer")
library("EBImage")
library("tiff")
library("reshape")
# install.packages('Rtsne') ### installing packages and libraries for tSNE plots in R
library(Rtsne)



# Set directory where files will be saved
# Set output_dir as path/outputs/
output_dir = "output/"


# Set filenm as path/input_celldata_complete.csv"
filenm = "Figures_input/input_celldata_complete.csv"
# Read in dataset
celldata = read.csv(filenm)



##############################################################################################
# Supplementary Figure 4 

# Add as Sup Fig 4a heatmap of all clusters
# Function to generate a heatmap with markers versus clusters, labeled with clustername
# markerlist = c("CD45", "aSMA", "MHCcII", "F480", "CD68", "EPCAM", "CD44", "LY6G", "CD3", "CD103", "CD8", "CD4", "B220", "CD11c", "NKp46", "PECAM", "PVR", "Foxp3")
# all markers (27):
markerlist = grep("MI_", names(celldata), value = TRUE)
markerlist = markerlist[c(2:12, 15:27, 30:32)]
markerlist = substr(markerlist, 4, 1000L)

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
pdf(file=paste(output_dir, "_heatmap_clustername_supplemental.pdf", sep = ""), width=14, height=12)

colors = c(seq(0,0.6,length=100),seq(0.61,4,length=100))
my_palette <- colorRampPalette(c("yellow", "blue"))(n = 199)

heatmap.2(heatmap_df, col= my_palette, breaks = colors,
          scale = "none", cexRow = 1.9, cexCol = 1.9, Rowv = FALSE,
          density.info="none",  # turns off density plot inside color legend
          trace="none", keysize=0.6,
          margins =c(8,25))
dev.off()




# Supplementary Figure 4b
# Counts in clusters coloured by treatment
cluster_order = c("CD4 T cells","Regulatory T cells","CD8 T cells","B cells", "NK cells",
                  "Dendritic cells cDC1", "Dendritic cells other", "Macrophages","Neutrophils",
                  "Endothelium","Epithelium","Fibroblasts","Tumour", "Unclassified")
percentages = celldata %>% group_by(clustername2) %>% tally()
percentages = percentages[match(cluster_order, percentages$clustername2),]
percentages$percentages = percentages$n/nrow(celldata)*100
percentages$percentages = paste(as.character(round(percentages$percentages, digits = 1)), "%", sep = "")

p = ggplot(celldata, aes(x = clustername2, fill = factor(treatment, levels = c("Vehicle", "MRTX")))) +
  geom_bar(stat = "count") +
  annotate("text", x = percentages$clustername2, y = percentages$n, label = percentages$percentages, vjust = -0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,size =12, hjust = 1),
        axis.text.y = element_text(size =12),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        plot.title = element_text(size = 12), 
        axis.title=element_text(size=12),
        )  +
  scale_x_discrete(limits = cluster_order) +
  scale_y_continuous(limits = c(0,120000)) +
  labs(y = "count") +
  guides(colour = guide_legend(override.aes = list(size=5)))
  

p
filename = "cluster counts treatment.pdf"
ggsave(plot = p, device = "pdf", width=9, height=4, dpi=300, path = output_dir, filename = filename)



# Supplementary Figure 4c
# Principal component analysis, per ROI
selectcolumns = c("ROI_name3", "treatment", "MouseID", "MI_CD45", "MI_aSMA", "MI_MHCcII", "MI_Vimentin", "MI_PECAM", "MI_F480", "MI_CD68", "MI_EPCAM", "MI_CD44",
                  "MI_LY6G", "MI_CD3", "MI_PDL1", "MI_CD103", "MI_Foxp3", "MI_TCRgd", "MI_PVR", "MI_CD86",
                  "MI_CD8", "MI_CD206", "MI_pS6", "MI_CD4", "MI_casp3", "MI_Ki67", "MI_B220", "MI_CD11c") #"CD11b",CD24
cd = celldata[,selectcolumns]
PCA_input = data.frame()
MI_Ch = grep("MI_", selectcolumns, value = TRUE)
for (n in MI_Ch){
  print(n)
  for (rn in unique(cd$ROI_name3)){
    md = cd[which(cd$ROI_name3 == rn),]
    PCA_input[rn,n] = mean(md[[n]])
    PCA_input[rn,"ROI_name3"] = unique(md$ROI_name3)
    PCA_input[rn,"MouseID"] = unique(md$MouseID)
    PCA_input[rn,"treatment"] = unique(md$treatment)
  }
}

PCA_results = prcomp(PCA_input[,MI_Ch])
PCA_results_ROI = as.data.frame(PCA_results$x)
PCA_results_ROI$treatment = PCA_input$treatment
PCA_results_ROI$ROI_name3 = PCA_input$ROI_name3
PCA_results_ROI$MouseID = PCA_input$MouseID

PCA_results_ROI = cbind(PCA_results_ROI, PCA_input[,MI_Ch])
PCA_results_ROI

p = ggplot(PCA_results_ROI, aes(x = PC1, y = PC2, colour = factor(treatment, levels = c("Vehicle", "MRTX")))) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  geom_point() +
  guides(colour = guide_legend(override.aes = list(size=5)))
p
filename = "PCA plot ROI.pdf"
ggsave(plot = p, device = "pdf", width=4, height=3, dpi=300, path = output_dir, filename = filename)


# PCA Principal component analysis, per mouse
selectcolumns = c("ROI_name3", "treatment", "MouseID", "MI_CD45", "MI_aSMA", "MI_MHCcII", "MI_Vimentin", "MI_PECAM", "MI_F480", "MI_CD68", "MI_EPCAM", "MI_CD44",
                  "MI_LY6G", "MI_CD3", "MI_PDL1", "MI_CD103", "MI_Foxp3", "MI_TCRgd", "MI_PVR", "MI_CD86",
                  "MI_CD8", "MI_CD206", "MI_pS6", "MI_CD4", "MI_casp3", "MI_Ki67", "MI_B220", "MI_CD11c")
cd = celldata[,selectcolumns]
PCA_input = data.frame()
MI_Ch = grep("MI_", selectcolumns, value = TRUE)
for (n in MI_Ch){
  print(n)
  for (rn in unique(cd$MouseID)){
    md = cd[which(cd$MouseID == rn),]
    PCA_input[rn,substr(n, start = 4, stop = nchar(n))] = mean(md[[n]])
    #PCA_input[rn,"ROI_name3"] = unique(md$ROI_name3)
    PCA_input[rn,"MouseID"] = unique(md$MouseID)
    PCA_input[rn,"treatment"] = unique(md$treatment)
  }
}
MI_Ch = sapply(X = MI_Ch, FUN = substr, start = 4, stop = 20)
PCA_results = prcomp(PCA_input[,MI_Ch])

PCA_results_mouse = as.data.frame(PCA_results$x)
PCA_results_mouse$treatment = PCA_input$treatment
# PCA_results_mouse$ROI_name3 = PCA_input$ROI_name3
PCA_results_mouse$MouseID = PCA_input$MouseID

PCA_results_mouse = cbind(PCA_results_mouse, PCA_input[,MI_Ch])
PCA_results_mouse

p = ggplot(PCA_results_mouse, aes(x = PC1, y = PC2, colour = factor(treatment, levels = c("Vehicle", "MRTX")))) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  geom_point() +
  guides(colour = guide_legend(override.aes = list(size=5)))
p
filename = "PCA plot MouseID.pdf"
ggsave(plot = p, device = "pdf", width=4, height=3, dpi=300, path = output_dir, filename = filename)



# 
#tSNE
# Sup Fig 4d
# Run tSNE, or skip this step and load the dataset with tSNE coordinates
# selectcolumns = c("MI_CD45", "MI_aSMA", "MI_MHCcII", "MI_PECAM", "MI_F480", "MI_CD68", "MI_EPCAM", 
#                   "MI_CD44", "MI_LY6G", "MI_CD3", "MI_CD103", "MI_PVR",
#                   "MI_CD8", "MI_CD4", "MI_B220", "MI_CD11c", "MI_NKp46") 
selectcolumns = grep("MI_", names(celldata), value = TRUE)
selectcolumns = selectcolumns[c(2:12, 15:27, 30:32)]

tSNE_results = Rtsne(celldata[,selectcolumns], num_threads = 4, perplexity = 100, verbose = TRUE, check_duplicates = FALSE, max_iter = 1500)
celldata$tSNE1b = tSNE_results$Y[,1]
celldata$tSNE2b = tSNE_results$Y[,2]

# Load celldata including tSNE coordinates
filenm = "~path/Figures_input/input_celldata_complete.csv"
celldata = read.csv(filenm)
p = ggplot(celldata[sample(nrow(celldata), 50000), ], aes(x=tSNE1, y=tSNE2, color = factor(clustername2, levels = 
                                                                                                c("CD4 T cells","Regulatory T cells","CD8 T cells", "B cells", "NK cells",
                                                                                                  "Dendritic cells cDC1", "Dendritic cells other", "Macrophages","Neutrophils",
                                                                                                  "Endothelium","Epithelium","Fibroblasts","Tumour", "Unclassified")))) +
  geom_point(size = 0.1) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(size=5)))
p
filename = "celldata_tSNE_cluster.png"
ggsave(plot = p, device = "png", width=6.8, height=5, dpi=300, path = output_dir, filename = filename)


# Sup Fig 4d
# tSNE coloured by treatment
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = append(c("Black", "White"), gg_color_hue(14))

p = ggplot(celldata[sample(nrow(celldata), 50000), ], aes(x=tSNE1, y=tSNE2, color = factor(treatment, levels = c("Vehicle", "MRTX")))) +
  geom_point(size = 0.1) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size=5)))
p
filename = "celldata_tSNEb_treatment.png"
ggsave(plot = p, device = "png", width=6.2, height=5, dpi=300, path = output_dir, filename = filename)


# Sup Fig 4d
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
filename = "celldata_tSNEb_domains.png"
ggsave(plot = p, device = "png", width=6.4, height=5, dpi=300, path = output_dir, filename = filename)




# Supplementary Figure 4e
# Cell types plotted using X-Y coordinates
p = ggplot(celldata, aes(x= Location_Center_X, y = -Location_Center_Y, color = factor(clustername2, levels = 
                                                                                        c("CD4 T cells","Regulatory T cells","CD8 T cells", "B cells", "NK cells",
                                                                                          "Dendritic cells cDC1", "Dendritic cells other", "Macrophages","Neutrophils",
                                                                                          "Endothelium","Epithelium","Fibroblasts","Tumour", "Unclassified")))) +
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
  facet_wrap(~ ROI_name3, nrow = 2, ncol = 6) +
  guides(colour = guide_legend(override.aes = list(size=5)))
# p
filename = "XY metaclusters.pdf"
ggsave(plot = p, device = "pdf", width=10, height=4, dpi=300, path = output_dir, filename = filename)



##############################################################################################
# Supplementary Figure 6 - T cells 



####################
# Supplementary Figure 6a
# Proportions of cells in total tissue
cd_prop = data.frame()
for (r in unique(celldata$ROI_name3)){
  print(r)
  cd = celldata[which(celldata$ROI_name3 == r),]
  cd_prop[r, "treatment"] = unique(cd$treatment)
  cd_prop[r, "ROI_name3"] = unique(cd$ROI_name3)
  cd_prop[r, "MouseID"] = unique(cd$MouseID)
  total_cc = nrow(cd)
  for (cl in unique(celldata$clustername3)){
    #     print(cl)
    cd_prop[r, cl] = nrow(cd[which(cd$clustername3 == cl),])/total_cc
  }
}  
cd_prop[,1:10]

for (i in c("CD4 T cells", "CD8 T cells", "Regulatory T cells")){
  p =   ggplot(cd_prop, aes(x = factor(treatment, levels = c("Vehicle", "MRTX")), y = !!sym(i))) +  
    geom_dotplot(binaxis='y', stackdir='center', aes(fill = factor(MouseID))) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="crossbar", color="grey", width=0.2) +
    theme_classic() +
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          legend.text = element_text(size = 14))  +
    ylab(i) #+
    # stat_compare_means(method = "wilcox.test", aes(label = paste("p =", ..p.format.., sep = " ")), 
    #                    label.x = 1.5, label.y = max(cd_prop[,i])+0.005, hjust = 0.5, size = 5 )  # (method = "t.test") default = "wilcox.test" (non-parametric) +
  
  filename = paste(i, "proportions in tissue.pdf", sep = " ")
  ggsave(plot = p, device = "pdf", width=5, height=4, dpi=300, path = output_dir, filename = filename)
  
}

# Supplementary Figure 6b
# Proportions of cells in tumour domain
cd_prop = data.frame() 
for (r in unique(celldata$ROI_name3)){
  print(r)
  cd = celldata[which(celldata$ROI_name3 == r & celldata$domain2 == "Tumour"),]
  cd_prop[r, "treatment"] = unique(cd$treatment)
  cd_prop[r, "ROI_name3"] = unique(cd$ROI_name3)
  cd_prop[r, "MouseID"] = unique(cd$MouseID)
  total_cc = nrow(cd)
  for (cl in unique(celldata$clustername3)){
    #     print(cl)
    cd_prop[r, cl] = nrow(cd[which(cd$clustername3 == cl),])/total_cc
  }
}  
cd_prop[,1:10]

for (i in c("CD4 T cells", "CD8 T cells", "Regulatory T cells")){
  p =   ggplot(cd_prop, aes(x = factor(treatment, levels = c("Vehicle", "MRTX")), y = !!sym(i))) +  
    geom_dotplot(binaxis='y', stackdir='center', aes(fill = factor(MouseID))) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="crossbar", color="grey", width=0.2) +
    theme_classic() +
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          legend.text = element_text(size = 14))  +
    ylab(i) #+
    # stat_compare_means(method = "wilcox.test", aes(label = paste("p =", ..p.format.., sep = " ")), 
    #                    label.x = 1.5, label.y = max(cd_prop[,i])+0.005, hjust = 0.5, size = 5 )  # (method = "t.test") default = "wilcox.test" (non-parametric) +
  
  filename = paste(i, "proportions in tumour.pdf", sep = " ")
  ggsave(plot = p, device = "pdf", width=5, height=4, dpi=300, path = output_dir, filename = filename)
  
}



# Supplementary Figure 6c
# Visualise T cells onto x-y coordinates
p = ggplot(celldata, aes(x= Location_Center_X, y = -Location_Center_Y)) +
  geom_point(size = 0.2, color = "black") +
  geom_point(data = celldata[which(celldata$cluster2 %in% c("06B")),], size = 0.2, color = "red") +
  geom_point(data = celldata[which(celldata$cluster2 %in% c("25")),], size = 0.2, color = "blue") +
  geom_point(data = celldata[which(celldata$cluster2 %in% c("06A")),], size = 0.2, color = "yellow") +
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

filename = "T cells XY in ROIs.png"
ggsave(plot = p, device = "png", width=12, height=5, dpi=300, path = output_dir, filename = filename)


# Supplementary Figure 6d
# Boxplots of distances
# CD4 T cells
cluster_order = c("Regulatory T cells","CD8 T cells", "B cells",
                  "Dendritic cells cDC1", "Dendritic cells other", "Macrophages type 1", "Macrophages type 2","Neutrophils",
                  "Endothelium","Epithelium","Fibroblasts","Tumour")
p = ggplot(celldata[which(!celldata$clustername3 %in% c("Macrophages", "Unclassified", "CD4 T cells", "NK cells")),], 
           aes(x = clustername3, y=`dist_clustername3_CD4 T cells`, fill = factor(treatment, levels = c("Vehicle", "MRTX")))) + 
  
  
  geom_boxplot() +
  scale_color_gradientn(colours = rainbow(5)) +
  scale_x_discrete(limits = cluster_order) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right") +
  ylab("Distance to nearest CD4 T cell (px)") +
  scale_y_continuous(trans='log2') 
p
filename = "distance_metacluster2_to_CD4.pdf"
ggsave(plot = p, device = "pdf", width=9, height=5, dpi=300, path = output_dir, filename = filename)


# Supplementary Figure 6e
# Regulator T cells
cluster_order = c("CD4 T cells","CD8 T cells", "B cells",
                  "Dendritic cells cDC1", "Dendritic cells other", "Macrophages type 1", "Macrophages type 2","Neutrophils",
                  "Endothelium","Epithelium","Fibroblasts","Tumour")
p = ggplot(celldata[which(!celldata$clustername3 %in% c("Macrophages", "Unclassified", "Regulatory T cells", "NK cells")),], 
           aes(x = clustername3, y=`dist_clustername3_Regulatory T cells`, fill = factor(treatment, levels = c("Vehicle", "MRTX")))) + 
  
  
  geom_boxplot() +
  scale_color_gradientn(colours = rainbow(5)) +
  scale_x_discrete(limits = cluster_order) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right") +
  ylab("Distance to nearest Regulatory T cell (px)") +
  scale_y_continuous(trans='log2') 
p
filename = "distance_metacluster2_to_Treg.pdf"
ggsave(plot = p, device = "pdf", width=9, height=5, dpi=300, path = output_dir, filename = filename)



# ####################################################################################################
# # Linear Regression analysis for Supplementary Table 4
####################################################################################################

selectcolumns = grep("MI_", names(celldata), value = TRUE)

lm_result = data.frame()
for (i in selectcolumns){
  c = lm(get(i) ~ treatment, data = celldata)
  lm_result[i,"estimate"] = summary(c)[[4]][2]
  lm_result[i,"Std_error"] = summary(c)[[4]][4]
  lm_result[i,"t_value"] = summary(c)[[4]][6]
  lm_result[i,"p-value"] = summary(c)[[4]][8]
}

lm_result = lm_result[order(lm_result$estimate),]
lm_result
write.csv(lm_result, file = paste(output_dir, date, "lm_result.csv", sep = ""))

################################ 


# Supplementary FIGURE 7


# `summarise()` regrouping output by 'MouseID', 'ROI_name3', 'treatment', 'domain2' 

stat_summary_means = celldata %>%
      group_by(MouseID, ROI_name3, treatment, domain, clustername3) %>%
      summarise(
          mean_vimentin = mean(MI_Vimentin, na.rm = TRUE),
          mean_distance_to_CD8 = mean(dist_clustername3_CD8.T.cells,  na.rm = TRUE)
        )
filenm = paste(output_dir, "means_table_domains",  ".csv", sep = "")
write.csv(stat_summary_means, filenm)

p = ggplot(stat_summary_means[which(stat_summary_means$clustername3 == "Endothelium" & stat_summary_means$domain != "n/a"),], aes(x=factor(domain, levels = c("Normal", "Interface", "Tumour", "Structural")), y=mean_vimentin, col=factor(treatment, levels = c("Vehicle", "MRTX")))) +
  geom_boxplot() +
  geom_point( position=position_dodge(0.9)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())
filename = "Vimentin_endothelium.pdf"
ggsave(plot = p, device = "pdf", width=5, height=3, dpi=300, path = output_dir, filename = filename)



