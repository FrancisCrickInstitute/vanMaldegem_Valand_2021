### Simple RE model for distances etc
library(nlme)
library(tidyverse)
library(emmeans)
library("readr")
library("ggplot2")
library("dplyr")

# Set output_dir as ~working_dir/outputs/
output_dir = "output/"
filenm = "Figures_input/input_celldata_complete.csv"
celldata = read.csv(filenm)
# Make table with mean intensities
stat_summary_means = celldata[which(celldata$clustername3 %in% c("CD4 T cells","Regulatory T cells","CD8 T cells","B cells",
                                                                "Dendritic cells cDC1", "Dendritic cells other", "Macrophages type 1", "Macrophages type 2","Neutrophils",
                                                                "Endothelium","Epithelium","Fibroblasts","Tumour")),] %>%
  group_by(MouseID, ROI_name3, treatment, clustername3) %>%
  summarise(
    mean_vimentin = mean(MI_Vimentin, na.rm = TRUE),
    dCD8 = mean(dist_clustername3_CD8.T.cells,  na.rm = TRUE),
    dCD4 = mean(dist_clustername3_CD4.T.cells,  na.rm = TRUE),
    dTreg = mean(dist_clustername3_Regulatory.T.cells,  na.rm = TRUE)
  )
filenm = paste(output_dir, "means_table",  ".csv", sep = "")
write.csv(stat_summary_means, filenm)


filenm = "Stats_input/means_table.csv"
dat <- read_csv(filenm) %>%
  rename(celltype=clustername3,
         roi = ROI_name3) %>%
  mutate(treatment=relevel(factor(treatment), "Vehicle"))


cluster_order = c("CD4 T cells","Regulatory T cells","CD8 T cells", "B cells",
                  "Dendritic cells cDC1", "Dendritic cells other", "Macrophages type 1", "Macrophages type 2","Neutrophils",
                  "Endothelium","Epithelium","Fibroblasts","Tumour")


fit_vimentin <- lme(
  log(mean_vimentin) ~ treatment*celltype , random= ~1|MouseID/roi,
  data=dat)
anova(fit_vimentin)
# 'celltype' main effect significant, so different celltypes have different vimentin signals
# 'treatment' is significant, so vimentin signal differs between treatment groups
# Interaction is significant, so there different treatments have different vimentin x celltype profiles
em_vimentin <- emmeans(fit_vimentin, ~celltype + treatment, type="response")
write.csv(as.data.frame(em_vimentin), file="output/vimentin.csv", quote=FALSE, row.names=FALSE)


#If you want p-values.  
as.data.frame(emmeans(fit_vimentin, revpairwise~treatment| celltype, type="response")$contrast) %>%
  mutate(padj= p.adjust(`p.value`, method="fdr"))



# Statistics on distances to nearest CD8 T cell
dat_noCD8 <- subset(dat, celltype!="CD8 T cells")
ggplot(dat_noCD8, aes(x=celltype, y=dCD8, col=treatment)) +
  geom_point() + geom_line(aes(group=roi)) + scale_y_log10()

fit_cd8 <- lme(
  log(dCD8) ~ treatment*celltype , random= ~1|MouseID/roi,
  data=dat_noCD8)
anova(fit_cd8)
# 'celltype' main effect significant, so different celltypes are different distances from CD8
# Interaction is significant, so there different treatments have different distance x celltype profiles
em_cd8 <- emmeans(fit_cd8, ~celltype+treatment, type="response")
write.csv(as.data.frame(em_cd8), file="output/cd8.csv", quote=FALSE, row.names=FALSE)

as.data.frame(emmeans(fit_cd8, revpairwise~treatment| celltype, type="response")$contrast) %>%
  mutate(padj= p.adjust(`p.value`, method="fdr"))


# Statistics on distances to nearest CD4 T cell
dat_noCD4 <- subset(dat, celltype!="CD4 T cells")
ggplot(dat_noCD4, aes(x=celltype, y=dCD4, col=treatment)) +
  geom_point() + geom_line(aes(group=roi)) + scale_y_log10()

fit_cd4 <- lme(
  log(dCD4) ~ treatment*celltype , random= ~1|MouseID/roi,
  data=dat_noCD4)
anova(fit_cd4)
# 'celltype' main effect significant, so different celltypes are different distances from CD4
# Interaction is significant, so there different treatments have different distance x celltype profiles
as.data.frame(emmeans(fit_cd4, revpairwise~treatment| celltype, type="response")$contrast) %>%
  mutate(padj= p.adjust(`p.value`, method="fdr"))



# Statistics on distances to nearest Regulatory T cell
dat_noTreg <- subset(dat, celltype!="Regulatory T cells")
ggplot(dat_noTreg, aes(x=celltype, y=dTreg, col=treatment)) +
  geom_point() + geom_line(aes(group=roi)) + scale_y_log10()

fit_Treg <- lme(
  log(dTreg) ~ treatment*celltype , random= ~1|MouseID/roi,
  data=dat_noTreg)
anova(fit_Treg)
# 'celltype' main effect significant, so different celltypes are different distances from Regulatory T cells
# Interaction is significant, so there different treatments have different distance x celltype profiles
as.data.frame(emmeans(fit_Treg, revpairwise~treatment| celltype, type="response")$contrast) %>%
  mutate(padj= p.adjust(`p.value`, method="fdr"))

