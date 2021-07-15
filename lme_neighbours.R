library(tidyverse)
library(lme4)
library(lattice)
library(emmeans)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# # BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
# BiocManager::install("optimx")
library(optimx)

# Set output_dir as ~working_dir/outputs/
output_dir = "output/"


dat <- read_csv("Stats_input/neighbouRhood_figure4d_lme_input.csv")
dat[which(dat$domain == "Inteface"), "domain"] = "Interface"
unique(dat$domain)
names(dat)
# dat = dat[,-1]
write.csv(dat, paste(output_dir,"neighbouRhood_figure4d_lme_input.csv", sep =""))

names(dat) <- c("log2fc","obj1", "obj2", "domain", "roi", "mouse")
dat$treatment <- ifelse(grepl("_MRTX", dat$mouse), "MRTX","Vehicle")
dat[which(dat$log2fc == -Inf),"log2fc"] = NA
dat = dat[!is.na(dat$log2fc),]


dat$roi <- as.factor(dat$roi)


# Loop over domains to save plots and stats
for (dom in unique(dat$domain)){
  
  fit <- lmer(log2fc ~ obj1*obj2 + treatment + (1|mouse)+(1|mouse:roi),
              data=subset(dat, domain==dom)
  )
  em <- emmeans(fit,  ~ obj1|obj2)
  cluster_order = c("CD4 T cells","Regulatory T cells","CD8 T cells","B cells", 
                    "Dendritic cells cDC1", "Dendritic cells other", "Macrophages type 1", "Macrophages type 2","Neutrophils",
                    "Endothelium","Epithelium","Fibroblasts","Tumour")
  p = ggplot(as.data.frame(em), aes(x=obj2, y=emmean, colour=obj1)) +
    geom_point() +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL)) +
    scale_x_discrete(limits = rev(cluster_order)) +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          axis.title.y = element_blank()) +
    ylab("estimated marginal means") +
    # ylim(-2.5,2.5) +
    coord_flip()
  filename = paste("Comparing_neighbourhoods_", dom, sep = "")
  ggsave(plot = p, device = "png", width=6, height=4, dpi=300, path = output_dir, filename = filename)
  
  
  em <- emmeans(fit, pairwise ~ obj1 | obj2)
  df <- as.data.frame(em$contrasts)
  filenm = paste(output_dir,"stats_emmeans_", dom, ".csv", sep = "")
  write.csv(df, filenm)
  print(paste("Stats and plot saved for", dom, sep = " "))
}

# Or for total tissue
fit <- lmer(log2fc ~ obj1*obj2 + treatment + (1|mouse)+(1|mouse:roi),
            data=dat
)
em <- emmeans(fit,  ~ obj1|obj2)
cluster_order = c("CD4 T cells","Regulatory T cells","CD8 T cells","B cells", 
                  "Dendritic cells cDC1", "Dendritic cells other", "Macrophages type 1", "Macrophages type 2","Neutrophils",
                  "Endothelium","Epithelium","Fibroblasts","Tumour")
p = ggplot(as.data.frame(em), aes(x=obj2, y=emmean, colour=obj1)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL)) +
  scale_x_discrete(limits = rev(cluster_order)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.y = element_blank()) +
  ylab("estimated marginal means") +
  coord_flip()
filename = "Comparing_neighbourhoods_Total.png"
ggsave(plot = p, device = "png", width=6, height=4, dpi=300, path = output_dir, filename = filename)


em <- emmeans(fit, pairwise ~ obj1 | obj2)
df <- as.data.frame(em$contrasts)
filenm = paste(output_dir, "stats_emmeans_Total.csv", sep = "")
write.csv(df, filenm)

