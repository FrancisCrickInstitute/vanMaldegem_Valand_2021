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

# Make table with counts
#  columns needed: MouseID, ROI_name, domain2, clustername2, count
counts_table = data.frame()
n = 1
for (r in unique(celldata$ROI_name3)){
  m = unique(celldata[which(celldata$ROI_name3 == r), "MouseID"])
  for (d in unique(celldata$domain2)){
    for (cl in unique(celldata$clustername2)){
      cd = celldata[which(  celldata$ROI_name3 == r &
                              celldata$domain2 == d &
                              celldata$clustername2 == cl),]
      counts_table[n, c("MouseID", "ROI_name", "domain2", "clustername2")] = c(m,r,d,cl)
      counts_table[n, "count"] = nrow(cd)
      n = n+1
      
    }
  }
}
counts_table[1:10, 1:5]
dim(counts_table)
sum(counts_table$count)
# clean up table (remove "n/a" assignments for domain and cell type "Unclassified"):
counts_table = counts_table[which(counts_table$domain2 != "n/a"),]
counts_table = counts_table[which(counts_table$clustername2 != "Unclassified"),]
dim(counts_table)
sum(counts_table$count)

filenm = paste(output_dir, "counts_table",  ".csv", sep = "")
write.csv(counts_table, filenm)
# dat <- read_csv(filenm)[,-1]
dat <- counts_table

names(dat) <- c("mouse","roi", "domain", "celltype", "count")
dat$treatment <- ifelse(grepl("_MRTX", dat$roi), "MRTX","Vehicle")


################################################################
## fig 3c
################################################################


# Simplest approach -  pool counts for rois within mice and mice within treatment.
# Does the way counts split differently between domains for a celltype depend on the treatment
# Exploratory just to get a feel for things.  Don't use.
fit0 <- glm(count ~ (treatment+celltype+domain)^2 , data=dat, family="poisson")
fit <- glm(count ~ (treatment+celltype+domain)^2 + treatment:celltype:domain, data=dat, family="poisson")
anova(fit, fit0, test="Chisq")


# Fit mouse and domain-within-mouse as having a random celltype effect.
# Not used - sanity check for convergence
if (FALSE) {
  fit0 <- glm(count ~ (treatment+celltype+domain)^2 +mouse:celltype+mouse:domain, data=dat, family="poisson")
  fit <- glm(count ~ (treatment+celltype+domain)^2 + treatment:celltype:domain + mouse + mouse:celltype + mouse:domain, data=dat, family="poisson")
  anova(fit, fit0, test="Chisq")
  
  fit0 <- glm(count ~ (treatment+celltype+domain)^2 +mouse + mouse:roi + mouse:roi:domain + mouse:roi:celltype, data=dat, family="poisson")
  fit <- glm(count ~ (treatment+celltype+domain)^2 + treatment:celltype:domain + mouse + mouse:roi + mouse:roi:domain + mouse:roi:celltype, data=dat, family="poisson")
  anova(fit, fit0, test="Chisq")
}

## Turn counts into proportions and pretend normal
## Again, not used, just sometimes quicker and more rigorous p-values that in the poisson/multionomial case
if (FALSE) {
  mouse_dat <- dat %>%
    group_by(mouse) %>%
    mutate(freq=count/sum(count))
  fit0 <- lmer(count ~ (treatment+celltype+domain)^2 +(celltype+domain|mouse)+(celltype+domain|mouse:roi),
               dat=mouse_dat)
  fit <- lmer(count ~ (treatment+celltype+domain)^2 + treatment:celltype:domain + (celltype+domain|mouse)+(celltype+domain|mouse:roi),
              dat=mouse_dat)
}

# More complicated - There's a treatment effect, around which
# mouse-within-treatment varies, and around that mouse's expected value
# the ROIs will vary
if (file.exists(fname <- "Stats_input/glmer_null_3c.rds")) {
  fit0 <- readRDS(file=fname)
} else {
  fit0 <- glmer(count ~ (treatment+celltype+domain)^2 +(celltype+domain|mouse)+(celltype+domain|mouse:roi),
                data=dat,
                family="poisson",
                control = glmerControl(optimizer = "nloptwrap")
  )
  saveRDS(fit0, file=fname)
}


################################################################
#### The main approach for the whole dataset
#### METHODS:
#### We use the `lme4` [@lme4] package within R [@R] to fit a mixed-effects
#### model to account for a fixed effects of domain, celltype and treatment,
#### whilst allowing for a per-mouse and ROI-within-mouse variation of distribution
#### of cells between celltypes.  We use a poisson model as a surrogate to fit the
#### multinomial distribution of cell counts across the celltypes. Individual
#### comparisons are carried out using a Wald test.

#### [lme4]  Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
#### [R] R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.



################################################################
## Supplementary Figure 4f
################################################################

# Takes a while to run, loading of .rds file is recommended.
if (file.exists(fname <- "Stats_input/glmer_full_3c.rds")) {
  fit <- readRDS(file=fname)
} else{
  fit <- glmer(count ~ (treatment+celltype+domain)^2 +treatment:celltype:domain+(celltype+domain|mouse)+(celltype+domain|mouse:roi),
              data=dat,
              family="poisson",
              control = glmerControl(optimizer = "nloptwrap")
              )
  saveRDS(fit, file=fname)
}

em <- emmeans(fit, pairwise~treatment|domain+celltype)
df <- as.data.frame(em$contrasts)

pdf(file="output/treatment_per_type_by_domain.pdf",  width=9, height=6)
ggplot(df, aes(x=celltype, y=estimate, fill=p.value<0.05)) +
  geom_col() +
  coord_flip() + 
  facet_wrap(~domain) + labs(fill="P<0.05",y="Log(MRTX/Vehicle)" , x="") +
  theme_bw(base_size=18)
dev.off()

################################################################
## Supplementary Figure 6a
################################################################

# p-values for Supplementary Figure 6a, top row plots, stats for T cells in whole tissue
em <- emmeans(fit, pairwise~treatment|celltype)
df <- as.data.frame(em$contrasts)
df[which(df$celltype == "CD4 T cells"),]
df[which(df$celltype == "Regulatory T cells"),]
df[which(df$celltype == "CD8 T cells"),]


# p-values for Supplementary Figure 6a, bottom row plots, stats for T cells in the tumour domain
em <- emmeans(fit, pairwise~treatment|domain+celltype)
df <- as.data.frame(em$contrasts)
df[which(df$domain == "Tumour" & df$celltype == "CD4 T cells"),]
df[which(df$domain == "Tumour" & df$celltype == "Regulatory T cells"),]
df[which(df$domain == "Tumour" & df$celltype == "CD8 T cells"),]



################################################################
## Stats related to Figure 3d
################################################################
f3d <- subset(dat, domain=="Tumour")

## The 'correct' way is to say that mice vary around their expected treatment group:
if (file.exists(fname <- "Stats_input/glmer_null_3d.rds")) {
  fit0 <- readRDS(file=fname)
} else {
  fit0 <- glmer(count ~ treatment+celltype + (celltype|mouse) + (celltype|mouse:roi),
               data=f3d,
               family="poisson",
               control = glmerControl(optimizer = "nloptwrap")
               #             control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
               )
  saveRDS(fit0, file=fname)
}

if (file.exists(fname <- "Stats_input/glmer_full_3d.rds")) {
  fit <- readRDS(file=fname)
} else {
fit <- glmer(count ~ treatment+celltype +treatment:celltype+(celltype|mouse)+(celltype|mouse:roi),
            data=f3d,
            family="poisson",
            control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
            )
  saveRDS(fit, file=fname)
}


em <- emmeans(fit, pairwise~treatment|celltype)
df <- as.data.frame(em$contrasts)

ggplot(df, aes(x=celltype, y=estimate, fill=p.value<0.05)) +
  geom_col() +
  coord_flip()


## 
