
# Make 3D plot of clusters arranged over the domains Normal & Tumour without Interface, split by treatment
library("plotly")

# Set filenm as path/input_celldata_complete.csv"
filenm = "~path/Figure_3_input/input_celldata_complete.csv"
# Read in dataset
celldata = read.csv(filenm)



####################################################################################
# 3D plot with CLUSTERS
####################################################################################
dt = data.frame()

for (t in unique(factor(celldata$treatment, levels = c("Vehicle", "MRTX")))){
  cd = celldata[which(celldata$treatment == t),] 
  for (cl in sort(unique(cd$cluster2))){
    print(cl)
    df = cd[which(cd$cluster2 == cl),]
    dt[paste(t, cl, sep = "_"),"Tumour"] = 100*nrow(df[which(df$domain2 == "Tumour"),])/dim(df)[1]
    dt[paste(t, cl, sep = "_"),"Normal"] = 100*nrow(df[which(df$domain2 == "Normal"),])/dim(df)[1]
    dt[paste(t, cl, sep = "_"),"Interface"] = 100*nrow(df[which(df$domain2 == "Interface"),])/dim(df)[1]
    dt[paste(t, cl, sep = "_"),"cl"] = cl
    dt[paste(t, cl, sep = "_"), "t"] = t
    dt[paste(t, cl, sep = "_"), "clustername"] = unique(df[which(df$cluster2 == cl),"clustername"])
  }
}

# To generate the traces
tt = list()
nn = list()
ii = list()

for (unicl in unique(dt$cl)){
  print(unicl)
  cl_name = as.character(dt[grep(unicl, dt$cl), "clustername"][1])
  tt[[cl_name]] = dt[grep(unicl, row.names(dt)), "Tumour"]
  nn[[cl_name]] = dt[grep(unicl, row.names(dt)), "Normal"]
  ii[[cl_name]] = dt[grep(unicl, row.names(dt)), "Interface"]
}


# Make plot
tf <- list(
  family = "sans serif",
  size = 12,
  color = toRGB("grey50"))
p = plot_ly(data = dt, x = ~Tumour, y = ~Normal, z = ~Interface, type="scatter3d", 
            mode="markers", text = ~clustername, color = ~factor(t, levels = c("Vehicle", "MRTX")), 
            colors = c("lightcoral", "mediumturquoise"))
p = p %>% add_text(textfont = tf, textposition = "top", showlegend = FALSE)
for (cl_name in unique(dt$clustername)){
  p = p %>% add_trace(x = tt[cl_name][[1]], y = nn[cl_name][[1]], z = ii[cl_name][[1]], type = "scatter3d", mode = "lines", 
                      name = "lines", showlegend = FALSE, inherit=FALSE, text = cl_name)
}
p


####################################################################################
# 3D plot with METACLUSTERS
####################################################################################
dt = data.frame()

for (t in unique(factor(celldata$treatment, levels = c("Vehicle", "MRTX")))){
  cd = celldata[which(celldata$treatment == t),] 
  for (cl in sort(unique(cd$metacluster2)[-11])){
    print(cl)
    df = cd[which(cd$metacluster2 == cl),]
    dt[paste(t, cl, sep = "_"),"Tumour"] = 100*nrow(df[which(df$domain2 == "Tumour"),])/dim(df)[1]
    dt[paste(t, cl, sep = "_"),"Normal"] = 100*nrow(df[which(df$domain2 == "Normal"),])/dim(df)[1]
    dt[paste(t, cl, sep = "_"),"Interface"] = 100*nrow(df[which(df$domain2 == "Interface"),])/dim(df)[1]
    dt[paste(t, cl, sep = "_"),"cl"] = cl
    dt[paste(t, cl, sep = "_"), "t"] = t
    # dt[paste(t, cl, sep = "_"), "clustername"] = unique(df[which(df$cluster2 == cl),"clustername"])
  }
}

# To generate the traces
tt = list()
nn = list()
ii = list()

for (unicl in unique(dt$cl)){
  print(unicl)
  cl_name = as.character(dt[grep(unicl, dt$cl), "cl"][1])
  tt[[cl_name]] = dt[grep(unicl, row.names(dt)), "Tumour"]
  nn[[cl_name]] = dt[grep(unicl, row.names(dt)), "Normal"]
  ii[[cl_name]] = dt[grep(unicl, row.names(dt)), "Interface"]
}


# Make plot
tf <- list(
  family = "sans serif",
  size = 12,
  color = toRGB("grey50"))
p = plot_ly(data = dt, x = ~Tumour, y = ~Normal, z = ~Interface, type="scatter3d", 
            mode="markers", text = ~cl, color = ~factor(t, levels = c("Vehicle", "MRTX")), 
            colors = c("lightcoral", "mediumturquoise"))
p = p %>% add_text(textfont = tf, textposition = "top", showlegend = FALSE)

for (cl_name in unique(dt$cl)){
  p = p %>% add_trace(x = tt[cl_name][[1]], y = nn[cl_name][[1]], z = ii[cl_name][[1]], type = "scatter3d", mode = "lines", 
                      name = "lines", showlegend = FALSE, inherit=FALSE, text = cl_name)
}
p 
