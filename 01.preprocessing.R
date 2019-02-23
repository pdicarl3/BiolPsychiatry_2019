######################################################################################################
#### EXPRESSION DATA PRE-PROCESSING - remove unwanted variation with "RUVcorr" R package #############
######################################################################################################
source("https://bioconductor.org/biocLite.R")
biocLite("RUVcorr")
library(RUVcorr)

### LOAD EXPRESSION DATA MATRIX [samples, genes] (e.g., log2(RPKM+1) expression values) 
### LOAD gene annotations
### LOAD demographics
dataset = "a"           # the name of your dataset
nameData =  "b"         # the name of your analysis
expr = expr             # expr = matrix with samples as rows and genes as columns
map = map               # gene annotations
demo = demo             # demographics
CPU = 6                 # CPUs for parallel computation

### LOAD SEPARATE VECTORS OF POSITIVE CONTROL GENE NAMES (known co-expressed pathways)
index_MHC = which(map$Gene%in% positive_MHC)
index_PRC2 = which(map$Gene %in% positive_PRC2)
index_sodium = which(map$Gene %in% positive_sodium)
posCon.index = list(MHC = index_MHC, PRC2 = index_PRC2, sodium = index_sodium)

### LOAD A VECTOR OF HOUSEKEEPING GENE NAMES
index_PGC = which(map$Gene %in% as.character(PGC_genes))

### EXCLUDE SCHIZOPRHENIA RISK GENES FROM NEGATIVE CONTROL GENES (RANDOM GENES)
exclude_index = index_PGC
index_nc = empNegativeControls(expr, exclude=exclude_index, nc=3000)

### MEAN/INTERQUARTILE PLOT (HIGHLIGHTS HOUSEKEEPING GENES) 
### FIGURE S01 A-B
pdf(paste0(dataset, ".", nameData, ".IQR.Mean.plot.pdf"))
genePlot(expr, index=index_hk_13, legend="HK genes", title="IQR-Mean Plot")
dev.off()

### SELECT RUV PARAMETER TO COMPARE
library(snowfall)
k <- c(1,2,3,4,5,6,7,8)
nu <- c(0,500,1000,5000)
k.nu.matrix <- cbind(rep(k, each=4), rep(nu, 4))
k.nu.matrix <- as.list(as.data.frame(t(k.nu.matrix)))

### PREPROCESS DATA USING DIFFERENT PARAMETER COMBINATIONS
sfInit(parallel=TRUE, cpus= CPU)
sfLibrary(RUVcorr)
sfExport("expr", "k.nu.matrix", "index_hk_13")
expr_AllRUV <- sfLapply(k.nu.matrix, function(x) RUVNaiveRidge(expr, center=TRUE, nc_index = index_hk_13, nu=x[2], kW=x[1]))
sfStop()

### DIAGNOSTIC PLOTS FOR POSITIVE CONTROL GENE SETS
### FIGURE S02 B-D
pdf("RUV.parameter.posControls.pdf")
for(i in names(posCon.index)) {
  cor_AllRUV <- lapply(expr_AllRUV, function(x) cor(x[,posCon.index[[i]]]))
  cor_Raw <- cor(expr[,posCon.index[[i]]])
  
  par(mfrow=c(2,2))
  lapply(1:4, function(i) histogramPlot(cor_AllRUV[seq(0,31,4)+i], cor_Raw,
                                        col.X = c("grey", "cyan", "greenyellow", "blue", "orange", "magenta", "darkgreen", "purple"),
                                        title=paste("nu=", nu[i]), 
                                        legend=c(paste("k=", k), "Raw")))
}
dev.off()

### DIAGNOSTIC PLOTS FOR NEGATIVE CONTROL GENE SETS
### FIGURE S02 A
index_random <- background(expr, nBG=500, exclude=exclude_index, nc_index=index_hk_13)
cor_AllRUV_random <- lapply(expr_AllRUV, function(x) cor(x[,index_random]))
cor_Raw_random <- cor(expr[,index_random])
pdf("RUV.parameter.RandomControls.pdf")
par(mfrow=c(2,2))
lapply(1:4, function(i) histogramPlot(cor_AllRUV_random[seq(0,31,4)+i], cor_Raw_random,
                                      col.X = c("grey", "cyan", "greenyellow", "blue", "orange", "magenta", "darkgreen", "purple"),
                                      title=paste("nu=", nu[i]), 
                                      legend=c(paste("k=", k), "Raw")))
dev.off()

### OTHER DIAGNOSTIC PLOTS (not shown, see RUVcorr package manual)
pdf("RUV.RLE.parameter.pdf")
for(i in seq(0,31,4)) {
  par(mfrow=c(2,2))
  lapply(1:4, function(x) RLEPlot(expr, expr_AllRUV[[i+x]], 
                                  name=c("Raw", "RUV"), title=paste("nu=", nu[x]),
                                  method="IQR.boxplots"))
}
dev.off()

### OTHER DIAGNOSTIC PLOTS (not shown, see RUVcorr package manual)
pdf("RUV.RLE.2.parameter.pdf")
for(i in 1:length(expr_AllRUV)) {
  par(mfrow=c(1,1))
  RLEPlot(expr, expr_AllRUV[[i]], name=c("Raw", "RUV"), 
          title=paste("Race", i), method="IQR.points", anno=demo, 
          Factor="Race", numeric=FALSE)
  RLEPlot(expr, expr_AllRUV[[i]], name=c("Raw", "RUV"), 
          title=paste("Sex", i), method="IQR.points", anno=demo, 
          Factor="Sex", numeric=FALSE)
  RLEPlot(expr, expr_AllRUV[[i]], name=c("Raw", "RUV"), 
          title=paste("RIN", i), method="IQR.points", anno=demo, 
          Factor="RIN", numeric=TRUE)
  RLEPlot(expr, expr_AllRUV[[i]], name=c("Raw", "RUV"), 
          title=paste("Age", i), method="IQR.points", anno=demo, 
          Factor="Age", numeric=TRUE)
  RLEPlot(expr, expr_AllRUV[[i]], name=c("Raw", "RUV"),
          title=paste("Diagnosis", i), method="IQR.points", anno=demo,
          Factor="Dx", numeric=FALSE)
  RLEPlot(expr, expr_AllRUV[[i]], name=c("Raw", "RUV"),
          title=paste("totalMapped", i), method="IQR.points", anno=demo,
          Factor="totalMapped", numeric=TRUE)
  RLEPlot(expr, expr_AllRUV[[i]], name=c("Raw", "RUV"),
          title=paste("mitoMapped", i), method="IQR.points", anno=demo,
          Factor="mitoMapped", numeric=TRUE)
  
}
dev.off()

#########################################################################################################
### PRE-PROCESS EXPRESSION DATA WITH RUV SURROGATE VARIABLES ############################################
#########################################################################################################    
### LOAD MODIFIED RUV FUNCTION TO EXTRACT SURROGATE VARIABLES
source("RUVNaiveRidge_mod.R")
### select the number of surrogate variables to use
k.var = 1:12   
### run
RUVexpr.list = list()
for(k in k.var) {
  RUV = RUVNaiveRidge.default(expr, center = TRUE, nc_index = index_hk_13, nu = nu, kW = k, check.input = TRUE)
  RUVexpr.list[[k]] = RUV$Yhat
}
names(RUVexpr.list) = paste0("RUV.k", 1:max(k.var))
RUVexpr = RUVexpr.list[[max(k.var)]]
if(identical(row.names(demo), row.names(RUVexpr))) {
  demo.add = cbind(demo, RUV$What) } else {"row.names do not match"}
colnames(demo.add)[(ncol(demo.add)-max(k.var)+1) :ncol(demo.add)] = paste0("RUV", 1:max(k.var), "_hk13")

### save
assign(paste("RUV", dataset, nameData, sep="."), RUVexpr.list)
assign(paste("demo", dataset, nameData, sep="."), demo.add)
assign(paste("geneMap", dataset, nameData, sep="."), map)
save(list = c(paste("RUV", dataset, nameData, sep="."), 
              paste("demo", dataset, nameData, sep="."),
              paste("geneMap", dataset, nameData, sep="."),
              "index_hk_13"),
     file = paste0(dataset, ".", nameData, ".", "RUV", ".k", max(k.var), ".nu", nu, ".", ".RData"))

