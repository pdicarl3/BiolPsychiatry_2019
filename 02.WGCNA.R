###################################################################################################################
##### WEIGHTED GENE CO-EXPRESSION NETWORK ANALYSIS (WGCNA) ########################################################
###################################################################################################################
library(WGCNA)
options(stringsAsFactors = FALSE)
nameData = ""   # the name of your analysis

### LOAD IMPUT EXPRESSION MATRIX - pre-processed expression data (rows as samples, columns as genes)
k = 5                             # number of RUV surrogate variables to use 
input_matrix = RUVexpr.list[k]    # from the output of the previous script "01.preprocessing.R" or any other expression data matrix
pd = demo.add                     # from the output of the previous script "01.preprocessing.R" 
CPU = 2                  

### SOFT-THRESHOLDING - obtain the parameter (BETA) that satisfy scale invariance of connectivity distribution in the network
enableWGCNAThreads(nThreads = CPU)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(input_matrix, powerVector = powers, 
                        RsquaredCut = 0.8, networkType = "unsigned", verbose = 5)
pdf(paste("sft_", nameData, ".pdf", sep = ""))
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red", lwd = 3)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
rm(cex1, powers)

### GET BETA
beta_value = sft$powerEstimate
beta_value

### OBTAIN CO-EXPRESSION MODULES
enableWGCNAThreads(nThreads = CPU)
net = blockwiseModules(input_matrix, power = beta_value, minModuleSize = 40, corType = "pearson",
                       reassignThreshold = 0, mergeCutHeight = 0.05, networkType = "unsigned",
                       numericLabels = TRUE, pamRespectsDendro = FALSE, TOMType = "signed",
                       saveTOMs = FALSE, saveTOMFileBase = "Signling network", verbose = 3)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
mergedColors = labels2colors(net$colors)

### NETWORK DENDROGRAM
geneTree = net$dendrograms[[1]]
pdf(paste("dendrogram.", nameData, ".pdf", sep = ""))
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", 
                    main = paste("Genes=", length(moduleColors),  
                                 "Beta=", beta_value, "\nModules=", 
                                 length(levels(as.factor(moduleColors)))-1, 
                                 "Grey=", table(net$colors)[1]),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, cex.main = 2)
dev.off()
rm(geneTree)

### CALCULATE MODULE EIGENGENES
MEs = moduleEigengenes(input_matrix, moduleColors, softPower = beta_value)
eigen0 = data.frame(MEs$eigengenes, row.names = row.names(input_matrix))
MEs$isPC
varexp = data.frame(MEs$varExplained)
row.names(varexp) = "varexp"
names(varexp) = names(eigen0)
size = data.frame(table(moduleColors))
eigenval = list(eigenval=eigen0, stats=rbind(varexp, moduleSize=size[,2]))
rm(eigen0, size, varexp)

### plot variance explained and module size
s = data.frame(t(eigenval$stats))
pdf(paste("modulestats.", nameData, ".pdf", sep=""))
hist(s[!row.names(s) %in% "MEgrey",1]*100, main="Variance Explained by Module Eigengenes", xlab="Variance Explained %")
hist(s[!row.names(s) %in% "MEgrey",2], main="Module Size", xlab="Number of Genes")
require(car)
scatterplot(varexp*100~moduleSize, data=s[!row.names(s) %in% "MEgrey",], 
            xlab="Module Size", ylab="Variance Explained %")
dev.off()

### CALCULATE ADJACENCY MATRIX
enableWGCNAThreads(nThreads = CPU)
adjacency <- adjacency(input_matrix, selectCols = NULL, 
                       type = "unsigned", power = beta_value,
                       corFnc = "cor", corOptions = "use = 'p'",
                       distFnc = "dist", distOptions = "method = 'euclidean'")

### ESITMATE GENE CONNECTIVITY 
IMk <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = FALSE)
IMk.module <- data.frame(IMk, moduleColors)
IMksc <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = TRUE)
IMk.module.sc <- data.frame(IMksc, moduleColors)
IMk.module <- data.frame(IMk.module, K.sc = IMk.module.sc$kWithin)
rm(adjacency, IMk, IMksc, IMk.module.sc)

### save
assign(paste0("eigenval.", nameData), eigenval)
assign(paste0( "IMk.module.", nameData), IMk.module) 
assign(paste0("net.", nameData), net)
assign(paste0("sft.", nameData), sft)
save(list = c(paste0("eigenval.", nameData),
              paste0( "IMk.module.", nameData),
              paste0("net.", nameData),
              paste0("sft.", nameData)),
     file = paste0("WGCNA.", nameData, ".RData"))
