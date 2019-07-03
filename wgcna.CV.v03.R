##########################################################################################
######## NETWORK ANALYSIS CROSS VALIDATION ###############################################
##########################################################################################
cv.round = 1
num.folds = 3
L = 1
### load RNA sequencing data
load("/users/pdicarlo/rnaseq-dlpfc/rse_gene_BrainSeq_Phase1_hg19_TopHat2_EnsemblV75.rda")
require(SummarizedExperiment)
## envinronment
nameData = 'LIBD'
CPU = 2
## network parameter
networkType = 'signed hybrid'
corType = 'bicor'       ## 'pearson' or 'bicor'
corType.sft = 'bicor'   ## 'cor' or 'bicor'
modSize = 40
beta_fix = 6
## covariates
retCov = c('Age', 'Sex')  ## variable to retain in removing variation
obsCov = c('snpPC1','snpPC2','snpPC3','snpPC4','snpPC5','snpPC6','snpPC7','snpPC8','snpPC9','snpPC10',
           'RIN','mitoRate','totalAssignedGene')             
## pre-processing parameters
preMethod = 'RUV'  ## default = 'RUV'  alternative 'qSV'
modMatrix = 'mod0'    ## default = 'mod0'  alternative 'mod' (model matrix with observed covariates)
obsCovMethod = FALSE  ## TRUE if use observed covariates are added
nameData = paste0(nameData, '.', num.folds, 'kCV.', modMatrix,'.', preMethod, '.obsCov', obsCovMethod)
nameFile = paste0(nameData, '.r', cv.round, '.fld', L)
### demographics
rpkm = .1
age.min = 17
age.max = 65
case = 'Control'
ethnicity = c('AA', 'CAUC')
rin = 6
##########################################################################################
### pre-processing #######################################################################
##########################################################################################
### subset subjects
se_gene = rse_gene[,rse_gene$Age >= age.min &
                     rse_gene$Age <= age.max & 
                     rse_gene$Dx %in% case &
                     rse_gene$Race %in% ethnicity &
                     rse_gene$RIN >= rin] 
require(recount)
p = getRPKM(se_gene)
pd = data.frame(colData(se_gene))
map = data.frame(rowRanges(se_gene))
row.names(map) = row.names(p)

### subset genes RPKM >= 0.1
map$rpkm_median = apply(p, 1, median)
map = map[map$rpkm_median >= rpkm,]
p = p[row.names(map),]
p = log2(p +1)

### subset degradatin matrix 
se_qsvs = rse_qsvs[,rse_qsvs$Age >= age.min &
                     rse_qsvs$Age <= age.max & 
                     rse_qsvs$Dx %in% case &
                     rse_qsvs$Race %in% ethnicity &
                     rse_qsvs$RIN >= rin] 
degrMat = assays(se_qsvs)$count
degrMat = log2(degrMat +1)

##########################################################################################
###### pre-processing ####################################################################
### study the number of qSV to extract   #################################################
require(sva)
dim(degrMat)  # SAMPLES must be in COLUMNS for sva
mod = model.matrix(~ snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + 
                     snpPC6 + snpPC7 + snpPC8 + snpPC9 + snpPC10 + 
                     RIN + mitoRate + totalAssignedGene,
                   data = pd)
mod0 = model.matrix(~1, data = pd)

### check the number of component to extract
n.qsv = num.sv(degrMat, mod = if(modMatrix == 'mod0') {mod0} else {mod}, method="be", B = 20)
n.qsv

### perform PCA on degradation matrix
degrMat = t(degrMat)  # transpose for PCA
qsvs = prcomp(degrMat, scale = T, center = T)$x[,1:n.qsv, drop=F]
if(n.qsv > 0) { 
  colnames(qsvs) = paste0('qSV', 1:n.qsv) } else { "zero qSVs selected" }
### add qSV to demographics
pd = cbind(pd, qsvs)

##########################################################################################
###### pre-processing ####################################################################
### study the number of RUV to extract   #################################################
### load housekeeping genes
load("/users/pdicarlo/rnaseq-dlpfc/hk_Eisenberg2013.RData")
require(RUVcorr)
require(sva)
dim(p)
#### obtain list of housekeeping genes (exclude genes of interest)
hk = setdiff(housekeeping_2013, Ripke2014)
index.hk13 = which(map$Symbol %in% hk)
index.Ripke2014 = which(map$Symbol %in% Ripke2014)
### obtain housekeeping genes matrix
hk.expr = p[index.hk13,]
dim(hk.expr)

############# model matrix  ########################################################
mod = model.matrix(~ snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + 
                     snpPC6 + snpPC7 + snpPC8 + snpPC9 + snpPC10 + 
                     RIN + mitoRate + totalAssignedGene,
                   data = pd)
mod0 = model.matrix(~1, data = pd)
############# COMPUTE SURROGATE VARIABLE ON THE RIGHT MATRIX  ######################
##### check the number of component to extract  # SAMPLES must be in COLUMNS for sva
n.ruv = num.sv(hk.expr, mod = if(modMatrix == 'mod0') {mod0} else {mod}, method="be", B = 20)
n.ruv

### OUTLIERS DETECTION ######################################################
### compute PCA #############################################################
pr = prcomp(t(p), center=TRUE, scale=TRUE)

#####################################################
### Cook's distance #################################
fit = lm(pr$x[,1]~ pr$x[,2])
cooksd = cooks.distance(fit)
influential = names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))]  # influential row numbers
cooksd.1 = names(cooksd)[(cooksd > 1)]  # influential row numbers
require(Rmisc)
ci95 = CI(cooksd, ci = .95)
influential.ci95 = names(cooksd)[(cooksd > ci95[1])]

#####################################################
### inter array correlation #########################
numSD = 2
numSD.3 = 3
# Calculating IACs for all pairs of samples and examining the distribution of IACs in the dataset:
IAC = cor(p,use="p") ## genes should be rows
## Performing hierachical clustering (average linkage) using 1-IAC as a distance metric:
require(cluster)
cluster1 = hclust(as.dist(1-IAC), method="average")
## Another way to visualize outliers is to calculate the mean IAC for each array and examine this distribution:
meanIAC = apply(IAC, 2, mean)
sdCorr = sd(meanIAC)
numbersd = (meanIAC-mean(meanIAC))/sdCorr
sdout = -numSD
sdout.3 = -numSD.3
outliers.2sd = names(numbersd)[numbersd < sdout]
outliers.3sd = names(numbersd)[numbersd < sdout.3]
#######################################################
### REMOVE OULIERS ####################################
removeOut = union(cooksd.1, outliers.3sd)
good.Samples = setdiff(row.names(pd), removeOut)
p = p[,good.Samples]
pd = pd[good.Samples,]
identical(row.names(t(p)), row.names(pd))

#############################################################################################
### WGCNA ###################################################################################
require(WGCNA)
options(stringsAsFactors = FALSE)
####### USE RUV TO GET RESIDUALS #############################################################
source("/users/pdicarlo/R/custom_functions/RUVNaiveRidge_mod.r")
RUV = RUVNaiveRidge.default(t(p), center = TRUE, nc_index = index.hk13, 
                            nu = 0, kW = n.ruv, check.input = TRUE)
RUVexpr = RUV$Yhat
RUVcomp = RUV$What

### add RUV components to demographics
if(n.ruv > 0) {
  colnames(RUVcomp) = paste0('RUV', 1:n.ruv) } else { "zero RUVs selected" }
if(identical(row.names(pd), row.names(RUVexpr))) {
  pd = cbind(pd, RUVcomp) } else { 'ERROR: row.names do not match'}

#############################################################################################
### RUV or qSVs
if(preMethod == 'RUV') {
  if(n.ruv > 0) {
    uwv = pd[,c(paste0('RUV', 1:n.ruv))]
  } else { uwv = NULL } 
} else { if(n.qsv > 0) { 
  uwv = pd[,c(paste0('qSV', 1:n.qsv))] } else { uwv = NULL } } 

### observed covariates or NOT
if(obsCovMethod == TRUE) { 
  uwv = cbind(uwv, pd[,obsCov]) 
} else { uwv = uwv }

### covariates to be retained (modeled in empirical bayes)
retCovMat = pd[,retCov]

#############################################################################################
### empirical bayes function ################################################################
eBLM = empiricalBayesLM(data = t(p), removedCovariates = uwv,  retainedCovariates = retCovMat,
                        weights = NULL, weightType = "empirical", stopOnSmallWeights = TRUE,
                        tol = 1e-4, maxIterations = 1000, scaleMeanToSamples = NULL,
                        robustPriors = TRUE, automaticWeights = "bicov", aw.maxPOutliers = 0.1,
                        verbose = 5, indent = 3)


##############################################################################################
### CREATE NETWORK  ##########################################################################
##############################################################################################
input_tot = eBLM$adjustedData
dim(input_tot)
### save #####################################################################################
save(input_tot, pd, map, cooksd.1, outliers.3sd, n.qsv, n.ruv,
     file = paste0("input.", nameData, ".RData"))

##############################################################################################
####### output lists #########################################################################
res.LOADcorr = list()
res.varExpTrain = list()
res.varExpTest = list()
res.varExpDiff = list()
res.coClust = list()
res.GOBPfdr = list()
res.ratioModules = list()
res.cohend = list()
##############################################################################################
### START CROSS VALIDATION  ##################################################################
##############################################################################################
require(caret)
### set seed 
set.seed(cv.round)
### create folds
flds <- createFolds(1:nrow(input_tot), k = num.folds, list = TRUE, returnTrain = FALSE)
print(flds$Fold1[1:20])

## train
input_matrix = input_tot[-flds[[L]],]
pd_cv = pd[-flds[[L]],]
## test
test_matrix = input_tot[flds[[L]],]
test_pd = pd[flds[[L]],]
###########################################################################################
# #### soft-thresholding ##################################################################
# enableWGCNAThreads(nThreads = CPU)
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# sft = pickSoftThreshold(input_matrix, dataIsExpr = TRUE,
#                         powerVector = powers, blockSize = NULL,
#                         RsquaredCut = 0.8, networkType = networkType, corFnc = corType.sft,
#                         verbose = 5)
# #### get beta power
# beta_value = sft$powerEstimate
beta_value = beta_fix

###########################################################################################
### NETWORK CONSTRUCTION ##################################################################
### TRAIN set  
enableWGCNAThreads(nThreads = CPU)
net = blockwiseModules(input_matrix, power = beta_value, minModuleSize = modSize, 
                       corType = corType, networkType = networkType,
                       reassignThreshold = 0, mergeCutHeight = 0.05, 
                       numericLabels = TRUE, pamRespectsDendro = FALSE, 
                       TOMType = "signed", saveTOMs = FALSE, maxBlockSize = 4000,
                       saveTOMFileBase = "Signling network",
                       verbose = 3)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
### TEST set
enableWGCNAThreads(nThreads = CPU)
netTest = blockwiseModules(test_matrix, power = beta_value, minModuleSize = modSize, 
                           corType = corType, networkType = networkType,
                           reassignThreshold = 0, mergeCutHeight = 0.05, 
                           numericLabels = TRUE, pamRespectsDendro = FALSE, 
                           TOMType = "signed", saveTOMs = FALSE, maxBlockSize = 4000,
                           saveTOMFileBase = "Signling network",
                           verbose = 3)
table(netTest$colors)
moduleLabelsTest = netTest$colors
moduleColorsTest = labels2colors(netTest$colors)

##########################################################################################
## Recalculate module eigengenes TRAIN
MEs = moduleEigengenes(input_matrix, moduleColors, softPower = beta_value)
## Recalculate module eigengenes TEST (based on TRAIN MODULES)
MEsTest = moduleEigengenes(test_matrix, moduleColors, softPower = beta_value)
moduleData = data.frame(table(moduleColors), 
                        varExpTrain = t(MEs$varExplained), 
                        varExpTest = t(MEsTest$varExplained))
row.names(moduleData) = moduleData$moduleColors
moduleData$diff = moduleData$varExpTrain - moduleData$varExpTest


##### get intramodular connectivity #################################################################
enableWGCNAThreads(nThreads = CPU)
IMk <- intramodularConnectivity.fromExpr(datExpr = input_matrix, colors =  moduleColors,
                                         corFnc = corType, corOptions = "use = 'p'",
                                         distFnc = "dist", distOptions = "method = 'euclidean'",
                                         networkType = networkType, power = if (networkType=="distance") 1 else 6,
                                         scaleByMax = FALSE, ignoreColors = if (is.numeric(colors)) 0 else "grey",
                                         getWholeNetworkConnectivity = TRUE)
IMksc <- intramodularConnectivity.fromExpr(datExpr = input_matrix, colors =  moduleColors,
                                           corFnc = corType, corOptions = "use = 'p'",
                                           distFnc = "dist", distOptions = "method = 'euclidean'",
                                           networkType = networkType, power = if (networkType=="distance") 1 else 6,
                                           scaleByMax = TRUE, ignoreColors = if (is.numeric(colors)) 0 else "grey",
                                           getWholeNetworkConnectivity = FALSE)
imkx <- data.frame(moduleColors, IMk, IMksc)
row.names(imkx) = colnames(input_tot)
################################################################################################
##### compute module eigengene through PCA #####################################################
### get gene sets indeces
modLabels = sort(unique(imkx$moduleColors))
modLabels = modLabels[!modLabels %in% "grey"]
###
list_index = lapply(modLabels, function(mod)  which(imkx$moduleColors %in% mod))
names(list_index) = modLabels
### train
list_modules = lapply(list_index, function(index) input_matrix[,index])
list_pca = lapply(list_modules, function(matr) prcomp(matr))
loadings1 = lapply(list_modules, function(matr) prcomp(matr)$rotation[,1])
eigen1 = sapply(list_modules, function(matr) prcomp(matr)$x[,1])
### test
test_modules = lapply(list_index, function(index) test_matrix[,index])
test_pca = lapply(test_modules, function(matr) prcomp(matr))
test_loadings1 = lapply(test_modules, function(matr) prcomp(matr)$rotation[,1])
test_eigen1 = sapply(test_modules, function(matr) prcomp(matr)$x[,1])

################################################################################################
##### OUTCOME MEASURES #########################################################################
##### 1. correlation between factor loadings 
res.LOADcorr[[L]] = as.vector(abs(sapply(modLabels, function(mod) cor(loadings1[[mod]], test_loadings1[[mod]], method = 'spearman'))))
##### 2. variance explained
moduleData = moduleData[modLabels,]
res.varExpTrain[[L]] = moduleData$varExpTrain*100
res.varExpTest[[L]] = moduleData$varExpTrain*100
##### 3. difference of variance explained
res.varExpDiff[[L]] = moduleData$diff*100
##### 4. co-clustering
res.coClust[[L]] = as.vector(coClustering(moduleColors, moduleColorsTest, tupletSize = 2, unassignedLabel = "grey"))

##################################################################################
###### Gene Ontology with gprofiler ##############################################
##################################################################################
require(gProfileR)
## set parameter
query.list = sapply(modLabels, function(x) row.names(imkx)[imkx$moduleColors == x])
onto = c('GO:BP', 'GO:MF', 'GO:CC')
background = row.names(imkx)
minSize = 20
maxSize = 2000
minOverlap = 0
maxPvalue = 1
correction = 'fdr'
## get results for all ontologies available
GOresults = list()
for(o in onto) {
  GOresults[[o]] =  gprofiler(query = query.list, organism = "hsapiens", 
                              sort_by_structure = T, ordered_query = F, significant = T, 
                              exclude_iea = F, underrep = F, evcodes = F, region_query = F, 
                              max_p_value = maxPvalue, 
                              min_set_size = minSize, max_set_size = maxSize, min_isect_size = minOverlap, 
                              correction_method = correction,
                              hier_filtering = "moderate", domain_size = "annotated", 
                              custom_bg = background,
                              numeric_ns = "", 
                              png_fn = NULL, include_graph = F, 
                              src_filter = o)
}
gobp = GOresults$`GO:BP`
goodgenes = sum(imkx$moduleColors != 'grey')

##### OUTCOME
##### 5. number of GO:BP per module   
res.GOBPfdr[[L]] = nrow(gobp)/length(modLabels)

########################################################################################
################# Johnson TOM replication script #######################################
########################################################################################
IMk_module = imkx
IMk_module$Gene = row.names(IMk_module)
### set permutations
nOfPerm = 1000

### calculate TOM
enableWGCNAThreads()
TOM <- TOMsimilarityFromExpr(test_matrix, corType = corType, networkType = networkType, 
                             power = beta_value, TOMType = 'signed', TOMDenom = "min",
                             maxPOutliers = 1, quickCor = 0, pearsonFallback = "individual",
                             cosineCorrelation = FALSE, nThreads = CPU,
                             verbose = 3, indent = 0)
dimnames(TOM) <- list(colnames(input_matrix), colnames(input_matrix))

########################################
### create empty list to collect results
annotation_list <- list()
module_tom_list <- list()
median_random_list <- list()
cohen_d_list <- list()
p_list <- list()
p_empirical <- list()
###
for(i in modLabels) {
  annotation <- IMk_module$Gene[IMk_module$moduleColors %in% i]
  annotation <- unique(annotation)
  annotation_list[[i]] <- annotation
  tomIndex = which(colnames(TOM) %in% annotation)
  module_tom <- TOM[tomIndex, tomIndex]
  diag(module_tom) <- NA
  module_tom_list[[i]] <- module_tom
  n_genes <- ncol(module_tom)
  ####
  cohen.d <- list()
  median_random_tom <- list()  
  for(r in 1:nOfPerm) {
    s <- sample(colnames(TOM), n_genes)
    s_Index <- which(colnames(TOM) %in% s)
    tom_m <- TOM[s,s]
    diag(tom_m) <- NA
    ## compute random median
    median_random_tom[[r]] <- median(tom_m, na.rm = TRUE)
    ## compute cohed d vs. random TOM
    require(effsize)
    cohen.d[[r]] = cohen.d(module_tom[upper.tri(module_tom, diag = FALSE)], 
                           tom_m[upper.tri(tom_m, diag = FALSE)])$estimate 
  }
  ### aggregate per module
  median_random_list[[i]] <- unlist(median_random_tom)
  cohen_d_list[[i]] <- unlist(cohen.d)
  
  ### compute empirical p-value
  p_list[[i]] <- signif(pnorm(median(module_tom_list[[i]], na.rm = TRUE), 
                              mean = mean(median_random_list[[i]]), sd=sd(median_random_list[[i]]), 
                              lower.tail = FALSE),4)
  z <- median(module_tom_list[[i]], na.rm = TRUE) <= median_random_list[[i]]
  p_empirical[[i]] <- sum(z)/nOfPerm
  print(i)
}
##
p_list <- unlist(p_list)
p_list <- data.frame(p_list)
p_list$bonferroni <- p.adjust(p_list$p_list, method="bonf")
p_empirical = unlist(p_empirical)
p_empirical = data.frame(p_empirical)
p_empirical$bonferroni = p.adjust(p_empirical$p_empirical, method="bonf")
###
p_empirical$emp_05 = ifelse(p_empirical$p_empirical < 0.05, "p<.05", "not")
p_empirical$bonf_05 = ifelse(p_empirical$bonferroni < 0.05, "p-corr<.05", "not")
p_empirical$emp_min = ifelse(p_empirical$p_empirical < 1/nOfPerm, paste0("p<", 1/nOfPerm), "not")
###
p_Johnson = if(identical(row.names(p_list), row.names(p_empirical))) {
  cbind(p_list, p_empirical) } else {"row.names do no match"}
p_Johnson$moduleSize.Ref = if(identical(row.names(p_Johnson), names(annotation_list))) {
  sapply(annotation_list, length)} else {"row.names do no match"}
p_Johnson$moduleSize.Test = if(identical(row.names(p_Johnson), names(module_tom_list))) {
  sapply(module_tom_list, function(i) ncol(i))} else {"row.names do no match"}

##### OUTCOME 
##### 6. % of replicated modules   
res.ratioModules[[L]] = sum(p_Johnson$p_empirical < .001)/length(modLabels)*100
##### 7. cohen d of module replication
res.cohend[[L]] = as.vector(sapply(cohen_d_list, function(x) median(x)) )

### unlist output
res.LOADcorr = unlist(res.LOADcorr)
res.varExpTrain = unlist(res.varExpTrain)
res.varExpTest = unlist(res.varExpTest)
res.varExpDiff = unlist(res.varExpDiff)
res.coClust = unlist(res.coClust)
res.GOBPfdr = unlist(res.GOBPfdr)
res.ratioModules = unlist(res.ratioModules)
res.cohend = unlist(res.cohend)

### save  
save(res.LOADcorr, res.varExp, res.varExpDiff, res.coClust, 
     res.GOBPfdr, res.ratioModules, res.cohend, flds,
     file = paste0(nameFile, '.RData'))


rm(list=ls(all=TRUE))