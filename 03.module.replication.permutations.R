#######################################################################################################################
################# MODULE REPLICATION (PERMUTATIONS) ###################################################################
#######################################################################################################################
library(WGCNA)
options(stringsAsFactors = FALSE)
networkType = "unsigned"    # network type
corType = "pearson"         # correlation type
nOfPerm = 10000             # number of permutations 
beta_value = 5              # network soft-threshold

##############################################################################
### load TEST set expression matrix [samples,genes]  #########################
test.Data = "test.Data"  
input_matrix = input_matrix

##############################################################################
### load WGCNA output or any module labels (see "02.WGCNA.R") ################
ref.Data = "ref.Data"       
IMk_module = IMk.module # output of "02.WGCNA.R" is a data.frame with rows as genes and a column with module labels for each gene.
IMk_module$Gene = row.names(IMk_module)

### obtain a vector of module labels and exclude grey
col = sort(as.character(unique(IMk_module$moduleColors)))
col = col[!col %in% "grey"]
col

###############################################################################
### calculate TOPOLOGICAL OVERLAP MATRIX (TOM) ################################
enableWGCNAThreads()
ptm <- proc.time()
TOM <- TOMsimilarityFromExpr(input_matrix, corType = corType, networkType = networkType, 
                             power = beta_value, TOMType = "signed", TOMDenom = "min",
                             maxPOutliers = 1, quickCor = 0, pearsonFallback = "individual",
                             cosineCorrelation = FALSE, nThreads = 2,
                             verbose = 3, indent = 0)
dimnames(TOM) <- list(colnames(input_matrix), colnames(input_matrix))
t = proc.time() - ptm

################################################################################
### create empty list to collect results #######################################
annotation_list <- list()
module_tom_list <- list()
median_random_list <- list()
p_list <- list()
p_empirical <- list()
ptm <- proc.time()
for(i in col) {   ## looping modules
  ## select genes within a module
  annotation <- IMk_module$Gene[IMk_module$moduleColors %in% i]
  annotation <- unique(annotation)
  annotation_list[[i]] <- annotation
  ## subset TOM with module genes
  tomIndex = which(colnames(TOM) %in% annotation)
  module_tom <- TOM[tomIndex, tomIndex]
  diag(module_tom) <- NA
  module_tom_list[[i]] <- module_tom
  n_genes <- ncol(module_tom)
  ## random subset of TOM with the same number of genes as reference module
  median_random_tom <- list()  
  for(r in 1:nOfPerm) {
    s <- sample(colnames(TOM), n_genes)
    s_Index <- which(colnames(TOM) %in% s)
    tom_m <- TOM[s,s]
    diag(tom_m) <- NA
    ## take the median of each random TOM subset
    median_random_tom[[r]] <- median(tom_m, na.rm = TRUE)
  }
  median_random_list[[i]] <- unlist(median_random_tom)
  ## compute p-value as Z-distribution comparing actual module TOM median with random distribution of medians
  p_list[[i]] <- signif(pnorm(median(module_tom_list[[i]], na.rm = TRUE), 
                              mean = mean(median_random_list[[i]]), sd=sd(median_random_list[[i]]), 
                              lower.tail = FALSE),4)
  ## compute empirical p-value comparing actual module TOM median with random distribution of medians
  z <- median(module_tom_list[[i]], na.rm = TRUE) <= median_random_list[[i]]
  p_empirical[[i]] <- sum(z)/nOfPerm
  print(i)
}
t1 = proc.time() - ptm
## collect results
p_list <- unlist(p_list)                                       ## list of Z-distribution p-values
p_list <- data.frame(p_list)
p_list$bonferroni <- p.adjust(p_list$p_list, method="bonf")    ## list of Z-distribution p-values (Bonferroni correction)
p_empirical = unlist(p_empirical)                              ## list of empirical p-value
p_empirical = data.frame(p_empirical)
p_empirical$bonferroni = p.adjust(p_empirical$p_empirical, method="bonf")  ## list of empirical p-value (Bonferroni correction)
p_empirical$emp_05 = ifelse(p_empirical$p_empirical < 0.05, "p<.05", "not")  ## modules with p-empirical < 0.05
p_empirical$bonf_05 = ifelse(p_empirical$bonferroni < 0.05, "p-corr<.05", "not") ## modules with p-empirical Bonferroni < 0.05
p_empirical$emp_min = ifelse(p_empirical$p_empirical < 1/nOfPerm, paste0("p<", 1/nOfPerm), "not") ## modules with p-empirical < 0.0001
## the "p_Johnson" data.frame contains the results combine "p_list" and "p_empirical" in a single data.frame
p_Johnson = if(identical(row.names(p_list), row.names(p_empirical))) {
  cbind(p_list, p_empirical) } else {"row.names do no match"}
p_Johnson$moduleSize.Ref = if(identical(row.names(p_Johnson), names(annotation_list))) {
  sapply(annotation_list, length)} else {"row.names do no match"}
p_Johnson$moduleSize.Test = if(identical(row.names(p_Johnson), names(module_tom_list))) {
  sapply(module_tom_list, function(i) ncol(i))} else {"row.names do no match"}

### save
save(module_tom_list, annotation_list, median_random_list, p_Johnson,
     file=paste("Johnson.replication", "test", test.Data, "ref", ref.Data, "RData", sep = "."))

############################################################################
### plot to visualize results  #############################################
pdf(paste("Johnson.replication.empirical", "test", test.Data, "ref", ref.Data, "pdf", sep = "."))

## histograms of null distribustion -log10(median[TOM])
pdf(paste("Johnson.replication.empirical.log", "test", test.Data, "ref", ref.Data, "pdf", sep = "."))
for(i in col) {
  hist(-log10(median_random_list[[i]]), main=paste("test: ", test.Data, " (ref: ", ref.Data, ")", "\n", i, "\np-empirical=",
                                                   signif(p_empirical[i,"p_empirical"],4), sep=""),  
       xlab = "-log10 (TOM median)", xlim = c(-log10(min(median_random_list[[i]])), -log10(max(median_random_list[[i]])*3)))
  abline(v=median(-log10(module_tom_list[[i]]), na.rm = TRUE), col="black", lwd = 3)
}
dev.off()

## barplot of replicated modules under certain threshold
pdf(paste("Johnson.bars", "test", test.Data, "ref", ref.Data, "pdf", sep="."))
par(mfrow=c(1,3), las = 1)
barplot(table(p_empirical$emp_05), col = c("lightgoldenrodyellow", "indianred2"),  ## modules with p-empirical < 0.05
        ylab = "N of gene sets", xlab = "empirical p-value", ylim = c(0,60), 
        cex.names = 1.5, cex.axis = 1.2, cex.lab = 1.5)
barplot(table(p_empirical$bonf_05), col = c("lightgoldenrodyellow", "indianred3"), ## modules with p-empirical Bonferroni < 0.05
        ylab = "N of gene sets", xlab = "empirical p-value", ylim = c(0,60), 
        cex.names = 1.5, cex.axis = 1.2, cex.lab = 1.5)
barplot(table(p_empirical$emp_min), col = c("lightgoldenrodyellow", "indianred4"), ## modules with p-empirical < 0.0001
        ylab = "N of gene sets", xlab = "empirical p-value", ylim = c(0,60), 
        cex.names = 1.5, cex.axis = 1.2, cex.lab = 1.5)
dev.off()
