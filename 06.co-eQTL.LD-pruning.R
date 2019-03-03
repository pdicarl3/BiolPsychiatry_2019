##################################################################################################
### co-eQTL LD-PRUNING ###########################################################################
##################################################################################################
### Based on the co-eQTL meta-analysis rank of p-values, perform linkage-disequilibrium SNP pruning
nameData <- ""
### load the results of a co-eQTL analysis (or meta-analysis) with SNP annotations (chromosome, position, etc.)
assoc = assoc  # data.frame [SNPs, annotations and summary statistics]
### the data.frame should have at least the following structure
colnames(assoc) = c("Marker", "Chromosome", "Position", "t_value", "p_value")
### rank SNPs based on p-value (low p-values at the top)
assoc = assoc[order(assoc$p_value),]
c <- assoc  

### load the SNP genotype matrix [samples, SNPs] in allelic dosage mode (0,1,2 as the number of reference alleles) 
geneGenotype = geneGenotype
a <- data.frame(apply(geneGenotype, 2, as.double))

### set parameter for LD-based pruning
f <- 250000  ## base pairs window around the SNP position
rsq <- 0.1   ## r-squared cut-off for LD-based pruning of correlalated SNPs

### load custom fucntion for SNP LD-pruning
source("SNP.LD-pruning.function.R")

### start iteration, start with the top associated SNPs and filter correlated SNPs (above r-squared threshold) in the selected window
### output a data.frame of independent co-eQTLs
for(i in 1:nrow(c)) { 
  if (i <= nrow(c)){
    pruning <- SNPr2(genotypes = a, 
                     flank = f,
                     association = c, 
                     snp = i, 
                     rsquared = rsq)
    a <- pruning$new.genotypes 
    c <- pruning$new.association
  } else { break }
}

### compute correction for multiple comparisons
c$FDR <- p.adjust(c[,ncol(c)], method = "fdr")
c$Bonferroni <- p.adjust(c[,(ncol(c)-1)], method = "bonferroni")
c$rank <- rank(c[,(ncol(c)-2)]) 

### save
assign(paste0("RP_", nameData), c)
save(list = c(paste0("RP_", nameData)), 
     file = paste0("RP_", nameData, ".RData"))
