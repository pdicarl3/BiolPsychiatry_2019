#########################################################################################
### POLYGENIC CO-EXPRESSION INDEX (PCI) computation #####################################
#########################################################################################
### load function to compute A-prime 
source("PCI_function.R")

### name data
database = ""
PCI = "PCI"
nameData = ""

### load the vector of co-expression values (module eigengene) of interest
## single vector of trait variable (expression, co-expression, etc.) ## one value for each subject in your sample
expr <- module_eigengene            
### load the SNP genotype matrix [samples, SNPs] in allelic dosage mode, plink format (0,1,2 as the number of reference alleles (minor allele)) 
### used the top 150 independent co-eQTL to calculate your weights
gen <- data.frame(gen)   

### PCI computation start
## computate weights for each SNP-genotype
## output a data.frame with three columns (SNP, genotype, weight), these weights are then used to calculete PCI
AA <- NULL
for(i in colnames(gen)) {
  for(l in levels(as.factor(gen[,i]))) {
    ## the homozygouse of the major allele homozygouse (coded as 0, zero copies of the minor allele)
    hit <- pnorm(mean(expr[as.character(gen[,i]) == "0"], na.rm=TRUE),  
                 mean = mean(expr[gen[,i]==l], na.rm=TRUE),
                 sd = sd(expr[gen[,i]==l], na.rm=TRUE)
    )
    fa <- pnorm(mean(expr[as.character(gen[,i]) == "0"], na.rm=TRUE),
                mean = mean(expr[as.character(gen[,i]) == "0"], na.rm=TRUE),
                sd = sd(expr[as.character(gen[,i]) == "0"], na.rm=TRUE)
    )
    
    a <- aprime(hit, fa)                 ## here, use the sourced function
    names(a) <- paste(i, l, sep=" ")
    AA <- c(AA, a)
  }
}
AA <- data.frame(AA)
y <- strsplit(row.names(AA), " ")   
AA$snp <- sapply(y, "[", 1)
AA$genotype <- sapply(y, "[", 2)

### REPLACE WEIGHTS TO GENOTYPES
### this could be used also to replace weights in an independent dataset (be careful, in an independent dataset the concordance of the reference allele should be checked)
GS_score <- data.frame(apply(gen, 2, function(d) as.factor(d)))
for(i in colnames(GS_score)) {
  levels(GS_score[,i]) <- AA[AA[,2]==i, 1]
}
### output a data.frame with SNP-genotype weights [samples,SNPs]
GS_score <- data.frame(apply(GS_score, 2, as.double))
row.names(GS_score) <- row.names(gen)

#########  !!!!!!!
### in the following, THIS SCRIPT COMPUTES THE PCI WITH ALL THE SNPs SELECTED IN THE INPUT (e.g., 150 top ranked co-eQTLs)
### TO COMPUTE PCIs WITH ANY NUMBER OF SNPs, JUST SUBSET THE "GS_score" data.frame AS YOUR CONVENIENCE (e.g., calculating the PCI with the top 2 co-eQTLs, top 3, top 4, ect.) 
#########

### calculate the PCI as the mean of SNP-genotype weights for each subject (rows)
### the PCI is defined by custom function as negatively correlated to the module eigengene (negatively indexing co-expression) 
GS_score <- transform(GS_score, GS_neg = rowMeans(GS_score, na.rm=TRUE))
### multiply by -1 to obtain a PCI positevely indexing co-expression
GS_score$GS_pos <- GS_score$GS_neg*(-1) 

### save
assign(paste(PCI, database, nameData, sep="."), GS_score)
assign(paste("geno", database, nameData, sep="."), gen)
assign(paste("AA", database, nameData, sep="."), AA)
save(list = c(paste(PCI, database, nameData, sep="."),
              paste("geno", database, nameData, sep="."),
              paste("AA", database, nameData, sep=".")),
     file = paste0(PCI, "_", database, "_", nameData, ".RData"))
