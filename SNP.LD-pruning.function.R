#########################################################################################
### SNP LD-PRUNING ######################################################################
#########################################################################################
## function to be loaded 
SNPr2 <- function(genotypes, flank, association, snp, rsquared) {
  index <- which(association[,2]==association[snp,2] &
                   association[,3]>= association[snp,3]-flank &
                   association[,3]<= association[snp,3]+flank)
  marker <- row.names(association)[index] 
  window <- genotypes[,(colnames(genotypes) %in% marker), drop=F]
  r2 <- cor(window, genotypes[,association[snp,1]], use="pairwise.complete.obs")^2 
  cut <- subset(r2, r2[,1] >= rsquared)
  snp_r2 <- row.names(cut)[!(row.names(cut) %in% association[snp,1])]
  new_genotypes <- genotypes[,!(colnames(genotypes) %in% snp_r2)]
  new_assocation <- association[!(association[,1] %in% snp_r2),]
  results <- list(cut=cut, 
                  new.genotypes=new_genotypes, 
                  new.association=new_assocation,
                  snp=snp_r2)
  return(results)
}
