#########################################################################################################
### co-eQTL - association between module eigengene and SNP genotypes ####################################
#########################################################################################################
nameModule = "darkgreen"
nameData = "LIBD"
f = "lmRob"
CPU = 4 # set number of CPU for parallel computing
### load a data.frame with module eigengenes (ME) [samples, ME] containing also demographics information and important covariates
pd = pd  # demo-eigengene data.frame
### laod SNP genotype data [samples, SNP]
GEN =  data.frame(GEN)  
                        
### parallel adjust for covariates - get the p-value of the ME-SNP association
### use robust linear model
require(snow)
require(robust)
require(car)
clus = makeCluster(CPU)
clusterExport(clus, c("pd", "GEN", "lmRob")) ## if you need to export a variable
ptm = proc.time()
###
QTL = parLapply(clus, GEN, function(G) 
  summary(lmRob(darkgreen ~ Age + RIN + mitoMapped + totalMapped + as.factor(Sex) +  ## your covariates
                  EV1 + EV2 + EV3 + EV4 + EV5 +                                      ## your covariates
                  EV6 + EV7 + EV8 + EV9 + EV10 +                                     ## your covariates
                  as.factor(Dx) +                                                    ## your covariates
                  G,                                                                 ## SNP genotypes
                data = pd))$coefficients) 

t1 = proc.time() - ptm
stopCluster(clus)

### aggregate results (get the last row of coefficients table)
clus = makeCluster(CPU)
clusterExport(clus, c("QTL"))
coeQTL = parSapply(clus, QTL, function(l) l[nrow(l),])
stopCluster(clus)
### give variables name
row.names(coeQTL) = c("beta", "st_err", "t_value", "p_value")
coeQTL = data.frame(t(coeQTL))

### multiple comparisons correction
coeQTL$FDR = p.adjust(coeQTL$p_value, method = "fdr")
coeQTL$Bonferroni = p.adjust(coeQTL$p_value, method = "bonferroni")

### assign
assign(paste0("QTL_", f, nameData), QTL)
assign(paste0("coeQTL_", f, nameData), coeQTL)

### save
save(list=c(paste0("QTL_", f, nameData), 
            paste0("coeQTL_", f, nameData)), 
     file=paste0("coeQTL_", f, nameData, ".RData"))