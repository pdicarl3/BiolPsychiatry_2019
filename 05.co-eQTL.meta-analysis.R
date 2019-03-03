##############################################################################################################
### META ANALYSIS of two co-eQTL analyses (different datasets) ###############################################
##############################################################################################################
### load two co-eQTL analyses (data.frames with SNPs as rows containing summary statitstics of co-eQTL analysis: t-value, standard error, etc.)
nameData = ""
results = cbind(coeQTL.1, coeQTL.2) # combine the two data.frames
CPU = 4
meth = "FE"               # model method 'FE' fixed-effect
nameMethod = "fixed_"

### get the effect size measure 
source("t2r.R")
results$r.1 = t2r(results$t.value.1, df1) ## trasform t-value to Pearson's r  ## df, degree of freedom
results$r.2 = t2r(results$t.value.2, df2) ## trasform t-value to Pearson's r
results$SEr.1 = rsq2se(results$r.1, n1)  ## standard error  ## n, sample size
results$SEr.2 = rsq2se(results$r.2, n2)  ## standard error

### perform meta-analysis in parallel
require(snow)
require(metafor)
clus = makeCluster(CPU)
clusterExport(clus, c("results", "rma.uni", "meth")) ## if you need to export a variable
ptm = proc.time()
MA = data.frame(QEp = parRapply(clus, results, function(x)
  rma.uni(yi = c(x["r.1"], x["r.2"]), 
          sei = c(x["SEr.1"], x["SEr.2"]), 
          method = meth, measure = "COR")$QEp),
  b = parRapply(clus, results, function(x)                ## meta-analysis effect size
    rma.uni(yi = c(x["r.1"], x["r.2"]), 
            sei = c(x["SEr.1"], x["SEr.2"]), 
            method = meth, measure = "COR")$b),
  ci.lb = parRapply(clus, results, function(x)            ## confidence interval lower bound 
    rma.uni(yi = c(x["r.1"], x["r.2"]), 
            sei = c(x["SEr.1"], x["SEr.2"]), 
            method = meth, measure = "COR")$ci.lb),
  ci.ub = parRapply(clus, results, function(x)            ## confidence interval upper bound 
    rma.uni(yi = c(x["r.1"], x["r.2"]), 
            sei = c(x["SEr.1"], x["SEr.2"]), 
            method = meth, measure = "COR")$ci.ub),
  zval = parRapply(clus, results, function(x)             ## meta-analysis Z-value
    rma.uni(yi = c(x["r.1"], x["r.2"]), 
            sei = c(x["SEr.1"], x["SEr.2"]), 
            method = meth, measure = "COR")$zval),
  pval = parRapply(clus, results, function(x)             ## meta-analysis p-value
    rma.uni(yi = c(x["r.1"], x["r.2"]), 
            sei = c(x["SEr.1"], x["SEr.2"]), 
            method = meth, measure = "COR")$pval)
)
t1 = proc.time() - ptm
stopCluster(clus)                

### save
if(identical(row.names(results), row.names(MA))){
  results_MA = cbind(results, MA)
} else { "row.names do not match"}
results_MA$FDR = p.adjust(results_MA$pval, method = "fdr")
results_MA$Bonferroni = p.adjust(results_MA$pval, method = "bonferroni")
assign(paste0("MA_", nameMethod, nameData), results_MA)
save(list=c(paste0("MA_", nameMethod, nameData)), 
     file=paste0("MA_", nameMethod, nameData, ".RData"))     


