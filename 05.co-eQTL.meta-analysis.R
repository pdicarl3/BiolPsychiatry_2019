##############################################################################################################
### META ANALYSIS of two co-eQTL analysis (different datasets) ###############################################
##############################################################################################################
### set results database
nameData = ""
results = cbind(coeQTL.1, coeQTL.2)
CPU = 4
meth = "FE"               # model method 'FE' fixed-effect
nameMethod = "fixed_"

### get the effect size measure 
source("t2r.R")
results$r.1 = t2r(results$t.value.1, df) ## follow t-value to Pearson's r
results$r.2 = t2r(results$t.value.2, df) ## follow t-value to Pearson's r
results$SEr.1 = rsq2se(results$r.1, n)  ## standard error
results$SEr.2 = rsq2se(results$r_2, n)  ## standard error

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
  b = parRapply(clus, results, function(x)
    rma.uni(yi = c(x["r.1"], x["r.2"]), 
            sei = c(x["SEr.1"], x["SEr.2"]), 
            method = meth, measure = "COR")$b),
  ci.lb = parRapply(clus, results, function(x)
    rma.uni(yi = c(x["r.1"], x["r.2"]), 
            sei = c(x["SEr.1"], x["SEr.2"]), 
            method = meth, measure = "COR")$ci.lb),
  ci.ub = parRapply(clus, results, function(x)
    rma.uni(yi = c(x["r.1"], x["r.2"]), 
            sei = c(x["SEr.1"], x["SEr.2"]), 
            method = meth, measure = "COR")$ci.ub),
  zval = parRapply(clus, results, function(x)
    rma.uni(yi = c(x["r.1"], x["r.2"]), 
            sei = c(x["SEr.1"], x["SEr.2"]), 
            method = meth, measure = "COR")$zval),
  pval = parRapply(clus, results, function(x)
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


