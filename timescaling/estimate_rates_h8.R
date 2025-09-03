library(evorates)
library(ape)
library(reticulate)
# library('rlist')
library(lme4)
library(parallel)

# trees<-lapply,(filename))
np <- import("numpy")

#add trend parameter

estimate_trend <- function(t,xx,chains=4){
trend.fit <- fit.evorates(
  tree = t,
  trait.data =xx,
  chains = chains,
  cores=chains,
  trend = TRUE,
  include.warmup=TRUE,
  report.devs=FALSE,
  remove.trend = FALSE)
trend.fit=combine.chains(trend.fit)
# r=rbind(trend.fit$quantiles[,"R_mu",],trend.fit$means[,"R_mu",])
# list(means=rowMeans(r), sds=apply(r,1,sd))
c(trend.fit$quantiles[,'R_mu'],trend.fit$means["R_mu"])
}

datadir='../data/m0.4_h8/'
print(datadir)
trees<- read.tree(paste(datadir,'parent_trees.phy',sep=''))
ngenes=500
npz1 <- np$load(paste(datadir, "samples_ngenes",ngenes,".npz",sep=""))
npz1$files

replicates = 10
filename=paste(datadir, "R_mu_estimates_ngenes", ngenes,"_r", replicates, ".RData", sep="")
results= list()
for (sp_tree in 1:length(trees)){
    arr=paste("arr_",sp_tree-1,sep="")
    x<-npz1$f[[arr]][1:replicates,]
    t<-trees[[sp_tree]]
    colnames(x)=sapply(1:dim(x)[2],function(i)paste("T",i,sep=""))

    r = apply(x, 1, function(xx) estimate_trend(t,xx))
    results = list.append(results,r)
    if (! sp_tree%%3) {
        save(results,file=filename)
    }
    print(paste('finished ',sp_tree))
}
save(results,file=filename)