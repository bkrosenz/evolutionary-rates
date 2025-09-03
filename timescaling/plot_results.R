require(ape)

options(mc.cores=4) ## used by fitContinuous

replicates = 50
ngenes = 50

args = commandArgs(trailingOnly=TRUE)
if (length(args)>0){
  datadir = args
} else {
  datadir = '../data/m0.8_h9_s5/'
}

trees<- read.tree(paste(datadir,'parent_trees_4plus.phy',sep=''))

load(paste(datadir,'EB_a_estimates_ngenes50_r50.RData'))
d=results[[1]]
d=as.data.frame(do.call(rbind, d))