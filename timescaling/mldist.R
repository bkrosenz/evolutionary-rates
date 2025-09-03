require(ape)

trees = read.tree( '/N/project/phyloML/data/phyloDNN/test/trees_60.tsv')

trees_k80 = c()
dists = c()
for ( tid in 1:length(trees)){
  s=paste('/N/project/phyloML/data/phyloDNN/test/hky_60/', tid-1, '_60_tips.phy', sep='')
  d=dist.dna(read.dna(s), model='K80',as.matrix = T)
  d[is.infinite(d)] <- NA #max(d,na.rm=T)
  d[is.na(d)] <- 9 # 2*max(d,na.rm=T)
  t = bionj(d)
  dists = c(dists,RF.dist(t,trees[[tid]],normalize = T))
  trees_k80 = c(trees_k80,t)
}
