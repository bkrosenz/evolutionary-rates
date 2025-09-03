library(ape)
library(reticulate)

# reticulate::use_python('/N/project/phyloML/app/Carb/anaconda/envs/r-reticulate')
# use_condaenv('r-reticulate')
library(rlist)
library(geiger)
# library(lme4)
library(parallel)

# find ../data -maxdepth 1|parallel -j2 --lb Rscript estimate_rates.R  {}

# pstree -hp bkrosenz |grep R


##### functions
height_bl_corr <- function(tree) {
  # TODO: debug.
  #tree=tree[1]
  # print(tree)
  ntaxa = length(tree$tip.label)
  # leaves occupy indices 1:ntaxa
  e = tree$edge
  internal_edges = e[order(e[, 2])[-c(1:ntaxa)], ]
  bl = tree$edge.length[internal_edges[, 1]]
  depths = node.depth.edgelength(tree)
  root <-  setdiff(tree$edge[, 1], tree$edge[, 2])
  cor(bl, depths[-c(1:ntaxa, root)], method = 'spearman')
}

plot_ntaxa <- function(trees) {
  ns = data.frame(lapply(trees, function(t) {
    length(t$tip.label)
  }))
  hist(t(ns), breaks = 20)
}


get_field <- function(x, field) {
  mapply(function(tab)
    mapply(function(row)
      row[[field]], tab), x)
}

#add trend parameter

# EB infers high a and negligible sigsq
# delta infers high delta and high sigsq
# rate trend infers high (at boundary) slope

geiger_control = list(
  method = c("subplex", "L-BFGS-B"),
  niter = 400,
  hessian = TRUE,
  CI = 0.95
)

estimate_trend_geiger <- function(tree,
                                  xx,
                                  ncores = 4,
                                  model = "EB") {
  if (model == "EB") {
    bounds = list(a = c(-60, 60))
  } else if (model == "rate_trend") {
    bounds = list(slope = c(min = -900, max = 900))
  } else if (model == "delta") {
    bounds = list(delta = c(exp(-600), 30))
  } else {
    bounds = list()
  }
  
  bounds = list.append(bounds,
                       sigsq = c(min = exp(-400), max = exp(40)))
  
  fitContinuous(
    tree,
    xx,
    model = model,
    bounds = bounds,
    control = geiger_control,
    ncores = ncores
  )$opt
}

# evorates: mu<0 + low sig -> Early Burst
#           mu>0 + low sig -> Late Burst

estimate_trend <- function(t, xx, chains = 4) {
  library(evorates)
  trend.fit <- fit.evorates(
    tree = t,
    trait.data = xx,
    chains = chains,
    cores = chains,
    trend = TRUE,
    iter = 1500,
    include.warmup = TRUE,
    report.devs = FALSE,
    remove.trend = FALSE
  )
  trend.fit = combine.chains(trend.fit)
  # r=rbind(trend.fit$quantiles[,"R_mu",],trend.fit$means[,"R_mu",])
  # list(means=rowMeans(r), sds=apply(r,1,sd))
  c(trend.fit$quantiles[, 'R_mu'], trend.fit$means["R_mu"])
}

run_geiger_sims <- function(model = 'EB') {
  # uses global variables datadir, ngenes, replicates
  filename = paste(datadir,
                   '/',
                   model,
                   "_estimates_ngenes",
                   ngenes,
                   "_r",
                   replicates,
                   ".RData",
                   sep = "")
  results = list()
  for (sp_tree in 1:min(length(npz1),length(trees)) ){
    arr = paste("arr_", sp_tree - 1, sep = "")
    x <- npz1$f[[arr]][1:replicates, ]
    t <- trees[[sp_tree]]
    colnames(x) = sapply(1:dim(x)[2],
                         function(i)
                           paste("T", i, sep = ""))
    r = mclapply(asplit(x, 1),
                 function(xx)
                   estimate_trend_geiger(t, xx, model = model),
                 mc.cores = 4)
    # r = apply(x, 1, function(xx) estimate_trend_geiger(t,xx,model = model))
    results = list.append(results, r)
    # if (! sp_tree%%50) {
    #   save(results,file=filename)
    # }
  }
  print(paste('finished ', sp_tree))
  save(results, file = filename)
  results
}


run_sims <- function() {
  filename = paste(datadir,
                   "/R_mu_estimates_ngenes",
                   ngenes,
                   "_r",
                   replicates,
                   ".RData",
                   sep = "")
  results = list()
  # try ( load(filename))
  for (sp_tree in 1:length(trees)) {
    arr = paste("arr_", sp_tree - 1, sep = "")
    x <- npz1$f[[arr]][1:replicates, ]
    t <- trees[[sp_tree]]
    colnames(x) = sapply(1:dim(x)[2],
                         function(i)
                           paste("T", i, sep = ""))
    
    r = apply(x, 1,
              function(xx)
                estimate_trend(t, xx))
    results = list.append(results, r)
    if (!sp_tree %% 3) {
      save(results, file = filename)
    }
    print(paste('finished ', sp_tree))
  }
  save(results, file = filename)
}

np <- import("numpy")
