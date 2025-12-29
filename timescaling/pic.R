pic.X <- pic(X, tree, scaled=T)
get_mbl <- function (t){
  p=pic.X[as.character(t$name)]
  mbl = mean(t$edge.length)
  c(p,mbl)
}

s = subtrees(tree)
r = mapply(get_mbl, s)


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


