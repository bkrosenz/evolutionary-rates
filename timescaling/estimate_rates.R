# invocation:
# parallel -j2 Rscript estimate_rates.R  {} 8 ::: `find /N/project/phyloML/rate_timescaling/data/tree_sims/ -name 'h*' -type d -maxdepth 3` ::: rate_trend EB delta BM

# parallel -j2 Rscript estimate_rates.R  {} 6 1 ::: `find /N/project/phyloML/rate_timescaling/data/tree_sims/ -name 'no_ils' -type d -maxdepth 4` ::: rate_trend EB delta BM

#### parameters

replicates = 20


args = commandArgs(trailingOnly = TRUE)

print(args)

if (length(args) > 0) {
  datadir = args[1]
  model = args[2]
  cores = args[3]
} else {
  datadir = '/N/project/phyloML/rate_timescaling/data/tree_sims/l11_m5.5_h1/'
  model = 'rate_trend'
}

if (length(args) > 3) {
  ngenes = args[4]
} else {
  ngenes = 10
}



filename = paste(datadir,
                 '/',
                 model,
                 "_estimates_ngenes",
                 ngenes,
                 "_r",
                 replicates,
                 ".RData",
                 sep = "")

if (file.exists(filename)) {
  print(
    cat('file', filename, 'exists, exiting\n', sep = ' '))
  quit()
}

options(mc.cores = cores) ## used by fitContinuous
print(paste(datadir, model, sep = ' '))

#quit()

setwd('/N/project/phyloML/rate_timescaling/evorates')
source("functions.R")
trees <- tryCatch(
  read.tree(paste(datadir, '/parent_trees_4plus_pos_edges.phy', sep = '')),
  error = function(e) {
    read.tree(paste(datadir, '/../parent_trees_4plus_pos_edges.phy', sep = ''))
  },
  warning = function(e) {
    read.tree(paste(datadir, '/../parent_trees_4plus_pos_edges.phy', sep = ''))
  }
)

print(trees)

sample_file_name = paste(datadir,
                         "/samples_ngenes", ngenes, ".npz", sep =
                           "")

if (!file.exists(sample_file_name)) {
  print(cat('no samples in', sample_file_name, ', exiting\n', sep = ''))
  quit()
}

npz1 <- np$load(sample_file_name)


# uncomment this line to run evorates
# run_sims()
# filename=paste(datadir,
# "R_mu_estimates_ngenes", ngenes,"_r", replicates, ".RData",
# sep="")



results = run_geiger_sims(model)

quit()

#aicc = sapply(rdf['aicc',],as.numeric)
#     for (var in c("a", "sigsq",)){
#       v=rdf[var,aicc_max]
#         v = mean(sapply(rdf[var,],as.numeric))
#       }
#     z <- as.numeric(lapply(results[[1]] , function(x){mean(x$z0)}))
#   }
#   m <- as.numeric(lapply(results , function(x){mean(x['R_mu',])}))
#
# }




th <- as.numeric(lapply(trees[1:length(m)], height_bl_corr))
print(
  'R_mu:',
  m,
  'height bl correlation:',
  th,
  'Cor(mu,th)',
  cor(m, th, method = 'spearman', use = 'complete.obs')
)
