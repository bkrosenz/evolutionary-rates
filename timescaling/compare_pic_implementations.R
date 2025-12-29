# module load r/4.4.1

# invocation:

# Rscript compare_pic_implementations.R  /N/project/phyloML/rate_timescaling/data/tree_sims/l7.0/m6.650/h5 8 1 
#  Rscript compare_pic_implementations.R  /N/project/phyloML/rate_timescaling/data/tree_sims/l7.5/m6.750/h5/no_ils/ 8 1 
# datadir cores genes
# parallel -j8 --shuf Rscript compare_pic_implementations.R {} 8 1 ::: `find /N/project/phyloML/rate_timescaling/data/tree_sims/ -name 'no_ils' -type d -maxdepth 4` 

#### parameters

replicates = 20


args = commandArgs(trailingOnly = TRUE)

print(args)

if (length(args) > 0) {
  datadir = args[1]
  cores = args[2]
} else {
  datadir = '/N/project/phyloML/rate_timescaling/data/tree_sims/l11_m5.5_h1'
}

if (length(args) > 3) {
  ngenes = args[3]
} else {
  ngenes = 1
}


filename <- paste(datadir,
                  "/BM_ests_ngenes",
                  ngenes,
                  "_r",
                  replicates,
                  ".RData",
                  sep = ""
)

if (file.exists(filename)) {
  print(
    cat('file', filename, 'exists, exiting\n', sep = ' '))
  quit()
}

options(mc.cores = cores) ## used by fitContinuous

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

# uses global trees and npz1 variable
results = run_ape_and_geiger(filename)
