find /N/project/phyloML/rate_timescaling/data/tree_sims/ -maxdepth 3 -mindepth 3 -type d -exec bash -c  "if diff {}/parent_trees_4plus_pos_edges.phy {}/parent_trees_4plus.phy >/dev/null; then true ; else echo \"{}\";fi" \; > need_to_redo_estimates.txt

parallel -j24 rm {}/no_ils/*RData :::: need_to_redo_estimates.txt

parallel --shuf -j2 Rscript estimate_rates.R  {1}/no_ils {2} 12 1 :::: need_to_redo_estimates.txt ::: EB delta rate_trend BM
parallel --shuf -j2 Rscript estimate_rates.R  {} 12 10 :::: need_to_redo_estimates.txt ::: EB delta rate_trend BM

