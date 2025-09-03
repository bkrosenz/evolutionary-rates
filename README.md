# evolutionary-rates
Most of the figures in the paper were generated in the two jupyter notebooks ``Bias and Variance of sigma estimators.ipynb`` and ``Birth Death Trees.ipynb`` 
simulation/ contains utils for simulating BD trees with dendropy
the timescaling/ folder contains R code for running the time-varying rates models in geiger and evorates
``calculate_pic_var.py`` contains symbolic algebra code for computing the bias and variance of the sigma estimators for a given phylogenetic covariance matrix, where all tree lengths are free variables.
