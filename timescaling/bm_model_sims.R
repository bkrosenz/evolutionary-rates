remove(list = ls())
library(ggplot2)
library(MASS)
library(tidyr)
library(reshape2)
library(rdist)
library(gridExtra)

#This script simulates quantitative traits on a three-taxon tree
#using the specified model parameters, and plots the results. Used to generate
#the covariance matrices and trait distances plotted in Figure 2.

plot_BM <- function(t1, t2, tm, delta2, delta3, sigma2, n_traits) {
   
   if(tm > t1 | t1 > t2) {
      stop("Invalid split times")
   }
   else if((delta2 + delta3) > 1){
      stop("Invalid rate of introgression")
   }
   
   ### Theory ###
   
   pt1_sort <- 1 - exp(-(t2 - t1))
   pt1_ILS <- (1/3)*exp(-(t2 - t1))
   pt2_sort <- 1 - exp(-(t2 - tm))
   pt2_ILS <- (1/3)*exp(-(t2 - tm))
   pt3_sort <- 1 - exp(-(t1 - tm))
   pt3_ILS <- (1/3)*exp(-(t1 - tm))
   no_intro <- 1 - (delta2 + delta3)
   
   cov_AB_1 <- (exp(t2)*(t2 - t1))/(exp(t2) - exp(t1))
   
   cov_AB <- sigma2*(no_intro*(pt1_sort*cov_AB_1 + pt1_ILS)
                     + delta2*pt2_ILS + delta3*pt3_ILS)
   
   cov_BC_1 <- (exp(t2)*(t2 - tm))/(exp(t2) - exp(tm))
   cov_BC_2 <- (exp(t1)*(t1 - tm))/(exp(t1) - exp(tm))
   
   cov_BC <- sigma2*(delta2*(pt2_sort*cov_BC_1 + pt2_ILS) + 
                     delta3*(pt3_sort*cov_BC_2 + pt3_ILS) + 
                     no_intro*pt1_ILS)
   
   cov_AC <- sigma2*(no_intro*pt1_ILS + delta2*pt2_ILS + delta3*pt3_ILS)
   
   var_all_1 <- no_intro*(pt1_sort*(t2 + 1) + exp(-(t2 - t1))*(t2 + 1 + (1/3)))
   var_all_2 <- delta2*(pt2_sort*(t2 + 1) + exp(-(t2 - tm))*(t2 + 1 + (1/3)))
   var_all_3 <- delta3*(pt3_sort*(t1 + 1) + exp(-(t1 - tm))*(t1 + 1 + (1/3)))
   
   var_all <- sigma2*(var_all_1 + var_all_2 + var_all_3)
   
   var_covar <- matrix(c(var_all, cov_AB, cov_AC,
                         cov_AB, var_all, cov_BC,
                         cov_AC, cov_BC, var_all), 
                       nrow = 3, ncol = 3)
   
   var_covar_plot <- var_covar
   
   var_covar_plot[1,1] <- NA
   var_covar_plot[2,2] <- NA
   var_covar_plot[3,3] <- NA
   # var_covar_plot[2,1] <- NA
   # var_covar_plot[3,1] <- NA
   # var_covar_plot[3,2] <- NA
   
   ### Trait value simulations ### 
   
   traits <- data.frame(mvrnorm(n = n_traits, mu = c(0, 0, 0), Sigma = var_covar))
   colnames(traits) <- c("A", "B", "C")
   long_traits <- gather(traits, species, trait_val)
   
   ### Plot variance/covariance matrix ### 
   
   long_varcovar <- melt(var_covar_plot)
   long_varcovar[long_varcovar=="1"] <- "A"
   long_varcovar[long_varcovar=="2"] <- "B"
   long_varcovar[long_varcovar=="3"] <- "C"
   
   long_varcovar$Var1 <- factor(long_varcovar$Var1)
   long_varcovar$Var2 <- factor(long_varcovar$Var2,
                                levels = c("C", "B", "A"))
   
   
   varcovar_plot <- ggplot(long_varcovar, aes(x =  Var1, y = Var2))
   varcovar_plot <- varcovar_plot + geom_raster(aes(fill=value))
   varcovar_plot <- varcovar_plot + scale_fill_gradient(low = "grey90",
                                                        high = "red", 
                                                        na.value = "white")
   varcovar_plot <- varcovar_plot + labs(x = "", y = "", 
                                         title = "Trait variance/covariance matrix")
   varcovar_plot <- varcovar_plot + geom_text(aes(label=round(value,
                                                              digits = 3)))
   varcovar_plot <- varcovar_plot + theme_minimal(base_size = 15)
   varcovar_plot <- varcovar_plot + theme(legend.position = "none",
                                          plot.title = element_text(hjust = 0.5))
   varcovar_plot <- varcovar_plot + scale_y_discrete(limits = rev(levels(long_varcovar$Var2)))

   ### Plot simulated trait differences ### 
   
   AB_dist <- apply(cbind(long_traits[long_traits$species=="A",2], 
                          long_traits[long_traits$species=="B",2]),
                    1, FUN = diff)
   
   BC_dist <- apply(cbind(long_traits[long_traits$species=="B",2], 
                          long_traits[long_traits$species=="C",2]),
                    1, FUN = diff)
   
   AC_dist <- apply(cbind(long_traits[long_traits$species=="A",2], 
                          long_traits[long_traits$species=="C",2]),
                    1, FUN = diff)
   
   AB_dist <- (abs(AB_dist))
   BC_dist <- (abs(BC_dist))
   AC_dist <- (abs(AC_dist))
   

   dists_df <- data.frame(cbind(AB_dist, BC_dist, AC_dist))
   colnames(dists_df) <- c("q1q2", "q2q3", "q1q3")
   write.csv(dists_df, "fig2_intro_sims.csv")
   long_dists <- gather(dists_df, pair, value)

   dists_plot <- ggplot(long_dists, aes(x = pair, y = value, group = 1))
   dists_plot <- dists_plot + stat_summary(fun.data=mean_se,
                                           geom="pointrange",
                                           size = 1)
   
   dists_plot <- dists_plot + scale_x_discrete(
      breaks = c("q1q2", "q1q3", "q2q3"),
      labels = c(expression(q[1]~q[2]), 
                 expression(q[1]~q[3]),
                 expression(q[2]~q[3])))
   
   dists_plot <- dists_plot + theme_bw(base_size = 15)
   dists_plot <- dists_plot + labs(x = "Trait pair",
                                   y = "Mean trait distance",
                                   title = "Simulated pairwise trait distances")
   dists_plot <- dists_plot + theme(plot.title = element_text(hjust = 0.5))
   
   final_plot <- grid.arrange(varcovar_plot, dists_plot, 
                              ncol = 2)
   
   return(final_plot)
}

test <- plot_BM(1, 1.3, 0.5, 0.1, 0, 1, 20000)
