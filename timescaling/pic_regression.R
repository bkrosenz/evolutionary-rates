library(dplyr)
library(tidyr)
library(xtable)
library(yaml)
library(ggplot2)
library(hrbrthemes)

regression <- function(formula, data) {
  g = glm(formula = formula,
          data = data,
          family = binomial())
  print(summary(g))
  print(paste('BIC', BIC(g), 'AIC', g$aic))
  g
}

df = read.csv(   
  "/N/project/phyloML/rate_timescaling/data/pic_predictions_no_ils.csv.gz"
)
df$X=NULL
# df['lam_inv'] = 1/df['lam']
# df['lam_minus_mu'] = 1/f['lam']

for (var in c('sigma','sigma_js','sigma_cherry','sigma_paired')){
  
  f=  as.formula(
    paste(
      var, 
      # '~ lam + mu + mu_over_lam + h_mu_over_lam + lam_minus_mu + lam_minus_mu + (n + log(n) ) * rho * tau * ( bl_max + bl_median + bl_min + log(bl_min) + bl_q1 + bl_q3 )')
      '~ lam + mu + mu_over_lam + h_mu_over_lam + h_lam_minus_mu + lam_minus_mu + n + log(n) + rho + tau + bl_max + bl_median + bl_min + log(bl_min) + bl_q1 + bl_q3 ')
    ) # add inverse lam term
  data = df%>% drop_na() %>% filter(n>4 & n< 1000 & get(var)>1)
  
  paste0('MAE_',var)
  data[var] = abs(data[var]-1)
  
  model = glm(
    as.formula(paste(var, '~ bl_min+n')),
    data = data,
    family = 'gaussian'
  )
  model_forward = step(model,
                       k = log(nrow(data)), # turns AIC into BIC
                       scope = f)
  
  print(xtable(model_forward), 
        file = paste0(
          "/N/project/phyloML/rate_timescaling/feature_importances/",
          var,
          '_overest_regression.tex')
    )
  save(model_forward, 
       file = paste0(
          "/N/project/phyloML/rate_timescaling/feature_importances/",
          "forward_model_overest_",
          var,
          ".rda")
       )

  
}

f = formula(sigmas ~ poly(taus, 2) * poly(lam, 2) * poly(mu, 2) * poly(rho, 2) *
              poly(h, 2))


predictions = predictions %>% mutate(across(c(taus, sigmas), matrix2list))

pred_long = predictions %>% unnest_longer(c(sigmas, taus)) # %>% mutate(across(!model, scale))
# 

intercept_only <-
  glm(sigmas ~ 1, data = pred_long %>% filter(model == 'BM'))
# g_eb = regression(accelerating~taus*(lam+mu+rho)*h,data=pac_long%>%filter(model=='EB'))
f = formula(sigmas ~ poly(taus, 2) * poly(lam, 2) * poly(mu, 2) * poly(rho, 2) *
              poly(h, 2))
# g_eb_fw = step(intercept_only,
#                scope=f,
#                direction='forward',
#                criteria='BIC')

ggpairs(
  predictions %>% filter(model == 'BM') %>% filter(taus < .5) %>% select(sigmas, taus, h),
  lower = list(
    continuous = "cor",
    combo = "box_no_facet",
    discrete = "count",
    na = "na"
  ),
  upper = list(
    continuous = "points",
    combo = "facethist",
    discrete = "facetbar",
    na = "na"
  )
)

pp = (predictions %>% filter(model == 'BM') %>% select(sigmas) %>% drop_na())$sigmas

eps = .05


f_mae = function(thresh,winsorize=F) {
  preds = (pred_long %>% filter(model == 'BM') %>% filter(
    (taus > thresh - eps) &
      (thresh + eps > taus) &
      !is.na(sigmas) ) %>% select(sigmas) %>% drop_na()
  )$sigmas
  if (!length(preds)) return(NaN)
  err = abs(preds-1)
  return( mean(if(winsorize) Winsorize(err) else err))
}

x = seq(-.5, 1, length.out = 100)
y = mapply(f_mae, x, MoreArgs = list(winsorize=T))

data = data.frame(x = x, MAE = y)

p <- ggplot(data, aes(x = x, y = MAE)) +
  geom_line(
    color = "#69b3a2",
    size = .7,
    alpha = 0.9,
    linetype = 1
  ) +
  labs(x = expression(tau)) +
  ggtitle("Estimation error (Winsorized)") +   scale_y_continuous(
    trans = 'log10',
    breaks = trans_breaks('log10', function(x)
      10 ^ x),
    labels = trans_format('log10', math_format(10 ^
                                                 .x))
  ) +
  theme(text = element_text(size = 20))
p
figdir = '/N/project/phyloML/rate_timescaling/figures/'


ggsave(paste(figdir, "mae_tau_winsorized_log.pdf", sep = ''),
       width = 10,
       height = 10)


f = function(thresh) {
  rmse(1,
       (
         pred_long %>% filter(model == 'BM') %>% filter(taus > thresh - eps &
                                                          thresh + eps > taus &
                                                          !is.na(sigmas)) %>% select(sigmas) %>% drop_na()
       )$sigmas)
}

f_r = function(thresh, eps=.1) {
  rmse(1,
       (
         pred_long %>% filter(model == 'BM') %>% filter(r > thresh - eps &
                                                          thresh + eps > r &
                                                          !is.na(sigmas)) %>% select(sigmas) %>% drop_na()
       )$sigmas)
}

library(scales)

x = seq(min(pred_long$r), max(pred_long$r), length.out = 100)
y = mapply(f_r, x)


x = seq(-.5, 1, length.out = 100)
y = mapply(f, x)

data = data.frame(x = x, RMSE = y)

p <- ggplot(data, aes(x = x, y = RMSE)) +
  geom_line(
    color = "#69b3a2",
    size = .7,
    alpha = 0.9,
    linetype = 1
  ) +
  labs(x = expression(tau)) +
  ggtitle("Estimation error") +   scale_y_continuous(
    trans = 'log10',
    breaks = trans_breaks('log10', function(x)
      10 ^ x),
    labels = trans_format('log10', math_format(10 ^
                                                 .x))
  ) +
  theme(text = element_text(size = 20))
p
figdir = '/N/project/phyloML/rate_timescaling/figures/'


ggsave(paste(figdir, "rmse_tau_log.pdf", sep = ''),
       width = 10,
       height = 10)