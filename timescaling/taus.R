source('analysis_functions.R')
library(data.table)
library(dplyr)
library(zeallot)
library(xtable)
library(Metrics)
library(tidyr)
library(DescTools)
library(crunch)
library(runner)

replicates = 20

get_full_results = function(globstr, ngenes, sigmas = F) {
  predictions = data.frame()
  for (model in c('EB', 'rate_trend', 'delta', 'BM')) {
    for (d in Sys.glob(globstr, dirmark = FALSE)) {
      filename = paste(d,
                       '/',
                       model,
                       "_estimates_ngenes",
                       ngenes,
                       "_r",
                       replicates,
                       ".RData",
                       sep = "")
      tau_file = paste(d, 'taus.npy', sep = '/')
      # print(paste(filename,tau_file))
      if (file.exists(tau_file)) {
        taus = np$load(tau_file)[, 1:replicates]
      } else {
        tau_file = paste(d, 'no_ils', 'taus.npy', sep = '/')
        if (file.exists(tau_file)) {
          taus = np$load(tau_file)[, 1:replicates]
        }
        
        else{
          next
        }
      }
      if (file.exists(filename)) {
        load(filename)
        # print(results[[1]][[1]]$CI)
        e = get_all_ml_estimates(results, model)
        c(lam, mu, h) %<-% sapply(tail(strsplit(gsub(
          '/no_ils', '', d
        ) , '/')[[1]], 3),
        function(x)
          as.numeric(substr(x, 2, 100)))
        if (sigmas) {
          a = list(sigmas = e$sigs)
        } else{
          a = get_accel(e, model, reduce = F) # if as_mean==F, will be nested list
        }
        
        a = list.append(
          a,
          lam = lam,
          mu = mu,
          h = h,
          model = model,
          taus = taus[1:(dim(a$sigmas)[1]), ]
        )
        predictions = bind_rows(predictions, a)
      }
    }
  }
  rdf = predictions
  rdf$r = rdf$lam - rdf$mu
  rdf$rho = rdf$mu / rdf$lam
  return(rdf)
}

l = '2.0'
m = '1.900'
h = '3'
path = paste(
  '/N/project/phyloML/rate_timescaling/data/tree_sims/l',
  l,
  '/m',
  m,
  '/h',
  h,
  '/no_ils',
  sep = ''
)

load(paste(path, 'EB_estimates_ngenes1_r20.RData', sep = '/'))
e = get_all_ml_estimates(results, model)
a = get_accel(e, model, reduce = F)
taus = np$load(paste(path, 'taus.npy', sep = '/'))[, 1:20]

df = data.frame('acc' = as.vector(a$accelerating), 'tau' = as.vector(taus))
df = drop_na(df)


ga = glm(formula = acc ~ tau,
         data = df,
         family = binomial)
summary(ga)
newdata <- data.frame('tau' = seq(min(df$tau), max(df$tau), len = 500))
newdata$vs = predict(ga, newdata, type = "response")
plot(
  acc ~ tau,
  data = df,
  col = "steelblue",
  xlim = c(0, max(df$tau))
)
lines(vs ~ tau, newdata, lwd = 2)

cor(df, method = 'spearman')
print(path)
print(length(taus))
print(colMeans(df))
print(dim(df))

cors = c()
for (i in 1:59) {
  cors = c(cors, cor(taus[i, ], a$accelerating[i, ], method = 'spearman'))
}
hist(cors)

colMeans(df %>% filter(acc == 0))
colMeans(df %>% filter(acc == 1))

# TODO: binomial glm for NO ILS with all sims, all conditions
g <- glm(
  formula = acc ~ tau + lam + mu + r * h + rho * h,
  data = df,
  family = binomial
)

######## full reg

matrix2list = function(foo) {
  lapply(seq_len(nrow(foo)), function(x)
    foo[x, ])
}
pac = get_full_results('/N/project/phyloML/rate_timescaling/data/tree_sims/l*/m*/h*/no_ils',
                       1)
save(pac, file = "/N/project/phyloML/rate_timescaling/data/tree_sims/geiger_res_no_ils.Rda")           # Save data frame
load(file = "/N/project/phyloML/rate_timescaling/data/tree_sims/geiger_res_no_ils.Rda")           # Save data frame

pac = pac %>% mutate(across(c(accelerating, CI, taus, value), matrix2list))
pac_long = pac %>% unnest_longer(c(accelerating, CI, value, taus))
# TODO: add tree index 1...20
write.csv.gz(
  pac_long,
  '/N/project/phyloML/rate_timescaling/data/tree_sims/geiger_res_no_ils.csv.gz'
)
pac_long = read.csv('/N/project/phyloML/rate_timescaling/data/tree_sims/geiger_res_no_ils.csv.gz')
pac_long = pac_long %>% mutate(across(!accelerating &
                                        !CI & !value & !model, scale))

regression <- function(formula, data) {
  g = glm(formula = formula,
          data = data,
          family = binomial())
  print(summary(g))
  print(paste('BIC', BIC(g), 'AIC', g$aic))
  g
}

g = regression(accelerating ~ taus * lam * mu * h * model + h * rho)
g_back = step(g, criteria = 'BIC')
print(xtable(gback)) #, file = "filename.tex", compress = FALSE)

intercept_only <-
  glm(accelerating ~ 1,
      data = pac_long %>% filter(model == 'EB'),
      family = binomial)
# g_eb = regression(accelerating~taus*(lam+mu+rho)*h,data=pac_long%>%filter(model=='EB'))
f = formula(accelerating ~ poly(taus, 2) * poly(lam, 2) * poly(mu, 2) *
              poly(rho, 2) * poly(h, 2))

g_eb_fw = step(
  intercept_only,
  scope = f,
  direction = 'forward',
  criteria = 'BIC'
)

print(xtable(g_eb_fw), file = "eb_model.tex", compress = FALSE)


g_rt = regression(accelerating ~ taus * (lam + mu + rho) * h,
                  data = pac_long %>% filter(model == 'rate_trend'))
g_rt_fw = step(g_rt, method = 'seqrep', criteria = 'BIC')



###### sigmas only
predictions = get_full_results('/N/project/phyloML/rate_timescaling/data/tree_sims/l*/m*/h*',
                               25,
                               sigmas = T)

predictions = get_full_results(
  '/N/project/phyloML/rate_timescaling/data/tree_sims/l*/m*/h*/no_ils',
  1,
  sigmas = T)
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
