# invocation:
# parallel -j2 Rscript estimate_rates.R  {} 8  ::: `find /N/project/phyloML/rate_timescaling/data/tree_sims/ -name 'h*' -type d -maxdepth 3` ::: rate_trend EB delta BM
#### parameters
require(tidyverse)
# require(crunch)
require(zeallot)
library(rlist)
library(dplyr)
library("GGally")                      # Load GGally package
p_ <- GGally::print_if_interactive
library('ggplot2')
replicates = 20

ngenes = 1 #0 

args = commandArgs(trailingOnly=TRUE)
if (length(args)>0){
  datadir = args[1]
  model=args[2]
  cores=args[3]
} else {
  datadir = '/N/project/phyloML/rate_timescaling/data/tree_sims/l2.5/m1.750/h7'
  model='rate_trend'
  
}
output_dir = 'output'
filename=paste(datadir,'/',
               model,"_estimates_ngenes", 
               ngenes,"_r", replicates, ".RData", 
               sep="")


load(filename)


############ functions ###########

source('analysis_functions.R')


# USE ALL THE DATAPOINTS:
predictions = get_full_results(
  data_dir,
  genes = 0 # 25
  )

# for NO ILS condition Pr(Accelerating) \propto H*(lambda-mu):

predictions<-  predictions %>% mutate(mae=abs(1-sig)) #hr=h*r,hrho=h*rho,
predictions %>% write_csv("../data/evorates_preds.csv.gz")

algnames =  c('EB','rate_trend','delta')
for (m in algnames ){
  cor_mat = predictions%>% filter(model==m) %>% drop_na(sig) %>% select(!model) %>% as.matrix %>%
    cor(method="spearman") %>%
    as.data.frame
  print(m)
  print(cor_mat)
  d = predictions %>% filter(model==m)
  print(mean(d$sig,na.rm=T))
  ggplot(predictions %>% filter(model==m), 
         aes(sig,fill=as.factor(model))) + geom_histogram() + scale_x_log10()
}


ggplot(
  predictions %>% mutate(model=as.factor(model)) %>% filter(hr<15),# %>%filter(model=='rate_trend'),
       aes(x=hr,y=mae,
           colour=model,
           shape=model)
           )  + 
  # geom_point(size=1) + 
  geom_smooth(level=0.50, method='gam')+ scale_y_log10()  +
  xlab(expression(italic(H)(lambda - mu))) +
  ylab(expression(abs(hat(sigma)[0]^2 - sigma[0]^2))) 
ggsave('../figures/timescaling/mae.png')



ggplot(
  predictions %>% mutate(model=as.factor(model)) %>% filter(hr<15),
  aes(x=hr,
      y=accelerating,
      colour=model,
      shape=model)
       )  +  
  # geom_point(size=1) +  
  geom_smooth(linetype="dashed",level=.5, method='gam')+ 
  xlab(expression(italic(H)(lambda - mu))) +
  ylab('Pr(Accelerating)')
ggsave('../figures/timescaling/accelerating.png')

predictions$h_inv=1/predictions$h

pac= predictions %>%select(accelerating, r,rho,h,h_inv,model)


pac$accelerating=as.data.frame(pac$accelerating)
pac=pac%>%unnest(accelerating) %>% pivot_longer(cols=starts_with('V'),values_to = 'accelerating',names_to = 'block') %>% drop_na()

g = glm(formula=accelerating ~ model+(r+rho)*h+h_inv+block,data = pac,family=binomial) 
# a = aov(as.formula('accelerating ~ model+r+rho+h+block'),data = pac) # not for binary data!
summary(g)

predictions$std = apply(predictions$value,1, function(r) sd(r,na.rm = T))

predictions$mean_value = apply(predictions$value,1, function(r) mean(r,na.rm = T))

save(predictions, 
     file = paste(output_dir,'preds_ils.Rdata', sep='/'))
predictions_no_ils = get_full_results(
  paste(datadir,'*/*/h*/no_ils',sep='/'),
  ngenes = 1)
save(predictions_no_ils, 
     file = paste(output_dir,'preds_no_ils.Rdata', sep='/'))

##################### use means
rdf = analyze_results(
  paste(datadir,'*/*/h*/no_ils',sep='/'),
  ngenes = 10)
# rdf=rdf[!is.na(rdf$CI),]

rdf$r_h = rdf$r * rdf$h
rdf$rho_h = rdf$rho* rdf$h
rdf$r_over_h = rdf$r / rdf$h
rdf$rho_over_h = rdf$rho / rdf$h
save(rdf, 
     file = paste(datadir,'rdf_ils.Rdata',sep='/'))

rdf_null = analyze_results(
  paste(datadir,'*/*/h*/no_ils',sep='/'),
  ngenes = 1)
rdf_null$r_h = rdf_null$r * rdf_null$h
rdf_null$rho_h = rdf_null$rho* rdf_null$h
rdf_null$r_over_h = rdf_null$r / rdf_null$h
rdf_null$rho_over_h = rdf_null$rho / rdf_null$h

save(rdf_null, 
     file = paste(datadir,'rdf_no_ils.Rdata',sep='/'))


figdir='/N/project/phyloML/rate_timescaling/figures/'

for (m in c('delta','rate_trend','EB')){
  # title=paste("Model: ",m, "(with ILS, scaled)")
  title=paste("Model: ",m, "(with ILS)")
  pm<-ggpairs(
    rdf %>% filter(model==!!m)%>%select(c(accelerating,CI,r,rho,h)),
    lower = list(continuous = "cor", combo = "box_no_facet", discrete = "count", na = "na"),
    upper = list(continuous = "points", combo = "facethist", discrete = "facetbar", na =
                   "na"),
    title=title  )
  p_(pm)
  ggsave(paste(figdir,title,".pdf"),width=11,height = 11)
  dev.off()
  
}

cor(rdf[rdf$model=='EB', -7], method = 'spearman')

figdir='/N/project/phyloML/rate_timescaling/figures/'

make_pair_plots = function(df,scaling,condition){
  for (m in c('delta','rate_trend','EB')){
    title= if(scaling) paste("Model: ",m, "(",condition,", scaled)") else paste("Model: ",m, "(",condition,")")
    if(scaling){
      pm <- ggpairs(
        
        df  %>% filter(model==!!m)%>%select(c(accelerating,CI,r_h,rho_h)),
        lower  = list(continuous = "cor",
                      combo = "box_no_facet", discrete = "count", na = "na"),
        upper = list(continuous = "points",
                     combo = "facethist", discrete = "facetbar", na =
                       "na"),
        title=title  )}
    else{
      pm <- ggpairs(
        df %>% filter(model==!!m)%>%select(c(accelerating,CI,r,rho,h)),
        lower  = list(continuous = "cor",
                      combo = "box_no_facet", discrete = "count", na = "na"),
        upper = list(continuous = "points",
                     combo = "facethist", discrete = "facetbar", na =
                       "na"),
        title=title  )
    }
    p_(pm)
    ggsave(paste(figdir,title,".pdf"),width=11,height = 11)
    dev.off()
    
  }
  
}

for (scaling in c(T,F)){
  make_pair_plots(rdf,scaling,'with ILS')
  make_pair_plots(rdf_null,scaling,'no ILS')
  }

cor(rdf_null[rdf_null$model=='EB', -7], method = 'spearman')

f=na.omit(
  inner_join(rdf, rdf_null,
             by=c('lam','mu','h','r','rho','model','r_h','rho_h'),
             suffix=c('.ils','.null')
             ))
mean(f %>% select(CI.ils) > f%>% select( CI.null,))
mean(f %>% select(mean_value.ils) > f%>% select( mean_value.null,))
mean(f %>% select(accelerating.ils) > f%>% select( accelerating.null,))

f$h_inv=1/f$h
for (m in c('delta','rate_trend','EB')){
  pairs(
    f %>% filter(model==!!m)%>%select(c(acc_greater,CI,r,rho,h)),
    main=paste("Model: ",m, "()")  )
}
f$acc_greater=f %>% select(accelerating.ils) > f%>% select( accelerating.null,)
f$CI_greater=f %>% select(CI.ils) > f%>% select( CI.null,)

acc_model = glm(formula=acc_greater ~ model+(r+rho+lam+mu)*h_inv,data = f,family=binomial)
summary(acc_model)
acc_model = glm(formula=acc_greater ~ model+(r+rho+lam+mu)*h,data = f,family=binomial)
summary(acc_model)

acc_model = glm(formula=acc_greater ~ model+rho+r*h,data = f,family=binomial)
ga = glm(formula=CI_greater ~ model+(rho+r+lam+mu)*h,data = f,family=binomial)
summary(ga)
ga = glm(formula=CI_greater ~ model+(r+rho+lam+mu)*h_inv,data = f,family=binomial)
summary(ga)
rh_model = glm(formula=acc_greater ~ rh,data = f,family=binomial)
newdata <- data.frame('rh'=seq(min(f$r_h), max(f$r_h)+2,len=500))
newdata$vs = predict(rh_model, newdata, type="response")
plot(acc_greater ~ r_h, data=f, col="steelblue",xlim=c(0,max(f$r_h)+2))
lines(vs ~ rh, newdata, lwd=2)

f %>% filter(CI.ils-CI.null>.7) %>% select(accelerating.ils,accelerating.null,CI.ils,CI.null,lam,mu,h,r,rho,r_h,rho_h)

compare_stretched_trees <- function(d,model,ngenes=10,title='Scaled Trees'){
  filename=paste(d,'/',
                 model, "_estimates_ngenes", 
                 ngenes, "_r", replicates, ".RData", 
                 sep="")
  
  load(filename)
  # print(results[[1]][[1]]$CI)
  e = get_all_ml_estimates(results,model)
  a = get_accel(e, model,reduce='mean')
  hist(a$accelerating,xlab='Pr(Acceleration)',breaks=10,xlim=c(0,1),ylim=c(0,95),main=title)
}

map <- c("delta" = "Delta", "EB" = "ACDC", "rate_trend" = "LT")

preds$model=map[preds$model]
cp = preds %>% filter(!is.na(CI)) %>% group_by(model,hr)%>%summarise(EI_CI=mean(CI))
ip = preds %>% filter(!is.na(accelerating)) %>% group_by(model,hr)%>%summarise(EI = mean(accelerating))
ggplot(ip)+aes(x=hr,y=EI,color=as.factor(model))+geom_point() + labs(x = 'Hr',y=expression(EI)) + theme(legend.position="none") #+      guides(color = guide_legend(title = "Model"))
ggsave('../figures/timescaling/accel.png')

ggplot(cp)+aes(x=hr,y=EI_CI,color=as.factor(model))+geom_point()+ labs(x = 'Hr',y=expression('EI'['CI'])) +      guides(color = guide_legend(title = "Model"))
ggsave('../figures/timescaling/ci.png')
