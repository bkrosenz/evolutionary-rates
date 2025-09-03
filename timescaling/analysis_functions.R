# TODO: fix this, since the replicates are actually different samples
# var_name = switch (model,
#                    'EB' = 'a',
#                    'rate_trend'='slope',
#                    'delta'='delta',
# )
#
# N=length(results)
# params =vector(length = N)
# sigs=vector(length = N)
# bounds = matrix(data=NA, nrow=N, ncol=2)
#
# for (tree_no in 1:N){
#   r=results[[tree_no]]
#   aicc = sapply(r, function(x){x['aicc'][[1]]})
#   aicc_min = which.min(aicc)
#
#   h= r[[aicc_min]]
#
#   v = h[var_name][[1]]
#   s=h$sigsq
#
#   b = if (all(is.na(h$CI))){
#     c(NA,NA)
#   } else {
#     unlist(h$CI[,var_name])
#   }
#
#   sigs[tree_no]=s
#   params[tree_no]=v
#   bounds[tree_no,]=b
#
#   # rdf=as.data.frame(
#   #   do.call(cbind,
#   #           lapply(r, function(x){x[c(var_name,'sigsq')]}))
#   #   )
# }
# list(sigs=sigs,params=params,bounds=bounds)
# }


get_all_ml_estimates = function(results, model) {
  # TODO: fix this, since the replicates are actually different samples
  var_name = switch (
    model,
    'EB' = 'a',
    'rate_trend' = 'slope',
    'delta' = 'delta',
    'BM' = NA,
  )
  get_bound = function(h, ix) {
    # print(h)
    if (all(is.na(h$CI))) {
      NA
    } else {
      # print(h$CI)
      # if (){NA}else{
      (h$CI[ix, var_name])
    }
    # }
  }
  
  N = length(results)
  params = matrix(data = NA,
                  nrow = N,
                  ncol = replicates)
  sigs = matrix(data = NA,
                nrow = N,
                ncol = replicates)
  lbounds = matrix(data = NA,
                   nrow = N,
                   ncol = replicates)
  ubounds = matrix(data = NA,
                   nrow = N,
                   ncol = replicates)
  
  for (tree_no in 1:N) {
    r = results[[tree_no]]
    r = r[sapply(r, typeof) != 'character'] # some Rdata files are corrupted
    M = length(r)
    #
    # if (typeof(r)=="character"){
    #   print('error in loading data',r)
    #   next
    # }
    sigs[tree_no, 1:M] = sapply(r, function(x) {
      x$sigsq[[1]]
    })
    if (is.na(var_name)) {
      params[tree_no, 1:M] =     lbounds[tree_no, 1:M] =     ubounds[tree_no, 1:M] =
        NA
    } else{
      params[tree_no, 1:M] = sapply(r, function(x) {
        x[var_name][[1]]
      })
      lbounds[tree_no, 1:M] = sapply(r, get_bound, ix = 'lb')
      ubounds[tree_no, 1:M] = sapply(r, get_bound, ix = 'ub')
    }
  }
  list(
    sigs = sigs,
    params = params,
    lbounds = lbounds,
    ubounds = ubounds
  )
}

# get_accel=function(ests,model){
#   thresh  = switch (
#     model,
#     'EB' = 0,
#     'rate_trend'=0,
#     'delta'=1,
#   )
#   value = mean(ests$params)
#   acc=mean(ests$params>thresh )
#   x=ests$bounds[,1]
#   ci=mean(x[!is.na(x)] > thresh)
#   list(accelerating=acc,CI=ci, mean_value=value)
# }



get_accel = function(ests, model, reduce = 'grand') {
  if (model == 'BM') {
    return (
      list(
        accelerating = NA,
        CI = NA ,
        value = NA
      ))}
  
  thresh = switch (model,
                   'EB' = 0,
                   'rate_trend' = 0,
                   'delta' = 1,)
  value = ests$params
  acc = ests$params > thresh
  # x=ests$bounds[,1]
  x = ests$lbounds
  ci = x > thresh
  if (reduce == 'grand') {
    return(list(
      accelerating = mean(acc),
      CI = mean(ci, na.rm = T),
      mean_value = mean(value),
      nas = mean(is.na(ci))
    ))
  } else if (reduce == 'mean') {
    return (list(
      accelerating = rowMeans(acc, na.rm = T),
      CI = rowMeans(ci, na.rm = T),
      value = rowMeans(value, na.rm = T)
    ))
  } else {
    return (list(
      accelerating = acc,
      CI = ci ,
      value = value
    ))
  }
}

get_full_results = function(globstr, ngenes) {
  if (ngenes==0){
    globstr = paste(globstr,'no_ils',sep='/')
    ngenes=1
  }
  print(globstr)
  predictions = data.frame()
  for (model in c('EB', 'rate_trend', 'delta')) {
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
      if (file.exists(filename)) {
        print(filename)
        load(filename)
        # print(results[[1]][[1]]$CI)
        e = get_all_ml_estimates(results, model)
        c(lam, mu, h) %<-% sapply(tail(strsplit(gsub(
          '/no_ils', '', d
        ) , '/')[[1]], 3),
        function(x)
          as.numeric(substr(x, 2, 100)))
        a = get_accel(e, model, reduce = 'NA') # if as_mean==F, will be nested list
        # a = get_accel(e, model, reduce = 'mean') # if as_mean==F, will be nested list
        a$sig = e$sigs
        a = as.data.frame(sapply(a, as.vector))
        a$lam=lam
        a$h=h
        a$mu=mu
        a$r =lam-mu
        a$rho=mu/lam
        a$hr=h*(lam-mu)
        a$hrho=h*mu/lam
        a$model=model
        # a = list.append(
        #   a,
        #   sig = e$sigs,
        #   lam = lam,
        #   mu = mu,
        #   h = h,
        #   model = model
        # )
        predictions = bind_rows(predictions, a)
      }
    }
  }
  return(predictions)
}

# TODO:
analyze_results = function(globstr, ngenes) {
  predictions = list()
  for (model in c('EB', 'rate_trend', 'delta')) {
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
      if (file.exists(filename)) {
        load(filename)
        # print(results[[1]][[1]]$CI)
        e = get_all_ml_estimates(results, model)
        c(lam, mu, h) %<-% sapply(tail(strsplit(gsub(
          '/no_ils', '', d
        ) , '/')[[1]], 3),
        function(x)
          as.numeric(substr(x, 2, 100)))
        a = get_accel(e, model, reduce = 'grand') # if as_mean==F, will be nested list
        a = list.append(
          a,
          lam = lam,
          mu = mu,
          h = h,
          model = model
        )
        
        predictions = list.append(predictions, a)
      }
    }
  }
  rdf = bind_rows(predictions)
  rdf$r = rdf$lam - rdf$mu
  rdf$rho = rdf$mu / rdf$lam
  
  return(rdf)
}



