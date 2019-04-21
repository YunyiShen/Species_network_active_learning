SampleZ = function(env,beta,graph,n,method = "CFTP"){
  nspp = ncol(beta)
  require(IsingSampler)
  thr = env%*%beta
  res = IsingSampler(thresholds = thr,graph = graph,responses = c(-1L,1L),n=n,method = method)
}

getGraph = function(graphPar,nspp){
  graph = matrix(0,nspp,nspp)
  for(i in 2:nspp-1){
    graph[i,(i+1):nspp]=graphPar[1:(nspp-i)]
    graphPar = graphPar[-c(1:(nspp-i))]
  }
  graph = graph + t(graph)
  return(graph)
}

getGraphpar = function(graph,nspp){
  graphPar = numeric()
  for(i in 2:nspp-1){
    graphPar=c(graphPar,graph[i,(i+1):nspp])
  }
  return(graphPar)
}

logLik = function(theta,s,env,nspp){
  require(IsingSampler)
  ncov = ncol(env)
  nsite = nrow(env)
  beta = matrix(theta[1:(ncov*nspp)],ncov,nspp)
  thr = env%*%beta
  graphpar = theta[(ncov*nspp+1):length(theta)]
  graph = getGraph(graphpar,nspp)
  
  logliki = sum(log(apply(matrix(1:nsite),1,function(k,s,thr,graph){
    IsingStateProb(s[k,],graph,thr[k,],1,c(-1L,1L))
  },s,thr,graph)))
  return(-logliki)
}

TrFI = function(env,theta,nspp,nsample = 5000,method = "MH"){
  ncov = length(env)
  betas = matrix(theta[1:(ncov*nspp)],ncov,nspp)
  graphpar = theta[(ncov*nspp+1):length(theta)]
  graph = getGraph(graphpar,nspp)
  Z_samples = SampleZ(env,betas,graph,nsample,method = method)
  Tr = matrix(0,nrow = nsample,ncol = .5*(nspp-1)*nspp)
  k=1
  
  for(i in 2:nspp-1){
    for(j in (i+1):nspp){
    Tr[,k] = (Z_samples[,i]*Z_samples[,j])
    k = k + 1
    }
  }
  Tr = cbind(Tr,Z_samples)
  FI = cov(Tr)
  #cat(Tr,"\n") # for test
  return(1/sum(1/eigen(FI,T,T)$value))
}

sufstat = function(env,theta,nspp,nsample = 500,method = "MH"){
  ncov = length(env)
  betas = matrix(theta[1:(ncov*nspp)],ncov,nspp)
  graphpar = theta[(ncov*nspp+1):length(theta)]
  graph = getGraph(graphpar,nspp)
  Z_samples = SampleZ(env,betas,graph,nsample,method = method)
  Tr = matrix(0,nrow = nsample,ncol = .5*(nspp-1)*nspp)
  k=1
  suf_beta = matrix(0,nrow = nsample,ncol = ncov*nspp)
  for(i in 2:nspp-1){
    for(j in (i+1):nspp){
      Tr[,k] = (Z_samples[,i]*Z_samples[,j])
      k = k + 1
    }
    suf_beta[,1:2 + (i-1)*2]=Z_samples[,i]%*%t(matrix(env))
  }
  i = nspp
  suf_beta[,1:2 + (i-1)*2]=Z_samples[,i]%*%t(matrix(env))
  #suf_beta = apply(Z_samples,2,function(w,k){w%*%k},t(as.matrix(env)))
  Tr = cbind(Tr,(suf_beta))
  return(Tr)
}

InvCR_group = function(env,theta,nspp,nsample = 500,method = "MH"){
  env_list = split(env,row(env)) # make it a list of sites to loop with
  suf_stat_list = lapply(env_list,sufstat,theta,nspp,nsample,method)
  
  suf_stat = Reduce('+',suf_stat_list) # sum the sufficient statistics
  FI = cov(suf_stat)
  return(1/sum(1/eigen(FI,T,T)$value))
}

