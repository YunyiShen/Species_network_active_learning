SampleZ = function(env,beta,graph,n){
  nspp = ncol(beta)
  require(IsingSampler)
  thr = env%*%beta
  res = IsingSampler(thresholds = thr,graph = graph,responses = c(-1L,1L),n=n,method = "CFTP")
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

TrFI = function(env,theta,nspp,nsample = 500){
  ncov = length(env)
  betas = matrix(theta[1:(ncov*nspp)],ncov,nspp)
  graphpar = theta[(ncov*nspp+1):length(theta)]
  graph = getGraph(graphpar,nspp)
  Z_samples = SampleZ(env,betas,graph,nsample)
  Tr = matrix(0,nrow = nsample,ncol = length(theta))
  k=1
  for(i in 1:ncov){
    for(j in 1:nspp){
      Tr[,k] = env[i]*Z_samples[,j]
      k = k + 1 
    }
  }
  
  for(i in 2:nspp-1){
    for(j in (i+1):nspp)
    Tr[,k] = (Z_sample[,i]*Z_sample[,j])
    k = k + 1
  }
  #cat(Tr,"\n") # for test
  return(sum(eigen(Tr,T,T)))
}