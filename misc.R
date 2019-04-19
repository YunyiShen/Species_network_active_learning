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
