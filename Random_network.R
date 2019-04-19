## generate a random graph
source("misc.R")
nspp = 4
set.seed(1996)
linking = matrix(runif(nspp^2)<0.7,nspp,nspp)
linking = linking * t(linking)
diag(linking) = 0
raster::plot(raster::raster(linking))
strength = matrix(runif(nspp^2,-.5,.5),nspp,nspp)
strength = .5*(strength+t(strength))

graph = strength * linking
raster::plot(raster::raster(graph))

## generate some random environment
nlat = 30 # nlat*nlat grid (not really too important)

landscape = (runif(nlat*nlat))
raster::plot(raster::raster(matrix(landscape,nlat,nlat)))

betas = matrix(runif(2*nspp,-1,1),2,nspp)
raster::plot(raster::raster(betas))
env = cbind(1,landscape)

## get some distribution
Z_sample = apply(env,1,SampleZ,betas,graph,1)
Z_sample = t(Z_sample)


## Alright, try full data
theta = c(as.vector(betas),getGraphpar(graph,nspp))
theta_ini = theta + 0.05*runif(length(theta),-1,1)
MLE_full = optim(theta_ini,logLik,s=Z_sample,env=env,nspp=nspp,method = "L-BFGS-B",control = list(maxit=1000))
graphpar_est = MLE_full$par[-(1:length(betas))]
relatgraph = graph
diag(relatgraph)=1
raster::plot(raster::raster(
  (sqrt(getGraph(graphpar_est,nspp)-graph)^2))

)
raster::plot(raster::raster(getGraph(graphpar_est,nspp)))
## start active leaning
### first round
data_using = 90
Iter_1 = sample(nlat^2,data_using)
env_1 = env[Iter_1,]
Z_sample_1 = Z_sample[Iter_1,]

MLE_1 = optim(theta_ini,logLik,s=Z_sample_1,env=env_1,nspp=nspp,method = "L-BFGS-B",control = list(maxit=1000))
graphpar_est_1 = MLE_1$par[-(1:length(betas))]
raster::plot(raster::raster(getGraph(graphpar_est_1,nspp)))
