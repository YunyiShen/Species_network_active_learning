## generate a random graph
source("misc.R")
nspp = 4
set.seed(1996)
linking = matrix(runif(nspp^2)<0.7,nspp,nspp)
linking = linking * t(linking)
diag(linking) = 0
raster::plot(raster::raster(linking))
strength = matrix(runif(nspp^2,-1,1),nspp,nspp)
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



## Numeric experiment
data_using = 90
Iter_1 = sample(nlat^2,data_using)
env_1 = env[Iter_1,]
Z_sample_1 = Z_sample[Iter_1,]
Sampled = matrix(0,1,nlat^2)
Sampled[,Iter_1]=1

MLE_1 = optim(theta_ini,logLik,s=Z_sample_1,env=env_1,nspp=nspp,method = "L-BFGS-B",control = list(maxit=1000))
graphpar_est_1 = MLE_1$par[-(1:length(betas))]
graphpar_real = getGraphpar(graph,nspp)
L2_dif = matrix(1,1,nlat^2/data_using)
L2_dif[1]=sqrt(sum((graphpar_est_1-graphpar_real)^2))
#raster::plot(raster::raster(getGraph(graphpar_est_1,nspp)))

for(i in 2:(nlat^2/data_using)){
  FImap = apply(env,1,TrFI,MLE_1$par,nspp)
  filename = paste0("./figs/Iter_",i,"_FImap.jpg")
  jpeg(filename)
  raster::plot(raster::raster(matrix(FImap,nlat,nlat)))
  dev.off()
  
  FI_unsurveied = FImap * (1-Sampled)
  Iter_2 = order(FI_unsurveied,decreasing=TRUE)[1:data_using]
  Iter_1 = c(Iter_1,Iter_2)
  
  env_1 = env[c(Iter_1),]
  Z_sample_1 = Z_sample[Iter_1,]
  Sampled = matrix(0,1,nlat^2)
  Sampled[,Iter_1]=1
  filename = paste0("./figs/Iter_",i,"_sampled.jpg")
  jpeg(filename)
  raster::plot(raster::raster(matrix(Sampled,nlat,nlat)))
  dev.off()
  
  MLE_1 = optim(theta_ini,logLik,s=Z_sample_1,env=env_1,nspp=nspp,method = "L-BFGS-B",control = list(maxit=1000))
  graphpar_est_1 = MLE_1$par[-(1:length(betas))]
  L2_dif[i]=sqrt(sum((graphpar_est_1-graphpar_real)^2))
  filename = paste0("./figs/Iter_",i,"_graph.jpg")
  jpeg(filename)
  raster::plot(raster::raster(getGraph(graphpar_est_2,nspp)))
  dev.off()
  
  cat("Iter_",i,"done with L2_dif",L2_dif[i],"\n")
  
}