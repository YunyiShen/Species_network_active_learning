## generate a random graph
source("misc.R")
require(parallel)
require(optimParallel)

nspp = 3
set.seed(12345)
#linking = matrix(runif(nspp^2)<0.7,nspp,nspp)
#linking = linking * t(linking)
linking = matrix(1,nspp,nspp)
diag(linking) = 0
raster::plot(raster::raster(linking))
strength = matrix(runif(nspp^2,-1,1),nspp,nspp)
strength = .5*(strength+t(strength))

graph = strength * linking
raster::plot(raster::raster(graph))

## generate some random environment
nlat = 30 # nlat*nlat grid (not really important)

landscape1 = (rnorm(nlat*nlat,1,1))
landscape2 = (rnorm(nlat*nlat,0,1))
raster::plot(raster::raster(matrix(landscape2,nlat,nlat)))

set.seed(42)
betas = rbind(rnorm(nspp,0,1),rnorm(nspp,0,1),rnorm(nspp,0,1))
raster::plot(raster::raster(betas))
env = cbind(1,landscape1,landscape2)

## get some distribution
Z_sample = apply(env,1,SampleZ,betas,graph,1)
Z_sample = t(Z_sample)


ncore = 6
cl = makeCluster(getOption("cl.cores", ncore))
clusterExport(cl,c("env","Shan_ent","nspp","getGraph","logLik"))

## Numeric experiment
data_using = 60
#Iter_1 = sample(nlat^2,data_using)
Iter_1 = order(env[,2],decreasing=TRUE)[(1:data_using-1)* (nlat^2/data_using) + 1]
env_1 = env[Iter_1,]
Z_sample_1 = Z_sample[Iter_1,]
Sampled = matrix(0,1,nlat^2)
Sampled[,Iter_1]=1

theta = c(as.vector(betas),getGraphpar(graph,nspp))
theta_ini = theta + 0.05*runif(length(theta),-1,1)

MLE_1 = optimParallel(theta_ini,logLik,s=Z_sample_1,env=env_1,nspp=nspp,method = "L-BFGS-B",control = list(maxit=1000),parallel = list(cl=cl))
graphpar_est_1 = MLE_1$par[-(1:length(betas))]
graphpar_real = getGraphpar(graph,nspp)
L2_dif = matrix(1,1,nlat^2/data_using)
#L2_dif[1]=sqrt(sum((graphpar_est_1-graphpar_real)^2))
real_value = c( as.vector( betas), graphpar_real)
L2_dif[1]=sqrt(sum((MLE_1$par-real_value)^2))
L2_dif_random = L2_dif
Iter_1_random = Iter_1
Sampled_random = Sampled
raster::plot(raster::raster(getGraph(graphpar_est_1,nspp)))


for(i in 2:floor(nlat^2/(data_using))){
  clusterExport(cl,"MLE_1")
  cat("Making FImap",i, "...\n\n")
  #FImap = apply(env,1,TrFI,MLE_1$par,nspp)
  #FImap = apply(env,1,TrFI,real_value,nspp)#cheating
  #FImap = abs(FImap)
  FImap = parApply(cl,env,1,Shan_ent,MLE_1$par,nspp)
  filename = paste0("./figs/Iter_",i,"_SNmap.jpg")
  jpeg(filename)
  raster::plot(raster::raster(matrix(FImap,nlat,nlat)))
  dev.off()
  cat("Done \n\n")
  
  FI_unsurveyed = FImap[Sampled==0]
  env_unsurveyed = env[Sampled==0,]
  
  site_unsuried = (1:(nlat^2))[Sampled==0]
  
  dataperlevel = floor( sum(1-Sampled)/data_using)
  max_FI = apply(matrix(1:data_using),1,function(k,FIm,env2,dataperlevel){
    order_env = order(env2,decreasing=TRUE)
    env_level = env2 [order_env[1:dataperlevel + (k-1)*dataperlevel]]
    FI_level = FIm[order_env[1:dataperlevel + (k-1)*dataperlevel]]
    return(which(env2==env_level[FI_level==max(FI_level)]))
  },FI_unsurveyed,env_unsurveyed[,2],dataperlevel)
  
  Iter_2 = order(FImap*(1-Sampled),decreasing=TRUE)[1:data_using]
  Iter_1 = c(Iter_1,Iter_2)
  
  
  env_1 = env[c(Iter_1),]
  Z_sample_1 = Z_sample[Iter_1,]
  Sampled = matrix(0,1,nlat^2)
  Sampled[,Iter_1]=1
  filename = paste0("./figs/Iter_",i,"_sampled.jpg")
  jpeg(filename)
  raster::plot(raster::raster(matrix(Sampled,nlat,nlat)))
  dev.off()
  
  cat("random sample\n\n")
  
  env_unsurveyed_rand = env[Sampled_random==0,]
  rand_sample = apply(matrix(1:data_using),1,function(dummy,dataperlevl){sample(dataperlevel,1)},dataperlevel)
  rand_sample = rand_sample + (1:data_using-1)*dataperlevel
  
  site_unsuried_rand = (1:(nlat^2))[Sampled_random==0]
  #Iter_2_random = sample(which(Sampled_random==0),data_using)
  order_env_rand = order(env_unsurveyed_rand[,2],decreasing=TRUE)
  Iter_2_random = site_unsuried_rand[order_env_rand[rand_sample]]
  Iter_1_random = c(Iter_1_random,Iter_2_random)
  
  env_1_rand = env[c(Iter_1_random),]
  Z_sample_1_rand = Z_sample[Iter_1_random,]
  Sampled_random = matrix(0,1,nlat^2)
  Sampled_random[,Iter_1_random]=1
  filename = paste0("./figs/Iter_",i,"_sampled_random.jpg")
  jpeg(filename)
  raster::plot(raster::raster(matrix(Sampled_random,nlat,nlat)))
  dev.off()
  
  cat("MLE",i, "\n\n")
  MLE_1 = optimParallel(theta_ini,logLik,s=Z_sample_1,env=env_1,nspp=nspp,method = "L-BFGS-B",control = list(maxit=1000),parallel = list(cl=cl))
  graphpar_est_1 = MLE_1$par[-(1:length(betas))]
  L2_dif[i]=sqrt(sum((MLE_1$par-real_value)^2))
  filename = paste0("./figs/Iter_",i,"_graph.jpg")
  jpeg(filename)
  raster::plot(raster::raster(getGraph(graphpar_est_1,nspp)))
  dev.off()
  cat("Done\n\n")
  
  cat("Iter_",i,"done with L2_dif",L2_dif[i],"\n\n\n")
  
  
  
  cat("MLE of random sample",i, "\n\n")
  MLE_1_rand = optimParallel(theta_ini,logLik,s=Z_sample_1_rand,env=env_1_rand,nspp=nspp,method = "L-BFGS-B",control = list(maxit=1000),parallel = list(cl=cl))
  graphpar_est_1_rand = MLE_1_rand$par[-(1:length(betas))]
  L2_dif_random[i]=sqrt(sum((MLE_1_rand$par-real_value)^2))
  filename = paste0("./figs/Iter_",i,"_graph_random.jpg")
  jpeg(filename)
  raster::plot(raster::raster(getGraph(graphpar_est_1_rand,nspp)))
  dev.off()
  cat("Done\n\n")
  
  cat("Iter_",i,"done with L2_dif_rand",L2_dif_random[i],"\n\n\n")
  
}

stopCluster(cl)

require(ggplot2)
data_temp = data.frame(Sample_Size = 1:9 * 60,L2_difference = as.numeric( L2_dif[1:9]),method = "Active learning")
temp = data.frame(Sample_Size = 1:9 * 60,L2_difference = as.numeric( L2_dif_random[1:9]),method = "Environmental Gradient")
data_temp = rbind(data_temp,temp)
rm(temp)

ggplot(data = data_temp,aes(x=Sample_Size,y=L2_difference))+
  geom_point(aes(color = method))+
  geom_line(aes(color = method))

require(export)
graph2ppt(file="3spp_activelearning_2_biasenv_mu1_sd1_mu0_sd1.pptx")
