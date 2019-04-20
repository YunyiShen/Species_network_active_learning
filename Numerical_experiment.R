## generate a random graph
source("misc.R")

nspp = 4
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
nlat = 30 # nlat*nlat grid (not really too important)

landscape = (rnorm(nlat*nlat))
raster::plot(raster::raster(matrix(landscape,nlat,nlat)))

set.seed(42)
betas = rbind(runif(nspp,-.3,.3),rnorm(nspp,0,3))
raster::plot(raster::raster(betas))
env = cbind(1,landscape)

## get some distribution
Z_sample = apply(env,1,SampleZ,betas,graph,1)
Z_sample = t(Z_sample)



## Numeric experiment
data_using = 50
#Iter_1 = sample(nlat^2,data_using)
Iter_1 = order(env[,2],decreasing=TRUE)[(1:data_using-1)* (nlat^2/data_using) + 1]
env_1 = env[Iter_1,]
Z_sample_1 = Z_sample[Iter_1,]
Sampled = matrix(0,1,nlat^2)
Sampled[,Iter_1]=1

theta = c(as.vector(betas),getGraphpar(graph,nspp))
theta_ini = theta + 0.05*runif(length(theta),-1,1)

MLE_1 = optim(theta_ini,logLik,s=Z_sample_1,env=env_1,nspp=nspp,method = "L-BFGS-B",control = list(maxit=1000))
graphpar_est_1 = MLE_1$par[-(1:length(betas))]
graphpar_real = getGraphpar(graph,nspp)
L2_dif = matrix(1,1,nlat^2/data_using)
L2_dif[1]=sqrt(sum((graphpar_est_1-graphpar_real)^2))
L2_dif_random = L2_dif
Iter_1_random = Iter_1
Sampled_random = Sampled
raster::plot(raster::raster(getGraph(graphpar_est_1,nspp)))


for(i in 2:floor(nlat^2/(data_using*2))){
  cat("Making FImap",i, "...\n\n")
  FImap = apply(env,1,TrFI,MLE_1$par,nspp)
  FImap = abs(FImap)
  filename = paste0("./figs/Iter_",i,"_FImap.jpg")
  jpeg(filename)
  raster::plot(raster::raster(matrix(FImap,nlat,nlat)))
  dev.off()
  cat("Done \n\n")
  
  FI_unsurveied = FImap[Sampled==0]
  env_unsurveied = env[Sampled==0,]
  
  site_unsuried = (1:(nlat^2))[Sampled==0]
  
  dataperlevel = floor( sum(1-Sampled)/data_using)
  max_FI = apply(matrix(1:data_using),1,function(k,FIm,env2,dataperlevel){
    order_env = order(env2,decreasing=TRUE)
    env_level = env2 [order_env[1:dataperlevel + (k-1)*dataperlevel]]
    FI_level = FIm[order_env[1:dataperlevel + (k-1)*dataperlevel]]
    return(which(env2==env_level[FI_level==max(FI_level)]))
  },FI_unsurveied,env_unsurveied[,2],dataperlevel)
  
  data_using_for_graph = ceiling(data_using/2)
  data_using_for_field = floor(data_using/2)
  
  
  
  Iter_2_graph = order(FImap*(1-Sampled),decreasing=TRUE)[1:data_using_for_graph]
  #Iter_2 = site_unsuried[max_FI]
  Iter_1_temp = c(Iter_1,Iter_2_graph)
  Iter_2_field = sample( (1:nlat^2)[-Iter_1_temp],data_using_for_field )
  Iter_1 = c(Iter_1_temp,Iter_2_field)
  
  
  env_1 = env[c(Iter_1),]
  Z_sample_1 = Z_sample[Iter_1,]
  Sampled = matrix(0,1,nlat^2)
  Sampled[,Iter_1]=1
  filename = paste0("./figs/Iter_",i,"_sampled.jpg")
  jpeg(filename)
  raster::plot(raster::raster(matrix(Sampled,nlat,nlat)))
  dev.off()
  
  cat("random sample\n\n")
  
  rand_sample = apply(matrix(1:data_using),1,function(dummy,dataperlevl){sample(dataperlevel,1)},dataperlevel)
  rand_sample = rand_sample + (1:data_using-1)*dataperlevel
  
  #Iter_2_random = sample(which(Sampled_random==0),data_using)
  order_env = order(env_unsurveied[,2],decreasing=TRUE)
  Iter_2_random = site_unsuried[order_env[rand_sample]]
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
  MLE_1 = optim(theta_ini,logLik,s=Z_sample_1,env=env_1,nspp=nspp,method = "L-BFGS-B",control = list(maxit=1000))
  graphpar_est_1 = MLE_1$par[-(1:length(betas))]
  L2_dif[i]=sqrt(sum((graphpar_est_1-graphpar_real)^2))
  filename = paste0("./figs/Iter_",i,"_graph.jpg")
  jpeg(filename)
  raster::plot(raster::raster(getGraph(graphpar_est_1,nspp)))
  dev.off()
  cat("Done\n\n")
  
  cat("Iter_",i,"done with L2_dif",L2_dif[i],"\n\n\n")
  
  
  
  cat("MLE of random sample",i, "\n\n")
  MLE_1_rand = optim(theta_ini,logLik,s=Z_sample_1_rand,env=env_1_rand,nspp=nspp,method = "L-BFGS-B",control = list(maxit=1000))
  graphpar_est_1_rand = MLE_1_rand$par[-(1:length(betas))]
  L2_dif_random[i]=sqrt(sum((graphpar_est_1_rand-graphpar_real)^2))
  filename = paste0("./figs/Iter_",i,"_graph_random.jpg")
  jpeg(filename)
  raster::plot(raster::raster(getGraph(graphpar_est_1_rand,nspp)))
  dev.off()
  cat("Done\n\n")
  
  cat("Iter_",i,"done with L2_dif_rand",L2_dif_random[i],"\n\n\n")
  
}

require(ggplot2)
data_temp = data.frame(Sample_Size = 1:10 * 45,L2_difference = L2_dif[1:10],method = "Active learning")
temp = data.frame(Sample_Size = 1:10 * 45,L2_difference = L2_dif_random[1:10],method = "Random sampling")
data_temp = rbind(data_temp,temp)
rm(temp)

ggplot(data = data_temp,aes(x=Sample_Size,y=L2_difference))+
  geom_point(aes(color = method))+
  geom_line(aes(color = method))

require(export)
graph2ppt(file="5spp_activelearning_45perstep.pptx")
