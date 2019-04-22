current_temperature = function(iter, s_curve_amplitude, s_curve_center, s_curve_width) {
  s_curve_amplitude * s_curve(iter, s_curve_center, s_curve_width)
}

s_curve = function(x, center, width) {
  1 / (1 + exp((x - center) / width))
} 
# SANN Modified from Todd. Do data selection based on FI

run_intermediate_annealing_process = function(env_full,already_sampled,n,theta,nspp,
                                              starting_iteration, number_of_iterations,
                                              s_curve_amplitude, s_curve_center, s_curve_width,nsample = 500,method = "MH",ncore = 4) {
  source("misc.R")
  require(parallel)
  cl = makeCluster(getOption("cl.cores", ncore))
  clusterExport(cl,c("sufstat","getGraph","SampleZ"))
  candidate = (1:nrow(env_full))[-already_sampled]
  sample_curr = candidate[sample(length(candidate),n)]
  FI_curr = InvCR_group(env_full[c(already_sampled,sample_curr),],theta,nspp,nsample,method,cl)
  still = (1:nrow(env_full))[-c(already_sampled,sample_curr)]
  best_sample = sample_curr
  best_FI = FI_curr
  plot(1,1,type="n",xlim = c(0,number_of_iterations),ylim = c(0,1.5*best_FI))
  for(i in 1:number_of_iterations) {
    iter = starting_iteration + i
    temp = current_temperature(iter, s_curve_amplitude, s_curve_center, s_curve_width)
    
    sample_prop = sample_curr
    still_prop = still
    swap1 = sample(n,1)
    swap2 = sample(length(still),1)
    
    sample_prop[swap1]=still[swap2]
    still_prop[swap2]=sample_curr[swap1]
    FI_prop = InvCR_group(env_full[c(already_sampled,sample_prop),],theta,nspp,nsample,method,cl)
    
    if (temp > 0) {
      ratio = exp((FI_curr - FI_prop) / temp)
    } else {
      ratio = as.numeric(FI_prop < FI_curr)
    }
    FI_prev = FI_curr
    if (runif(1) < ratio) {
      sample_curr = sample_prop
      FI_curr = FI_prop
      still=still_prop
      if (FI_curr < best_FI) {
        best_sample = sample_curr
        best_FI = FI_curr
      }
    }
    #if(i%%1==0){cat(i,"\n FI_curr",FI_curr,"\n\n")}
    svMisc::progress(((i-1)/number_of_iterations)*100,progress.bar = T)
    lines(c(i-1,i),c(FI_prev,FI_curr))
  }
  
  return(list(sample=sample_curr, FI=FI_curr, best_sample=best_sample, best_FI=best_FI))
}
