current_temperature = function(iter, s_curve_amplitude, s_curve_center, s_curve_width) {
  s_curve_amplitude * s_curve(iter, s_curve_center, s_curve_width)
}

s_curve = function(x, center, width) {
  1 / (1 + exp((x - center) / width))
} 
# SANN Modified from Todd. Do data selection based on FI

run_intermediate_annealing_process = function(env_full,already_sampled,n,theta,nspp,
                                              starting_iteration, number_of_iterations,
                                              s_curve_amplitude, s_curve_center, s_curve_width,nsample = 500,method = "MH") {
  
  candidate = (1:nrow(env_full))[-already_sampled]
  sample_curr = candidate[sample(length(candidate),n)]
  FI_curr = InvCR_group(env_full[c(already_sampled,sample_curr),],theta,nspp,nsample,method)
  still = (1:nrow(env_full))[-c(already_sampled,sample_curr)]
  best_sample = sample_curr
  best_FI = FI_curr
  for(i in 1:number_of_iterations) {
    iter = starting_iteration + i
    temp = current_temperature(iter, s_curve_amplitude, s_curve_center, s_curve_width)
    
    sample_prop = sample_curr
    sample_prop[sample(n,1)]=still[sample(length(still),1)]
    FI_prop = InvCR_group(env_full[c(already_sampled,sample_prop),],theta,nspp,nsample,method)
    if(i%%1000==0){cat(i,"\n\n")}
    if (temp > 0) {
      ratio = exp((FI_curr - FI_prop) / temp)
    } else {
      ratio = as.numeric(FI_prop < FI_curr)
    }
    
    if (runif(1) < ratio) {
      sample_curr = sample_prop
      FI_curr = FI_prop
      
      if (FI_curr < best_distance) {
        best_sample = sample_curr
        best_FI = FI_curr
      }
    }
  }
  
  return(list(sample=sample_curr, FI=FI_curr, best_sample=best_sample, best_FI=best_FI))
}
