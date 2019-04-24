source("SAN.R")
set.seed(42)
SANN = run_intermediate_annealing_process(env,already_sampled=Iter_1,n=50,theta=MLE_1$par,nspp=4,nswap = 3,
                                          starting_iteration=1, number_of_iterations=15000,
                                          s_curve_amplitude=4, s_curve_center=0, 
                                          s_curve_width=2000,nsample = 300,method = "MH",ncore = 8) 
