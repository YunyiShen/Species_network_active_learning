source("SAN.R")
SANN = run_intermediate_annealing_process(env,already_sampled=Iter_1,n=50,theta=MLE_1$par,nspp=4,
                                          starting_iteration=100, number_of_iterations=5000,
                                          s_curve_amplitude=500, s_curve_center=0, 
                                          s_curve_width=3000,nsample = 300,method = "MH",ncore = 2) 
