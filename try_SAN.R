source("SAN.R")
SANN = run_intermediate_annealing_process(env,already_sampled=Iter_1,n=50,theta=MLE_1$par,nspp=4,
                                          starting_iteration=1, number_of_iterations=5000,
                                          s_curve_amplitude=4000, s_curve_center=0, 
                                          s_curve_width=3000,nsample = 200,method = "MH",ncore = 4) 
