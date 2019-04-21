source("SAN.R")
SANN = run_intermediate_annealing_process(env,already_sampled=Iter_1,n=100,theta,nspp=4,
                                          starting_iteration=1, number_of_iterations=15000,
                                          s_curve_amplitude=4000, s_curve_center=0, 
                                          s_curve_width=3000,nsample = 500,method = "MH") 
