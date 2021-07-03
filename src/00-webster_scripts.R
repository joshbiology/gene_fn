
#Designed for writing shell scripts that operate Webster runs

#The core algorithm behind Webster (DGRDL) is a Matlab script. I wrote a wrapper
#for that script to accept parameters used in Webster.

#There are machine specific dependencies to think about in the following, including:
#1. Where is your matlab path?
#2. Where is your install path for the Webster matlab dependencies?
#3. Does matlab have the correct folder access from the location in which the batch script is run?
#4. Does it have access to your data?

library(ProjectTemplate); load.project()

# Define environments ---------------------------------------------------

local <- list(matlab_path = "/Applications/MATLAB_R2019b.app/bin/matlab",
              install_path = "")

cluster <- list(matlab_path = "matlab",
                install_path = "")

#define your own environment here!

script_out_path <- file.path(".", "cluster_scripts")

# Define parameters -------------------------------------------------------


generate_webster_params <- function(input_data_path,
                                    out_name,
                                    out_path = script_out_path,
                                    K_param = 10,
                                    T_param = 2,
                                    alpha = 0.2,
                                    beta = 0.6,
                                    iternum = 20,
                                    seed = 1,
                                    num_neighbor_gene = 5,
                                    num_neighbor_cl = 5,
                                    knn_weight_mode = 'Cosine') {
  
  
  #The following inputs are fixed strings.
  #input_path
  #out_path
  #knn_weight_mode
  
  #The following inputs are numeric scalars or vectors.
  #K
  #T
  #alpha
  #beta
  #iternum
  #seed
  #num_neighbor_gene
  #num_neighbor_cl

  param_grid <- expand.grid(K_param = K_param,
                            T_param = T_param,
                            alpha = alpha,
                            beta = beta,
                            iternum = iternum,
                            seed = seed,
                            num_neighbor_gene = num_neighbor_gene,
                            num_neighbor_cl = num_neighbor_cl)
  

  out_strings <- map(1:nrow(param_grid), function(x) {
    tmp <- slice(param_grid,x)
    param_string <- paste(names(tmp),tmp,sep = "=",collapse = ";")
    filename_string <- paste(Sys.Date(),"-",out_name, ";", param_string, ".mat", sep = "")
    file.path(out_path, filename_string)
  })
  param_grid %>% 
    mutate(knn_weight_mode = knn_weight_mode,
           input_data_path = input_data_path,
           out =out_strings)
  
}



# Scripting function ------------------------------------------------------

webster_batch_script <- function(params, environment, script_name) {
  
  matlab_commands <- map_chr(1:nrow(params), function(index) {
    x <- params %>% slice(index)
    
    sprintf("%s -batch \"%sDGRDL_wrapper('%s', 'K', %d, 'T', %d, 'alpha', %f, 'beta', %f, 'iternum', %d, 'seed', %d, 'num_neighbor_gene', %d, 'num_neighbor_cl', %d, 'knn_weight_mode', '%s', 'filename_out', '%s');exit\"", 
                                                              environment$matlab_path,
                                                              environment$install_path,
                                                              x$input_data_path, 
                                                              x$K_param, 
                                                              x$T_param,
                                                              x$alpha,
                                                              x$beta,
                                                              x$iternum,
                                                              x$seed,
                                                              x$num_neighbor_gene,
                                                              x$num_neighbor_cl,
                                                              x$knn_weight_mode,
                                                              x$out)
  })
  
  cat(matlab_commands, sep = '\n',file=file.path(script_out_path, paste(script_name,"txt", sep = ".")))
}



# Export scripts -----------------------------------------------------------

#Genotoxic data
generate_webster_params(input_data_path = "./output/02-genotoxic/durocher_matlab.csv", out_name = "genotoxic", K_param = 1:31, T_param = 1:2) %>% 
  webster_batch_script(local,"genotoxic_batch")

#DepMap grid sweep
generate_webster_params(input_data_path = "/broad/hahnlab2/joshpan/cache/03-depmap_preprocess/avana_19q4_webster.csv", 
                        out_name = "depmap_grid", 
                        out_path = "/broad/hahnlab2/joshpan/output",
                        K_param = seq(25,675, by = 25), 
                        T_param = 1:8) %>% 
  webster_batch_script(cluster,"depmap_grid_batch")


#DepMap wide sweep
generate_webster_params(input_data_path = "/broad/hahnlab2/joshpan/cache/03-depmap_preprocess/avana_19q4_webster.csv", 
                        out_name = "depmap_wide", 
                        out_path = "/broad/hahnlab2/joshpan/output",
                        K_param = seq(25,675, by = 5), 
                        T_param = 4) %>% 
  webster_batch_script(cluster,"depmap_wide_batch")

#DepMap deep sweep
generate_webster_params(input_data_path = "/broad/hahnlab2/joshpan/cache/03-depmap_preprocess/avana_19q4_webster.csv",  
                        out_name = "depmap_deep", 
                        out_path = "/broad/hahnlab2/joshpan/output",
                        K_param = seq(150,350, by = 5), 
                        T_param = 4,
                        seed = 1:5) %>% 
  webster_batch_script(cluster,"depmap_deep_batch")


#DepMap trained on noisy data.
list("0", "25", "50", "100", "150") %>% 
  map_dfr(~generate_webster_params(input_data_path = sprintf("./output/04-depmap_postprocess/avana_19q4_webster_train_noise_%s.csv", .), 
                                   out_name = sprintf("depmap_denoise_%s", .), out_path = "./data/interim/matlab/depmap_denoise", K_param = 220, T_param = 4, seed = 1:5)) %>% 
  webster_batch_script(local,"denoise")



