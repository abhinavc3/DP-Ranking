

#Paralalising Code

list.of.packages <- c("foreach",
                      "doParallel",
                      "tidyverse",
                      "ggplot2",
                      "snow",
                      "Rmpi")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE, repos = 'http://cran.us.r-project.org')
}

#loading packages
for (package.i in list.of.packages) {
  suppressPackageStartupMessages(library(package.i,
                                         character.only = TRUE))
}

cores <- strtoi(Sys.getenv("NSLOTS"))-1
#create the cluster
my.cluster <- makeCluster(cores, methods = FALSE, type = "MPI")
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

## Don't forget to parallalise

### import helper functions
source("./Helper-Script.R")



simulator <- function(p_seq, epsilon_seq, n_seq, no_of_iter = 15) {
  empty_data_frame <-
    data.frame(matrix(
      nrow = length(p_seq) * length(epsilon_seq) * length(n_seq),
      ncol = no_of_iter
    ))
  
  empty_data_frame <-
    cbind(expand.grid(n_seq, p_seq, epsilon_seq), empty_data_frame)
  colnames(empty_data_frame)[1:3] <- c("n", "p", "epsilon")
  L_infty_error_priv <- empty_data_frame
  L_2_error_priv <- empty_data_frame
  hamming_error_est_priv <- empty_data_frame
  hamming_error_count_priv <- empty_data_frame
  
  writeLines(c(""), "log.txt")
  
  for (i in 1:length(n_seq)) {
    n <- n_seq[i]
    k <- round(n / 4)
    fixed_theta <-
      log(c(runif(n - k, min =  0.2 , max = 0.7), rep(1, k)))
    #fixed_theta <- exp(runif(n, min = 0.5, max = 1))
    #fixed_theta <- log(runif(n, min = 0.5, max = 1))
    #fixed_theta <- log(seq(0.5, 1, length.out = n))
    
    fixed_centered_theta <- fixed_theta - mean(fixed_theta)
    S_k = recover_top_k_set(fixed_centered_theta , k = k)
    
    for (ip in 1:length(p_seq)) {
      p <- p_seq[ip]
      #hamming_error <- numeric(no_of_iter)
      #hamming_error_priv <- numeric(no_of_iter)
      sim_list <-
        foreach(iter = 1:no_of_iter, .combine = "rbind", .export = ls(globalenv())) %dopar% {
          L_infty_error_priv_inside <- numeric(length(epsilon_seq))
          L_2_error_priv_inside <- numeric(length(epsilon_seq))
          hamming_error_est_priv_inside <- numeric(length(epsilon_seq))
          hamming_error_count_priv_inside <- numeric(length(epsilon_seq))
          
          data <- data_gen(n = n,
                           theta = fixed_centered_theta,
                           p = p)
          #hat_S_k <- top_k_borda_count(win_mat = data, k = k, lambda = 0)
          #hat_S_k_priv <- top_k_borda_count(win_mat = data, k = k)
          for (ieps in 1:length(epsilon_seq)) {
            epsilon <- epsilon_seq[ieps]
            if (epsilon > 500) {
              lambda_mle <- 0
              lambda_borda <- 0
            } else{
              lambda_mle <- 8 / epsilon
              lambda_borda = 2 / epsilon
            }
            optim_theta_priv <- optim(
              par = fixed_centered_theta,
              fn = obj_func,
              gr = grad_obj_func,
              win_mat = data,
              n = nrow(data),
              p = p,
              H = H_logistic,
              grad_H = grad_H_logistic,
              lambda = lambda_mle,
              method = "BFGS"
            )
            L_infty_error_priv_inside[ieps] <-
              log(
                max(abs(
                  optim_theta_priv$par  - fixed_centered_theta
                ) /
                  max(abs(
                    fixed_centered_theta
                  )))
              )
            L_2_error_priv_inside[ieps] <-
              log(
                sqrt(sum(
                  (optim_theta_priv$par -fixed_centered_theta)^2
                )) /
                sqrt(sum(
                    (fixed_centered_theta)^2
                  ))
                 )
            hat_S_k_est_priv = recover_top_k_set(optim_theta_priv$par, k = k)
            hamming_error_est_priv_inside[ieps] <-
              1 - length(intersect(S_k, hat_S_k_est_priv)) / k
            
            #Private borda counting
            
            hat_S_k_count_priv = top_k_borda_count(win_mat = data,
                                                   k = k,
                                                   lambda = lambda_borda)
            hamming_error_count_priv_inside[ieps] <-
              1 - length(intersect(S_k, hat_S_k_count_priv)) / k
            cat(
              paste(
                "\n",
                "Iteration ",
                iter,
                " for n, p, epsilon = ",
                n,
                p,
                epsilon
              ),
              file = "log.txt",
              append = TRUE
            )
          }
          return(
            list(
              L_infty_error_priv = L_infty_error_priv_inside,
              L_2_error_priv = L_2_error_priv_inside,
              hamming_error_est_priv = hamming_error_est_priv_inside,
              hamming_error_count_priv = hamming_error_count_priv_inside
            )
          )
        }
      sim_df <- as.data.frame(sim_list)
      
      L_infty_error_priv[L_infty_error_priv$n == n_seq[i] &
                           L_infty_error_priv$p == p_seq[ip], -(1:3)] <-
        matrix(unlist(sim_df$L_infty_error_priv), nrow = 2)
      
      L_2_error_priv[L_2_error_priv$n == n_seq[i] &
                           L_2_error_priv$p == p_seq[ip], -(1:3)] <-
        matrix(unlist(sim_df$L_2_error_priv), nrow = 2)
      
      hamming_error_est_priv[hamming_error_est_priv$n == n_seq[i] &
                               hamming_error_est_priv$p == p_seq[ip], -(1:3)] <-
        matrix(unlist(sim_df$hamming_error_est_priv), nrow = 2)
      
      hamming_error_count_priv[hamming_error_count_priv$n == n_seq[i] &
                                 hamming_error_count_priv$p == p_seq[ip], -(1:3)] <-
        matrix(unlist(sim_df$hamming_error_count_priv), nrow = 2)
      print(i)
    }
  }
  list(
    L_infty_error_priv = L_infty_error_priv,
    L_2_error_priv = L_2_error_priv,
    hamming_error_est_priv = hamming_error_est_priv,
    hamming_error_count_priv = hamming_error_count_priv
  )
}

######## Tester #########
simulator(p_seq = c(1), epsilon_seq = c(0.5), n_seq = c(100), no_of_iter = 2)


###### For Experiment 1 #######

no_of_iter <- 90

p_seq <- c( 1)
epsilon_seq <- c(0.5, 1, 2.5, 1000)
##n_seq <- c(40, 50)
n_seq <- seq(100, 700, 75)

simul_lists = simulator(p_seq = p_seq, epsilon_seq = epsilon_seq, n_seq = n_seq, no_of_iter = no_of_iter)

##############Saving


write.csv(simul_lists$hamming_error_est_priv, "./Simulation-Results/accuracy-vs-n/hamming-error-est-priv.csv" , row.names = F)

write.csv(simul_lists$hamming_error_count_priv, "./Simulation-Results/accuracy-vs-n/hamming-error-count-priv.csv", row.names = F)

write.csv(simul_lists$L_infty_error_priv, "./Simulation-Results/accuracy-vs-n/L-infty-error-priv.csv", row.names = F)

write.csv(simul_lists$ L_2_error_priv, "./Simulation-Results/accuracy-vs-n/L-2-error-priv.csv", row.names = F)

###### For Experiment 2 #######

no_of_iter <- 90

p_seq <- seq(0.1, 1, 0.075)
epsilon_seq <- c(0.5, 1, 2.5, 1000)
n_seq <- c(300)

simul_lists = simulator(p_seq = p_seq, epsilon_seq = epsilon_seq, n_seq = n_seq, no_of_iter = no_of_iter)

##############Saving


write.csv(simul_lists$hamming_error_est_priv, "./Simulation-Results/accuracy-vs-p/hamming-error-est-priv.csv" , row.names = F)

write.csv(simul_lists$hamming_error_count_priv, "./Simulation-Results/accuracy-vs-p/hamming-error-count-priv.csv", row.names = F)

write.csv(simul_lists$L_infty_error_priv, "./Simulation-Results/accuracy-vs-p/L-infty-error-priv.csv", row.names = F)

write.csv(simul_lists$L_2_error_priv, "./Simulation-Results/accuracy-vs-p/L-2-error-priv.csv", row.names = F)

###### For Experiment 3 #######

no_of_iter <- 90

p_seq <- c(0.25, 0.5, 0.75, 1)
epsilon_seq <- c(1/5, 1/4, 1/3, 1/2, 1, 3/2, 2, 5/2, 5, 1000)
n_seq <- c(300)

simul_lists = simulator(p_seq = p_seq, epsilon_seq = epsilon_seq, n_seq = n_seq, no_of_iter = no_of_iter)

##############Saving



write.csv(simul_lists$hamming_error_est_priv, "./Simulation-Results/accuracy-vs-epsilon/hamming-error-est-priv.csv" , row.names = F)

write.csv(simul_lists$hamming_error_count_priv, "./Simulation-Results/accuracy-vs-epsilon/hamming-error-count-priv.csv", row.names = F)

write.csv(simul_lists$L_infty_error_priv, "./Simulation-Results/accuracy-vs-epsilon/L-infty-error-priv.csv", row.names = F)

write.csv(simul_lists$L_2_error_priv, "./Simulation-Results/accuracy-vs-epsilon/L-2-error-priv.csv", row.names = F)
