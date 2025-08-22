#Paralalising Code

list.of.packages <- c("foreach",
                      "doParallel",
                      "ranger",
                      "palmerpenguins",
                      "tidyverse",
                      "kableExtra")

#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#if(length(new.packages) > 0){
#  install.packages(new.packages, dep=TRUE)
#}

#loading packages
#for (package.i in list.of.packages) {
 # suppressPackageStartupMessages(library(package.i,
#                                         character.only = TRUE))
#}
parallel::detectCores()
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(n.cores,
                                    type = "PSOCK")
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

#########################################################
######################### Importing Helper functions
#########################################################
source("./Helper-Script.R")

# Experiment

#n_seq <- seq(300, 450, 500)
epsilon <- 2.5
n <- 300
k <- round(n / 4)

param_seq <- seq(0.1, 1, 0.075)
#param_seq = c(0.5,1)
no_of_iter <- 15

#initializing theta and top k set
fixed_theta <- log(c(runif(n - k, min =  0.2 , max = 0.7), rep(1, k)))
fixed_centered_theta <- fixed_theta - mean(fixed_theta)
S_k = recover_top_k_set(fixed_centered_theta , k = k)


L_infty_error <- matrix(nrow = length(param_seq), ncol = no_of_iter)
L_infty_error_priv <-
  matrix(nrow = length(param_seq), ncol = no_of_iter)

hamming_error_est <- matrix(nrow = length(param_seq), ncol = no_of_iter)
hamming_error_est_priv <-
  matrix(nrow = length(param_seq), ncol = no_of_iter)

hamming_error_count <-
  matrix(nrow = length(param_seq), ncol = no_of_iter)
hamming_error_count_priv <-
  matrix(nrow = length(param_seq), ncol = no_of_iter)

writeLines(c(""), "log.txt")
for (i in 1:length(param_seq)) {
  p <- param_seq[i]
  
  #fixed_theta <- exp(runif(n, min = 0.5, max = 1))
  #fixed_theta <- log(runif(n, min = 0.5, max = 1))
  #fixed_theta <- log(seq(0.5, 1, length.out = n))
  
  

  
  #hamming_error <- numeric(no_of_iter)
  #hamming_error_priv <- numeric(no_of_iter)
  sim_list <-
    foreach(iter = 1:no_of_iter, .combine = "rbind") %dopar% {
      data <- data_gen(n = n, theta = fixed_centered_theta, p = p)
      #hat_S_k <- top_k_borda_count(win_mat = data, k = k, lambda = 0)
      #hat_S_k_priv <- top_k_borda_count(win_mat = data, k = k)
      optim_theta <- optim(
        par = fixed_centered_theta,
        fn = obj_func,
        gr = grad_obj_func,
        win_mat = data,
        n = nrow(data),
        p = p,
        H = H_logistic,
        grad_H = grad_H_logistic,
        lambda = 0,
        method = "BFGS"
      )
      L_infty_error <-
        max(abs( exp(optim_theta$par)  - exp(fixed_centered_theta)) /
              max(abs( exp(fixed_centered_theta) )))
      #S_k = recover_top_k_set(fixed_centered_theta, k = k)
      hat_S_k_est = recover_top_k_set(optim_theta$par, k = k)
      hamming_error_est <-
        1 - length(intersect(S_k, hat_S_k_est)) / k
      optim_theta_priv <- optim(
        par = fixed_centered_theta,
        fn = obj_func,
        gr = grad_obj_func,
        win_mat = data,
        n = nrow(data),
        p = p,
        H = H_logistic,
        grad_H = grad_H_logistic,
        eps = epsilon,
        method = "BFGS"
      )
      L_infty_error_priv <-
        max(abs( exp(optim_theta_priv$par)  - exp(fixed_centered_theta)) /
              max( abs(exp(fixed_centered_theta))) )
      hat_S_k_est_priv = recover_top_k_set(optim_theta_priv$par, k = k)
      hamming_error_est_priv <-
        1 - length(intersect(S_k, hat_S_k_est_priv)) / k
      
      #Private borda counting
      hat_S_k_count = top_k_borda_count(win_mat = data,
                                        k = k,
                                        lambda = 0)
      hat_S_k_count_priv = top_k_borda_count(win_mat = data,
                                             k = k,
                                             eps = epsilon)
      hamming_error_count <-
        1 - length(intersect(S_k, hat_S_k_count)) / k
      hamming_error_count_priv <-
        1 - length(intersect(S_k, hat_S_k_count_priv)) / k
      cat(paste("\n", "Iteration ", iter, " for p = ", p),
          file = "log.txt",
          append = TRUE)
      return(
        list(
          L_infty_error = L_infty_error,
          L_infty_error_priv = L_infty_error_priv,
          hamming_error_est = hamming_error_est,
          hamming_error_est_priv = hamming_error_est_priv,
          hamming_error_count = hamming_error_count,
          hamming_error_count_priv = hamming_error_count_priv
        )
      )
    }
  sim_df <- as.data.frame(sim_list)
  L_infty_error[i,] <- unlist(sim_df$L_infty_error)
  L_infty_error_priv[i,] <- unlist(sim_df$L_infty_error_priv)
  hamming_error_est[i,] <- unlist(sim_df$hamming_error_est)
  hamming_error_est_priv[i,] <-
    unlist(sim_df$hamming_error_est_priv)
  hamming_error_count[i,] <- unlist(sim_df$hamming_error_count)
  hamming_error_count_priv[i,] <-
    unlist(sim_df$hamming_error_count_priv)
  print(i)
  
}

############# Saving Simulation Results #########

setwd("C:/Users/abch/Dropbox (Penn)/Pairwise Comparison (Abhinav & Yichen)/Numerical-Study/")

write.csv(param_seq, "./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/p_seq.csv")

write.csv(hamming_error_est, "./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/hamming-error-est.csv")
write.csv(hamming_error_est_priv, "./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/hamming-error-est-priv.csv")

write.csv(hamming_error_count, "./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/hamming-error-count.csv")
write.csv(hamming_error_count_priv, "./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/hamming-error-count-priv.csv")


write.csv(L_infty_error, "./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/L-infty-error.csv")
write.csv(L_infty_error_priv, "./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/L-infty-error-priv.csv")

############ Reading simulation results ########

param_seq <- read.csv("./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/p_seq.csv")

hamming_error_est <- read.csv("./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/hamming-error-est.csv")
hamming_error_est_priv <- read.csv("./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/hamming-error-est-priv.csv")

hamming_error_count <- read.csv("./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/hamming-error-count.csv")
hamming_error_count_priv <- read.csv("./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/hamming-error-count-priv.csv")

L_infty_error <- read.csv("./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/L-infty-error.csv")
L_infty_error_priv <- read.csv("./Simulation-Results/accuracy-vs-p/epsilon=2.5 n=300/L-infty-error-priv.csv")


############ Plot for For Hamming Loss

#require(ggplot2)
#source("./Helper-Script.R")

ggsave(
  "./Plots/accuracy-vs-p/epsilon=2.5 n=300/hamming loss via estimation.pdf",
  plot_sim_results(
    hamming_error_est[ ,-1],
    hamming_error_est_priv[ ,-1],
    param_seq[, -1],
    "hamming loss via estimation",
    "p"
  ),
  height = 15,
  width = 12
)

ggsave(
  "./Plots/accuracy-vs-p/epsilon=2.5 n=300/hamming loss via counting.pdf",
  plot_sim_results(
    hamming_error_count[, -1],
    hamming_error_count_priv[, -1],
    param_seq[, -1],
    "hamming loss via counting",
    "p"
  ),
  height = 15,
  width = 12
)

########### For L-infty loss

ggsave(
  "./Plots/accuracy-vs-p/epsilon=2.5 n=300/l infinity loss.pdf",
  plot_sim_results(log(L_infty_error[, -1]),
                   log(L_infty_error_priv[, -1]),
                   param_seq[, -1],
                   latex2exp::TeX(r'(Logarithm of Relative $L_\infty$ loss)'),
                   "p"),
  height = 15,
  width = 12
)

