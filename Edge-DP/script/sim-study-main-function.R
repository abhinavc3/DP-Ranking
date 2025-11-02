#Comparing the three methods for top k selction
#Method 1 : perturbed MLE
#Method 2 : noisy borda count
#Method 3 : one-shot-top-k with spectral method


library(foreach)
library(doParallel)

source("Helper-Script.R")

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

# simulation study





do_sim_all <-
  function(repetition = 15,
           p_seq = c(1),
           n_seq =  c(50, 150),
           epsilon_seq = c(0.05, 100000)) {
    #repetition = 15
    #p_seq = c(1)
    #n_seq = c(100, 400)
    #k_seq = round(n_seq/4)
    #epsilon_seq = c(0.05, 100000)
    
    # First we define a data frame containing the parameter settings
    
    empty_data_frame <- expand.grid(p_seq, n_seq, epsilon_seq)
    
    theta = numeric(nrow(empty_data_frame))
    for (row in 1:nrow(empty_data_frame)) {
      n = empty_data_frame[row, 2]
      k = round(n / 4)
      fixed_theta <-
        log(c(runif(n - k, min =  0.2 , max = 0.7), rep(1, k)))
      fixed_centered_theta <- fixed_theta - mean(fixed_theta)
      theta[row] = list(fixed_centered_theta)
    }
    param_data_frame = cbind(empty_data_frame, I(theta), data.frame(matrix(
      nrow = length(p_seq) * length(n_seq) * length(epsilon_seq),
      ncol = repetition
    )))
    
    colnames(param_data_frame)[1:4] <-
      c("p", "n", "epsilon", "theta")
    
    
    noisy_borda_count_ham = param_data_frame
    noisy_one_shot_spectral = param_data_frame
    noisy_pert_mle = param_data_frame
    
    for (row in 1:nrow(empty_data_frame)) {
      n = param_data_frame$n[row]
      p = param_data_frame$p[row]
      k = round(n / 4)
      epsilon = param_data_frame$epsilon[row]
      fixed_centered_theta = param_data_frame$theta[[row]]
      S_k = recover_top_k_set(fixed_centered_theta , k = k)
      sim_list <-
        foreach(
          iter = 1:repetition,
          .combine = "rbind",
          .export = ls(globalenv())
        ) %dopar% {
          data = data_gen(n = n, theta = fixed_centered_theta, p = p)
          hat_S_k_count_borda = top_k_borda_count(win_mat = data,
                                                  k = k,
                                                  eps = epsilon)
          hat_S_k_count_spectral = top_k_one_shot_spectral(data,
                                                           p = p,
                                                           k = k,
                                                           epsilon = epsilon)
          hat_S_k_count_pert_mle = top_k_purturbed_mle(
            data,
            k = k,
            p = p,
            epsilon =  epsilon,
            fixed_centered_theta = fixed_centered_theta
          )
          hamming_error_borda <-
            1 - length(intersect(S_k, hat_S_k_count_borda)) / k
          hamming_error_spectral <-
            1 - length(intersect(S_k, hat_S_k_count_spectral)) / k
          hamming_error_pert_mle <-
            1 - length(intersect(S_k, hat_S_k_count_pert_mle)) / k
          
          list(
            hamming_error_borda = hamming_error_borda,
            hamming_error_spectral = hamming_error_spectral,
            hamming_error_pert_mle = hamming_error_pert_mle
          )
        }
      noisy_borda_count_ham[row,-(1:4)] = unlist(sim_list[, 1])
      noisy_one_shot_spectral[row,-(1:4)] = unlist(sim_list[, 2])
      noisy_pert_mle[row,-(1:4)] = unlist(sim_list[, 3])
      print(param_data_frame[row, 1:4])
    }
    return (
      list(
        borda = noisy_borda_count_ham,
        one_shot = noisy_one_shot_spectral,
        pert_mle = noisy_pert_mle
      )
    )
  }
