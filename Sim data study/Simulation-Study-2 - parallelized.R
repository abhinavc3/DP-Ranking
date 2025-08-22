#automatic install of packages if they are not installed already
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
#  suppressPackageStartupMessages(library(package.i,
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



####################################################################################
####################################################################################
#########################   Helper Functions #######################################
####################################################################################

H_logistic <- function(x) {
  return(1 / (1 + exp(-x)))
}
grad_H_logistic <- function(x) {
  return(exp(-x) / ((1 + exp(-x)) ^ 2))
}



#' Title : data generating function
#'
#' @param n : number of players
#' @param p : probability of two players competing
#' @param theta : strength vector of the players
#' @param H : function determining win probabilities based on relative strength
#'
#' @return : matrix of wins from pairwise comparison (-1 indicates no matches were played)
#' @export
#'
#' @examples data_gen(H = H_logistic)
data_gen <-
  function(n = 10,
           p = 0.5,
           theta = runif(n),
           H = H_logistic) {
    win_matrix <- matrix(rep(-1, n * n), nrow = n, ncol = n)
    for (i in 1:n) {
      j <- i + 1
      while (j <= n) {
        if (rbinom(1, 1, p)) {
          win_matrix[i, j] <- rbinom(1, 1, H(theta[i] - theta[j]))
          win_matrix[j, i] <-  1 - win_matrix[i, j]
        }
        j <- j + 1
      }
    }
    return(win_matrix)
  }

#' Title : log-likelihood under the $H$ model
#'
#' @param n : number of players
#' @param theta : strength of players
#' @param win_mat : the win matrix
#' @param H : function determining win probabilities based on relative strength
#'
#' @return log-likelihood of the observed data at theta
#' @export
#'
#' @examples log_likelihood(win_mat = data_gen())
log_likelihood <-
  function(n = 10,
           theta = runif(n),
           win_mat,
           H = H_logistic) {
    log_hood = 0
    for (i in 1:n) {
      j <- i + 1
      while (j <= n) {
        log_hood = log_hood +
          (win_mat[i, j] >= 0) * (-win_mat[i, j] * log(H(theta[i] - theta[j])) - win_mat[j, i] * log(1 - H(theta[i] - theta[j])))
        j <- j + 1
      }
    }
    return(log_hood)
  }


#' Title
#'
#' @param n : number of players
#' @param theta : strength of players
#' @param win_mat : the win matrix
#' @param grad_H : gradient of H
#'
#' @return gradient of the log likelihood
#' @export
#'
#' @examples grad_log_lik(win_mat = data_gen())
grad_log_lik <-
  function(n = 10,
           theta = runif(n),
           win_mat,
           H = H_logistic,
           grad_H = grad_H_logistic) {
    grad_log_hood = 0
    for (i in 1:n) {
      j <- i + 1
      while (j <= n) {
        grad_log_hood = grad_log_hood +
          (win_mat[i, j] >= 0) *
          (H(theta[i] - theta[j]) - win_mat[i, j]) *
          (grad_H(theta[i] - theta[j]) /
             (H(theta[i] - theta[j]) * (1 - H(
               theta[i] - theta[j]
             )))) *
          (replace(numeric(n), i, 1) - replace(numeric(n), j, 1))
        j <- j + 1
      }
    }
    return(grad_log_hood)
  }

obj_func <- function(theta,
                     win_mat,
                     n,
                     H = H_logistic,
                     grad_H = grad_H_logistic,
                     p = 0.5,
                     eps = 0.3,
                     gamma = 2 * sqrt(n * p * log(n)),
                     lambda = 8 / eps) {
  set.seed(13)
  unif = runif(n)
  w = sign(unif - 0.5) * log(1 - 2 * abs(unif - 0.5))
  return(
    log_likelihood(
      n = n,
      theta = theta,
      win_mat = win_mat,
      H = H
    ) + 0.5 * gamma * sum(theta * theta) + lambda *
      sum(w * theta)
  )
}

grad_obj_func <- function(theta,
                          n = length(theta),
                          win_mat,
                          H = H_logistic,
                          grad_H = grad_H_logistic,
                          p = 0.5,
                          eps = 0.3,
                          gamma = 2 * sqrt(n * p * log(n)),
                          lambda = 8 / eps) {
  set.seed(13)
  unif = runif(n)
  w = sign(unif - 0.5) * log(1 - 2 * abs(unif - 0.5))
  return(
    grad_log_lik(
      n = n,
      theta = theta,
      win_mat = win_mat,
      H = H,
      grad_H = grad_H
    ) + gamma * theta + w * lambda
  )
}

################################################################################
# Non - Parametric
################################################################################
#' Title : True Top k set corresponding to theta
#'
#' @param theta : strength vector of players
#' @param k : the k of the top k
#'
#' @return : index of the top-k players in decresing order
#' @export
#'
#' @examples recover_top_k_set(c(2,3,1,4), 2)
recover_top_k_set <- function(theta, k = round(length(theta) / 3)) {
  return(order(theta, decreasing = T)[1:k])
}

#' Title borda counting algorithm
#'
#' @param win_mat : matrix of wins and losses
#' @param k : k of the top-k
#'
#' @return : the top k players ranked by decresing number of wins
#' @export
#'
#' @examples top_k_borda_count(win_mat = data_gen())
top_k_borda_count <-
  function(win_mat,
           k = round(nrow(win_mat) / 3),
           eps = 0.1,
           lambda = 2 / eps) {
    unif = runif(nrow(win_mat))
    w = sign(unif - 0.5) * log(1 - 2 * abs(unif - 0.5))
    return(order(rowSums(win_mat == 1) + w * lambda, decreasing = T)[1:k])
  }

################################################################################
########################## Plotting Helpers ####################################


function_summary <- function(x) {
  mean = mean(x)
  sd = sd(x) 
  return(c(
    mean = mean,
    sd = sd,
    lower_quantile = mean - 1.96 * sd / sqrt(length(x)),
    upper_quantile = mean + 1.96 * sd / sqrt(length(x))
  ))
}


#' Title : Plots the loss value in private v/s non private setting with increasing values of n
#'
#' @param sim_results : a matrix containing the loss value of for each value of n and iteration
#' @param sim_results_priv : same as last but private counterpart
#' @param label_of_y : label of the y axis
#'
#' @return ggplot
#' @export
#'
#' @examples
plot_sim_results <-
  function(sim_results,
           sim_results_priv,
           label_of_y = "Loss",
           label_of_x) {
    # mean and sd for non private sim study
    sim_summary <- t(apply(sim_results, 1, function_summary))
    sim_summary <-
      data.frame(type = rep("NP", length(n_seq)), n_seq, sim_summary)
    
    
    # mean and sd for private sim study
    sim_summary_priv <-
      t(apply(sim_results_priv, 1, function_summary))
    sim_summary_priv <-
      data.frame(type = rep("P", length(n_seq)), n_seq, sim_summary_priv)
    
    merged_sim_summary <- rbind(sim_summary, sim_summary_priv)
    
    ggplot(merged_sim_summary, aes(x = n_seq)) +
      geom_line(aes(y = mean, colour = type)) +
      geom_point(aes(y = mean, colour = type)) +
      geom_errorbar(aes(
        y = mean,
        ymin = lower_quantile,
        ymax = upper_quantile,
        colour = type
      )) +
      xlab(paste0(label_of_x)) +
      ylab(paste0(label_of_y)) +
      theme_bw()
    
  }


#Experiment - 1

#n_seq <- seq(300, 450, 500)
epsilon <- 0.5
n_seq <- c(200, 300, 450, 700)
no_of_iter <- 30
L_infty_error <- matrix(nrow = length(n_seq), ncol = no_of_iter)
L_infty_error_priv <-
  matrix(nrow = length(n_seq), ncol = no_of_iter)

hamming_error_est <- matrix(nrow = length(n_seq), ncol = no_of_iter)
hamming_error_est_priv <-
  matrix(nrow = length(n_seq), ncol = no_of_iter)

hamming_error_count <-
  matrix(nrow = length(n_seq), ncol = no_of_iter)
hamming_error_count_priv <-
  matrix(nrow = length(n_seq), ncol = no_of_iter)

writeLines(c(""), "log.txt")
for (i in 1:length(n_seq)) {
  n <- n_seq[i]
  k <- round(n / 4)
  fixed_theta <- log(c(runif(n - k, min =  0.2 , max = 0.7), rep(1, k)))
  #fixed_theta <- exp(runif(n, min = 0.5, max = 1))
  #fixed_theta <- log(runif(n, min = 0.5, max = 1))
  #fixed_theta <- log(seq(0.5, 1, length.out = n))
  
  
  fixed_centered_theta <- fixed_theta - mean(fixed_theta)
  S_k = recover_top_k_set(fixed_centered_theta , k = k)
  #hamming_error <- numeric(no_of_iter)
  #hamming_error_priv <- numeric(no_of_iter)
  sim_list <-
    foreach(iter = 1:no_of_iter, .combine = "rbind") %dopar% {
      data <- data_gen(n = n, theta = fixed_centered_theta, p = 1)
      #hat_S_k <- top_k_borda_count(win_mat = data, k = k, lambda = 0)
      #hat_S_k_priv <- top_k_borda_count(win_mat = data, k = k)
      optim_theta <- optim(
        par = fixed_centered_theta,
        fn = obj_func,
        gr = grad_obj_func,
        win_mat = data,
        n = nrow(data),
        p = 1,
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
        p = 1,
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
      cat(paste("\n", "Iteration ", iter, " for n = ", n),
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

############ Plot for For Hamming Loss

#require(ggplot2)

ggsave(
  "./Plots/hamming loss via estimation.pdf",
  plot_sim_results(
    hamming_error_est,
    hamming_error_est_priv,
    "hamming loss via estimation",
    "n"
  )
)

ggsave(
  "./Plots/hamming loss via counting.pdf",
  plot_sim_results(hamming_error_count, hamming_error_count_priv, "hamming loss via counting", "n")
)

########### For L-infty loss

ggsave(
  "./Plots/l infinity loss.pdf",
  plot_sim_results(L_infty_error, L_infty_error_priv, "l infty loss", "n")
)
