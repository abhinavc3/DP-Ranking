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
  set.seed(NULL)
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
  set.seed(NULL)
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


top_k_purturbed_mle <- function(data, p, epsilon, fixed_centered_theta){
  k = round(n/4)
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
  return(recover_top_k_set(optim_theta_priv$par, k = k))
}


################################################################################
########################### One Shot Method ####################################

#' Title
#'
#' @param data the win matrix
#' @param p sampling probability
#' @param k 
#' @param epsilon 
#' @param n 
#' @param delta 
#'
#' @return indices of the top k set
#' @export
#'
#' @examples
top_k_one_shot_spectral <-
  function(data,
           p,
           k,
           epsilon,
           n = nrow(data),
           delta = 1 / n ^ 2) {
    # normalization constant
    d = n * p
    # Defining the Transition matrix
    
    # creating the trans matrix
    # the data matrix is transposed because y_{ij} defined in their paper is equal
    # to y_{ji} for us.
    data_mod = t(data)
    
    # Replace all occurrences of -1 with 0
    data_mod[data_mod == -1] = 0  # makes the diagonal and missing comparisons 0
    
    # Calculate the sum of each row in the modified data matrix
    diagonal = d - rowSums(data_mod)
    
    # Add the diagonal values to the modified data matrix
    data_mod = data_mod + diag(diagonal)
    
    # Normalize the modified data matrix to create the transition matrix
    transition_matrix = data_mod / d
    
    
    # finding the stationary dist
    
    # Compute the eigenvalues and eigenvectors of the transpose of the transition matrix
    eigen_info <- eigen(t(transition_matrix))
    
    # Find the index of the eigenvalue closest to 1 (within a small tolerance)
    index <- which(abs(eigen_info$values - 1) < 1e-8)
    
    # Extract the corresponding eigenvector
    stationary_distribution <- Re(eigen_info$vectors[, index])
    
    # Normalize the stationary distribution
    stationary_distribution <-
      stationary_distribution / sum(stationary_distribution)
    erg_coef = ergodicity_coef(transition_matrix)
    
    top_k_set = top_k_one_shot(
      stationary_distribution,
      k = k,
      d = d,
      rho = erg_coef,
      eps = epsilon
    )
    return(top_k_set)
  }

# defining a function to recover the top-k employing their method

#' Title
#'
#' @param prob_vec (stationary distribution) vector for which top-k is to extracted
#' @param d normalization for transition matrix is set at n*p
#' @param rho an upper bound for ergodicity coefficient
#' @param k 
#' @param eps 
#' @param delta 
#' @param L number of pairwise comparisons
#'
#' @return DP top-k set
#' @export
#'
#' @examples
top_k_one_shot <-
  function(prob_vec,
           d,
           rho,
           k = round(length(prob_vec) / 3),
           eps = 0.1,
           delta = 1 / length(prob_vec) ^ 2,
           L = 1
  ) {
    n = length(prob_vec)
    s = 2 / (d * L * (1 - rho))
    lambda = (8 * s * sqrt(k * log(n / delta)))/eps 
    unif = runif(n)
    w = sign(unif - 0.5) * log(1 - 2 * abs(unif - 0.5))
    return(order(prob_vec + w * lambda, decreasing = T)[1:k])
  }

#alternate expression for erg coef

#ergodicity_coef <- function(S) {
#  result <- 1 - min(apply(S, 1, function(row_i) {
#    apply(S, 1, function(row_j) {
#      sum(pmin(row_i, row_j))
#    })
#  }))
#  return(result)
#}

#' Title
#'
#' @param S a square matrix
#'
#' @return the ergodocity coefficient
#' @export
#'
#' @examples
ergodicity_coef <- function(S) {
  max_sum <- max(apply(S, 1, function(row_i) {
    apply(S, 1, function(row_j) {
      sum(abs(row_i - row_j))
    })
  }))
  result <- 0.5 * max_sum
  return(result)
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

library(ggplot2)

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
           param_seq,
           label_of_y = "Loss",
           label_of_x) {
    # mean and sd for non private sim study
    sim_summary <- t(apply(sim_results, 1, function_summary))
    sim_summary <-
      data.frame(type = rep("Non Private", length(param_seq)), param_seq, sim_summary)
    
    
    # mean and sd for private sim study
    sim_summary_priv <-
      t(apply(sim_results_priv, 1, function_summary))
    sim_summary_priv <-
      data.frame(type = rep("Private", length(param_seq)), param_seq, sim_summary_priv)
    
    merged_sim_summary <- rbind(sim_summary, sim_summary_priv)
    
    ggplot(merged_sim_summary, aes(x = param_seq)) +
      geom_line(aes(y = mean, colour = type)) +
      geom_point(aes(y = mean, colour = type)) +
      geom_errorbar(aes(
        y = mean,
        ymin = lower_quantile,
        ymax = upper_quantile,
        colour = type
      )) +
      xlab(label_of_x) +
      ylab(label_of_y) +
      scale_colour_discrete(name='Method') +
      theme_bw() +
      theme(legend.position = c(0.7,0.7),
            text = element_text(size = 40)
            ) 
    
  }

plot_sim_results_single <-
  function(sim_results,
           param_seq,
           label_of_y = "Loss",
           label_of_x) {
    # mean and sd for non private sim study
    sim_summary <- t(apply(sim_results, 1, function_summary))
    sim_summary <-
      data.frame(param_seq, sim_summary)
    
    
   
    
    ggplot(sim_summary, aes(x = param_seq)) +
      geom_line(aes(y = mean)) +
      geom_point(aes(y = mean)) +
      geom_errorbar(aes(
        y = mean,
        ymin = lower_quantile,
        ymax = upper_quantile
      )) +
      xlab(label_of_x) +
      ylab(label_of_y) +
      theme_bw() +
      theme(text = element_text(size = 40)) +
      scale_x_reverse()
  }

sim_summary <- function(sim_result) {
  cbind(sim_result[, 1:3], t(apply(sim_result[,-(1:3)], 1, function_summary)))
}