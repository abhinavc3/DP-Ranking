H_logistic <- function(x) {
  return(1 / (1 + exp(-x)))
}
grad_H_logistic <- function(x) {
  return(exp(-x) / ((1 + exp(-x)) ^ 2))
}

#' Title : negative log-likelihood under the $H$ model with multiple comparison
#'
#' @param n : number of players
#' @param theta : strength of players
#' @param win_mat : number of comaparisons where i wins over j
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
          (-win_mat[i, j] * log(H(theta[i] - theta[j])) - win_mat[j, i] * log(1 - H(theta[i] - theta[j])))
        j <- j + 1
      }
    }
    return(log_hood)
  }

#' Title gradient of -ve log likelihood
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
          (
            -win_mat[i, j] * grad_H(theta[i] - theta[j]) / H(theta[i] - theta[j]) +
              win_mat[j, i] * grad_H(theta[i] - theta[j]) / (1 - H(theta[i] - theta[j]))
          ) *
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
                     L = 1,
                     gamma = 2 * sqrt(n * p * log(n) / L),
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
                          L = 1,
                          gamma = 2 * sqrt(n * p * log(n) / L),
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

#' Title
#'
#' @param p : proportion of items compared
#' @param win_mat : complete win matrix
#'
#' @return sub_win_mat : sub sampled win matrix
#' @export
#'
#' @examples
data_gen_sub_sampled <- function(p = 1,
                     win_mat = diag(3)) {
  random_matrix = matrix(rbinom(nrow(win_mat) ^ 2, 1, p), nrow = nrow(win_mat))
  random_matrix[upper.tri(random_matrix, diag = T)] = 0
  sub_sample = random_matrix + t(random_matrix)
  sub_win_mat = win_mat * sub_sample
  return(sub_win_mat)
}



#### plotting helper ####
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

sim_summary <- function(sim_result) {
  cbind(sim_result[, 1:2], t(apply(sim_result[,-(1:2)], 1, function_summary)))
}
