##### One-shot-DP-top-k with spectral method########

source("./Helper-Script.R")


#Initialise parameters
p <- 1
k = 10
n = 1000
fixed_theta <- log(c(runif(n - k, min =  0.2 , max = 0.7), rep(1, k)))
fixed_centered_theta <- fixed_theta - mean(fixed_theta)
S_k = recover_top_k_set(fixed_centered_theta , k = k)
epsilon <- 2.5
delta = 1/n^2

data <- data_gen(n = n, theta = fixed_centered_theta, p = p)

data

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

ergodicity_coef <- function(S) {
  max_sum <- max(apply(S, 1, function(row_i) {
    apply(S, 1, function(row_j) {
      sum(abs(row_i - row_j))
    })
  }))
  result <- 0.5 * max_sum
  return(result)
}

top_k_one_shot_spectral(data, p= p, k = k, eps = 1000 )
top_k_borda_count(data, k = k, eps = epsilon)

S_k
