source(
  "../Real Data/helper-functions.R"
)
eps_grid = c(0.1, 0.5, 1, 2, 4, 8)
p_seq = c(0.2, 0.4, 0.6, 0.8, 0.9, 1)

param = eps_grid
n_param = length(param)

n_rep = 500

empty_data_frame <-
  data.frame(matrix(nrow = length(p_seq) * length(param),
                    ncol = n_rep))

empty_data_frame <-
  cbind(expand.grid(p_seq, param), empty_data_frame)
colnames(empty_data_frame)[1:2] <- c("p", "epsilon")



L_2_error = empty_data_frame
L_infinity_error = empty_data_frame
param_hamming_error = empty_data_frame
non_param_hamming_error = empty_data_frame


for (p_index in 1:length(p_seq)) {
  # non private estimates of theta
  p = p_seq[p_index]
  cat(p,"\n")
  for (irep in 1:n_rep) {
    pref_mat_sub = data_gen_sub_sampled(p = p, win_mat = pref_mat)
    if (sum(pref_mat_sub != 0) == 0)
      next
    optim_theta <- optim(
      par = rep(1, nrow(pref_mat_sub)),
      fn = obj_func,
      gr = grad_obj_func,
      win_mat = pref_mat_sub,
      n = nrow(pref_mat_sub),
      H = H_logistic,
      grad_H = grad_H_logistic,
      p = p,
      #proportion of paired items observed
      L = L,
      lambda = 0,
      method = "BFGS"
    )
    theta_non_priv = optim_theta$par - mean(optim_theta$par)
    rank_param_non_private <-
      items[order(optim_theta$par, decreasing = T)]
    
    prop_win = numeric(nrow(pref_mat_sub))
    prop_win_priv = numeric(nrow(pref_mat_sub))
    
    for (i in 1:nrow(pref_mat_sub)) {
      total_comparison = sum(pref_mat_sub[i, ]) + sum(pref_mat_sub[, i])
      total_win = sum(pref_mat_sub[i, ])
      prop_win[i] = total_win / total_comparison
    }
    rank_non_param_non_private <-
      items[order(prop_win, decreasing = T)]
    
    for (index in 1:n_param) {
      eps = param[index]
      #parametric
      optim_theta_priv <- optim(
        par = rep(1, nrow(pref_mat_sub)),
        fn = obj_func,
        gr = grad_obj_func,
        win_mat = pref_mat_sub,
        n = nrow(pref_mat_sub),
        H = H_logistic,
        grad_H = grad_H_logistic,
        p = p,
        L = L,
        lambda = 8 / eps,
        method = "BFGS"
      )
      theta_priv = optim_theta_priv$par - mean(optim_theta_priv$par)
      rank_param_priv = items[order(optim_theta_priv$par, decreasing = T)]
      
      #non-parametric ranking
      for (i in 1:nrow(pref_mat_sub)) {
        total_comparison = sum(pref_mat_sub[i, ]) + sum(pref_mat_sub[, i])
        total_win = sum(pref_mat_sub[i, ])
        prop_win[i] = total_win / total_comparison
        unif = runif(1)
        w = sign(unif - 0.5) * log(1 - 2 * abs(unif - 0.5))
        prop_win_priv[i] = prop_win[i] + w / (total_comparison * eps)
      }
      
      rank_non_param_priv <-
        items[order(prop_win_priv, decreasing = T)]
      
      
      # calculating the evaluation metrics
      L_2_error[L_2_error$p == p &
                  L_2_error$epsilon == eps, 2 + irep] = sqrt(sum((theta_priv - theta_non_priv) ^ 2))
      L_infinity_error[L_infinity_error$p == p &
                        L_infinity_error$epsilon == eps, 2 + irep] = max(abs(theta_priv - theta_non_priv))
      param_hamming_error[param_hamming_error$p == p &
                            param_hamming_error$epsilon == eps, 2 + irep] = 1 - mean(rank_param_priv == rank_param_non_private)
      non_param_hamming_error[non_param_hamming_error$p == p &
                                non_param_hamming_error$epsilon == eps, 2 + irep] = 1 - mean(rank_non_param_priv == rank_non_param_non_private)
    }
  }
}

# clipping the evaluation metrics (using the fact that theta is bounded by 1 in each coordinate)

L_2_error[-(1:2)][is.na(L_2_error[-(1:2)]) | L_2_error[-(1:2)] > 4*(nrow(pref_mat)-1)] = 4*(nrow(pref_mat)-1)
L_infinity_error[-(1:2)][is.na(L_infinity_error[-(1:2)]) | L_infinity_error[-(1:2)] > 2] = 2
param_hamming_error[is.na(param_hamming_error)] = 1
non_param_hamming_error[is.na(non_param_hamming_error)] = 1

#### Saving simulation results ####
write.csv(L_2_error, "./simulation-results-2/L_2_error.csv", row.names = F)
write.csv(L_infinity_error, "./simulation-results-2/L_infinity_error.csv", row.names = F)
write.csv(param_hamming_error, "./simulation-results-2/param_hamming_error.csv", row.names = F)
write.csv(non_param_hamming_error, "./simulation-results-2/non_param_hamming_error.csv", row.names = F)
