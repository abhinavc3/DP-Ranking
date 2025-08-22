load("Sim-results/numerical-results/rel_ham_vs_n.RData")
sim_output_df = dplyr::bind_rows(rel_ham_vs_n_NP, .id = "Method")
sim_output_df

source('Helper-Script.R')

sim_summary_data <- sim_summary(sim_output_df)
#sim_summary_data$epsilon[sim_summary_data$epsilon == 1000] = "Non-Private"
library(ggplot2)


# Define custom colors for each Method
method_colors <- c("borda" = "#619CFF", "one_shot" = "#F8766D", "pert_mle" = "#00BA38")

ggplot(sim_summary_data, aes(x = n)) +
  geom_point(aes(y = mean, colour = as.factor(Method)), size = 3) +
  geom_line(aes(
    y = mean,
    colour = as.factor(Method),
    linetype = as.factor(epsilon)
  ),
  linewidth = 1) +
  geom_errorbar(aes(
    y = mean,
    ymin = lower_quantile,
    ymax = upper_quantile,
    colour = as.factor(Method)
  )) +
  ylab("Relative Hamming Distance") +
  labs(linetype = latex2exp::TeX(r'($\epsilon$)')) +
  scale_color_manual(values = method_colors, name = "Method", labels = c("Copeland Counting", "One Shot (baseline)", "Perturbed MLE")) +
  theme_bw()

  #scale_linetype(guide = F) +
  #theme_bw() +
  #theme(text = element_text(size = 40), legend.position = c(0.8,0.8))


## relative hamming loss vs epsilon
load("Sim-results/numerical-results/rel_ham_vs_eps.RData")
sim_output_df = dplyr::bind_rows(rel_ham_vs_eps, .id = "Method")
sim_output_df

sim_summary_data <- sim_summary(sim_output_df)
#sim_summary_data$epsilon[sim_summary_data$epsilon == 1000] = "Non-Private"
ggplot(sim_summary_data, aes(x = epsilon)) +
  geom_point(aes(y = mean, colour = as.factor(Method)), size = 3) +
  geom_line(aes(
    y = mean,
    colour = as.factor(Method),
    linetype = as.factor(n)
  ),
  linewidth = 1) +
  geom_errorbar(aes(
    y = mean,
    ymin = lower_quantile,
    ymax = upper_quantile,
    colour = as.factor(Method)
  )) +
  ylab("Relative Hamming Distance") +
  xlab(latex2exp::TeX(r'($\epsilon$)')) +
  labs(linetype = 'n') +
  scale_color_manual(values = method_colors, name = "Method", labels = c("Copeland Counting", "One Shot (baseline)", "Perturbed MLE")) +
  theme_bw()
