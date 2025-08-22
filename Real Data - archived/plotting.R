library(ggplot2)
# loading data

L_2_error <- read.csv("./simulation-results-2/L_2_error.csv")
L_infinity_error <- read.csv("./simulation-results-2/L_infinity_error.csv")
param_hamming_error <- read.csv("./simulation-results-2/param_hamming_error.csv")
non_param_hamming_error <- read.csv("./simulation-results-2/non_param_hamming_error.csv")


plot_vs_p <- function(data, label_of_y) {
  ggplot(sim_summary(data), aes(x = p)) +
    geom_point(aes(y = mean, shape = as.factor(epsilon),color = as.factor(epsilon)), size = 3) +
    geom_line(aes(
      y = mean,
      linetype = as.factor(epsilon),
      color = as.factor(epsilon)
    ),
    linewidth = 1) +
    geom_errorbar(aes(
      y = mean,
      ymin = lower_quantile,
      ymax = upper_quantile,
      colour = as.factor(epsilon)
    )) +
    ylab(label_of_y) +
    #labs(color =latex2exp::TeX(r'($\epsilon$)')) +
    #labs(linetype = "n") +
    theme_bw() +
    theme(text = element_text(size = 40))
}

plot_vs_epsilon <- function(data, label_of_y) {
  ggplot(sim_summary(data), aes(x = epsilon)) +
    geom_point(aes(y = mean, colour = as.factor(p)), size = 3) +
    geom_line(aes(
      y = mean,
      colour = as.factor(p)
    ),
    linewidth = 1) +
    geom_errorbar(aes(
      y = mean,
      ymin = lower_quantile,
      ymax = upper_quantile,
      colour = as.factor(p)
    )) +
    ylab(label_of_y) +
    xlab(latex2exp::TeX(r'($\epsilon$)')) +
    labs(colour = latex2exp::TeX(r'($p$)')) +
    labs(linetype = "n") +
    theme_bw() +
    theme(text = element_text(size = 40))
}

ggsave(
  "./Plots/accuracy-vs-p/L_2_loss.pdf",
  plot_vs_p(
    L_2_error,
    latex2exp::TeX(r'($L_2$ loss)')
  ),
  height = 15,
  width = 14
)

ggsave(
  "./Plots/accuracy-vs-p/L_infinity_loss.pdf",
  plot_vs_p(
    L_infinity_error,
    latex2exp::TeX(r'($L_\infty$ loss)')
  ),
  height = 15,
  width = 14
)
ggsave(
  "./Plots/accuracy-vs-p/param_hamming_error.pdf",
  plot_vs_p(
    param_hamming_error,
    latex2exp::TeX(r'(hamming loss via estimation)')
  ),
  height = 15,
  width = 14
)

ggsave(
  "./Plots/accuracy-vs-p/non_param_hamming_error.pdf",
  plot_vs_p(
    non_param_hamming_error,
    latex2exp::TeX(r'(hamming loss via counting)')
  ),
  height = 15,
  width = 14
)


ggsave(
  "./Plots/accuracy-vs-epsilon/L_2_loss.pdf",
  plot_vs_epsilon(
    L_2_error,
    latex2exp::TeX(r'($L_2$ loss)')
  ),
  height = 15,
  width = 14
)

ggsave(
  "./Plots/accuracy-vs-epsilon/L_infinity_loss.pdf",
  plot_vs_epsilon(
    L_infinity_error,
    latex2exp::TeX(r'($L_\infty$ loss)')
  ),
  height = 15,
  width = 14
)
ggsave(
  "./Plots/accuracy-vs-epsilon/param_hamming_error.pdf",
  plot_vs_epsilon(
    param_hamming_error,
    latex2exp::TeX(r'(hamming loss via estimation)')
  ),
  height = 15,
  width = 14
)

ggsave(
  "./Plots/accuracy-vs-epsilon/non_param_hamming_error.pdf",
  plot_vs_epsilon(
    non_param_hamming_error,
    latex2exp::TeX(r'(hamming loss via counting)')
  ),
  height = 15,
  width = 14
)

