#Paralalising Code

list.of.packages <- c("foreach",
                      "doParallel",
                      "tidyverse",
                      "ggplot2",
                      "snow",
                      "Rmpi",
                      "latex2exp")

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


################## Experiment 1 Plots ##############################

##### Reading data ######

hamming_error_est_priv <-
  read.csv("./Simulation-Results/accuracy-vs-n/hamming-error-est-priv.csv")

hamming_error_count_priv <-
  read.csv("./Simulation-Results/accuracy-vs-n/hamming-error-count-priv.csv")

L_infty_error_priv <-
  read.csv("./Simulation-Results/accuracy-vs-n/L-infty-error-priv.csv")

L_2_error_priv <-
  read.csv("./Simulation-Results/accuracy-vs-n/L-2-error-priv.csv")

sim_summary <- function(sim_result) {
  cbind(sim_result[, 1:3], t(apply(sim_result[,-(1:3)], 1, function_summary)))
}

plot_vs_n <- function(data, label_of_y) {
  sim_summary_data <- sim_summary(data)
  sim_summary_data$epsilon[sim_summary_data$epsilon == 1000] = "Non-Private"
  ggplot(sim_summary_data, aes(x = n)) +
    geom_point(aes(y = mean, colour = as.factor(epsilon)), size = 3) +
    geom_line(aes(
      y = mean,
      colour = as.factor(epsilon),
      linetype = as.factor(p)
    ),
    linewidth = 1) +
    geom_errorbar(aes(
      y = mean,
      ymin = lower_quantile,
      ymax = upper_quantile,
      colour = as.factor(epsilon)
    )) +
    ylab(label_of_y) +
    labs(colour = latex2exp::TeX(r'($\epsilon$)')) +
    scale_linetype(guide = F) +
    theme_bw() +
    theme(text = element_text(size = 40), legend.position = c(0.8,0.8))
  
}

ggsave(
  "./Plots/accuracy-vs-n/l infinity loss.pdf",
  plot_vs_n(
    L_infty_error_priv,
    latex2exp::TeX(r'(logarithm of relative $L_\infty$ loss)')
  ),
  height = 12,
  width = 12
)
ggsave(
  "./Plots/accuracy-vs-n/l 2 loss.pdf",
  plot_vs_n(
    L_2_error_priv,
    latex2exp::TeX(r'(logarithm of relative $L_2$ loss)')
  ),
  height = 12,
  width = 12
)
ggsave(
  "./Plots/accuracy-vs-n/hamming loss via estimation.pdf",
  plot_vs_n(hamming_error_est_priv, "relative hamming loss via estimation"),
  height = 12,
  width = 12
)

ggsave(
  "./Plots/accuracy-vs-n/hamming loss via counting.pdf",
  plot_vs_n(hamming_error_count_priv, "relative hamming loss via counting"),
  height = 12,
  width = 12
)

################## Experiment 2 Plots ##############################

##### Reading data ######

hamming_error_est_priv <-
  read.csv("./Simulation-Results/accuracy-vs-p/hamming-error-est-priv.csv")

hamming_error_count_priv <-
  read.csv("./Simulation-Results/accuracy-vs-p/hamming-error-count-priv.csv")

L_infty_error_priv <-
  read.csv("./Simulation-Results/accuracy-vs-p/L-infty-error-priv.csv")

L_2_error_priv <-
  read.csv("./Simulation-Results/accuracy-vs-p/L-2-error-priv.csv")

sim_summary <- function(sim_result) {
  cbind(sim_result[, 1:3], t(apply(sim_result[,-(1:3)], 1, function_summary)))
}

plot_vs_p <- function(data, label_of_y) {
  sim_summary_data <- sim_summary(data)
  sim_summary_data$epsilon[sim_summary_data$epsilon == 1000] = "Non-Private"
  ggplot(sim_summary_data, aes(x = p)) +
    geom_point(aes(y = mean, colour = as.factor(epsilon)), size = 3) +
    geom_line(aes(
      y = mean,
      colour = as.factor(epsilon),
      linetype = as.factor(n)
    ),
    linewidth = 1) +
    geom_errorbar(aes(
      y = mean,
      ymin = lower_quantile,
      ymax = upper_quantile,
      colour = as.factor(epsilon)
    )) +
    ylab(label_of_y) +
    labs(colour = latex2exp::TeX(r'($\epsilon$)')) +
    scale_linetype(guide = F) +
    theme_bw() +
    theme(text = element_text(size = 40), legend.position = c(0.8,0.8))
  
}

ggsave(
  "./Plots/accuracy-vs-p/l infinity loss.pdf",
  plot_vs_p(
    L_infty_error_priv,
    latex2exp::TeX(r'(logarithm of relative $L_\infty$ loss)')
  ),
  height = 12,
  width = 12
)

ggsave(
  "./Plots/accuracy-vs-p/l 2 loss.pdf",
  plot_vs_p(
    L_2_error_priv,
    latex2exp::TeX(r'(logarithm of relative $L_2$ loss)')
  ),
  height = 12,
  width = 12
)

ggsave(
  "./Plots/accuracy-vs-p/hamming loss via estimation.pdf",
  plot_vs_p(hamming_error_est_priv, "relative hamming loss via estimation"),
  height = 12,
  width = 12
)

ggsave(
  "./Plots/accuracy-vs-p/hamming loss via counting.pdf",
  plot_vs_p(hamming_error_count_priv, "relative hamming loss via counting"),
  height = 12,
  width = 12
)

################## Experiment 3 Plots ##############################

##### Reading data ######

hamming_error_est_priv <-
  read.csv("./Simulation-Results/accuracy-vs-epsilon/hamming-error-est-priv.csv")

hamming_error_count_priv <-
  read.csv("./Simulation-Results/accuracy-vs-epsilon/hamming-error-count-priv.csv")

L_infty_error_priv <-
  read.csv("./Simulation-Results/accuracy-vs-epsilon/L-infty-error-priv.csv")

L_2_error_priv <-
  read.csv("./Simulation-Results/accuracy-vs-epsilon/L-2-error-priv.csv")
  
sim_summary <- function(sim_result) {
  cbind(sim_result[, 1:3], t(apply(sim_result[,-(1:3)], 1, function_summary)))
}

plot_vs_epsilon <- function(data, label_of_y) {
  ggplot(sim_summary(data), aes(x = 1/epsilon)) +
    geom_point(aes(y = mean, colour = as.factor(p)), size = 3) +
    geom_line(aes(
      y = mean,
      colour = as.factor(p),
      linetype = as.factor(n)
    ),
    linewidth = 1) +
    geom_errorbar(aes(
      y = mean,
      ymin = lower_quantile,
      ymax = upper_quantile,
      colour = as.factor(p)
    )) +
    ylab(label_of_y) +
    xlab(latex2exp::TeX(r'($1/\epsilon$)')) +
    labs(colour = "p") +
    scale_linetype(guide = F) +
    theme_bw() +
    theme(text = element_text(size = 40), legend.position = c(0.8,0.8)) +
    scale_x_reverse()
  
}

ggsave(
  "./Plots/accuracy-vs-epsilon/l infinity loss.pdf",
  plot_vs_epsilon(
    L_infty_error_priv,
    latex2exp::TeX(r'(logarithm of relative $L_\infty$ loss)')
  ),
  height = 12,
  width = 12
)

ggsave(
  "./Plots/accuracy-vs-epsilon/l 2 loss.pdf",
  plot_vs_epsilon(
    L_2_error_priv,
    latex2exp::TeX(r'(logarithm of relative $L_2$ loss)')
  ),
  height = 12,
  width = 12
)


ggsave(
  "./Plots/accuracy-vs-epsilon/hamming loss via estimation.pdf",
  plot_vs_epsilon(hamming_error_est_priv, "relative hamming loss via estimation"),
  height = 12,
  width = 12
)

ggsave(
  "./Plots/accuracy-vs-epsilon/hamming loss via counting.pdf",
  plot_vs_epsilon(hamming_error_count_priv, "relative hamming loss via counting"),
  height = 12,
  width = 12
)
