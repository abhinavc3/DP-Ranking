source('sim-study-main-function.R')

### relative hamming loss vs n ###
n_seq = c(150, 250, 350, 450, 550, 650)
eps_seq = c(0.1, 2)
p_seq = c(1)
reps = 45
rel_ham_vs_n = do_sim_all(repetition = reps, n_seq = n_seq, epsilon_seq = eps_seq,p_seq = p_seq)
save(rel_ham_vs_n, file = "Sim-results/numerical-results/rel_ham_vs_n.RData")

### relative hamming loss vs n Non_private ######

n_seq = c(150, 250, 350, 450, 550, 650)
eps_seq = c(100000)
p_seq = c(1)
reps = 30
rel_ham_vs_n_NP = do_sim_all(repetition = reps, n_seq = n_seq, epsilon_seq = eps_seq,p_seq = p_seq)
save(rel_ham_vs_n_NP, file = "Sim-results/numerical-results/rel_ham_vs_n_NP.RData")


### relative hamming loss vs epsilon ###
n_seq = c(350)
eps_seq = c(0.1, 0.5, 1, 2, 5)
p_seq = c(1)
reps = 45
rel_ham_vs_eps = do_sim_all(repetition = reps, n_seq = n_seq, epsilon_seq = eps_seq,p_seq = p_seq)
save(rel_ham_vs_eps, file = "Sim-results/numerical-results/rel_ham_vs_eps.RData")
