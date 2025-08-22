########loading CEMS data ############

library(BradleyTerry2)
data("CEMS")
items = c("London","Paris","Milano","St.Gallen","Barcelona","Stockholm")

df.pref = CEMS$preferences
library(dplyr)
#df.pref %>% filter(school1 == "London" & school2 == "Milano") %>% select(win1.adj) %>% sum()

L = nrow(CEMS$students)
### defining the win matrix #####
pref_mat = matrix(0,nrow = 6, ncol = 6)

for(i in 1:6) {
  j = i + 1
  while (j <= 6) {
    pref_mat[i, j] = df.pref %>% filter(school1 == items[i] &
                                          school2 == items[j]) %>% select(win1.adj) %>% sum(na.rm = T)
    pref_mat[j, i] = df.pref %>% filter(school1 == items[i] &
                                          school2 == items[j]) %>% select(win2.adj) %>% sum(na.rm = T)
    j = j + 1
  }
}

########## loading immigration data ##########
library(prefmod)
data("immig")
items = c("crimRate", "position", "socBurd", "culture")
L = nrow(immig)
pref_mat = matrix(0,nrow = 4, ncol = 4)
pref_mat[1,2] = sum(immig$V12==1, na.rm = T) 
pref_mat[2,1] = sum(immig$V12==-1, na.rm = T)
pref_mat[1,3] = sum(immig$V13==1, na.rm = T) 
pref_mat[3,1] = sum(immig$V13==-1, na.rm = T)
pref_mat[1,4] = sum(immig$V14==1, na.rm = T) 
pref_mat[4,1] = sum(immig$V14==-1, na.rm = T)
pref_mat[2,3] = sum(immig$V23==1, na.rm = T) 
pref_mat[3,2] = sum(immig$V23==-1, na.rm = T)
pref_mat[2,4] = sum(immig$V24==1, na.rm = T) 
pref_mat[4,2] = sum(immig$V24==-1, na.rm = T)
pref_mat[3,4] = sum(immig$V34==1, na.rm = T) 
pref_mat[4,3] = sum(immig$V34==-1, na.rm = T)

############## paramteric ############
############## running optim ########

# non-private #
optim_theta <- optim(
  par = rep(1,nrow(pref_mat)),
  fn = obj_func,
  gr = grad_obj_func,
  win_mat = pref_mat,
  n = nrow(pref_mat),
  H = H_logistic,
  grad_H = grad_H_logistic,
  p = 1,
  L = L,
  lambda = 0,
  method = "BFGS"
)
items[order(optim_theta$par, decreasing = T)]
# private #
eps = 1
optim_theta_priv <- optim(
  par = rep(1,4),
  fn = obj_func,
  gr = grad_obj_func,
  win_mat = pref_mat,
  n = nrow(pref_mat),
  H = H_logistic,
  grad_H = grad_H_logistic,
  p = 1, 
  L = L,
  lambda = 8 / eps,
  method = "BFGS"
)
items[order(optim_theta_priv$par, decreasing = T)]

############### non parametric ##########

## we rank using the proportion of pairwise comparison where that item wins ##

prop_win = numeric(nrow(pref_mat))
prop_win_priv = numeric(nrow(pref_mat))

for(i in 1:nrow(pref_mat)){
  total_comparison = sum(pref_mat[i,]) + sum(pref_mat[,i])
  total_win = sum(pref_mat[i,])
  prop_win[i] = total_win/total_comparison
  unif = runif(1)
  w = sign(unif - 0.5) * log(1 - 2 * abs(unif - 0.5))
  prop_win_priv[i] = prop_win[i] + w/(total_comparison*eps)
}
prop_win
items[order(prop_win, decreasing = T)]
prop_win_priv
items[order(prop_win_priv, decreasing = T)]

