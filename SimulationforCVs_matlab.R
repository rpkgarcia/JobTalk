# Simulation for CVS direct recreation of Matlab
library(MASS)

big_T = 1000
nrep = 50
m_max = 12

kfrac = seq(0, .99, by = 0.005)
kfrac = matrix(kfrac)        # b values 
kvec = floor(big_T*kfrac)    # M values
n_kvec = length(kvec)

wght = matrix(0, n_kvec, big_T)
wght[,1] = 1

# Make a weights matrix 
for(ii in 1:n_kvec){
  k = kvec[ii, ]
  
  if(k >0){
    wght[ii, 2:(k+1)] = k:1/(k+1)
  }
}

F_save = array(NA, dim = c(nrep, n_kvec, m_max))


R = function(h, the_sim_data){
  big_T = nrow(the_sim_data)
  index = 1:(big_T -h)
  
  # Already centered
  est = lapply(index, function(i){
    est = the_sim_data[i, ]%*%t(the_sim_data[(i+h), ])
  })
  
  # Sum together
  autocov_s = Reduce('+', est)/big_T
  
  # Because of symmetry 
  if(h!=0){
    autocov_s = autocov_s + t(autocov_s)
  }
  
  return(autocov_s)
}


u_matlab <- read.csv("~/Documents/GitHub/2022NISS_GradConference/Recreations/u_data_from_matlab.csv", 
                               header=FALSE)
u_matlab = as.matrix(u_matlab)
u = u_matlab

for(irep in 1:nrep){
  u = mvrnorm(big_T, mu = rep(0, m_max), Sigma = diag(m_max))
  ubar = colMeans(u)
  
  u_center = apply(u, 1, function(row) row - ubar)
  u_center = t(u_center)
  all_autocovariances <-sapply(0:(big_T-1), R, the_sim_data= u_center)
  # First 12 
  all_autocovariances <- t(all_autocovariances)
  
  acv = array(NA, dim = c(big_T, m_max, m_max))
  for(i in 1:big_T){
    acv[i, , ] = matrix(all_autocovariances[i, ], nrow=m_max, ncol = m_max, 
                        byrow = T)
  }
  
  omega_hat_mat = array(NA, dim = c(n_kvec, m_max, m_max))
  
  for( ii in 1:m_max){
    for(jj in ii:m_max){
      omega_hat_mat[, ii, jj] = c(wght%*%acv[, ii, jj])/big_T
      omega_hat_mat[,jj, ii] = omega_hat_mat[,ii, jj]
    }
  }
  
  for(ik in 1:n_kvec){
    for(m in 1:m_max){
      omega_hat = omega_hat_mat[ik, 1:m, 1:m]
      omega_inv = solve(omega_hat)
      F_save[irep, ik, m]  = t(ubar[1:m])%*%omega_inv%*%matrix(ubar[1:m])/m
    }
  }
  # F_save[, , ]: 
      # first dimension has replication (tables)
      # second dimension has b value (rows)
      # third dimension has m value (dimensions of data set)
}