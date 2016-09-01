###############################################################################
###############################################################################
###GIBBS SAMPLER FOR Bayesian Network-Response MODEL###########################
###############################################################################
###############################################################################
# logit(\pi^i_{vu} = Z_{vu} + sum_{r=1}^R X^i_{vr}*X^i_{ur}
# Z_{vu} ~ N(mu_z,phi^{-1}_z)
# X^i_{vr} ~ N(mu_{vr}(y),phi^{-1}_x)
# mu_{vr}(y) = sum_{k=1}^K G_{vk}W_{kr}(y)
# assume GP prior W_{kr}(y) ~ GP(0,c_w)
# G_{vk} ~ N(0,tau^{-1}_k)
# N subjects and n unique IQ level

###############################################################################
###############################################################################
#CLEAR WORK SPACE AND LOAD USEFUL LIBRARIES####################################
###############################################################################
###############################################################################

rm(list=ls())

library(BayesLogit)
library(MASS) #mvrnorm
library(gdata) #lowerTriangle
library(coda) # HPDinterval
library(ROCR)
library(parallel) # mclapply
detectCores()
library(ggplot2)

#################################################################################
#################################################################################
#LOAD DATA (available via VxVxN array with N the number of networks and V the nodes)
#################################################################################
#################################################################################
load("brain_binary_IQ.RData")

#Define Model Dimensions
V<-dim(A)[1] # number of regions
N<-dim(A)[3] # number of subjects

Y=FSIQ
Y_s = sort(unique(Y))
n = length(Y_s)
c_num = as.integer(table(Y))

D = diag(c_num)
map = matrix(0, nr=n, nc=max(c_num) ) # store the indices of Y that correspond to each value of Y_s
for (j in 1:n){
  map[j,1:c_num[j]] = which(Y == Y_s[j])
}

index_y = match(Y,Y_s) # sometimes we need to expand W

###############################################################################
###############################################################################
#NUMBER OF GIBBS SAMPLES AND TRUNCATION LEVELS#################################
###############################################################################
###############################################################################

N_sampl<-5000
burnin = 1000
thin = 4
N_s = (N_sampl - burnin)/thin

#-----------------------------------------------------------------------------#
#Upper bound latent space dimensions
R=5
#-----------------------------------------------------------------------------#
#Upper bound number of basis functions in each dimension (columns of G)
K=5

###############################################################################
###############################################################################
#DEFINE HYPERPARAMETERS SETTINGS###############################################
###############################################################################
###############################################################################

#-----------------------------------------------------------------------------#
# hyperparamters for Tau, the column-shrinkage variance of G
# \Tau_k ~ Ga( a*q^3(k-1), q^2(k-1))
a=2
q=2
#-----------------------------------------------------------------------------#
# W_{kr}^(i) (y_i) ~ GP(0,c_w)
# c_w(y_1,y_1) = exp( -kappa_w*(y_1-y_2)^2)
###------ kappa_w ------###
kappa_w = 0.001
cov_W = matrix(0,nc=n,nr=n)
for (j in 1:n){
  cov_W[j,]= exp( - kappa_w*(Y_s[j] - Y_s)^2 )
}
inv_cov_W = solve(cov_W + diag(1e-6,n)) # add a small random term to c_w to make it invertible
#-----------------------------------------------------------------------------#
# precision for X
phi_x = 1 # var_x = 1/phi_x
#-----------------------------------------------------------------------------#
# mean and variance of Z
mu_z = 0
phi_z = 0.1 # var_z = 1/phi_z = 10

################################################################################
################################################################################
#CREATE ALLOCATION MATRICES#####################################################
################################################################################
################################################################################

#-----------------------------------------------------------------------------#
# Basic log odds Z in matrix form
Z_samples = array(0, c(V,V,N_s))
#-----------------------------------------------------------------------------#
# X[v,r,i,]=X_{v,r}^{(i)}
X_samples <- array(0, c(V,R,N,N_s))
#-----------------------------------------------------------------------------#
# Shared behavior G
G_samples = array(0, c(V,K,N_s))
#-----------------------------------------------------------------------------#
# W with W[i,k,r]= W_{k,r}^{(i)}
W_samples = array(0, c(n,K,R,N_s))
#-----------------------------------------------------------------------------#
# Tau, the vector storing the column-shrinkage variance of G
Tau = numeric(K)
#-----------------------------------------------------------------------------#
# Other useful matrices for implementing the Gibbs
# kron_c_w = kronecker(diag(K),cov_W)
kron_inv_c_w = kronecker(diag(K), inv_cov_W)
W_td = array(0,c(N,K,R))
XX = array(0,c(V,V,N))
hat_X = array(0,c(V,R,n))

################################################################################
################################################################################
##INITIALIZE QUANTITIES#########################################################
################################################################################
################################################################################
set.seed(50)

#-----------------------------------------------------------------------------#
G = matrix(rnorm(V*K),nr=V, nc=K)
Tau = rep(1,K)
#-----------------------------------------------------------------------------#
# W[,k,r] ~ N(0,c_w)
W = array(0, c(n,K,R))
W = apply(W,c(2,3),function(x){ mvrnorm(1, mu=numeric(n),Sigma= cov_W) })
#-----------------------------------------------------------------------------#
X = array(rnorm(V*R*N,mean=0,sd=1/sqrt(phi_x)),c(V,R,N))
#-----------------------------------------------------------------------------#
# Z initialized using sample log-odds
Mean_A = apply(A,c(1,2),sum)/N
Mean_A[which(Mean_A==1,arr.ind=T)] = 0.99
Mean_A[which(Mean_A==0,arr.ind=T)] = 0.01
logit_Mean_A = log(Mean_A/(1-Mean_A))
Z = matrix(0,nr=V,nc=V)
lowerTriangle(Z) = rnorm(0.5*V*(V-1),mean=lowerTriangle(logit_Mean_A),sd=1)
Z = Z + t(Z)

################################################################################
################################################################################
#ALGORITHM Gibbs Sampler #######################################################
################################################################################
################################################################################
ptm=proc.time() 

for (t in 2:burnin){
  ################################################################################
  ################################################################################
  #SAMPLE AUGMENTED POLYA-GAMMA DATA (STORED IN MATRIX OMEGA)
  Omega<-array(0,c(V,V,N))
  
  for (i in 1:N){
    Omega_temp<-rpg.devroye(0.5*V*(V-1), n=1, lowerTriangle(Z + X[,,i] %*% t(X[,,i])))  
    lowerTriangle(Omega[,,i])<-Omega_temp
    Omega[,,i]<-Omega[,,i]+t(Omega[,,i])      
  } 
  
  ################################################################################
  ################################################################################
  #SAMPLE X 
  
  # expand W
  # W_td = array(0,c(N,K,R))
  for (i in 1:N){
    W_td[i,,] = W[index_y[i],,]  
  }
  
  for (i in 1:N){
    
    #-----------------------------------------------------------------------------#
    #Sample first row of X given the others
    
    Sigma_v<-solve( phi_x * diag(R) + t(X[-1,,i]) %*% diag(Omega[-1,1,i],V-1,V-1) %*% X[-1,,i] )
    mu_v<-Sigma_v %*% ( phi_x * t(W_td[i,,]) %*% G[1,] + t(X[-1,,i]) %*% (A[-1,1,i] - 0.5 - Omega[-1,1,i] * Z[-1,1]) )  
    X[1,,i]<-mvrnorm(1,mu_v,Sigma_v)
    
    #-----------------------------------------------------------------------------#
    #Sample other rows of X given the others
    for (v in 2:V){
      Sigma_v <- solve( phi_x * diag(R) + t(X[-v,,i]) %*% diag(Omega[-v,v,i],V-1,V-1) %*% X[-v,,i])
      mu_v<- Sigma_v %*% ( phi_x * t(W_td[i,,]) %*% G[v,] + t(X[-v,,i]) %*% (A[-v,v,i] - 0.5 - Omega[-v,v,i]*Z[-v,v]) )  
      X[v,,i]<-mvrnorm(1,mu_v,Sigma_v)
    }
    
  }
  
  ################################################################################
  ################################################################################
  # SAMPLE Z
  XX = apply(X,3,function(x){x %*% t(x)})
  XX = array(XX,c(V,V,N))
  Var_Z = 1/( lowerTriangle( apply(Omega,c(1,2),sum) ) + phi_z )
  mean_Z = Var_Z * ( lowerTriangle( apply(A-Omega*XX,c(1,2),sum) ) - N/2 + phi_z*mu_z)
  
  Z = matrix(0,nr=V,nc=V)
  lowerTriangle(Z) = rnorm(0.5*V*(V-1), mean = mean_Z, sd=sqrt(Var_Z))
  Z = Z + t(Z)
  
  ################################################################################
  ################################################################################
  # SAMPLE G
  Sigma_g_temp = apply(W_td, 1, function(x){x %*% t(x)})
  Sigma_g_temp = apply(Sigma_g_temp,1,sum)
  Sigma_g_temp = matrix(Sigma_g_temp, nr=K, nc=K)
  for (u in 1:V){
    Sigma_g = solve( diag(Tau,K,K) + phi_x * Sigma_g_temp )
    mu_g = 0
    for (i in 1:N){
      mu_g = mu_g + W_td[i,,] %*% X[u,,i]
    }
    mu_g = phi_x * mu_g
    mu_g = Sigma_g %*% mu_g
    G[u,] = mvrnorm(1,mu_g,Sigma_g)
  }
  
  ################################################################################
  ################################################################################
  # SAMPLE Tau
  
  for (k in 1:K){
    a_k = a*q^(3*(k-1))
    b_k = q^(2*(k-1))
    shape_tau = a_k + V/2
    rate_tau = b_k + 0.5 * t(G[,k]) %*% G[,k]
    Tau[k] = rgamma(1,shape=shape_tau, rate=rate_tau)
  }
  
  ################################################################################
  ################################################################################
  # SAMPLE W
  
  kron_GG = phi_x * kronecker( (t(G) %*% G), D)
  Sigma_w = solve(kron_inv_c_w + kron_GG) 
  
  # aggregated X
  # hat_X = array(0,c(V,R,n))
  for (j in 1:n){
    hat_X[,,j] = apply(X[,,map[j,1:c_num[j]] ],c(1,2),sum)  
  }
  hat_X_mat = apply(hat_X,2,function(x){as.vector(t(x))}) # nV x R
  
  for (r in 1:R){
    mu_w = Sigma_w %*%  (phi_x * kronecker(t(G[,]),diag(n)) ) %*% hat_X_mat[,r]
    W[,,r]= matrix(mvrnorm(1,mu_w,Sigma_w), nr=n, nc=K)   
  }
  
}

for (t in 1:N_s){
  for (s in 1:thin){
    ################################################################################
    ################################################################################
    #SAMPLE AUGMENTED POLYA-GAMMA DATA (STORED IN MATRIX OMEGA)
    Omega<-array(0,c(V,V,N))
    
    for (i in 1:N){
      Omega_temp<-rpg.devroye(0.5*V*(V-1), n=1, lowerTriangle(Z + X[,,i] %*% t(X[,,i])))  
      lowerTriangle(Omega[,,i])<-Omega_temp
      Omega[,,i]<-Omega[,,i]+t(Omega[,,i])      
    } 
    
    ################################################################################
    ################################################################################
    #SAMPLE X 
    
    # expand W
    # W_td = array(0,c(N,K,R))
    for (i in 1:N){
      W_td[i,,] = W[index_y[i],,]  
    }
    
    for (i in 1:N){
      
      #-----------------------------------------------------------------------------#
      #Sample first row of X given the others
      
      Sigma_v<-solve( phi_x * diag(R) + t(X[-1,,i]) %*% diag(Omega[-1,1,i],V-1,V-1) %*% X[-1,,i] )
      mu_v<-Sigma_v %*% ( phi_x * t(W_td[i,,]) %*% G[1,] + t(X[-1,,i]) %*% (A[-1,1,i] - 0.5 - Omega[-1,1,i] * Z[-1,1]) )  
      X[1,,i]<-mvrnorm(1,mu_v,Sigma_v)
      
      #-----------------------------------------------------------------------------#
      #Sample other rows of bar_X given the others
      for (v in 2:V){
        Sigma_v <- solve( phi_x * diag(R) + t(X[-v,,i]) %*% diag(Omega[-v,v,i],V-1,V-1) %*% X[-v,,i])
        mu_v<- Sigma_v %*% ( phi_x * t(W_td[i,,]) %*% G[v,] + t(X[-v,,i]) %*% (A[-v,v,i] - 0.5 - Omega[-v,v,i]*Z[-v,v]) )  
        X[v,,i]<-mvrnorm(1,mu_v,Sigma_v)
      }
      
    }
    
    ################################################################################
    ################################################################################
    # SAMPLE Z
    XX = apply(X,3,function(x){x %*% t(x)})
    XX = array(XX,c(V,V,N))
    Var_Z = 1/( lowerTriangle( apply(Omega,c(1,2),sum) ) + phi_z )
    mean_Z = Var_Z * ( lowerTriangle( apply(A-Omega*XX,c(1,2),sum) ) - N/2 + phi_z*mu_z)
    
    Z = matrix(0,nr=V,nc=V)
    lowerTriangle(Z) = rnorm(0.5*V*(V-1), mean = mean_Z, sd=sqrt(Var_Z))
    Z = Z + t(Z)
    
    ################################################################################
    ################################################################################
    # SAMPLE G
    Sigma_g_temp = apply(W_td, 1, function(x){x %*% t(x)})
    Sigma_g_temp = apply(Sigma_g_temp,1,sum)
    Sigma_g_temp = matrix(Sigma_g_temp, nr=K, nc=K)
    for (u in 1:V){
      Sigma_g = solve( diag(Tau,K,K) + phi_x * Sigma_g_temp )
      mu_g = 0
      for (i in 1:N){
        mu_g = mu_g + W_td[i,,] %*% X[u,,i]
      }
      mu_g = phi_x * mu_g
      mu_g = Sigma_g %*% mu_g
      G[u,] = mvrnorm(1,mu_g,Sigma_g)
    }
    
    ################################################################################
    ################################################################################
    # SAMPLE Tau
    
    for (k in 1:K){
      a_k = a*q^(3*(k-1))
      b_k = q^(2*(k-1))
      shape_tau = a_k + V/2
      rate_tau = b_k + 0.5 * t(G[,k]) %*% G[,k]
      Tau[k] = rgamma(1,shape=shape_tau, rate=rate_tau)
    }
    
    ################################################################################
    ################################################################################
    # SAMPLE W
    
    kron_GG = phi_x * kronecker( (t(G) %*% G), D)
    Sigma_w = solve(kron_inv_c_w + kron_GG) 
    
    # aggregated X
    # hat_X = array(0,c(V,R,n))
    for (j in 1:n){
      hat_X[,,j] = apply(X[,,map[j,1:c_num[j]] ],c(1,2),sum)  
    }
    hat_X_mat = apply(hat_X,2,function(x){as.vector(t(x))}) # nV x R
    
    for (r in 1:R){
      mu_w = Sigma_w %*%  (phi_x * kronecker(t(G[,]),diag(n)) ) %*% hat_X_mat[,r]
      W[,,r]= matrix(mvrnorm(1,mu_w,Sigma_w), nr=n, nc=K)   
    }
  }
    
  X_samples[,,,t] = X
  Z_samples[,,t] = Z
  G_samples[,,t] = G
  W_samples[,,,t] = W
  print(t)
 
}
proc.time()-ptm

######################################################
# individual probability #############################
######################################################
# store posterior samples of edge probabilities for each subject
Pi_sub = array(0,c(V,V,N,N_s)) 
mean_Pi_sub = array(0, c(V,V,N))# poseterior mean of individual edge prob

for (i in 1:N){
  XX = apply(X_samples[,,i,], 3, function(x){x %*% t(x)})
  XX = array(XX, c(V,V,N_s))
  Pi_sub[,,i,] = 1/(1 + exp(-Z_samples -XX) )
  mean_Pi_sub[,,i] = apply(Pi_sub[,,i,],c(1,2),sum)/N_s
}

#--------- ROC curve ---------#
pmean_vec=as.vector(apply(mean_Pi_sub,3,lowerTriangle))
A_vec = as.vector(apply(A,3,lowerTriangle))
pred = prediction(pmean_vec, A_vec)
perf.auc = performance(pred,"auc")
perf.auc@y.values
AUC=round(as.numeric(perf.auc@y.values),digits=3)
perf = performance(pred, "tpr","fpr")

########################################################
# expected edge probability for each value of y ########
########################################################
# number of samples for marginalizing out X
N_m = 1000

ptm=proc.time()

# store posterior samples for conditional expection of edge prob given each value of y
Pi_y = unlist( mclapply(1:N_s,function(t){
  mu_X_temp = apply(W_samples[,,,t],1,function(x){G_samples[,,t] %*% x})
  mu_X = array(mu_X_temp,c(V,R,n))
  mean_Pi_t = array(0,c(V,V,n))
  for(j in 1:n){
    X_sample = apply(mu_X[,,j],c(1,2),function(x){rnorm(N_m, mean=x, sd=sqrt(1/phi_x))}) # N_m x V x R
    mean_Pi_temp = apply(X_sample,1,function(x){1/(1+exp(-Z_samples[,,t] - x %*% t(x)))})
    mean_Pi_temp = array(mean_Pi_temp,c(V,V,N_m))
    mean_Pi_t[,,j] = apply(mean_Pi_temp,c(1,2),sum)/N_m
  }
  return(mean_Pi_t)
}) )

Pi_y = array(Pi_y,c(V,V,n,N_s))
proc.time()-ptm 

mean_Pi_y = apply(Pi_y,c(1,2,3),sum)/N_s  # posterior mean of conditional expection


############################################################################################
### check average connection probability for different anatomical areas in the brain #######
############################################################################################
#---- whole brain/network density ----#
pmean_network = apply(mean_Pi_y,3,lowerTriangle)
pmean_network = apply(pmean_network,2,mean)
## quantile
mean_network_q = apply(Pi_y,c(3,4),function(x){mean(lowerTriangle(x))})
mean_network_q = apply(mean_network_q,1,function(x){HPDinterval(as.mcmc(x),prob=0.95)})
#-------- ggplot ----------------#
net_den_data = data.frame(Y_s, pmean_network, lb=mean_network_q[1,], ub=mean_network_q[2,])
ggplot(data=net_den_data, aes(Y_s, pmean_network))+
  geom_line(col='red')+ xlab("IQ") + ylab("") +
  geom_ribbon(aes(ymin=lb,ymax=ub),fill="lightcoral",alpha=0.3) + theme_bw() +
  ggtitle('network density')

#-----------left hemisphere------------#
pmean_network = apply(mean_Pi_y[1:34,1:34,],3,lowerTriangle)
pmean_network = apply(pmean_network,2,mean)

mean_network_q = apply(Pi_y[1:34,1:34,,],c(3,4),function(x){mean(lowerTriangle(x))})
mean_network_q = apply(mean_network_q,1,function(x){HPDinterval(as.mcmc(x),prob=0.95)})
#-------- ggplot ----------------#
net_den_data = data.frame(Y_s, pmean_network, lb=mean_network_q[1,], ub=mean_network_q[2,])
ggplot(data=net_den_data, aes(Y_s, pmean_network))+
  geom_line(col='darkblue')+ xlab("IQ") + ylab("") +
  geom_ribbon(aes(ymin=lb,ymax=ub),fill="skyblue",alpha=0.3) + theme_bw() +
  ggtitle('left hemisphere network density')

#-----------right hemisphere------------#
pmean_network = apply(mean_Pi_y[35:68,35:68,],3,lowerTriangle)
pmean_network = apply(pmean_network,2,mean)

mean_network_q = apply(Pi_y[35:68,35:68,,],c(3,4),function(x){mean(lowerTriangle(x))})
mean_network_q = apply(mean_network_q,1,function(x){HPDinterval(as.mcmc(x),prob=0.95)})
#-------- ggplot ----------------#
net_den_data = data.frame(Y_s, pmean_network, lb=mean_network_q[1,], ub=mean_network_q[2,])
ggplot(data=net_den_data, aes(Y_s, pmean_network))+
  geom_line(col='darkblue')+ xlab("IQ") + ylab("") +
  geom_ribbon(aes(ymin=lb,ymax=ub),fill="steelblue",alpha=0.3) + theme_bw() +
  ggtitle('right hemisphere network density')

#-----------cross hemisphere------------#
pmean_network = apply(mean_Pi_y[1:34,35:68,],3,as.vector)
pmean_network = apply(pmean_network,2,mean)
#pdf("cross brain network density.pdf",width=5,height=5)

mean_network_q = apply(Pi_y[1:34,35:68,,],c(3,4),function(x){mean(as.vector(x))})
mean_network_q = apply(mean_network_q,1,function(x){HPDinterval(as.mcmc(x),prob=0.95)})
#-------- ggplot ----------------#
net_den_data = data.frame(Y_s, pmean_network, lb=mean_network_q[1,], ub=mean_network_q[2,])
ggplot(data=net_den_data, aes(Y_s, pmean_network))+
  geom_line(col='darkgreen')+ xlab("IQ") + ylab("") +
  geom_ribbon(aes(ymin=lb,ymax=ub),fill="springgreen",alpha=0.3) + theme_bw() +
  ggtitle('inter-hemisphere network density')
#dev.off()

#--------------check each lobe----------------#
frontal_lobe = c(3,37,32,66,12,46,14,48,17,51,18,52,19,53,20,54,24,58,27,61,28,62)
limbic_cortex = c(2,36,4,38,10,44,16,50,23,57,26,60)
occipital_lobe = c(5,39,11,45,13,47,21,55)
parietal_lobe = c(8,42,22,56,25,59,29,63,31,65)
temporal_lobe = c(1,35,6,40,7,41,9,43,15,49,30,64,33,67,34,68)

left_frontal_lobe = c(3,32,12,14,17,18,19,20,24,27,28)
right_frontal_lobe = c(37,66,46,48,51,52,53,54,58,61,62)

left_limbic_cortex = c(2,4,10,16,23,26)
right_limbic_cortex = c(36,38,44,50,57,60)

#-frontal_lobe-#
pmean_network = apply(mean_Pi_y[frontal_lobe,frontal_lobe,],3,lowerTriangle)
pmean_network = apply(pmean_network,2,mean)
#pdf("frontal lobe density.pdf",width=5,height=5)
mean_network_q = apply(Pi_y[frontal_lobe,frontal_lobe,,],c(3,4),function(x){mean(lowerTriangle(x))})
mean_network_q = apply(mean_network_q,1,function(x){HPDinterval(as.mcmc(x),prob=0.95)})

#-------- ggplot ----------------#
net_den_data = data.frame(Y_s, pmean_network, lb=mean_network_q[1,], ub=mean_network_q[2,])
ggplot(data=net_den_data, aes(Y_s, pmean_network))+
  geom_line(col='darkgreen')+ xlab("IQ") + ylab("") +
  geom_ribbon(aes(ymin=lb,ymax=ub),fill="springgreen",alpha=0.2) + theme_bw() +
  ggtitle('frontal lobe density')
#dev.off()

#-- cross frontal lobe --#
# check the average connection probability of cross-frontal-lobe links
pmean_network = apply(mean_Pi_y[left_frontal_lobe,right_frontal_lobe,],3,as.vector)
pmean_network = apply(pmean_network,2,mean)
#pdf("cross frontal lobe density.pdf",width=5,height=5)

mean_network_q = apply(Pi_y[left_frontal_lobe,right_frontal_lobe,,],c(3,4),function(x){mean(as.vector(x))})
mean_network_q = apply(mean_network_q,1,function(x){HPDinterval(as.mcmc(x),prob=0.95)})
#-------- ggplot ----------------#
net_den_data = data.frame(Y_s, pmean_network, lb=mean_network_q[1,], ub=mean_network_q[2,])
ggplot(data=net_den_data, aes(Y_s, pmean_network))+
  geom_line(col='coral4')+ xlab("IQ") + ylab("") +
  geom_ribbon(aes(ymin=lb,ymax=ub),fill="coral",alpha=0.2) + theme_bw() +
  ggtitle('inter frontal lobe density')
#dev.off()

#-Limbic Cortex-#
pmean_network = apply(mean_Pi_y[limbic_cortex,limbic_cortex,],3,lowerTriangle)
pmean_network = apply(pmean_network,2,mean)

mean_network_q = apply(Pi_y[limbic_cortex,limbic_cortex,,],c(3,4),function(x){mean(lowerTriangle(x))})
mean_network_q = apply(mean_network_q,1,function(x){HPDinterval(as.mcmc(x),prob=0.95)})

#-------- ggplot ----------------#
net_den_data = data.frame(Y_s, pmean_network, lb=mean_network_q[1,], ub=mean_network_q[2,])
ggplot(data=net_den_data, aes(Y_s, pmean_network))+
  geom_line(col='darkgreen')+ xlab("IQ") + ylab("") +
  geom_ribbon(aes(ymin=lb,ymax=ub),fill="springgreen",alpha=0.3) + theme_bw() +
  ggtitle('limbic cortex density')

#-- cross limbic cortex --#
# check the average connection probability of cross-limbic-cortex connections
pmean_network = apply(mean_Pi_y[left_limbic_cortex,right_limbic_cortex,],3,as.vector)
pmean_network = apply(pmean_network,2,mean)

mean_network_q = apply(Pi_y[left_limbic_cortex,right_limbic_cortex,,],c(3,4),function(x){mean(as.vector(x))})
mean_network_q = apply(mean_network_q,1,function(x){HPDinterval(as.mcmc(x),prob=0.95)})
#-------- ggplot ----------------#
net_den_data = data.frame(Y_s, pmean_network, lb=mean_network_q[1,], ub=mean_network_q[2,])
ggplot(data=net_den_data, aes(Y_s, pmean_network))+
  geom_line(col='darkgreen')+ xlab("IQ") + ylab("") +
  geom_ribbon(aes(ymin=lb,ymax=ub),fill="springgreen",alpha=0.3) + theme_bw() +
  ggtitle('inter limbic cortex density')

#-FL <-> LC-#
pmean_network = apply(mean_Pi_y[frontal_lobe,limbic_cortex,],3,as.vector)
pmean_network = apply(pmean_network,2,mean)
# pdf("FL_LC density.pdf",width=5,height=5)

mean_network_q = apply(Pi_y[frontal_lobe,limbic_cortex,,],c(3,4),function(x){mean(as.vector(x))})
mean_network_q = apply(mean_network_q,1,function(x){HPDinterval(as.mcmc(x),prob=0.95)})

#-------- ggplot ----------------#
net_den_data = data.frame(Y_s, pmean_network, lb=mean_network_q[1,], ub=mean_network_q[2,])
ggplot(data=net_den_data, aes(Y_s, pmean_network))+
  geom_line(col='darkgreen')+ xlab("IQ") + ylab("") +
  geom_ribbon(aes(ymin=lb,ymax=ub),fill="springgreen",alpha=0.2) + theme_bw() +
  ggtitle('FL_LC density')
# dev.off()

