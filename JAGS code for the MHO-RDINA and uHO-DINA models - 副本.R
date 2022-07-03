####################JAGs code for the MHO-RDINA model########
MHORDINA <-  function(){
  for(n in 1:N){
    for(i in 1:I){
      for (k in 1:K){w[n, i, k] <- pow(alpha[n, k], Q[i, k])}
      logit(p[n, i]) <- lamda0[i] + lamda1[i] *prod(w[n,i, ])
      Y[n, i] ~ dbern(p[n, i])}}
  for(n in 1:N){
    for(k in 1:K){
      logit(prob.a[n, k]) <- xi[k] * theta[n] - beta[k] 
      alpha[n, k] ~ dbern(prob.a[n, k])}
    theta[n] ~ dnorm(0,1)
  }
  for(k in 1:K){ 
    beta[k] ~ dnorm(0, 0.25) 
    xi[k] ~ dnorm(0, 0.25) %_%T(0,)}
  
  for(i in 1:I) {
    item_parameter[i,1:2] ~ dmnorm(item_mu[1:2], item_den[1:2, 1:2])
    lamda0[i] <- item_parameter[i,1]
    lamda1[i] <- item_parameter[i,2]
    logit(g[i]) <- lamda0[i]
    logit(ns[i]) <- lamda0[i] + lamda1[i]
    s[i] <- 1 - ns[i]
  }
  
  item_mu[1] ~ dnorm(-2.197,0.5)
  item_mu[2] ~ dnorm(4.394,0.5) %_%T(0,)
  
  R[1, 1] <- 1
  R[2, 2] <- 1
  R[1, 2] <- 0
  R[2, 1] <- 0
  item_den[1:2, 1:2] ~ dwish(R[1:2, 1:2], 2)
  Sigma_item[1:2, 1:2] <- inverse(item_den[1:2, 1:2])
  
  
}

####################JAGs code for the uHO-RDINA model########

UHODINA <-  function(){
  for(n in 1:N){
    for(i in 1:I){
      for (k in 1:K){w[n, i, k] <- pow(alpha[n, k], Q[i, k])}
      logit(p[n, i]) <- (exp(ga[n])) * (lamda0[i] + lamda1[i] *prod(w[n,i, ]))
      Y[n, i] ~ dbern(p[n, i])}}
  for(n in 1:N){
    for(k in 1:K){
      logit(prob.a[n, k]) <- xi[k] * theta[n] - beta[k] 
      alpha[n, k] ~ dbern(prob.a[n, k])}
    
  }
  
  for (n in 1:N) {
    theta[n] ~ dnorm(0,1)
    ga[n] ~ dnorm(0,0.25)
    gamma[n] <- exp(ga[n])
  }
  
  for(k in 1:K){ 
    beta[k] ~ dnorm(0, 0.25) 
    xi[k] ~ dnorm(0, 0.25) %_%T(0,)}
  
  for(i in 1:I) {
    item_parameter[i,1:2] ~ dmnorm(item_mu[1:2], item_den[1:2, 1:2])
    lamda0[i] <- item_parameter[i,1]
    lamda1[i] <- item_parameter[i,2]
    logit(g[i]) <- lamda0[i]
    logit(ns[i]) <- lamda0[i] + lamda1[i]
    s[i] <- 1 - ns[i]
  }
  
  item_mu[1] ~ dnorm(-2.197,0.5)
  item_mu[2] ~ dnorm(4.394,0.5) %_%T(0,)
  R[1, 1] <- 1
  R[2, 2] <- 1
  R[1, 2] <- 0
  R[2, 1] <- 0
  item_den[1:2, 1:2] ~ dwish(R[1:2, 1:2], 2)
  Sigma_item[1:2, 1:2] <- inverse(item_den[1:2, 1:2])
  
  person_mu[1] <- 0
  person_mu[2] <- 0
  L_theta[1, 1] <- 1
  L_theta[2, 2] ~ dgamma(1, 1)
  L_theta[2, 1] ~ dnorm(0, 1)
  L_theta[1, 2] <- 0
  Sigma_person <- L_theta %*% t(L_theta)
  person_den[1:2, 1:2] <- inverse(Sigma_person[1:2, 1:2])
  
  
}