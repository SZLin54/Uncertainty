#############JAGs code for the JRT-DINA model############

JRTDINA <-  function(){
  for(n in 1:N){
    for(i in 1:I){
      for (k in 1:K){w[n, i, k] <- pow(alpha[n, k], Q[i, k])}
      logit(p[n, i]) <- lamda0[i] + lamda1[i] *prod(w[n,i, ])
      Y[n, i] ~ dbern(p[n, i])
      logT[n, i] ~ dnorm(zeta[i] - tau[n], den_epsilon[i])
      }}
  for(n in 1:N){
    for(k in 1:K){
      logit(prob.a[n, k]) <- xi[k] * theta[n] - beta[k]
      alpha[n, k] ~ dbern(prob.a[n, k])}
  }
  
  for (n in 1:N) {
    person_parameter[n, 1:2] ~ dmnorm(person_mu[1:2], person_den[1:2, 1:2])
    theta[n] <- person_parameter[n, 1]
    tau[n] <- person_parameter[n, 2]
  }
  
  for(k in 1:K){ 
    beta[k] ~ dnorm(0, 0.25) 
    xi[k] ~ dnorm(0, 0.25) %_%T(0,)}
  for(i in 1:I) {
    item_parameter[i,1:3] ~ dmnorm(item_mu[1:3], item_den[1:3, 1:3])
    lamda0[i] <- item_parameter[i,1]
    lamda1[i] <- item_parameter[i,2]
    zeta[i] <- item_parameter[i,3]
    logit(g[i]) <- lamda0[i]
    logit(ns[i]) <- lamda0[i] + lamda1[i]
    s[i] <- 1 - ns[i]
    den_epsilon[i] ~ dgamma(1, 1)
  }
  
  item_mu[1] ~ dnorm(-2.197,0.5)
  item_mu[2] ~ dnorm(4.394,0.5) %_%T(0,)
  item_mu[3] ~ dnorm(3,0.5)
  R[1, 1] <- 1
  R[2, 2] <- 1
  R[3, 3] <- 1
  R[1, 2] <- 0
  R[1, 3] <- 0
  R[2, 1] <- 0
  R[2, 3] <- 0
  R[3, 1] <- 0
  R[3, 2] <- 0
  item_den[1:3, 1:3] ~ dwish(R[1:3, 1:3], 3)
  Sigma_item[1:3, 1:3] <- inverse(item_den[1:3, 1:3])
  
  person_mu[1] <- 0
  person_mu[2] <- 0
  L_theta[1, 1] <- 1
  L_theta[2, 2] ~ dgamma(1, 1)
  L_theta[2, 1] ~ dnorm(0,1)
  L_theta[1, 2] <- 0
  Sigma_theta <- L_theta %*% t(L_theta)
  person_den[1:2, 1:2] <- inverse(Sigma_theta[1:2, 1:2])
  
  
}

#############JAGs code for the uJRT-DINA model############

UJRTDINA <-  function(){
  for(n in 1:N){
    for(i in 1:I){
      for (k in 1:K){w[n, i, k] <- pow(alpha[n, k], Q[i, k])}
      logit(p[n, i]) <- exp(ga[n])*(lamda0[i] + lamda1[i] *prod(w[n,i, ]))
      Y[n, i] ~ dbern(p[n, i])
      logT[n, i] ~ dnorm(zeta[i] - tau[n], den_epsilon[i])
    }}
  
  for(n in 1:N){
    for(k in 1:K){
      logit(prob.a[n, k]) <- xi[k] * theta[n] - beta[k]
      alpha[n, k] ~ dbern(prob.a[n, k])}
    person_parameter[n, 1:3] ~ dmnorm(person_mu[1:3], person_den[1:3, 1:3])
    theta[n] <- person_parameter[n, 1]
    tau[n] <- person_parameter[n, 2]
    ga[n] <- person_parameter[n, 3]
  }
  
  for(k in 1:K){ 
    beta[k] ~ dnorm(0, 0.25) 
    xi[k] ~ dnorm(0, 0.25) %_%T(0,)}
  for(i in 1:I) {
    item_parameter[i,1:3] ~ dmnorm(item_mu[1:3], item_den[1:3, 1:3])
    lamda0[i] <- item_parameter[i,1]
    lamda1[i] <- item_parameter[i,2]
    zeta[i] <- item_parameter[i,3]
    logit(g[i]) <- lamda0[i]
    logit(ns[i]) <- lamda0[i] + lamda1[i]
    s[i] <- 1 - ns[i]
    den_epsilon[i] ~ dgamma(1, 1)
    Sigma_epsilon[i] <- 1/den_epsilon[i]
  }
  
  item_mu[1] ~ dnorm(-2.197,0.5)
  item_mu[2] ~ dnorm(4.394,0.5) %_%T(0,)
  item_mu[3] ~ dnorm(3,0.5)
  R[1, 1] <- 1
  R[2, 2] <- 1
  R[3, 3] <- 1
  R[1, 2] <- 0
  R[1, 3] <- 0
  R[2, 1] <- 0
  R[2, 3] <- 0
  R[3, 1] <- 0
  R[3, 2] <- 0
  item_den[1:3, 1:3] ~ dwish(R[1:3, 1:3], 3)
  Sigma_item[1:3, 1:3] <- inverse(item_den[1:3, 1:3])
  
  person_mu[1] <- 0
  person_mu[2] <- 0
  person_mu[3] <- 0
  L_theta[1, 1] <- 1
  L_theta[1, 2] <- 0
  L_theta[1, 3] <- 0
  L_theta[2, 1] ~ dnorm(0,1)
  L_theta[2, 2] ~ dgamma(1, 1)
  L_theta[2, 3] <- 0
  L_theta[3, 1] ~ dnorm(0,1)
  L_theta[3, 2] ~ dnorm(0,1)
  L_theta[3, 3] ~ dgamma(1, 1)
  Sigma_theta <- L_theta %*% t(L_theta)
  person_den[1:3, 1:3] <- inverse(Sigma_theta[1:3, 1:3])
  
}
