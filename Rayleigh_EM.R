E_M_Ray <- function(starting_values, times, thresh, T_max){
  omega <- starting_values$alpha/starting_values$delta
  delta <- starting_values$delta
  mu <- starting_values$A
  shift <- ifelse(is.null(starting_values$delay),0,starting_values$delay)
  t <- sort(times)
  n <- length(t)
  k <- 0
  
  while(T){
    k <- k+1
    p_i_j <- matrix(0, n, n)
    p_i_i <- rep(0, n)
    L <- 0
    M <- matrix(0, n, n)
    N <- 0
    p_i_i[1] <- 1
    
    #expectation step
    for(i in 2:n){
      O <- total_kern_ray(omega, delta, t, x=(i-1), shift)
      p_i_i[i] <- mu/(mu + O)
      for(j in 1:(i-1)){
        p_i_j[i, j] <- ((omega*delta*(t[i]-(t[j] + shift))*exp(-(delta/2)*(t[i]-(t[j] + shift))^2))/(mu + O))*(t[i] > (t[j] + shift))
        M[i, j] <- (((t[i]-(t[j]+shift))^2)/2)*p_i_j[i,j]*(t[i] > (t[j] + shift))
      }
      L <- L + exp(-(delta/2)*(T_max-(t[i] + shift))^2)
      N <- N + ((T_max-(t[i] + shift))^2)/2*exp(-(delta/2)*(T_max-(t[i] + shift))^2)
    }
    
    
    #maximisation step
    mu_new <- sum(p_i_i, na.rm=TRUE)/T_max
    omega_new <- sum(p_i_j, na.rm=TRUE)/(n-L)
    delta_new <- sum(p_i_j, na.rm=TRUE)/(sum(M, na.rm=TRUE)+omega*N)
    
    if( (abs(mu_new-mu) <= thresh) & (abs(omega_new-omega) <= thresh) & (abs(delta_new-delta) <= thresh) ){
      break
      
    } else {
      mu <- mu_new
      omega <- omega_new
      delta <- delta_new
    }
    
  }
  alpha <- omega*delta
  return(c(mu, alpha, delta))
}


total_kern_ray <- function(omega, delta, times, x, shift){
  current <- times[x+1]
  sum <- 0
  for(j in 1:x){
    sum <- sum + omega*delta*(current-(times[j] + shift))*(exp(-(delta/2)*(current-(times[j] + shift))^2))*(current > (times[j] + shift))
  }
  return(sum)
} 
