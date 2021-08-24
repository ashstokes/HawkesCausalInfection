constant_mu <- function(parameters, T_max){
  k <- rpois(1, parameters$A*T_max)
  time <- runif(k, 0, T_max)
  return(sort(time))
}

sinusoidal_mu <- function(parameters, T_max){
  return(hawkes_simulation(c(0), T_max = 100, N_max = Inf,
                           kernel=exp_kernel, parameters= parameters,
                           mu_fn = mu_sinusoidal_linear, mu_fn_diff = mu_diff_sinusoidal_linear, mu_t_max = mu_t_max,
                           imported = T, print_level = 1))
}

