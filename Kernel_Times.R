kernel_exp <- function(no_offspring, parameters){
  return(rexp(no_offspring, rate = parameters$delta))
}

kernel_ray <- function(no_offspring, parameters){
  return(rweibull(no_offspring, shape = 2, scale = (sqrt(2)/sqrt(parameters$delta))))
}
