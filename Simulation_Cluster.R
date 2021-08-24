Cluster_Simulation <- function(events, T_max, parameters, N_max, mu_fn, phi_fn, branching_factor){
  
  #generate times of each background event
  time_ex_events <- mu_fn(parameters, T_max)
  
  #collect the times of generation 0
  events <- c(events, time_ex_events)
  
  #allocating family numbers to each event of generation 0
  child_ID <- 1:length(events)
  
  
  #indexing of generations
  gen_number <- 1
  gen_size <- vector(mode="numeric")
  gen_size[gen_number] <- length(events)
  
  previous_events <- 0
  
  parents <- rep(NA, length(child_ID))
  children <- child_ID
  gen <- rep(gen_number, length(child_ID))
  
  while(T){
    
    offspring_time <- c()
    no_offspring <- c()
    
    for(i in 1:gen_size[gen_number]){
      no_offspring[i] <- rpois(1, branching_factor(parameters))
      offspring_time <- c(offspring_time, phi_fn(no_offspring[i], parameters) + events[previous_events + i] + parameters$delay)
    }
    
    #collecting family indexes and removing any offspring outside time period
    parent_ID <- rep(child_ID, no_offspring)
    parent_ID <- subset(parent_ID, offspring_time < T_max)
    child_ID <- seq(sum(gen_size)+1, sum(gen_size)+length(parent_ID))
    offspring_time <- subset(offspring_time, offspring_time < T_max)
    
    #stopping code if reaches extinction or maximum number of events allowed
    if( length(offspring_time) == 0 | length(events) >= N_max){
      break
      
    } else {
      
      events <- c(events, offspring_time)
      
      gen_number <- gen_number + 1
      gen_size[gen_number] <- length(offspring_time)
      
      previous_events <- sum(gen_size[-gen_number])
      
      children <- c(children, child_ID)
      parents <- c(parents, parent_ID)
      gen <- c(gen, rep(gen_number, length(child_ID)))
      
      
    }
  }
  output <- matrix(NA, nrow = length(events), ncol=4)
  output[,1] <- events
  output[,3] <- children
  output[,2] <- parents
  output[,4] <- gen
  output <- output[order(output[,1]),]
  return(output)
}




