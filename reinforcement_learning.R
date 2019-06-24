###### Generate transition matrices ##############

generate_transition_matrix <- function(nb_states, seed){
  set.seed(seed)
  
  transition_matrix <- matrix(runif(nb_states^2, min=0,max=1), ncol=nb_states)
  transition_matrix[nb_states,] <- c(rep(0,nb_states-1),1)
  
  
  for (i in 1:(nrow(transition_matrix)-1)){
    for (j in i:ncol(transition_matrix)){
      transition_matrix[i,j] <- j
      transition_matrix[i,1+floor(i/4)] <- 3*j
    }
  }
  
  for (ii in 1:(nrow(transition_matrix)-1)){
    transition_matrix[ii,] <- transition_matrix[ii,]/sum(transition_matrix[ii,])
  }
  
  
  return(transition_matrix)
}

##### Example #####

nb_states <- 25
tolerance <- 0.02
cost_vector <- c(rep(1,nb_states-1),0)
transition_matrix_list <- list()

for (j in 1:5){
  transition_matrix_list[[j]] <- generate_transition_matrix(nb_states=nb_states,seed=j)
}

TML <- transition_matrix_list


# This function outputs a vector containing in entry i the position of the minimum entry of a row i of matrix

minperrow <- function(matrix){
  nrow <- nrow(matrix)
  result <- numeric(nrow)
  for (k in 1:nrow){
    result[k] <- which.min(matrix[k,])
  }
  return(result)
}



q.learn <- function(N, TML1, TML2, TML3, TML4, TML5, costVector, target_state, nb_dice=5) {
  
  #costVector <- c(rep(1,nb_states-1),0)
  Q <- matrix(rep(0,nb_states*5), nrow=nb_states)
  
  stage <- 0
  
  for (i in 1:N){
    alpha <- 1/i
    # for each episode, choose an initial state at random
    cs <- sample(nrow(Q),1)
    ### iterate until we get to the tgt.state
    stage <- stage + 1
    if(stage%%500==0){
      print(paste("episode",stage)) # to keep track of how advanced the loop is
    }
    
    while(1){
      ## choose next state from possible actions at current state
      
      r <- rmultinom(1,1,rep(1/nb_states,nb_states))
      if (r[1,1]==1){ # choose next state from TML1 transition
        die <- 1
        p <- TML1[cs,]
      } else if (r[2,1]==1){ # choose next state from TML2 transition
        die <- 2
        p <- TML2[cs,]
      } else if (r[3,1]==1){ # choose next state from TML3 transition
        die <- 3
        p <-TML3[cs,]
      } else if (r[4,1]==1){  # choose next state from TML4 transition
        die <- 4
        p <-TML4[cs,]
      }
        
        else { # choose next state from TML5 transition
        die <- 5
        p <- TML5[cs,]
      }
      
      ns <- which(rmultinom(1,1,p)==1) 
      
      
      Q[cs,die] <- Q[cs,die] + alpha*(costVector[cs] + min(Q[ns,]) - Q[cs,die])
      
      ## break out of while loop if target state is reached
      ## otherwise, set next.state as current.state and repeat      
      if (ns == target_state){break}
      cs <- ns
    }
  }
  policy <- minperrow(Q)
  return(policy)
}

### Need a function that determines the average number of steps to absorption for a Markov chain
### P: matrix with one absorbing state (last row/column)
### s: index from which we start

expected_number_Of_steps_to_absorption <- function(s,P){
  m <- dim(P)[1] - 1 
  
  ones <- rep(1,m) 
  Q <- matrix(rep(0,m**2), ncol=m) #initialise matrix of 0s of dimensions mxm
  for (i in 1:m){
    for (j in 1:m){
      Q[i,j] <- P[i,j] #fill Q with entries of initial transition matrix : Q = P disregarding last column and last row
    }
  }
  N <- matrix(rep(0,m**2), ncol = m) 
  I_t <- matrix(rep(0,m**2), ncol = m)
  diag(I_t) <- 1 # I_t is now the identity matrix
  
  N <- solve(I_t-Q) # N = (I-Q)^-1
  t <- N%*%ones #N = Ne
  t <- c(t,0) 
  return(t[s]) #return number of steps to absorbtion from the starting square
}


opt_policy <- q.learn(5000, TML[[1]], TML[[2]], TML[[3]], TML[[4]], TML[[5]], costVector= cost_vector, target_state=nb_states, nb_dice=5)

#### Markov chain induced by the optimal policy
#### opt_TML is the resulting one-step transition matrix

opt_TML <- matrix(rep(NA,nb_states^2), nrow=nb_states)

for (l in 1:nb_states){
  opt_TML[l,] <- TML[[opt_policy[l]]][l,]
}

# Expected number of steps from each state given that one follows the optimal policy
expected_number_Of_steps_to_absorption(1:(nb_states-1),opt_TML)
