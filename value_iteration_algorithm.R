euclidean_distance <- function(vec1,vec2){
  sum <- 0
  for (i in 1:length(vec1)){
    sum <- sum + (vec1[i]-vec2[i])**2
  }
  result <- as.numeric(sqrt(sum))
  return(result)
}

##### Suppose we have a game that is modelled by a Markov decision process with 20 states.
##### The state 20 is the one absorbing state.
##### An action is the choice of one die one wants to roll among a collection of five dice 
##### that induce various transition one-step transition probabilities between the states.
##### The goal of the game is to reach the absorbing state (20) in as few moves as possible.
##### We use value iteration to determine the optimal policy of the game, i.e.
##### the mapping of every state to an action (the choice of a die) that minimises the expected 
##### number of steps to complete the game
##### The Ti's are the one-step transition probability matrices associated with the dice

valueIterationAlgorithm <- function(T1, T2, T3, T4, T5, costVector, tolerance){ 
 
  vector <- numeric(length(costVector))
  policy <- numeric(length(costVector))
  
  convergence <- FALSE
  
  for(t  in 1 : 1000){ # perform 1000 iterations; break as soon as there is convergence
    pVector <- vector # Will be useful to test convergence
    if (convergence==TRUE){break}
    
    for (i in 1 : length(costVector)){
      die1 <- costVector[i]+T1[i,]%*%vector
      die2 <- costVector[i]+T2[i,]%*%vector
      die3 <- costVector[i]+T3[i,]%*%vector
      die4 <- costVector[i]+T4[i,]%*%vector
      die5 <- costVector[i]+T5[i,]%*%vector
      
      vector[i] <- min(die1,die2,die3,die4,die5)
      
      die <- c(die1,die2,die3,die4,die5)
      
      for (k in 1:5){
          if (die[k]==vector[i]){
             policy[i] <- k
          }
      }
     
      
      if (i==length(costVector) && (euclidean_distance(pVector,vector)<tolerance)){
        convergence <- TRUE
      }
    }
  }
  return(policy)
}

######### Generate Markov chain with nb_states states s.t. only state n° nb_states is absorbing   #########

generate_transition_matrix <- function(nb_states, seed){
  set.seed(seed)
  
  transition_matrix <- matrix(runif(nb_states^2, min=0,max=1), ncol=nb_states)
  transition_matrix[nb_states,] <- c(rep(0,nb_states-1),1)
  
  for (q in 1:(ncol(transition_matrix)-1)){
     transition_matrix[,q] <- sqrt((nb_states-q))*transition_matrix[,q] # so that backward transition have high probability enough
     # that the dynamics of the game is not straightforward
  }
  
  for (i in 1:(nrow(transition_matrix)-1)){
     transition_matrix[i,] <- transition_matrix[i,]/sum(transition_matrix[i,])
     if (transition_matrix[i,i]==1){
         a <- runif(1,min=0,max=1)
         b <- floor(runif(1,min=1, max=nb_states))
         transition_matrix[i,i] <- a
         transition_matrix[i,b] <- 1 - a
     }
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

optimal_policy <- valueIterationAlgorithm(TML[[1]], TML[[2]], TML[[3]], TML[[4]], TML[[5]], cost_vector, tolerance)

#### Markov chain induced by the optimal policy
#### opt_TML is the resulting one-step transition matrix

opt_TML <- matrix(rep(NA,nb_states^2), nrow=nb_states)

for (l in 1:nb_states){
  opt_TML[l,] <- TML[[optimal_policy[l]]][l,]
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

# Expected number of steps from each state given that one follows the optimal policy
expected_number_Of_steps_to_absorption(1:(nb_states-1),opt_TML)
