
#Making the deck
  
  hand <- rep(c(1:13), times=4)
  shape <- rep(c(1:4),each=13)
  deck <- data.frame(shape, hand, stringsAsFactors = FALSE)



# Checks if the given hand is a tris
check_tris <- function(df){
  tris <- FALSE
  count <- c(replicate(13,0))
  for(i in 1:5){
      count[df$hand[i]] = count[df$hand[i]] + 1
  }
  for(i in 1:13){
    if(count[i]==2){
      return(0)
    }
    if(count[i]==3){
      tris = TRUE 
    }
  }
  if(tris == TRUE){
    return(1)
  }
  return(0)
}


#Checks if the given hand is a poker
check_four <- function(df){
  c <- 0 
  for(i in 1:2){
    count <- 0 
    for(j in (i+1):5){
      if(df$hand[i]==df$hand[j]){
        c <- 1 
        count <- count + c
      }
      if(count == 3){
        return(1)
      }
    }
  }
  return(0)
}

#Checks if the given hand is a full house
check_full <- function(df){
  full <- FALSE
  check_tris <- FALSE
  check_pair <- FALSE
  count <- c(replicate(13,0))
  for(i in 1:5){
    count[df$hand[i]] = count[df$hand[i]] + 1
  }
  for(i in 1:13){
    if(count[i]==2){
      check_pair <- TRUE
    }
    if(count[i]==3){
      check_tris <- TRUE
    }
    
  }
  if((check_tris & check_pair) == TRUE){
    return(1)
  }
  else{
    return(0)
  }
}



#Checks the case A K Q J 10 (in our model 1 10 11 12 13). This function supposes that the 
#hand is already sorted by value.

check_top_straight <- function(hand){
  if(hand[1] == 1 & hand[2] == 10 & hand[3] == 11 & hand[4] == 12 & hand[5] == 13)
    return(1)
  else{
    return(0)
  }
}


check_double_pair <- function(df){
  count <- c(replicate(13,0))
  for(i in 1:5){
    count[df$hand[i]] = count[df$hand[i]] + 1
  }
  count_pairs <- 0
  for(i in 1:13)
  {
    if(count[i] >= 3)
    {
      return(0)
    }
    if(count[i] == 2)
    {
      count_pairs = count_pairs + 1
    }
  }
  if(count_pairs==2)
  {
    return(1)
  }
  else{
    return(0)
  }
}


# Function that implements the montecarlo simulation

montecarlo <- function(f){
  prob <- c()
  count_f <- 0
  prob <- c()
  x <- c()
  j <- 1
  for(i in 1:nsim){
    
    df <- deck[sample(nrow(deck),5,replace = FALSE),]
    
    count_f <- count_f + f(df)
    if(i %% 100 == 0){
      prob[j] <- count_f / i 
      x[j] <- i
      j <- j +1 
    }
  }
  prob_f <- count_f / nsim
  print(prob_f)
  result <- data.frame(prob,x)
  return(result)
}


# Here we call all the montecarlo routines that we want to execute

# Set the number of trials
nsim <- 1000000


# Perform simulations

df <- montecarlo(check_tris)
#Plot tris probability
plot(log(df$x),df$prob)
title("Tris probability")
abline(h = 0.021128, col = "red")


df <- montecarlo(check_double_pair)
# Plot double pair probability
plot(log(df$x),df$prob)
title("Double pair probability")
abline(h = 0.047539, col = "red")


df <- montecarlo(check_four)
# Plot four probability
plot(log(df$x),df$prob) 
title("Poker probability")
abline(h = 0.0002401, col = "red")

df <- montecarlo(check_full)
# Plot full house probability
plot(log(df$x),df$prob) 
title("Full house probability")
abline(h = 0.001441, col = "red")


