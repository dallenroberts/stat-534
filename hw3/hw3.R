################################################################################
## Allen Roberts
## April 16, 2020
## Stat 534
## Homework 3
################################################################################
require(rcdd)

## Return the AIC value of a logistic regression. Also checks for convergence
## which might be redundant given that the algorithm is going to use the Geyer
## method
getLogisticAIC <- function(y, vars, data)
{
  # check if the regression has no explanatory variables
  if(0 == length(vars))
  {
    out = glm(data[,y] ~ 1, family = binomial);
  }
  
  # regression with at least one explanatory variable. We suppress
  # convergence warnings to reduce screen clutter. 
  else
  {
    out = suppressWarnings(glm(data[,y] ~ as.matrix(data[,as.numeric(vars)]),
                               family = binomial));
  }
  
  # compute the AIC
  AIC = out$deviance+2*(1+length(vars))
  
  # we must check whether the glm algorithm properly converged,
  # as judged by the algorithm.
  
  converged = out$converged
  
  # when glm fails to converge return NA, otherwise return AIC
  return(ifelse(converged, AIC, NA))
}

## Checks if logistic regression coefficient estimates are valid according to
## Geyer method
isValidLogisticRCDD <- function(response,explanatory,data)
{
  if(0==length(explanatory))
  {
    #we assume that the empty logistic regresion is valid
    return(TRUE);
  }
  logisticreg = suppressWarnings(glm(data[,response] ~ as.matrix(data[,as.numeric(explanatory)]),
                                     family=binomial(link=logit),
                                     x=TRUE));
  tanv = logisticreg$x;
  tanv[data[,response] == 1, ] <- (-tanv[data[,response] == 1, ]);
  vrep = cbind(0, 0, tanv);
  #with exact arithmetic; takes a long time
  #lout = linearity(d2q(vrep), rep = "V");
  
  lout = linearity(vrep, rep = "V");
  return(length(lout)==nrow(data));
}

## Generate random starting logistic regression model with k predictors. 
## Iter controls how many times new predictors will be sampled before 
## giving up.
getRandomModel <- function(data, y, iter = 100) {
  
  # Counter of number of tries to find valid model (max = iter)
  t <- 0
  
  while(t < iter) {
    
    ## Sample k random predictors, where k is uniformly sampled from the set of
    ## integers from 1 to ncol(data).
    vars <- 1:ncol(data)[-y]
    k <- sample(vars, size = 1)
    
    randomVars <- sample(vars, size = k)
    
    ## If regression is valid, return model
    if(isValidLogisticRCDD(y, randomVars, data)) {
      
      return(randomVars)
      
    } else {
      t <- t + 1
    }
  }
  
  cat("\n\nError generating initial random model - No valid function found after",
      iter, "iterations\n\n")
  break
}

## Identify neighbors: This function returns a matrix of the column indices of
## data corresponding to the variables included in each element of the set of 
## neighbors. 

identifyAllNeighbors <- function(data, y, currentMod) {
  p <- ncol(data) - 1
  
  ## Subtract one variable
  if(length(currentMod) > 0) {
    
    neighbors_sub <- matrix(0,
                            nrow = length(currentMod),
                            ncol = ncol(data))
    
    for(v in 1:length(currentMod)) {
      
      ## Subtract one of the variables in the current model
      newSet <- currentMod[-v]
      
      ## Add entry into neighbors matrix
      for(c in 1:ncol(neighbors_sub)) {
        neighbors_sub[v, c] <- ifelse(c %in% newSet, 1, 0)
      }
    }
  }
  
  ## Add one variable
  if(length(currentMod) < p) {
    
    neighbors_add <- matrix(0, 
                            nrow = p - length(currentMod),
                            ncol = ncol(data))
    ## Add one variable
    allVars <- 1:ncol(data)
    addVars <- allVars[-c(y, currentMod)]
    
    for(v in 1:length(addVars)) {
      ## Add variable
      newSet <- c(currentMod, addVars[v])
      
      for(c in 1:ncol(neighbors_add)) {
        
        neighbors_add[v,c] <- ifelse(c %in% newSet, 1, 0)
        
      }
    }
  }
  
  ## Combined set of neighbors
  if(length(currentMod) == 0) {
    if(sum(neighbors_add[, y]) > 0) {
      cat("\n\nError in identify neighbors: 
            trying to add the response variable!\n\n")
      break
    } else {
      return(neighbors_add)
    }
  } else if(length(currentMod) == p) {
    return(neighbors_sub)
  } else {
    if(sum(neighbors_add[, y]) > 0) {
      cat("\n\nError in identify neighbors:
          trying to add the response variable!\n\n")
      break
    } else {
      neighbors <- rbind(neighbors_sub, neighbors_add)
      return(neighbors)
    }
  }
}

## Identify valid neighbors. This function takes a matrix of neighbor indices
## and checks whether each set results in a valid logistic regression model.
## It then returns a subset of the neighbor marix  for which the models are 
## valid.
identifyValidNeighbors <- function(data, y, neighbors) {
  
  validNeighborSets <- NULL
  
  for(ii in 1:nrow(neighbors)) {
    
    vars <- which(neighbors[ii, ] == 1)
    
    if(isValidLogisticRCDD(y, vars, data)) {
      validNeighborSets <- c(validNeighborSets, ii)
    }
  }
  
  if(is.null(validNeighborSets)) {
    
    cat("\n\n No valid neighbors found!\n\n")
    break
    
  } else {
    validNeighbors <- neighbors[validNeighborSets, ]
    return(validNeighbors)
  }
  
}
## Implement Markov Chain Monte Carlo search for model that returns the best 
## AIC value.
MC3search <- function(response, data, n_iter, verbose = FALSE) {
  
  bestMod <- currentMod <- NULL
  bestAIC <- currentAIC <- NULL
  
  ## Iteration 0: Start with a random model
  iter <- 0
  
  vars <- getRandomModel(data, response)
  bestAIC <- currentAIC <- getLogisticAIC(response, vars, data)
  bestMod <- currentMod <- vars
  
  while(iter <= n_iter) {
    if(verbose == TRUE) {
      cat("\n\nIteration", iter, ":\n")
      cat("Best AIC:", bestAIC, "\n")
      cat("Current model:", currentMod, "\n")
      cat("Current AIC:", currentAIC, "\n")
    }
    ## Subsequent iterations
    ## Step 1: Identify neighbors of current logistic regression model
    neighbors <- identifyAllNeighbors(data, response, currentMod)
    
    ## Step 2: Eliminate models that don't form valid logistic regressions from
    ## the set of neighbors
    validNeighbors <- identifyValidNeighbors(data, response, neighbors)
  
    ## Step 3: Uniformly sample a model from the set of valid neighbors
    sampledRow <- validNeighbors[sample(1:nrow(validNeighbors), 1), ]
    newMod <- which(sampledRow == 1)
    
    ## Step 4: Uniformly sample a model from the set of valid neighbors of the 
    ## newly sampled model
    neighborsNewMod <- identifyAllNeighbors(data, response, newMod)
    validNeighborsNewMod <- identifyValidNeighbors(data, response, neighborsNewMod)
    
    ## Step 5: Calculate p_A'
    aicNew <- getLogisticAIC(response, newMod, data)
    pANew <- -aicNew - log(nrow(validNeighborsNewMod))

    
    ## Step 6: Calculate p_A
    pACurrent <- -currentAIC - log(nrow(validNeighbors))
    
    if(verbose == TRUE) {
      cat("Candidate model:", newMod, "\n")
      cat("Candidate AIC:", aicNew, "\n")
      cat("pA_r:", pACurrent, "\n")
      cat("pA':", pANew, "\n")
    }

    ## Step 7:
    if(pANew > pACurrent) {
      
      if(verbose == TRUE) {
        cat("pA' > pA_r. New model accepted.\n")
      }
      currentMod <- newMod
      currentAIC <- aicNew
      
      if(currentAIC < bestAIC) {
        if(verbose == TRUE) {
          cat("AIC improved.\n")
        }
        bestAIC <- currentAIC
        bestMod <- currentMod
      } else {
        if(verbose == TRUE) {
          cat("AIC not improved.\n")
        }
      }
    } else {
      ## Step 8: If pAPrime < pA, sample u from the uniform distribution on (0,1).
      u <- runif(n = 1, min = 0, max = 1)
      
      ## If log(u) < paPrime - pA, MaPrime becomes the current model. 
      ## Otherwise, Ma remains the current model.  
      if(log(u) < (pANew - pACurrent)) {
        
        currentMod <- newMod
        currentAIC <- aicNew
        if(verbose == TRUE) {
          cat("pA' <= pA_r, but new model accepted.\n")
        }

        if(currentAIC < bestAIC) {
          if(verbose == TRUE) {
            cat("AIC improved.\n")
          }

          bestAIC <- currentAIC
          bestMod <- currentMod
        } else {
          if(verbose == TRUE) {
            cat("AIC not improved.\n")
          }
        }
      } else {
        if(verbose == TRUE) {
          cat("pA' <= pA_r, and new model NOT accepted.\n")
        }
      }
    }
    
    iter <- iter + 1
  }
  if(verbose == TRUE) {
    cat("MC3 search finished\n\n")
  }
  return(list("bestAIC" = bestAIC, "bestAICvars" = sort(bestMod)))
  
}

main <- function(datafile) {
  
  ## For reproducibility
  set.seed(123)
  
  #read the data
  data = read.table(datafile,header=FALSE);
  response = ncol(data)
  
  ## Problem 1 - example with 5 iterations
  p1 <- MC3search(response, data, n_iter = 5, verbose = FALSE)
  cat("\n\nProblem 1: Demonstration of MC3search with 5 iterations\n")
  cat("Best model explanatory variable indices:", p1$bestAICvars, "\n")
  cat("Best AIC:", p1$bestAIC, "\n")

  ## Problem 2
  cat("\n\nProblem 2: 10 instances of MC3search, each with 25 iterations\n")
  for(ii in 1:10) {

    cat("\n\nInstance", ii, ":\n")
    bestMod <- MC3search(response, data, n_iter = 25, verbose = FALSE)
    # cat("Final model results:\n")
  cat("Best model explanatory variable indices:", bestMod$bestAICvars, "\n")
  cat("Best AIC:", bestMod$bestAIC, "\n")

  }
  
  
}

main('534binarydata.txt')

###############################################################################
## Comments
###############################################################################
## In problem 2, MC3search runs 10 times, each with 25 iterations. Each instance
## resulted in a different final model with a different best AIC value. The best
## AIC value ranged from 67.00736 to 90.93704, which is quite wide. Thes results
## suggest that 25 iterations is not nearly sufficient for this algorithm to 
## converge to the optimal model for these data. I suspect the algorithm is also
## sensitive to the value of k that is chosen, at least when the max number of 
## iterations is small.

## A few notes on the functions. I included a verbose flag on MC3search() that, 
## if set to TRUE, allows the results from each iteration to be printed.This 
## was useful in troubleshooting my function. 

## Second, the getRandomModel() function presented an interesting challenge, 
## because if k is large, then it's possible that there won't be any models that
## are valid logistic regression models. To address this potential issue, if 
## getRandomModel() gets an invalid model on the first try, it resamples k before
## it resamples another model. I suppose this makes it even more likely to start
## with a model with fewer indices. In the off-chance that it can't find a valid
## model at all, it stops after 100 tries and throws an error.
##
## Also, my model doesn't run very fast, and the slow step is the function that
## identifies neighbors (or valid neighbors). I'm sure there is a more efficient
## method than the for loop that I implemented, and I'd love to get some 
## other ideas.




