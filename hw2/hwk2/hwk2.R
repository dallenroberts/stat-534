rm(list = ls())

getLogisticAIC <- function(response, explanatory = NULL, data) {
  
  ## Check if the regression has no explanatory variables
  if(0 == length(explanatory)) {
    
    ## Regression with no explanatory variables
    deviance <- glm(data[, response] ~ 1, 
                   family = binomial(link = logit))$deviance
    
  } else {
    
    ## Regression with at least one explanatory variable
    deviance <- glm(data[, response] ~ as.matrix(data[, as.numeric(explanatory)]),
                   family = binomial(link = logit))$deviance
    
  }
  
  ## Calculate AIC
  return(deviance + 2 * (1 + length(explanatory)))
  
}

getLogisticBIC <- function(response, explanatory = NULL, data) {
  
  ## Check if the regression has no explanatory variables
  if(0 == length(explanatory)) {
    
    ## Regression with no explanatory variables
    deviance <- glm(data[, response] ~ 1, 
                    family = binomial(link = logit))$deviance
    
  } else {
    
    ## Regression with at least one explanatory variable
    deviance <- glm(data[, response] ~ as.matrix(data[, as.numeric(explanatory)]),
                    family = binomial(link = logit))$deviance
    
  }
  
  ## Calculate BIC
  return(deviance + log(nrow(data))*(1 + length(explanatory)))
  
}


forwardSearch <- function(response, data, lastPredictor, method) {
  
  stopifnot(method %in% c("AIC", "BIC"))
  
  #Start with the empty regression with no predictors
  bestRegression <- NULL
  
  #Calculate the AIC/BIC of the empty regression
  bestRegressionValue <- ifelse(method == "AIC", getLogisticAIC(response,bestRegression,data), getLogisticBIC(response,bestRegression,data))
  cat('\n\n\n\nforwardSearch :: The empty logistic regression has ', method, ' = ', bestRegression,'\n\n\n')
  
  #Vector that keeps track of all the variables that are currently NOT in the model
  VariablesNotInModel <- 1:lastPredictor
  
  #Add variables one at a time to the current regression and retain that variable that gives the smallest values of AIC/BIC associated
  #Make the model that includes that variable to be the current model if it gives a smaller AIC/BIC than the AIC/BIC of the current regression
  
  #stop when there are no variables that can be included in the model
  stepNumber = 0;
  while(length(VariablesNotInModel)>=1)
  {
    #record the number of steps performed
    stepNumber <- stepNumber + 1
    
    cat("Step Number: ", stepNumber, "\n")
    
    #create a vector that records the AIC/BIC values of the regressions
    #we are examining; the number of these regressions is equal
    #with the number of variables that are not in the model
    reg <- vector('numeric',length(VariablesNotInModel))
    
    # cat("Variables Not in Model: ", VariablesNotInModel, "\n")
    
    #take each variable that is not in the model and include it in the model
    for(ii in 1:length(VariablesNotInModel)) {
    
      # cat("\n", ii, "\n")
      
      addedPred <- VariablesNotInModel[ii]

      reg[ii] <- ifelse(method == "AIC", getLogisticAIC(response, c(bestRegression, addedPred), data), getLogisticBIC(response, c(bestRegression, addedPred), data))
      
    }
    
    # cat("\n", reg, "\n")

    if(min(reg) < bestRegressionValue) {
      
      bestRegressionValue <- min(reg)
      bestPredictor <- VariablesNotInModel[which.min(reg)]
      bestRegression <- c(bestRegression, bestPredictor)
      VariablesNotInModel <- VariablesNotInModel[VariablesNotInModel != bestPredictor]
      
      cat("\nforwardSearch :: Logistic regression with ", stepNumber, " predictors :: \nPredictorAdded: ", bestPredictor, "\nPredictor set: ", bestRegression, "\nBest ", method, " : ", bestRegressionValue, "\n\n\n\n")
      
    } else {
      cat("\nforwardSearch :: Logistic regression with ", stepNumber, " predictors did not improve ", method, ":: \n")
      
      break
    }
    
  }
  
  return(list(method = method, value=bestRegressionValue,reg=bestRegression));
}

backwardSearch <- function(response, data, lastPredictor, method) {
  
  stopifnot(method %in% c("AIC", "BIC"))
  
  #Start with the regression with all predictors
  bestRegression <- 1:lastPredictor
  
  #Calculate the AIC/BIC of the empty regression
  bestRegressionValue <- ifelse(method == "AIC", getLogisticAIC(response,bestRegression,data), getLogisticBIC(response,bestRegression,data))
  cat('\n\n\n\nbackwardSearch :: The full logistic regression has ', method, ' = ',bestRegressionValue,'\n\n\n')
  
  #Vector that keeps track of all the variables that are currently in the model
  VariablesInModel <- 1:lastPredictor
  
  # Subtract variables one at a time from the current regression and identify the variable that, when removed, gives the smallest values of AIC/BIC
  #Make the model that subtracts that variable to be the current model if it gives a smaller AIC/BIC than the AIC/BIC of the current regression
  
  #stop when there are no variables that can be removed from the model
  stepNumber = 0;
  while(length(VariablesInModel)>=1)
  {
    #record the number of steps performed
    stepNumber <- stepNumber + 1
    
    cat("Step Number: ", stepNumber, "\n")
    
    #create a vector that records the AIC/BIC values of the regressions
    #we are examining; the number of these regressions is equal
    #with the number of variables that are in the model
    reg <- vector('numeric',length(VariablesInModel))
    
    # cat("Variables In Model: ", VariablesInModel, "\n")
    
    #take each variable that is in the model and remove it from the model
    for(ii in 1:length(VariablesInModel)) {
      
      toDrop <- VariablesInModel[ii]
      newSet <- bestRegression[bestRegression != toDrop]
      reg[ii] <- ifelse(method == "AIC", getLogisticAIC(response, newSet, data), getLogisticBIC(response, newSet, data))
      
    }
    
    # cat("\n", reg, "\n")
    
    if(min(reg) < bestRegressionValue) {
      
      bestRegressionValue <- min(reg)
      removedPredictor <- VariablesInModel[which.min(reg)]
      bestRegression <- bestRegression[bestRegression != removedPredictor]
      VariablesInModel <- VariablesInModel[VariablesInModel != removedPredictor]
      
      cat("\nbackwardSearch :: Logistic regression with ", stepNumber, " predictors removed :: \nPredictorRemoved: ", removedPredictor, "\nPredictor set: ", bestRegression, "\nBest ", method, ": ", bestRegressionValue, "\n\n\n\n")
      
    } else {
      cat("\nbackwardSearch :: Logistic regression with ", stepNumber, " predictors removed did not improve ", method, ":: \n")
      
      break
    }
    
  }
  
  return(list(value=bestRegressionValue,reg=bestRegression));
}


main <- function(datafile)
{
  #read the data
  data <- read.table(datafile,header=FALSE)
  
  #the sample size is 148 (number of rows)
  #the explanatory variables are the first p-1 columns
  #the last column is the binary response
  response  <- ncol(data)
  
  lastPredictor <- ncol(data)-1
  
  ## Problem 2
  #perform a forward "greedy" search for the best logistic regression
  #i.e., the logistic regression with the smallest AIC
  forwardResultsAIC <- forwardSearch(response,data,lastPredictor, method = "AIC")
  
  ## Problem 3
  #perform a backward "greedy" search for the best logistic regression
  backwardResultsAIC <- backwardSearch(response,data,lastPredictor, method = "AIC");
  
  ## Problem 4
  #perform a forward "greedy" search for the best logistic regression
  #i.e., the logistic regression with the smallest BIC
  forwardResultsBIC <- forwardSearch(response,data,lastPredictor, method = "BIC")
  
  #perform a backward "greedy" search for the best logistic regression
  backwardResultsBIC <- backwardSearch(response,data,lastPredictor, method = "BIC");
  
  # output the results of our searches for AIC
   cat('\n\nForward search using AIC gives regression with ',length(forwardResultsAIC$reg),'explanatory variables [')
   if(length(forwardResultsAIC$reg)>=1)
   {
     for(i in 1:length(forwardResultsAIC$reg)) cat(' ',forwardResultsAIC$reg[i])
   }
   cat('] with AIC = ',forwardResultsAIC$value,'\n')
   
   cat('\n\nBackward search using AIC gives regression with ',length(backwardResultsAIC$reg),'explanatory variables [');
   if(length(backwardResultsAIC$reg)>=1)
   {
     for(i in 1:length(backwardResultsAIC$reg)) cat(' ',backwardResultsAIC$reg[i]);
   }
   cat('] with AIC = ',backwardResultsAIC$value,'\n');
   
   # output the results of our searches for BIC
   cat('\n\nForward search using BIC gives regression with ',length(forwardResultsBIC$reg),'explanatory variables [')
   if(length(forwardResultsBIC$reg)>=1)
   {
     for(i in 1:length(forwardResultsBIC$reg)) cat(' ',forwardResultsBIC$reg[i])
   }
   cat('] with BIC = ',forwardResultsBIC$value,'\n')
   
   cat('\n\nBackward search using BIC gives regression with ',length(backwardResultsBIC$reg),'explanatory variables [');
   if(length(backwardResultsBIC$reg)>=1)
   {
     for(i in 1:length(backwardResultsBIC$reg)) cat(' ',backwardResultsBIC$reg[i]);
   }
   cat('] with BIC = ',backwardResultsBIC$value,'\n');
}

main('534binarydata.txt')

## Comments for Problem 4
## Using AIC, forward and backward selection both resulted in logistic regression models with 10 explanatory variables that resulted in AIC of 22. However, the included variables differend between the two selection procedures. Of the 10 included variables, only four variables were shared between the forward and backward selection results. Similarly, forward and backward selection methods using BIC both resulted in logistic regression models with 10 explanatory variables and BIC of 54.96934. Again, the specific variables included differed between the forward and backward selection results, with only five variables shared between the two sets. The specific variables included in each regression model are printed when this code is run.
