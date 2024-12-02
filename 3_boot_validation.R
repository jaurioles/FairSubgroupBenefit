### Before we have validated the models using split-sampling.
### This has caused low precision, we will try to recuperate it
### using bootstrapping.

# Define Source Directory
setwd("/mnt/bmh01-rds/Sperrin_UKBB_Fairness/Code/study_1_benefit/Manuscript_Ready_Code/FairSubgroupBenefit/")


# Set random seed
RANDOM_SEED <- 42
set.seed(RANDOM_SEED)

#------------------------#
#------------------------#
####  Import libraries ####
#------------------------#
#------------------------#

library(probably)  # For calibration plot
library(tidyverse)  # For data managing
library(pROC)       # For calculating performance
library(glmnet)     # First model: Elastic-Net
library(xgboost)    # Second model: XGBoost
library(rBayesianOptimization)  # For HP tuning
library(caret)      # For confusion matrix
library(pracma)     # For trapezoidal integration
library(doParallel) # For parallel computing
library(KraljicMatrix) # For identifying Pareto
library(mice)       # Data imputation
library(glue)
library(boot)       # For bootstrapping

#------------------#
#------------------#
####  Parameters ####
#------------------#
#------------------#

# Multiple imputation parameters
# MICE_MAXIT <- 30 # Leave iterations to 30 for convergence
MICE_MAXIT <- 5 # To make faster for code snippet
MICE_M <- 1     # We keep multiple imputations to 1, as missing data is not the focus of this work

# Set parallel clusters
registerDoParallel(cores=5)

# Number of bootstraps
# NUM_BOOT <- 5 # In the manuscript we used 500
NUM_BOOT <- 5 

#------------------------------#
####  NET BENEFIT PARAMETERS ####
#------------------------------#

# For diabetes, we need thresholds over which to integrate, and their weight
diabetes_threshold <- 0.15
diabetes_effect <- 0.58
lung_cancer_threshold <- 0.015
lung_cancer_effect <- 0.20

# For dca, what is the range?
dcax_diabetes <- seq(0,0.995,0.005)
dcax_lungcancer <- seq(0,0.2,0.002)

#-----------------------#
#-----------------------#
####  Misc. Functions ####
#-----------------------#
#-----------------------#

# Calculates the net benefit curve (Corrected as per our work)
calcNetBenefit <- function(clinical_thresholds,model_thresholds,y_true,s_pred,new_version=TRUE) {
  
  # Calculate policies per threshold
  N <- length(y_true)
  pred_matrix <- sapply(model_thresholds, function(thresh) as.numeric(s_pred >= thresh))
  
  # Calculate true positives and false positives
  tp <- colSums(pred_matrix * y_true)
  fp <- colSums(pred_matrix * (1 - y_true))
  
  if (new_version==TRUE) {
    nb_vector_to_add <- (tp / clinical_thresholds - (1 / (1 - clinical_thresholds)) * fp)/N
  } else {
    nb_vector_to_add <- (tp - (clinical_thresholds / (1 - clinical_thresholds)) * fp)/N
  }
  
  return(nb_vector_to_add)
}

# Define diabetes Sensitive attribute from ethnicity
define_diabetes_sa <- function(ethnicity) {
  
  # Define 3 groups: White, Asian, Black and and Other
  sensitive_attribute <- as.character(ethnicity)
  
  # For fairness evaluations
  sensitive_attribute <- factor(sensitive_attribute,levels=c("White","Asian","Black","Other"))
  
  return(sensitive_attribute)
}

# Define lung cancer Sensitive attribute from IMD
define_lungcancer_sa <- function(imd,quintiles_imd) {
  
  # get the intervals of each individual according to their imd
  imd_q <- as.factor(findInterval(imd,quintiles_imd))
  
  return(imd_q)
}

# Turns dataframe into x and y variables ready for glmnet
prepareForGLMNET <- function(df) {
  
  # Data clean
  x <- df %>% select(-y,-sa)
  x <- model.matrix( ~ .,x)
  y <- df %>% select(y) %>% data.matrix
  sa <- df %>% select(sa)
  
  # Return data as named list
  return(list(x=x,y=y,sa=sa))
}

# Fit a CV LASSO regression on data to find which set of variables is important
find_important_variables_lasso <- function(x_matrix,y_matrix,s="lambda.min") {
  
  # Create an object `fit` that takes is a cv.glmnet lasso fit of y~x
  model_fit <- cv.glmnet(
    x = x_matrix,
    y = y_matrix,
    family = "binomial",
    type.measure = "deviance",
    nfolds = 5,
    trace.it = FALSE,
    parallel = TRUE
  )
  
  # Extract the coefficients index used  in the lambda.min
  model_coefs <- predict(model_fit,s=s,type="nonzero")$lambda.min
  
  return(model_coefs)
}

# For each group as defined by sa_levels, calculate the propensity score to identify that group
calculate_sensitive_propensity_score <- function(predictor_df,sensitive_attribute,sa_levels,s="lambda.min") {
  
  # We choose the variables that are most important for estimating Y|X
  # We can only do this if we have Y
  # Prepare predictor matrix
  if ("y" %in% colnames(predictor_df)) {
    
    # Predictors to include
    predictor_x <- predictor_df %>% select(-y)
    predictor_x <- model.matrix( ~ .,predictor_x)
    predictor_y <- predictor_df %>% select(y) %>% data.matrix
    
    # Slim down predictor_x by choosing only those that predict y in LASSO
    which_x <- find_important_variables_lasso(predictor_x,predictor_y,s=s)
    predictor_x <- predictor_x[,which_x]
    
    # Predictor matrix is X
    predictor_matrix <- predictor_x
    
  } else {
    
    predictor_x <- predictor_df %>% select(-y)
    predictor_x <- model.matrix( ~ .,predictor_x)
    
    # Without y, its just x values without cutting given y|x
    predictor_matrix <- predictor_x
    
  }
  
  
  # The list to return will be named according to each sa_level
  # And have a list of length dim(predictors)[1] with the propensity scores for each one
  propensity_list <- list()
  # We also save the models
  model_list <- list()
  
  # For each group, create correct model
  for (each_sa in sa_levels) {
    
    # What is our goal of prediction?
    sa_to_predict <- as.numeric(sensitive_attribute==each_sa)
    
    # Create a crossvalidated model
    propensity_model <- cv.glmnet(
      x = predictor_matrix,
      y = sa_to_predict,
      family = "binomial",
      type.measure = "deviance",
      nfolds = 5,
      trace.it = FALSE,
      parallel = TRUE
    )
    
    # Report AUC out of curiosity
    pred <- predict(propensity_model, newx = predictor_matrix, type = "response",s = "lambda.min")[,1]
    auroc <- (roc(c(sa_to_predict),c(pred),ci=FALSE, levels = c(0, 1), direction = "<"))$auc
    print(
      paste(
        "The AUC of the propensity score to discriminate ",as.character(each_sa),
        " is ",
        auroc
      )
    )
    
    # Create linear predictor for prop score
    propensity_scores <- predict(propensity_model, newx = predictor_matrix, type = "response",s = "lambda.min")[,1]
    propensity_scores <- propensity_scores/(1-propensity_scores)
    # You need to multiply by overall odds prevalence of group so og group distribution
    # Can have a weight of one
    propensity_scores <- propensity_scores*(1-mean(sa_to_predict))/mean(sa_to_predict)
    
    # Cap at 1
    # WE HAVE COMMENTED THIS OUT
    # The reason is that we apply this later, and we want the capping to happen
    # AFTER we apply the forgetting factor
    # propensity_scores <- pmin(propensity_scores, 1)
    
    # If original data is in sa, keep at 1
    propensity_scores[sa_to_predict==1] <- 1
    
    # Assign scores to corresponding group
    propensity_list[[each_sa]] <- propensity_scores
    
    # Also save models
    model_list[[each_sa]] <- propensity_model
  }
  
  # Create a combination list of variable names, propensity score list, and model list
  return_list <- list()
  return_list$propensity <- propensity_list
  return_list$model <- model_list
  
  # Return full list
  return(return_list)
}

# Cross-Entropy Loss
logLoss <- function(y_true,y_prob) {
  -1*sum((1-y_true)*log(1-y_prob) + y_true*log(y_prob))
}

# Get bootstrap info
get_bootstrap_vector <- function(lst) {
  result <- list()
  
  for (ii in seq_along(lst)) {
    for (name in names(lst[[ii]])) {
      # Recursively restructure nested lists
      if (is.list(lst[[ii]][[name]])) {
        result[[name]][[ii]] <- get_bootstrap_vector(lst[[ii]][[name]])
      } else {
        # Directly assign non-list elements
        result[[name]][[ii]] <-  c(result[[name]], lst[[ii]][[name]])
      }
    }
  }
  
  return(result)
}

#----------------#
#----------------#
####  DIABETES ####
#----------------#
#----------------#

#----------------------------#
####    Open data and impute ####
#----------------------------#

df <- read.csv(file = paste0("Data/diabetes_df_randomsplit.csv"))

# X is rownames
rownames(df) <- df$X
df <- df %>% select(-X)

# Turn factor values into factors
df$Ethnicity <- as.factor(df$Ethnicity)
df$Education <- as.factor(df$Education)
df$Smoker <- as.factor(df$Smoker)

# Remove values not used in prediction or y or sa
df[c("Censored","Started","Death","Centre","is_val")] <- NULL

# Change to factors strings and binary variables
df[, sapply(df, is.character)] <- lapply(df[, sapply(df, is.character)], as.factor)

# Get quintiles of Townsend
quintiles_townsend <- quantile(df$Townsend,probs = c(0.2,0.4,0.6,0.8),na.rm=TRUE)

# Exclude 'y' from being used as a predictor in MICE
predictorMatrix <- mice::make.predictorMatrix(df)
predictorMatrix[, "y"] <- 0   # Do not use 'y' as a predictor for any other variable
predictorMatrix["y", ] <- 0   # Do not impute 'y' (it has no missing values)

# Impute data using mice
df_mice <- mice(df,maxit = MICE_MAXIT, m = MICE_M, predictorMatrix = predictorMatrix)
df <- complete(df_mice, action = "stacked")
df <- na.omit(df)

# Define sa according to ethnicity
df$sa <- define_diabetes_sa(df$Ethnicity)

# Save and load
saveRDS(df, "Data/diabetes_df_randomsplit_imputed.rds")
df <- readRDS("Data/diabetes_df_randomsplit_imputed.rds")

#--------------------#
####    Prepare data ####
#--------------------#

# Turn validation and training data into right format
x_and_y <- prepareForGLMNET(df)
x <- x_and_y[["x"]]
y <- x_and_y[["y"]]

# What are the sa levels?
sa_levels <- levels(df$sa)

#-----------------------------------#
####    Calculate propensity scores ####
#-----------------------------------#

# Ideally, we would do the propensity score inside the bootstrap, but this is
# too computationall intensive.

# Find out propensity_score
propensity_train_df <- df %>% select(-sa,-Ethnicity)
propensity_score_result <- calculate_sensitive_propensity_score(
  propensity_train_df,df$sa,sa_levels, s = "lambda.min"
)


#------------------------------#
#------------------------------#
####    Calculate the Optimism ####
#------------------------------#
#------------------------------#

# Here is the dictionary we will be recording in the metrics per bootstrap
bootstrap_optimism <- list()

# Here, our validation data is x and y
x_val <- x
y_val <- y
sa_val <- df$sa 

set.seed(RANDOM_SEED)
bootstrap_optimism <- foreach(bb=1:NUM_BOOT) %dopar% {
  
  # Here we save everything
  per_bootstrap_metrics <- list()
  
  print(glue("Bootstrap Number {bb}"))
  
  # Select the new data
  bootstrap_index <- sample(nrow(df), replace = TRUE)
  
  # Training data
  x_train <- x_val[bootstrap_index,]
  y_train <- y_val[bootstrap_index]
  sa_train <- sa_val[bootstrap_index]
  
  for (model_name in c("LogNoSA", "LogSingleSA", "LogMultiSA")) {
    
    # Train the model
    if (model_name=="LogNoSA") {
      
      fit <- cv.glmnet(
        x = x_train[, !(colnames(x) %in% c("EthnicityBlack", "EthnicityOther", "EthnicityWhite"))],
        y = y_train,
        family = "binomial",
        type.measure = "deviance",
        nfolds = 5,
        trace.it = 0,
        parallel = FALSE
      )
      # What is the predictions of the model?
      pred_train <- predict(fit,
                            newx = x_train[, !(colnames(x) %in% c("EthnicityBlack", "EthnicityOther", "EthnicityWhite"))],
                            type = "response",s = "lambda.min")[,1]
      pred_val <- predict(fit,
                          newx = x[, !(colnames(x) %in% c("EthnicityBlack", "EthnicityOther", "EthnicityWhite"))],
                          type = "response",s = "lambda.min")[,1]
      
    } else if (model_name=="LogSingleSA") {
      
      fit <- cv.glmnet(
        x = x_train,
        y = y_train,
        family = "binomial",
        type.measure = "deviance",
        nfolds = 5,
        trace.it = 0,
        parallel = FALSE
      )
      # What is the predictions of the model?
      pred_train <- predict(fit,
                            newx = x_train,
                            type = "response",s = "lambda.min")[,1]
      pred_val <- predict(fit,
                          newx = x,
                          type = "response",s = "lambda.min")[,1]
      
    } else if (model_name=="LogMultiSA") {
      
      # For each group we fit a separate model
      model_per_group <- list()
      
      # For each group, train a model
      for (each_sa in sa_levels) {
        
        # Get best weights for SA
        training_weights <- propensity_score_result$propensity[[each_sa]][bootstrap_index]
        training_weights[sa_train != each_sa] <- pmin(1,training_weights[sa_train != each_sa])
        
        # Train model with weights
        trained_model <- cv.glmnet(
          x = x_train,
          y = y_train,
          family = "binomial",
          weights = training_weights,
          parallel = FALSE,
          type.measure = "deviance",
          nfolds = 5,
          trace.it = 0
        )
        
        # Add model to full list
        model_per_group[[each_sa]] <- trained_model
        
      }
      
      # What is the predictions of the model?
      pred_train <- predict(model_per_group[["White"]],
                            newx = x_train,
                            type = "response",
                            s = "lambda.min")[,1]
      pred_val <- predict(model_per_group[["White"]],
                          newx = x,
                          type = "response",
                          s = "lambda.min")[,1]
      for (each_sa in sa_levels) {
        pred_train[sa_train==each_sa] <- predict(model_per_group[[each_sa]],
                                        newx = x_train[sa_train==each_sa,],
                                        type = "response",
                                        s = "lambda.min")[,1]
        pred_val[sa_val==each_sa] <- predict(model_per_group[[each_sa]],
                                        newx = x_val[sa_val==each_sa,],
                                        type = "response",
                                        s = "lambda.min")[,1]
      }
      
    }
    
    # We are done fitting models! #
    # From here, we use pred_train and pred_val to calculate val and train performance.
    
    # From here, we can calculate the performance for each group
    for (each_sa in c(sa_levels,"Overall")) {
      
      # Get the selection index of the each_sa
      if (each_sa=="Overall") {
        select_train <- rep(TRUE, length(sa_train))
        select_val <- rep(TRUE, length(sa_val))
      } else {
        select_train <- sa_train == each_sa
        select_val <- sa_val == each_sa
      }
      x_train_sa <- x_train[select_train,]
      y_train_sa <- y_train[select_train]
      pred_train_sa <- pred_train[select_train]
      x_val_sa <- x_val[select_val,]
      y_val_sa <- y_val[select_val]
      pred_val_sa <- pred_val[select_val]
      
      
      # Here we start saving metrics #
      
      # ROC-AUC
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["AUROC"]] <- c(roc(y_train_sa, pred_train_sa, levels = c(0, 1), direction = "<")$auc)
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["AUROC"]] <- c(roc(y_val_sa, pred_val_sa, levels = c(0, 1), direction = "<")$auc)
      
      # Calibration slope/intercept
      lin_pred_train_sa <- -1*log((1-pred_train_sa)/pred_train_sa)
      lin_pred_val_sa <- -1*log((1-pred_val_sa)/pred_val_sa)
      # Slope
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["Slope"]] <- glm(
        y_train_sa~lin_pred_train_sa,family="binomial"
      )$coefficients[2]
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["Slope"]] <- glm(
        y_val_sa~lin_pred_val_sa,family="binomial"
      )$coefficients[2]
      # Intercept
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["Intercept"]] <- glm(
        y_train_sa~offset(lin_pred_train_sa)+1,family="binomial"
      )$coefficients[1]
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["Intercept"]] <- glm(
        y_val_sa~offset(lin_pred_val_sa)+1,family="binomial"
      )$coefficients[1]
      
      # OER
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["OER"]] <- mean(y_train_sa)/
        mean(pred_train_sa)
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["OER"]] <- mean(y_val_sa)/
        mean(pred_val_sa)
      
      # NB
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["NB"]] <- calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        y_train_sa,pred_train_sa,new_version = FALSE
      )
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["NB"]] <- calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        y_val_sa,pred_val_sa,new_version = FALSE
      )
      
      # sNB
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["sNB"]] <- 1 - mean(y_train_sa) +diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        y_train_sa,pred_train_sa,new_version = FALSE
      )
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["sNB"]] <- 1 - mean(y_val_sa) +diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        y_val_sa,pred_val_sa,new_version = FALSE
      )
      
      # DCA
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["DCA_x"]] <- dcax_diabetes
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["DCA_x"]] <- dcax_diabetes
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["DCA_y"]] <- calcNetBenefit(
        dcax_diabetes,dcax_diabetes,
        y_train_sa,pred_train_sa,new_version = FALSE
      )
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["DCA_y"]] <- calcNetBenefit(
        dcax_diabetes,dcax_diabetes,
        y_val_sa,pred_val_sa,new_version = FALSE
      )
      
    }
    
    # END each_sa in c()
    
  }
  
  # Write progress to log file
  print(paste("Completed:", bb))
  
  # At the end, return bootstrap metrics
  return(per_bootstrap_metrics)
  
}

# Save bootstrap
saveRDS(bootstrap_optimism,"Output/optimism_diabetes.rds")










#----------------#
#----------------#
####    Train final models ####
#-----------------------------#
#-----------------------------#

# Now our x_train is our x, and our x_val is our x
x_train <- x
x_val <- x
y_train <- y
y_val <- y
sa_train <- df$sa
sa_val <- df$sa

models_dict <- list()

for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  print(model_name)
  
  if (model_name=="LogNoSA") {
    
    fit <- cv.glmnet(
      x = x_train[, !grepl("Ethnicity", colnames(x_train))],
      y = y_train,
      family = "binomial",
      type.measure = "deviance",
      nfolds = 5,
      trace.it = 1,
      parallel = TRUE
    )
    
  } else if (model_name=="LogSingleSA") {
    
    fit <- cv.glmnet(
      x = x_train,
      y = y_train,
      family = "binomial",
      type.measure = "deviance",
      nfolds = 5,
      trace.it = 1,
      parallel = TRUE
    )
    
  } else if (model_name=="LogMultiSA") {
    
    # For each group we fit a separate model
    fit <- list()
    
    # Loop over Townsend
    for (each_sa in sa_levels) {
      
      # Get best weights for SA
      training_weights <- propensity_score_result$propensity[[each_sa]]
      training_weights[sa_train != each_sa] <- pmin(1,training_weights[sa_train != each_sa])
      
      # Train model with weights
      trained_model <- cv.glmnet(
        x = x_train,
        y = y_train,
        family = "binomial",
        weights = training_weights,
        parallel = TRUE,
        type.measure = "deviance",
        nfolds = 5,
        trace.it = 1
      )
      
      # Add model to full list
      fit[[each_sa]] <- trained_model
      
    }
  }
  
  # Save fit to respective model in dict
  models_dict[[model_name]] <- fit
  
}

# Now that they're in the dict, save and reload
saveRDS(models_dict,"Output//models_diabetes")
# Load
models_dict <- readRDS("Output/models_diabetes")

# Now generate predictions
preds_train <- list()
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  print(model_name)
  
  # Load model in
  fit <- models_dict[[model_name]]
  
  if (model_name=="LogNoSA") {
    
    # What is the predictions of the model?
    pred_train <- predict(fit,
                          newx = x_train[, !grepl("Ethnicity", colnames(x_train))],
                          type = "response",s = "lambda.min")[,1]
    
  } else if (model_name=="LogSingleSA") {
    
    # What is the predictions of the model?
    pred_train <- predict(fit,
                          newx = x_train,
                          type = "response",s = "lambda.min")[,1]
    
  } else if (model_name=="LogMultiSA") {
    
    # What is the predictions of the model?
    pred_train <- predict(fit[["Other"]],
                          newx = x_train,
                          type = "response",
                          s = "lambda.min")[,1]
    for (each_sa in sa_levels) {
      pred_train[sa_train==each_sa] <- predict(fit[[each_sa]],
                                               newx = x_train[sa_train==each_sa,],
                                               type = "response",
                                               s = "lambda.min")[,1]
    }
    
  }
  
  # Save predictions
  preds_train[[model_name]] <- pred_train
  
}


#--------------------------------#
#--------------------------------#
####    Calculate performance ####
#--------------------------------#
#--------------------------------#

# Load bootstrap optimism
bootstrap_optimism <- readRDS("Output/optimism_diabetes.rds")

# Save performance dictionnary for all the models
performance_dict <- list()

# For each model, save performance
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  # Load predictions
  pred_train <- preds_train[[model_name]]
  
  print(model_name)
  
  # For each group, including overall
  for (each_sa in c(sa_levels,"Overall")) {
    
    print(each_sa)
    
    # Get the selection index of the each_sa
    if (each_sa=="Overall") {
      select_train <- rep(TRUE, length(sa_train))
    } else {
      select_train <- sa_train == each_sa
    }
    x_train_sa <- x_train[select_train,]
    y_train_sa <- y_train[select_train]
    pred_train_sa <- pred_train[select_train]
    
    
    
    # First, roc
    auc.ci <- ci.auc(roc(y_train_sa, pred_train_sa))
    # correct for optimism
    auc.apparent <- sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Train"]][["AUROC"]]
    )
    auc.reality <- sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Val"]][["AUROC"]]
    )
    auc.optimism <- mean(auc.apparent-auc.reality)
    # Correct auc.ci with auc.optimism
    performance_dict[[model_name]][[each_sa]][["AUROC"]] <- auc.ci - auc.optimism
    
    
    
    # Calibration slope/intercept
    lin_pred_train_sa <- -1*log((1-pred_train_sa)/pred_train_sa)
    
    # Slope
    slope_model <- glm(
      y_train_sa~lin_pred_train_sa,family="binomial"
    )
    slope.ci <- c(confint(slope_model)[2,1],
                  slope_model$coefficients[2],
                  confint(slope_model)[2,2])
    # correct for optimism
    slope.apparent <- sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Train"]][["Slope"]]
    )
    slope.reality <- sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Val"]][["Slope"]]
    )
    slope.optimism <- mean(slope.apparent-slope.reality)
    # Correct
    performance_dict[[model_name]][[each_sa]][["Slope"]] <- slope.ci - slope.optimism
    
    # Intercept
    intercept.model <- glm(
      y_train_sa~offset(lin_pred_train_sa)+1,family="binomial"
    )
    intercept.ci <- c(
      confint(intercept.model)[1],
      intercept.model$coefficients[1],
      confint(intercept.model)[2]
    )
    # correct for optimism
    intercept.apparent <- sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Train"]][["Intercept"]]
    )
    intercept.reality <- sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Val"]][["Intercept"]]
    )
    intercept.optimism <- mean(intercept.apparent-intercept.reality)
    # Correct
    performance_dict[[model_name]][[each_sa]][["Intercept"]] <- intercept.ci - intercept.optimism
    
    
    
    # Now OER
    print("     Calculating OER...")
    oer.centre <- mean(y_train_sa)/mean(pred_train_sa)
    # Get confidence intervals through bootstrapping
    ratio_function <- function(data, indices) {
      d <- data[indices, ]
      mean(d$y_train_sa) / mean(d$pred_train_sa)
    }
    data <- data.frame(y_train_sa, pred_train_sa)
    results <- boot(data, statistic = ratio_function, R = 500,
                    parallel = "multicore", ncpus = 5)
    oer.boot.ci <- boot.ci(results, type = "perc")
    oer.ci <- c(
      oer.boot.ci$percent[length(oer.boot.ci$percent)-1],
      oer.centre,
      oer.boot.ci$percent[length(oer.boot.ci$percent)]
    )
    # correct for optimism
    oer.apparent <- sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Train"]][["OER"]]
    )
    oer.reality <- sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Val"]][["OER"]]
    )
    oer.optimism <- mean(oer.apparent-oer.reality)
    # Correct
    performance_dict[[model_name]][[each_sa]][["OER"]] <- oer.ci - oer.optimism
    
    
    
    # NB
    print("     Calculating NB...")
    nb.centre <- calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y_train_sa,pred_train_sa,new_version = FALSE
    )
    # Get confidence intervals through bootstrapping
    ratio_function <- function(data, indices) {
      d <- data[indices, ]
      calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        d$y_train_sa,d$pred_train_sa,new_version = FALSE
      )
    }
    data <- data.frame(y_train_sa, pred_train_sa)
    results <- boot(data, statistic = ratio_function, R = 500,
                    parallel = "multicore", ncpus = 5)
    nb.boot.ci <- boot.ci(results, type = "perc")
    nb.ci <- c(
      nb.boot.ci$percent[length(nb.boot.ci$percent)-1],
      nb.centre,
      nb.boot.ci$percent[length(nb.boot.ci$percent)]
    )
    # correct for optimism
    nb.apparent <- sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Train"]][["NB"]]
    )
    nb.reality <- sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Val"]][["NB"]]
    )
    nb.optimism <- mean(nb.apparent-nb.reality)
    # Correct
    performance_dict[[model_name]][[each_sa]][["NB"]] <- nb.ci - nb.optimism
    
    
    
    # sNB
    print("     Calculating sNB...")
    snb.centre <- 1 - mean(y_train_sa) + diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y_train_sa,pred_train_sa,new_version = FALSE
    )
    # Get confidence intervals through bootstrapping
    ratio_function <- function(data, indices) {
      d <- data[indices, ]
      1 - mean(d$y_train_sa) + diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        d$y_train_sa,d$pred_train_sa,new_version = FALSE
      )
    }
    data <- data.frame(y_train_sa, pred_train_sa)
    results <- boot(data, statistic = ratio_function, R = 500,
                    parallel = "multicore", ncpus = 5)
    snb.boot.ci <- boot.ci(results, type = "perc")
    snb.ci <- c(
      snb.boot.ci$percent[length(snb.boot.ci$percent)-1],
      snb.centre,
      snb.boot.ci$percent[length(snb.boot.ci$percent)]
    )
    # correct for optimism
    snb.apparent <- sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Train"]][["sNB"]]
    )
    snb.reality <- sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Val"]][["sNB"]]
    )
    snb.optimism <- mean(snb.apparent-snb.reality)
    # Correct
    performance_dict[[model_name]][[each_sa]][["sNB"]] <- snb.ci - snb.optimism
    
    
    
    
    # DCA
    print("     Calculating DCA...")
    performance_dict[[model_name]][[each_sa]][["DCA_x"]] <- dcax_diabetes
    dca.centre <- calcNetBenefit(
      dcax_diabetes,dcax_diabetes,
      y_train_sa,pred_train_sa,new_version = FALSE
    )
    ratio_function <- function(data, indices) {
      d <- data[indices, ]
      calcNetBenefit(
        dcax_diabetes,dcax_diabetes,
        d$y_train_sa,d$pred_train_sa,new_version = FALSE
      )
    }
    data <- data.frame(y_train_sa, pred_train_sa)
    results <- boot(data, statistic = ratio_function, R = 500,
                    parallel = "multicore", ncpus = 5)
    
    dca.boot.ci <-apply(results$t, 2, function(x) quantile(x, probs = c(0.05, 0.95)))
    dca.ci <- data.frame(
      DCA.Low = dca.boot.ci[1,],
      DCA = dca.centre,
      DCA.High = dca.boot.ci[2,]
    )
    # Calculate optimism
    dca.optimism <- rowMeans(sapply(
      bootstrap_optimism,
      function(x)
        x[[model_name]][[each_sa]][["Train"]][["DCA_y"]] - 
        x[[model_name]][[each_sa]][["Val"]][["DCA_y"]]
    ))
    # Correct
    performance_dict[[model_name]][[each_sa]][["DCA"]] <- sweep(dca.ci, 1, dca.optimism, FUN="-")
    
    
    
  }
  
}

# Save performance
saveRDS(performance_dict,"Output/performance_diabetes.rds")


#------------------------------------------------------------------#
#------------------------------------------------------------------#
####    Calculate performance CI a la carte ####

# Open optimism
bootstrap_optimism <- readRDS("Output/optimism_diabetes.rds")

# Open Ethnicity info
sa <- df$sa

####        Difference between Asian and White from treat.no.one to LogSingleSA ####
print("Calculating difference White Asian from Treat.No.One to LogSingleSA...")
snb.treat.no.one.centre <- 10000*((
  1 - mean(y[sa=="Asian"]) + 
    diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y[sa=="Asian"],0*preds_train[["LogSingleSA"]][sa=="Asian"],
      new_version = FALSE
    )
) - (
  1 - mean(y[sa=="White"]) + 
    diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y[sa=="White"],0*preds_train[["LogSingleSA"]][sa=="White"],
      new_version = FALSE
    )
))
snb.LogSingleSA.centre <- 10000*((
  1 - mean(y[sa=="Asian"]) + 
    diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y[sa=="Asian"],preds_train[["LogSingleSA"]][sa=="Asian"],
      new_version = FALSE
    )
) - (
  1 - mean(y[sa=="White"]) + 
    diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y[sa=="White"],preds_train[["LogSingleSA"]][sa=="White"],
      new_version = FALSE
    )
))
snb.centre <- snb.treat.no.one.centre - snb.LogSingleSA.centre

# For CI, define bootstrap function
ratio_function <- function(data, indices) {
  d <- data[indices, ]
  snb.treat.no.one.centre <- 10000*((
    1 - mean(d$y[d$sa=="Asian"]) + 
      diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        d$y[d$sa=="Asian"],0*d$LogSingleSA[d$sa=="Asian"],
        new_version = FALSE
      )
  ) - (
    1 - mean(d$y[d$sa=="White"]) + 
      diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        d$y[d$sa=="White"],0*d$LogSingleSA[d$sa=="White"],
        new_version = FALSE
      )
  ))
  snb.LogSingleSA.centre <- 10000*((
    1 - mean(d$y[d$sa=="Asian"]) + 
      diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        d$y[d$sa=="Asian"],d$LogSingleSA[d$sa=="Asian"],
        new_version = FALSE
      )
  ) - (
    1 - mean(d$y[d$sa=="White"]) + 
      diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        d$y[d$sa=="White"],d$LogSingleSA[d$sa=="White"],
        new_version = FALSE
      )
  ))
  snb.centre <- snb.treat.no.one.centre - snb.LogSingleSA.centre
  return(snb.centre)
}
data <- data.frame(y=y,sa=sa,LogSingleSA=preds_train[["LogSingleSA"]])
results <- boot(data, statistic = ratio_function, R = 500,
                parallel = "multicore", ncpus = 5)
snb.boot.ci <- boot.ci(results, type = "perc")
snb.ci <- c(
  snb.boot.ci$percent[length(snb.boot.ci$percent)-1],
  snb.centre,
  snb.boot.ci$percent[length(snb.boot.ci$percent)]
)

# correct for optimis
snb.LogSingleSA.apparent <- - 10000*(sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["Asian"]][["Train"]][["sNB"]]
) - sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["White"]][["Train"]][["sNB"]]
))
snb.LogSingleSA.reality <- - 10000*(sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["Asian"]][["Val"]][["sNB"]]
) - sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["White"]][["Val"]][["sNB"]]
))
snb.optimism <- mean(snb.LogSingleSA.apparent-snb.LogSingleSA.reality)
# Do correction
snb.ci <- snb.ci - snb.optimism 
# Show the difference:
print("Change in Asian-White from TreatNoOne to LogSingleSA")
print(snb.ci)







####        Difference between AsianBlackOther and White from treat.no.one to LogSingleSA ####
print("Calculating difference White Asian from Treat.No.One to LogSingleSA...")
snb.treat.no.one.centre <- ((
  1 - mean(y[sa=="Asian"]) + 
    diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y[sa=="Asian"],0*preds_train[["LogSingleSA"]][sa=="Asian"],
      new_version = FALSE
    )
  ) + (
  1 - mean(y[sa=="Black"]) + 
    diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y[sa=="Black"],0*preds_train[["LogSingleSA"]][sa=="Black"],
      new_version = FALSE
    )
  ) +
  (
  1 - mean(y[sa=="Other"]) + 
    diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y[sa=="Other"],0*preds_train[["LogSingleSA"]][sa=="Other"],
      new_version = FALSE
    )
  ) - 3*(
  1 - mean(y[sa=="White"]) + 
    diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y[sa=="White"],0*preds_train[["LogSingleSA"]][sa=="White"],
      new_version = FALSE
    )
))/3*10000


snb.LogSingleSA.centre <- ((
  1 - mean(y[sa=="Asian"]) + 
    diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y[sa=="Asian"],preds_train[["LogSingleSA"]][sa=="Asian"],
      new_version = FALSE
    )
) + (
  1 - mean(y[sa=="Black"]) + 
    diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y[sa=="Black"],preds_train[["LogSingleSA"]][sa=="Black"],
      new_version = FALSE
    )
) +
  (
    1 - mean(y[sa=="Other"]) + 
      diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        y[sa=="Other"],preds_train[["LogSingleSA"]][sa=="Other"],
        new_version = FALSE
      )
  ) - 3*(
    1 - mean(y[sa=="White"]) + 
      diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        y[sa=="White"],preds_train[["LogSingleSA"]][sa=="White"],
        new_version = FALSE
      )
  ))/3*10000
snb.centre <- snb.treat.no.one.centre - snb.LogSingleSA.centre

# For CI, define bootstrap function
ratio_function <- function(data, indices) {
  d <- data[indices, ]
  snb.treat.no.one.centre <- ((
    1 - mean(d$y[d$sa=="Asian"]) + 
      diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        d$y[d$sa=="Asian"],0*d$LogSingleSA[d$sa=="Asian"],
        new_version = FALSE
      )
  ) + (
    1 - mean(d$y[d$sa=="Black"]) + 
      diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        d$y[d$sa=="Black"],0*d$LogSingleSA[d$sa=="Black"],
        new_version = FALSE
      )
  ) +
    (
      1 - mean(d$y[d$sa=="Other"]) + 
        diabetes_effect * calcNetBenefit(
          diabetes_threshold,diabetes_threshold,
          d$y[d$sa=="Other"],0*d$LogSingleSA[d$sa=="Other"],
          new_version = FALSE
        )
    ) - 3*(
      1 - mean(d$y[d$sa=="White"]) + 
        diabetes_effect * calcNetBenefit(
          diabetes_threshold,diabetes_threshold,
          d$y[d$sa=="White"],0*d$LogSingleSA[d$sa=="White"],
          new_version = FALSE
        )
    ))/3*10000
  
  
  snb.LogSingleSA.centre <- ((
    1 - mean(d$y[d$sa=="Asian"]) + 
      diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        d$y[d$sa=="Asian"],d$LogSingleSA[d$sa=="Asian"],
        new_version = FALSE
      )
  ) + (
    1 - mean(d$y[d$sa=="Black"]) + 
      diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        d$y[d$sa=="Black"],d$LogSingleSA[d$sa=="Black"],
        new_version = FALSE
      )
  ) +
    (
      1 - mean(d$y[d$sa=="Other"]) + 
        diabetes_effect * calcNetBenefit(
          diabetes_threshold,diabetes_threshold,
          d$y[d$sa=="Other"],d$LogSingleSA[d$sa=="Other"],
          new_version = FALSE
        )
    ) - 3*(
      1 - mean(d$y[d$sa=="White"]) + 
        diabetes_effect * calcNetBenefit(
          diabetes_threshold,diabetes_threshold,
          d$y[d$sa=="White"],d$LogSingleSA[d$sa=="White"],
          new_version = FALSE
        )
    ))/3*10000
  snb.centre <- snb.treat.no.one.centre - snb.LogSingleSA.centre  
  return(snb.centre)
}
data <- data.frame(y=y,sa=sa,LogSingleSA=preds_train[["LogSingleSA"]])
results <- boot(data, statistic = ratio_function, R = 500,
                parallel = "multicore", ncpus = 5)
snb.boot.ci <- boot.ci(results, type = "perc")
snb.ci <- c(
  snb.boot.ci$percent[length(snb.boot.ci$percent)-1],
  snb.centre,
  snb.boot.ci$percent[length(snb.boot.ci$percent)]
)

# correct for optimis
snb.LogSingleSA.apparent <-  - 10000*(sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["Asian"]][["Train"]][["sNB"]]
) + sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["Black"]][["Train"]][["sNB"]]
) + sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["Other"]][["Train"]][["sNB"]]
) - 3*sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["White"]][["Train"]][["sNB"]]
))/3
snb.LogSingleSA.reality <- - 10000*(sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["Asian"]][["Val"]][["sNB"]]
) + sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["Black"]][["Val"]][["sNB"]]
) + sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["Other"]][["Val"]][["sNB"]]
) - 3*sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["White"]][["Val"]][["sNB"]]
))/3
snb.optimism <- mean(snb.LogSingleSA.apparent-snb.LogSingleSA.reality)
# Do correction
snb.ci <- snb.ci - snb.optimism 
# Show the difference:
print("Change in Asian-White from TreatNoOne to LogSingleSA")
print(snb.ci)



####        Difference in Asian between LogNoSA to LogSingleSA ####
print("Calculating difference  Asian from LogSingleSA to LogMultiSA...")
snb.treat.no.one.centre <- 10000*((
  1 - mean(y[sa=="Asian"]) + 
    diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y[sa=="Asian"],preds_train[["LogNoSA"]][sa=="Asian"],
      new_version = FALSE
    )
))
snb.LogSingleSA.centre <- 10000*((
  1 - mean(y[sa=="Asian"]) + 
    diabetes_effect * calcNetBenefit(
      diabetes_threshold,diabetes_threshold,
      y[sa=="Asian"],preds_train[["LogSingleSA"]][sa=="Asian"],
      new_version = FALSE
    )
))
snb.centre <- snb.treat.no.one.centre - snb.LogSingleSA.centre

# For CI, define bootstrap function
ratio_function <- function(data, indices) {
  d <- data[indices, ]
  snb.treat.no.one.centre <- 10000*((
    1 - mean(d$y[d$sa=="Asian"]) + 
      diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        d$y[d$sa=="Asian"],d$LogNoSA[d$sa=="Asian"],
        new_version = FALSE
      )
  ))
  snb.LogSingleSA.centre <- 10000*((
    1 - mean(d$y[d$sa=="Asian"]) + 
      diabetes_effect * calcNetBenefit(
        diabetes_threshold,diabetes_threshold,
        d$y[d$sa=="Asian"],d$LogSingleSA[d$sa=="Asian"],
        new_version = FALSE
      )
  ))
  snb.centre <- snb.treat.no.one.centre - snb.LogSingleSA.centre
  return(snb.centre)
}
data <- data.frame(y=y,sa=sa,
                   LogNoSA=preds_train[["LogNoSA"]],
                   LogSingleSA=preds_train[["LogSingleSA"]])
results <- boot(data, statistic = ratio_function, R = 500,
                parallel = "multicore", ncpus = 5)
snb.boot.ci <- boot.ci(results, type = "perc")
snb.ci <- c(
  snb.boot.ci$percent[length(snb.boot.ci$percent)-1],
  snb.centre,
  snb.boot.ci$percent[length(snb.boot.ci$percent)]
)

# correct for optimis
snb.LogSingleSA.apparent <-  10000*(sapply(
  bootstrap_optimism,
  function(x)
    x[["LogNoSA"]][["Asian"]][["Train"]][["sNB"]]
)) - 10000*(sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["Asian"]][["Train"]][["sNB"]]
))
snb.LogSingleSA.reality <- 10000*(sapply(
  bootstrap_optimism,
  function(x)
    x[["LogNoSA"]][["Asian"]][["Val"]][["sNB"]]
)) - 10000*(sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["Asian"]][["Val"]][["sNB"]]
))
snb.optimism <- mean(snb.LogSingleSA.apparent-snb.LogSingleSA.reality)
# Do correction
snb.ci <- snb.ci - snb.optimism 
# Show the difference:
print("Change in Asian from LogNoSA to LogSingleSA")
print(snb.ci)


####  LUNG CANCER ####
#----------------#
#----------------#

#----------------------------#
####    Open data and impute ####
#----------------------------#

df <- read.csv(file = paste0("Data/lungcancer_df_randomsplit.csv"))

# X is rownames
rownames(df) <- df$X
df <- df %>% select(-X)

# Turn factor values into factors
df$Ethnicity <- as.factor(df$Ethnicity)
df$Education <- as.factor(df$Education)
df$Smoker <- as.factor(df$Smoker)

# Remove values not used in prediction or y or sa
df[c("Censored","Started","Death","Centre","is_val")] <- NULL

# Change to factors strings and binary variables
df[, sapply(df, is.character)] <- lapply(df[, sapply(df, is.character)], as.factor)

# Get quintiles of Townsend
quintiles_townsend <- quantile(df$Townsend,probs = c(0.2,0.4,0.6,0.8),na.rm=TRUE)

# Exclude 'y' from being used as a predictor in MICE
predictorMatrix <- mice::make.predictorMatrix(df)
predictorMatrix[, "y"] <- 0   # Do not use 'y' as a predictor for any other variable
predictorMatrix["y", ] <- 0   # Do not impute 'y' (it has no missing values)

# Impute data using mice
df_mice <- mice(df,maxit = MICE_MAXIT, m = MICE_M, predictorMatrix = predictorMatrix)
df <- complete(df_mice, action = "stacked")
df <- na.omit(df)

# Define sa according to ethnicity
df$sa <- define_lungcancer_sa(df$Townsend,quintiles_townsend)

# Save and load
saveRDS(df, "Data/lungcancer_df_randomsplit_imputed.rds")
df <- readRDS("Data/lungcancer_df_randomsplit_imputed.rds")

#--------------------#
####    Prepare data ####
#--------------------#

# Turn validation and training data into right format
x_and_y <- prepareForGLMNET(df)
x <- x_and_y[["x"]]
y <- x_and_y[["y"]]

# What are the sa levels?
sa_levels <- levels(df$sa)

#----------------------------#
####    Calculate propensity ####
#----------------------------#

# Find out propensity_score
propensity_train_df <- df %>% select(-sa,-Townsend)
propensity_score_result <- calculate_sensitive_propensity_score(
  propensity_train_df,df$sa,sa_levels, s = "lambda.min"
)


#------------------------------#
#------------------------------#
####    Calculate the Optimism ####
#------------------------------#
#------------------------------#

# Here is the dictionary we will be recording in the metrics per bootstrap
bootstrap_optimism <- list()

# Here, our validation data is x and y
x_val <- x
y_val <- y
sa_val <- df$sa 

set.seed(RANDOM_SEED)
bootstrap_optimism <- foreach(bb=1:NUM_BOOT) %dopar% {
  
  # Here we save everything
  per_bootstrap_metrics <- list()
  
  print(glue("Bootstrap Number {bb}"))
  
  # Select the new data
  bootstrap_index <- sample(nrow(df), replace = TRUE)
  
  # Training data
  x_train <- x_val[bootstrap_index,]
  y_train <- y_val[bootstrap_index]
  sa_train <- sa_val[bootstrap_index]
  
  for (model_name in c("LogNoSA", "LogSingleSA", "LogMultiSA")) {
    
    # Train the model
    if (model_name=="LogNoSA") {
      
      fit <- cv.glmnet(
        x = x_train[, !(colnames(x) %in% c("Townsend"))],
        y = y_train,
        family = "binomial",
        type.measure = "deviance",
        nfolds = 5,
        trace.it = 0,
        parallel = FALSE
      )
      # What is the predictions of the model?
      pred_train <- predict(fit,
                            newx = x_train[, !(colnames(x) %in% c("Townsend"))],
                            type = "response",s = "lambda.min")[,1]
      pred_val <- predict(fit,
                          newx = x[, !(colnames(x) %in% c("Townsend"))],
                          type = "response",s = "lambda.min")[,1]
      
    } else if (model_name=="LogSingleSA") {
      
      fit <- cv.glmnet(
        x = x_train,
        y = y_train,
        family = "binomial",
        type.measure = "deviance",
        nfolds = 5,
        trace.it = 0,
        parallel = FALSE
      )
      # What is the predictions of the model?
      pred_train <- predict(fit,
                            newx = x_train,
                            type = "response",s = "lambda.min")[,1]
      pred_val <- predict(fit,
                          newx = x,
                          type = "response",s = "lambda.min")[,1]
      
    } else if (model_name=="LogMultiSA") {
      
      # For each group we fit a separate model
      model_per_group <- list()
      
      # Loop over Townsend
      for (each_sa in sa_levels) {
        
        # Get best weights for SA
        training_weights <- propensity_score_result$propensity[[each_sa]][bootstrap_index]
        training_weights[sa_train != each_sa] <- pmin(1,training_weights[sa_train != each_sa])
        
        # Train model with weights
        trained_model <- cv.glmnet(
          x = x_train,
          y = y_train,
          family = "binomial",
          weights = training_weights,
          parallel = FALSE,
          type.measure = "deviance",
          nfolds = 5,
          trace.it = 0
        )
        
        # Add model to full list
        model_per_group[[each_sa]] <- trained_model
        
      }
      
      # What is the predictions of the model?
      pred_train <- predict(model_per_group[["1"]],
                            newx = x_train,
                            type = "response",
                            s = "lambda.min")[,1]
      pred_val <- predict(model_per_group[["1"]],
                          newx = x,
                          type = "response",
                          s = "lambda.min")[,1]
      for (each_sa in sa_levels) {
        pred_train[sa_train==each_sa] <- predict(model_per_group[[each_sa]],
                                                 newx = x_train[sa_train==each_sa,],
                                                 type = "response",
                                                 s = "lambda.min")[,1]
        pred_val[sa_val==each_sa] <- predict(model_per_group[[each_sa]],
                                             newx = x_val[sa_val==each_sa,],
                                             type = "response",
                                             s = "lambda.min")[,1]
      }
      
    }
    
    # We are done fitting models! #
    # From here, we use pred_train and pred_val to calculate val and train performance.
    
    # From here, we can calculate the performance for each group
    for (each_sa in c(sa_levels,"Overall")) {
      
      # Get the selection index of the each_sa
      if (each_sa=="Overall") {
        select_train <- rep(TRUE, length(sa_train))
        select_val <- rep(TRUE, length(sa_val))
      } else {
        select_train <- sa_train == each_sa
        select_val <- sa_val == each_sa
      }
      x_train_sa <- x_train[select_train,]
      y_train_sa <- y_train[select_train]
      pred_train_sa <- pred_train[select_train]
      x_val_sa <- x_val[select_val,]
      y_val_sa <- y_val[select_val]
      pred_val_sa <- pred_val[select_val]
      
      
      # Here we start saving metrics #
      
      # ROC-AUC
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["AUROC"]] <- c(roc(y_train_sa, pred_train_sa, levels = c(0, 1), direction = "<")$auc)
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["AUROC"]] <- c(roc(y_val_sa, pred_val_sa, levels = c(0, 1), direction = "<")$auc)
      
      # Calibration slope/intercept
      lin_pred_train_sa <- -1*log((1-pred_train_sa)/pred_train_sa)
      lin_pred_val_sa <- -1*log((1-pred_val_sa)/pred_val_sa)
      # Slope
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["Slope"]] <- glm(
        y_train_sa~lin_pred_train_sa,family="binomial"
      )$coefficients[2]
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["Slope"]] <- glm(
        y_val_sa~lin_pred_val_sa,family="binomial"
      )$coefficients[2]
      # Intercept
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["Intercept"]] <- glm(
        y_train_sa~offset(lin_pred_train_sa)+1,family="binomial"
      )$coefficients[1]
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["Intercept"]] <- glm(
        y_val_sa~offset(lin_pred_val_sa)+1,family="binomial"
      )$coefficients[1]
      
      # OER
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["OER"]] <- mean(y_train_sa)/
        mean(pred_train_sa)
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["OER"]] <- mean(y_val_sa)/
        mean(pred_val_sa)
      
      # NB
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["NB"]] <- calcNetBenefit(
        lung_cancer_threshold,lung_cancer_threshold,
        y_train_sa,pred_train_sa,new_version = FALSE
      )
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["NB"]] <- calcNetBenefit(
        lung_cancer_threshold,lung_cancer_threshold,
        y_val_sa,pred_val_sa,new_version = FALSE
      )
      
      # sNB
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["sNB"]] <- 1 - mean(y_train_sa) +lung_cancer_effect * calcNetBenefit(
        lung_cancer_threshold,lung_cancer_threshold,
        y_train_sa,pred_train_sa,new_version = FALSE
      )
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["sNB"]] <- 1 - mean(y_val_sa) +lung_cancer_effect * calcNetBenefit(
        lung_cancer_threshold,lung_cancer_threshold,
        y_val_sa,pred_val_sa,new_version = FALSE
      )
      
      # DCA
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["DCA_x"]] <- dcax_lungcancer
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["DCA_x"]] <- dcax_lungcancer
      per_bootstrap_metrics[[model_name]][[each_sa]][["Train"]][["DCA_y"]] <- calcNetBenefit(
        dcax_lungcancer,dcax_lungcancer,
        y_train_sa,pred_train_sa,new_version = FALSE
      )
      per_bootstrap_metrics[[model_name]][[each_sa]][["Val"]][["DCA_y"]] <- calcNetBenefit(
        dcax_lungcancer,dcax_lungcancer,
        y_val_sa,pred_val_sa,new_version = FALSE
      )
      
    }
    
    # END each_sa in c()
    
  }
  
  print(glue("Completed {bb}"))
  
  # At the end, return bootstrap metrics
  return(per_bootstrap_metrics)
  
}

print("Optimism Calculation Done")

# Save bootstrap
saveRDS(bootstrap_optimism,"Output/optimism_lung_cancer.rds")


#-----------------------------#
#-----------------------------#
####    Train final models ####
#-----------------------------#
#-----------------------------#

# Now our x_train is our x, and our x_val is our x
x_train <- x
x_val <- x
y_train <- y
y_val <- y
sa_train <- df$sa
sa_val <- df$sa

models_dict <- list()

for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  print(model_name)
  
  if (model_name=="LogNoSA") {
    
    fit <- cv.glmnet(
      x = x_train[, !(colnames(x) %in% c("Townsend"))],
      y = y_train,
      family = "binomial",
      type.measure = "deviance",
      nfolds = 5,
      trace.it = 1,
      parallel = TRUE
    )
    
  } else if (model_name=="LogSingleSA") {
    
    fit <- cv.glmnet(
      x = x_train,
      y = y_train,
      family = "binomial",
      type.measure = "deviance",
      nfolds = 5,
      trace.it = 1,
      parallel = TRUE
    )
    
  } else if (model_name=="LogMultiSA") {
    
    # For each group we fit a separate model
    fit <- list()
    
    # Loop over Townsend
    for (each_sa in sa_levels) {
      
      # Get best weights for SA
      training_weights <- propensity_score_result$propensity[[each_sa]]
      training_weights[sa_train != each_sa] <- pmin(1,training_weights[sa_train != each_sa])
      
      # Train model with weights
      trained_model <- cv.glmnet(
        x = x_train,
        y = y_train,
        family = "binomial",
        weights = training_weights,
        parallel = TRUE,
        type.measure = "deviance",
        nfolds = 5,
        trace.it = 1
      )
      
      # Add model to full list
      fit[[each_sa]] <- trained_model
      
    }
  }
  
  # Save fit to respective model in dict
  models_dict[[model_name]] <- fit
    
  }

# Now that they're in the dict, save and reload
saveRDS(models_dict,"Output/models_lungcancer")

# Load
models_dict <- readRDS("Output//models_lungcancer")
x_train <- x
x_val <- x
y_train <- y
y_val <- y
sa_train <- df$sa
sa_val <- df$sa

# Now generate predictions
preds_train <- list()
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  print(model_name)
  
  # Load model in
  fit <- models_dict[[model_name]]
  
  if (model_name=="LogNoSA") {

    # What is the predictions of the model?
    pred_train <- predict(fit,
                          newx = x_train[, !(colnames(x) %in% c("Townsend"))],
                          type = "response",s = "lambda.min")[,1]
    
  } else if (model_name=="LogSingleSA") {
    
    # What is the predictions of the model?
    pred_train <- predict(fit,
                          newx = x_train,
                          type = "response",s = "lambda.min")[,1]
    
  } else if (model_name=="LogMultiSA") {
    
    # For each group we fit a separate model
    model_per_group <- list()
    
    # Loop over Townsend
    
    # What is the predictions of the model?
    pred_train <- predict(fit[["1"]],
                          newx = x_train,
                          type = "response",
                          s = "lambda.min")[,1]
    for (each_sa in sa_levels) {
      pred_train[sa_train==each_sa] <- predict(fit[[each_sa]],
                                               newx = x_train[sa_train==each_sa,],
                                               type = "response",
                                               s = "lambda.min")[,1]
    }
    
  }
  
  # Save predictions
  preds_train[[model_name]] <- pred_train
  
}

#--------------------------------#
#--------------------------------#
####      Calculate performance CI a la carte ####

# Open optimism
bootstrap_optimism <- readRDS("Output/optimism_lung_cancer.rds")

# Open Townsend info
sa <- df$sa

####        Difference in Q4 between Treat.No.One to LogSingleSA ####
print("Calculating difference Q4 from TreatNoOne to LogSingleSA")
snb.treat.no.one.centre <- 10000*((
  1 - mean(y[sa=="3"]) + 
    lung_cancer_effect * calcNetBenefit(
      lung_cancer_threshold,lung_cancer_threshold,
      y[sa=="3"],0*preds_train[["LogNoSA"]][sa=="3"],
      new_version = FALSE
    )
))
snb.LogSingleSA.centre <- 10000*((
  1 - mean(y[sa=="3"]) + 
    lung_cancer_effect * calcNetBenefit(
      lung_cancer_threshold,lung_cancer_threshold,
      y[sa=="3"],preds_train[["LogSingleSA"]][sa=="3"],
      new_version = FALSE
    )
))
snb.centre <- snb.treat.no.one.centre - snb.LogSingleSA.centre

# For CI, define bootstrap function
ratio_function <- function(data, indices) {
  d <- data[indices, ]
  snb.treat.no.one.centre <- 10000*((
    1 - mean(d$y[d$sa=="3"]) + 
      lung_cancer_effect * calcNetBenefit(
        lung_cancer_threshold,lung_cancer_threshold,
        d$y[d$sa=="3"],0*d$LogNoSA[d$sa=="3"],
        new_version = FALSE
      )
  ))
  snb.LogSingleSA.centre <- 10000*((
    1 - mean(d$y[d$sa=="3"]) + 
      lung_cancer_effect * calcNetBenefit(
        lung_cancer_threshold,lung_cancer_threshold,
        d$y[d$sa=="3"],d$LogSingleSA[d$sa=="3"],
        new_version = FALSE
      )
  ))
  snb.centre <- snb.treat.no.one.centre - snb.LogSingleSA.centre
  return(snb.centre)
}
data <- data.frame(y=y,sa=sa,
                   LogNoSA=preds_train[["LogNoSA"]],
                   LogSingleSA=preds_train[["LogSingleSA"]])
results <- boot(data, statistic = ratio_function, R = 500,
                parallel = "multicore", ncpus = 5)
snb.boot.ci <- boot.ci(results, type = "perc")
snb.ci <- c(
  snb.boot.ci$percent[length(snb.boot.ci$percent)-1],
  snb.centre,
  snb.boot.ci$percent[length(snb.boot.ci$percent)]
)

# correct for optimis
snb.LogSingleSA.apparent <- - 10000*(sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["3"]][["Train"]][["sNB"]]
))
snb.LogSingleSA.reality <- - 10000*(sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["3"]][["Val"]][["sNB"]]
))
snb.optimism <- mean(snb.LogSingleSA.apparent-snb.LogSingleSA.reality)
# Do correction
snb.ci <- snb.ci - snb.optimism 
# Show the difference:
print("Change in 3 from TreatNoOne to LogSingleSA")
print(snb.ci)



####        Difference in Q5 between Treat.No.One to LogSingleSA ####
print("Calculating difference Q5 from TreatNoOne to LogSingleSA")
snb.treat.no.one.centre <- 10000*((
  1 - mean(y[sa=="4"]) + 
    lung_cancer_effect * calcNetBenefit(
      lung_cancer_threshold,lung_cancer_threshold,
      y[sa=="4"],0*preds_train[["LogNoSA"]][sa=="4"],
      new_version = FALSE
    )
))
snb.LogSingleSA.centre <- 10000*((
  1 - mean(y[sa=="4"]) + 
    lung_cancer_effect * calcNetBenefit(
      lung_cancer_threshold,lung_cancer_threshold,
      y[sa=="4"],preds_train[["LogSingleSA"]][sa=="4"],
      new_version = FALSE
    )
))
snb.centre <- snb.treat.no.one.centre - snb.LogSingleSA.centre

# For CI, define bootstrap function
ratio_function <- function(data, indices) {
  d <- data[indices, ]
  snb.treat.no.one.centre <- 10000*((
    1 - mean(d$y[d$sa=="4"]) + 
      lung_cancer_effect * calcNetBenefit(
        lung_cancer_threshold,lung_cancer_threshold,
        d$y[d$sa=="4"],0*d$LogNoSA[d$sa=="4"],
        new_version = FALSE
      )
  ))
  snb.LogSingleSA.centre <- 10000*((
    1 - mean(d$y[d$sa=="4"]) + 
      lung_cancer_effect * calcNetBenefit(
        lung_cancer_threshold,lung_cancer_threshold,
        d$y[d$sa=="4"],d$LogSingleSA[d$sa=="4"],
        new_version = FALSE
      )
  ))
  snb.centre <- snb.treat.no.one.centre - snb.LogSingleSA.centre
  return(snb.centre)
}
data <- data.frame(y=y,sa=sa,
                   LogNoSA=preds_train[["LogNoSA"]],
                   LogSingleSA=preds_train[["LogSingleSA"]])
results <- boot(data, statistic = ratio_function, R = 500,
                parallel = "multicore", ncpus = 5)
snb.boot.ci <- boot.ci(results, type = "perc")
snb.ci <- c(
  snb.boot.ci$percent[length(snb.boot.ci$percent)-1],
  snb.centre,
  snb.boot.ci$percent[length(snb.boot.ci$percent)]
)

# correct for optimis
snb.LogSingleSA.apparent <- - 10000*(sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["4"]][["Train"]][["sNB"]]
))
snb.LogSingleSA.reality <- - 10000*(sapply(
  bootstrap_optimism,
  function(x)
    x[["LogSingleSA"]][["4"]][["Val"]][["sNB"]]
))
snb.optimism <- mean(snb.LogSingleSA.apparent-snb.LogSingleSA.reality)
# Do correction
snb.ci <- snb.ci - snb.optimism 
# Show the difference:
print("Change in 4 from TreatNoOne to LogSingleSA")
print(snb.ci)



####    Get Pareto Front corrected for Optimism ####

# Currently this only 'works' on real data.
# As it is specific to the incidence and predictive nature of the
# real lung cancer example.


# What are the possible thresholds explored?
pareto_thresholds <- seq(0.015,0.1,0.0005)

# Define the two groups
top_quintile <- sa_train=="4"

# First, for each model, we identify the Pareto optimal for each constrain
pareto_dict <- list()

# For each model
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  # For each constraint
  for (constraint in c("NoCon","3Perc","1Perc")) {
    
    # How many people can we at most screen?
    con <- case_when(
      constraint == "NoCon" ~ 1*nrow(x_train),
      constraint == "3Perc" ~ 0.03*nrow(x_train),
      constraint == "1Perc" ~ 0.01*nrow(x_train)
    )
    
    # Two dataframes, TopDepriv and Rest
    topdepriv_df <- data.frame(TopDeprivThreshold = pareto_thresholds)
    rest_df <- data.frame(RestThreshold = pareto_thresholds)
    
    # Calculate NB
    topdepriv_df$TopDeprivNB <- 1 - mean(y_train[top_quintile==1]) + 
      lung_cancer_effect * calcNetBenefit(
      rep(0.015,length(pareto_thresholds)),
      pareto_thresholds,
      y_train[top_quintile==1],
      preds_train[[model_name]][top_quintile==1],
      new_version = FALSE
    )
    rest_df$RestNB <- 1 - mean(y_train[top_quintile!=-7]) + 
      lung_cancer_effect * calcNetBenefit(
        rep(0.015,length(pareto_thresholds)),
        pareto_thresholds,
        y_train[top_quintile!=-7],
        preds_train[[model_name]][top_quintile!=-7],
        new_version = FALSE
      )

    # Get number Screened
    topdepriv_df$TopDeprivScreened <- sapply(
      pareto_thresholds,
      function(threshold) 
        sum(preds_train[[model_name]][top_quintile == 1] > threshold))
    rest_df$RestScreened <- sapply(
      pareto_thresholds,
      function(threshold) 
        sum(preds_train[[model_name]][top_quintile!=7] > threshold))
    
    # Combine two dfs
    pareto_model_constraint_df <- merge(topdepriv_df, rest_df, by = NULL)
    
    # How many people are there in the data?
    pareto_model_constraint_df$TotalNum <- length(top_quintile)
    
    # How many people have we screened
    pareto_model_constraint_df$TotalScreened <- 
      pareto_model_constraint_df$TopDeprivScreened +
      
      pareto_model_constraint_df$RestScreened
    
    # Remove Rows for which we exceed allowed number
    pareto_model_constraint_df <- pareto_model_constraint_df[pareto_model_constraint_df$TotalScreened<=con,]
  
    # Limit to Pareto and save
    pareto_front <- rownames(get_frontier(pareto_model_constraint_df,
                                    x=TopDeprivNB,y=RestNB,quadrant = "top.right"))
    pareto_model_constraint_df <- pareto_model_constraint_df[pareto_front,]
    
    # Save
    pareto_dict[[model_name]][[constraint]] <- pareto_model_constraint_df
    }
  
}


####      Optimism correction for Pareto ####
# Now we have to calculate the optimism for 500 bootstraps per model

# Here, our validation data is x and y
x_val <- x
y_val <- y
sa_val <- df$sa 

set.seed(RANDOM_SEED)
pareto_optimism <- foreach(bb=1:NUM_BOOT) %dopar% {
  
  # Here we save everything
  per_bootstrap_metrics <- list()
  
  print(glue("Bootstrap Number {bb}"))
  
  # Select the new data
  bootstrap_index <- sample(nrow(df), replace = TRUE)
  
  # Training data
  x_train <- x_val[bootstrap_index,]
  y_train <- y_val[bootstrap_index]
  sa_train <- sa_val[bootstrap_index]
  
  for (model_name in c("LogNoSA", "LogSingleSA", "LogMultiSA")) {
    
    # Train the model
    if (model_name=="LogNoSA") {
      
      fit <- cv.glmnet(
        x = x_train[, !(colnames(x) %in% c("Townsend"))],
        y = y_train,
        family = "binomial",
        type.measure = "deviance",
        nfolds = 5,
        trace.it = 0,
        parallel = FALSE
      )
      # What is the predictions of the model?
      pred_train <- predict(fit,
                            newx = x_train[, !(colnames(x) %in% c("Townsend"))],
                            type = "response",s = "lambda.min")[,1]
      pred_val <- predict(fit,
                          newx = x[, !(colnames(x) %in% c("Townsend"))],
                          type = "response",s = "lambda.min")[,1]
      
    } else if (model_name=="LogSingleSA") {
      
      fit <- cv.glmnet(
        x = x_train,
        y = y_train,
        family = "binomial",
        type.measure = "deviance",
        nfolds = 5,
        trace.it = 0,
        parallel = FALSE
      )
      # What is the predictions of the model?
      pred_train <- predict(fit,
                            newx = x_train,
                            type = "response",s = "lambda.min")[,1]
      pred_val <- predict(fit,
                          newx = x,
                          type = "response",s = "lambda.min")[,1]
      
    } else if (model_name=="LogMultiSA") {
      
      # For each group we fit a separate model
      model_per_group <- list()
      
      # Loop over Townsend
      for (each_sa in sa_levels) {
        
        # Get best weights for SA
        training_weights <- propensity_score_result$propensity[[each_sa]][bootstrap_index]
        training_weights[sa_train != each_sa] <- pmin(1,training_weights[sa_train != each_sa])
        
        # Train model with weights
        trained_model <- cv.glmnet(
          x = x_train,
          y = y_train,
          family = "binomial",
          weights = training_weights,
          parallel = FALSE,
          type.measure = "deviance",
          nfolds = 5,
          trace.it = 0
        )
        
        # Add model to full list
        model_per_group[[each_sa]] <- trained_model
        
      }
      
      # What is the predictions of the model?
      pred_train <- predict(model_per_group[["1"]],
                            newx = x_train,
                            type = "response",
                            s = "lambda.min")[,1]
      pred_val <- predict(model_per_group[["1"]],
                          newx = x,
                          type = "response",
                          s = "lambda.min")[,1]
      for (each_sa in sa_levels) {
        pred_train[sa_train==each_sa] <- predict(model_per_group[[each_sa]],
                                                 newx = x_train[sa_train==each_sa,],
                                                 type = "response",
                                                 s = "lambda.min")[,1]
        pred_val[sa_val==each_sa] <- predict(model_per_group[[each_sa]],
                                             newx = x_val[sa_val==each_sa,],
                                             type = "response",
                                             s = "lambda.min")[,1]
      }
      
    }
    
    # For each constraint
    for (constraint in c("NoCon","3Perc","1Perc")) {
      
      # How many people can we at most screen?
      con <- case_when(
        constraint == "NoCon" ~ 1*nrow(x_train),
        constraint == "3Perc" ~ 0.03*nrow(x_train),
        constraint == "1Perc" ~ 0.01*nrow(x_train)
      )
      
      # We are done fitting models! #
      
      
      
      # From here, we use pred_train and pred_val to 
      # calculate val and train benefit
      
      # TopDeprivNB
      per_bootstrap_metrics[[model_name]][[constraint]][["Train"]][["TopDeprivNB"]] <- 1 - 
        mean(y_train[top_quintile[bootstrap_index]==1]) + 
        lung_cancer_effect * 
        calcNetBenefit(
          rep(0.015,nrow(pareto_dict[[model_name]][[constraint]])),
          pareto_dict[[model_name]][[constraint]]$TopDeprivThreshold,
          y_train[top_quintile[bootstrap_index]==1],
          pred_train[top_quintile[bootstrap_index]==1],
          new_version = FALSE
        )
      per_bootstrap_metrics[[model_name]][[constraint]][["Val"]][["TopDeprivNB"]] <- 1 - 
        mean(y_val[top_quintile==1]) + 
        lung_cancer_effect * 
        calcNetBenefit(
          rep(0.015,nrow(pareto_dict[[model_name]][[constraint]])),
          pareto_dict[[model_name]][[constraint]]$TopDeprivThreshold,
          y_val[top_quintile==1],
          pred_val[top_quintile==1],
          new_version = FALSE
        )
      
      # RestNB
      per_bootstrap_metrics[[model_name]][[constraint]][["Train"]][["RestNB"]] <- 1 - 
        mean(y_train[top_quintile[bootstrap_index]!=-7]) + 
        lung_cancer_effect * 
        calcNetBenefit(
          rep(0.015,nrow(pareto_dict[[model_name]][[constraint]])),
          pareto_dict[[model_name]][[constraint]]$RestThreshold,
          y_train[top_quintile[bootstrap_index]!=-7],
          pred_train[top_quintile[bootstrap_index]!=-7],
          new_version = FALSE
        )
      per_bootstrap_metrics[[model_name]][[constraint]][["Val"]][["RestNB"]] <- 1 - 
        mean(y_val[top_quintile!=-7]) + 
        lung_cancer_effect * 
        calcNetBenefit(
          rep(0.015,nrow(pareto_dict[[model_name]][[constraint]])),
          pareto_dict[[model_name]][[constraint]]$RestThreshold,
          y_val[top_quintile!=-7],
          pred_val[top_quintile!=-7],
          new_version = FALSE
        )
      
      
    }
  }
  
  print(glue("Completed {bb}"))
  
  # At the end, return bootstrap metrics
  return(per_bootstrap_metrics)
}

# Now that we have the optimism we can correct
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  print(model_name)
  
  # For each group, including overall
  for (constraint in c("NoCon","3Perc","1Perc")) {
    
    # Get the optimism mean
    if (constraint=="NoCon") {
      optimism_topdepriv <- mean(sapply(
        pareto_optimism,
        function(x)
          x[[model_name]][[constraint]][["Train"]][["TopDeprivNB"]] - 
          x[[model_name]][[constraint]][["Val"]][["TopDeprivNB"]]
      ))
      optimism_rest <- mean(sapply(
        pareto_optimism,
        function(x)
          x[[model_name]][[constraint]][["Train"]][["RestNB"]] - 
          x[[model_name]][[constraint]][["Val"]][["RestNB"]]
      ))
    } else {
      optimism_topdepriv <- rowMeans(sapply(
        pareto_optimism,
        function(x)
          x[[model_name]][[constraint]][["Train"]][["TopDeprivNB"]] - 
          x[[model_name]][[constraint]][["Val"]][["TopDeprivNB"]]
      ))
      optimism_rest <- rowMeans(sapply(
        pareto_optimism,
        function(x)
          x[[model_name]][[constraint]][["Train"]][["RestNB"]] - 
          x[[model_name]][[constraint]][["Val"]][["RestNB"]]
      ))
    }
    
    
    # Save the corrected
    pareto_model_constraint_df <- pareto_dict[[model_name]][[constraint]]
    pareto_model_constraint_df$TopDeprivOCNB <- pareto_model_constraint_df$TopDeprivNB -
      optimism_topdepriv
    pareto_model_constraint_df$RestOCNB <- pareto_model_constraint_df$RestNB -
      optimism_rest
    pareto_dict[[model_name]][[constraint]] <- pareto_model_constraint_df 
    
  }
}

# Now we save pareto_dict for later plotting
saveRDS(pareto_dict,"Output/pareto_dict.rds")











####      Confidence intervals for Pareto line ####

# Add to the pareto the confidence intervals
pareto_dict_with_ci <- readRDS("Output/pareto_dict.rds")

# We do boostraps to calculate NB at each threshold
pareto_bootstraps <- list()

# Here, our validation data is x and y
x_val <- x
y_val <- y
sa_val <- df$sa 

set.seed(RANDOM_SEED)
pareto_bootstraps <- foreach(bb=1:NUM_BOOT) %dopar% {
  
  # Here we save everything
  per_bootstrap_metrics <- list()
  
  print(glue("Bootstrap Number {bb}"))
  
  # Select the new data
  bootstrap_index <- sample(nrow(df), replace = TRUE)
  
  # Training data
  x_train <- x_val[bootstrap_index,]
  y_train <- y_val[bootstrap_index]
  sa_train <- sa_val[bootstrap_index]
  
  for (model_name in c("LogNoSA", "LogSingleSA", "LogMultiSA")) {
    
    # Get prediction
    pred_train <- preds_train[[model_name]][bootstrap_index]
    
    
    # From here, we use pred_train and pred_val to calculate
    # The Overall NB and the NB in the Top Quintile of Deprivation.
    # For each constraint
    for (constraint in c("NoCon","3Perc","1Perc")) {
      
      # How many people can we at most screen?
      con <- case_when(
        constraint == "NoCon" ~ 1*nrow(x_train),
        constraint == "3Perc" ~ 0.03*nrow(x_train),
        constraint == "1Perc" ~ 0.01*nrow(x_train)
      )
      
      # Calculate NB
      per_bootstrap_metrics[[model_name]][[constraint]][["TopDeprivNB"]] <- 1 - mean(y_train[top_quintile[bootstrap_index]==1]) + 
        lung_cancer_effect * calcNetBenefit(
          rep(0.015,length(pareto_dict_with_ci[[model_name]][[constraint]]$TopDeprivThreshold)),
          pareto_dict_with_ci[[model_name]][[constraint]]$TopDeprivThreshold,
          y_train[top_quintile[bootstrap_index]==1],
          preds_train[[model_name]][top_quintile[bootstrap_index]==1],
          new_version = FALSE
        )
      per_bootstrap_metrics[[model_name]][[constraint]][["RestNB"]] <- 1 - mean(y_train[top_quintile[bootstrap_index]!=7]) + 
        lung_cancer_effect * calcNetBenefit(
          rep(0.015,length(pareto_dict_with_ci[[model_name]][[constraint]]$RestThreshold)),
          pareto_dict_with_ci[[model_name]][[constraint]]$RestThreshold,
          y_train[top_quintile[bootstrap_index]!=7],
          preds_train[[model_name]][top_quintile[bootstrap_index]!=7],
          new_version = FALSE
        )
      
    }
      
    }
  
  print(glue("Completed {bb}"))
  
  # At the end, return bootstrap metrics
  return(per_bootstrap_metrics)
  
}

# Now that we have the optimism we can correct
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  print(model_name)
  
  # For each group, including overall
  for (constraint in c("NoCon","3Perc","1Perc")) {
    
    # Get the optimism mean
    if (constraint=="NoCon") {
      
      
      # No Constraints
      # Low CI
      pareto_topdepriv_low <- quantile(sapply(
        pareto_bootstraps,
        function(x)
          x[[model_name]][[constraint]][["TopDeprivNB"]]
      ), 0.05)[[1]]
      pareto_rest_low <- quantile(sapply(
        pareto_bootstraps,
        function(x)
          x[[model_name]][[constraint]][["RestNB"]]
      ), 0.05)[[1]]
      # High CI
      pareto_topdepriv_high <- quantile(sapply(
        pareto_bootstraps,
        function(x)
          x[[model_name]][[constraint]][["TopDeprivNB"]]
      ), 0.95)[[1]]
      pareto_rest_high <- quantile(sapply(
        pareto_bootstraps,
        function(x)
          x[[model_name]][[constraint]][["RestNB"]]
      ), 0.95)[[1]]
      
      # What is the optimism?
      topdepriv_opt <- pareto_dict_with_ci[[model_name]][[constraint]][["TopDeprivNB"]] -
        pareto_dict_with_ci[[model_name]][[constraint]][["TopDeprivOCNB"]]
      pareto_topdepriv_low <- pareto_topdepriv_low - topdepriv_opt
      pareto_topdepriv_high <- pareto_topdepriv_high - topdepriv_opt
      rest_opt <- pareto_dict_with_ci[[model_name]][[constraint]][["RestNB"]] -
        pareto_dict_with_ci[[model_name]][[constraint]][["RestOCNB"]]
      pareto_rest_low <- pareto_rest_low - rest_opt
      pareto_rest_high <- pareto_rest_high - rest_opt
      
    } else {
      
      # Some Constraints
      # Low CI
      pareto_topdepriv_low <- apply(sapply(
        pareto_bootstraps,
        function(x) x[[model_name]][[constraint]][["TopDeprivNB"]]
      ), 1, quantile, probs = 0.05)
      pareto_rest_low <- apply(sapply(
        pareto_bootstraps,
        function(x) x[[model_name]][[constraint]][["RestNB"]]
      ), 1, quantile, probs = 0.05)
      # High CI
      pareto_topdepriv_high <- apply(sapply(
        pareto_bootstraps,
        function(x) x[[model_name]][[constraint]][["TopDeprivNB"]]
      ), 1, quantile, probs = 0.95)
      pareto_rest_high <- apply(sapply(
        pareto_bootstraps,
        function(x) x[[model_name]][[constraint]][["RestNB"]]
      ), 1, quantile, probs = 0.95)
      
      # What is the optimism?
      topdepriv_opt <- pareto_dict_with_ci[[model_name]][[constraint]][["TopDeprivNB"]] -
        pareto_dict_with_ci[[model_name]][[constraint]][["TopDeprivOCNB"]]
      pareto_topdepriv_low <- pareto_topdepriv_low - topdepriv_opt
      pareto_topdepriv_high <- pareto_topdepriv_high - topdepriv_opt
      rest_opt <- pareto_dict_with_ci[[model_name]][[constraint]][["RestNB"]] -
        pareto_dict_with_ci[[model_name]][[constraint]][["RestOCNB"]]
      pareto_rest_low <- pareto_rest_low - rest_opt
      pareto_rest_high <- pareto_rest_high - rest_opt
    }
    
    
    # Save the corrected
    pareto_model_constraint_df <- pareto_dict_with_ci[[model_name]][[constraint]]
    pareto_model_constraint_df$TopDeprivOCNB.Low <- pareto_topdepriv_low
    pareto_model_constraint_df$RestOCNB.Low <- pareto_rest_low
    pareto_model_constraint_df$TopDeprivOCNB.High <- pareto_topdepriv_high
    pareto_model_constraint_df$RestOCNB.High <- pareto_rest_high
    pareto_dict_with_ci[[model_name]][[constraint]] <- pareto_model_constraint_df 
    
  }
}

# Save new dataframe
saveRDS(pareto_dict_with_ci,"Output/pareto_with_ci_dict.rds")