# Define Source Directory
setwd("/mnt/bmh01-rds/Sperrin_UKBB_Fairness/Code/study_1_benefit/Manuscript_Ready_Code/FairSubgroupBenefit/")

# Set random seed
RANDOM_SEED <- 42
set.seed(RANDOM_SEED)

#------------------------#
#------------------------#
#### Import libraries ####
#------------------------#
#------------------------#

library(tidyverse)  # For data managing
library(pROC)       # For calculating performance
library(glmnet)     # First model: Elastic-Net
library(xgboost)    # Second model: XGBoost
library(rBayesianOptimization)  # For HP tuning
library(caret)      # For confusion matrix
library(pracma)     # For trapezoidal integration
library(probably)  # For calibration plot
library(doParallel) # For parallel computing
library(pmsampsize) # For sample size calculations
# For MICE, you need to import nlopt and cmake
library(mice)       # Data imputation

#------------------#
#------------------#
#### Parameters ####
#------------------#
#------------------#

# How many cores do we have?
registerDoParallel(cores=5)

# Number of optimisation runs per parameter
#NRUNS_PER_PARAM <- 10 # original number for manuscript is 10
NRUNS_PER_PARAM <- 1 # We choose 1 here instead to make code snippet run faster


# Multiple imputation parameters
#MICE_MAXIT <- 30 # original iterations were sent to 30
MICE_MAXIT <- 5 # We choose 5 here instead to make the code snippet run quicker
MICE_M <- 1     # We keep multiple imputations to 1, as missing data is not the focus of this work

# Forgetting factors explored
#forget_factor_range <- seq(0,1,0.1)  # original number
forget_factor_range <- seq(0,1,0.5)  # faster to run

### To load data later
which_output_path <- "_randomsplit"

#-----------------------#
#-----------------------#
#### Misc. Functions ####
#-----------------------#
#-----------------------#

# Define diabetes Sensitive attribute from ethnicity
define_diabetes_sa <- function(ethnicity) {
  
  # Define 3 groups: White, Asian, Black and and Other
  sensitive_attribute <- as.character(ethnicity)
  
  # For fairness evaluations
  sensitive_attribute <- factor(sensitive_attribute,levels=c("White","Asian","Black","Other"))
  
  return(sensitive_attribute)
}

# Define lung cancer Sensitive attribute from Townsend
define_lungcancer_sa <- function(twnsd,quintiles_twnsd) {
  
  # get the intervals of each individual according to their Townsend
  twnsd_q <- as.factor(findInterval(twnsd,quintiles_twnsd))
  
  return(twnsd_q)
}

# Cross-Entropy Loss
logLoss <- function(y_true,y_prob) {
  -1*sum((1-y_true)*log(1-y_prob) + y_true*log(y_prob))
}

# Turns dataframe into x and y variables ready for glmnet or xgboost
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
    trace.it = FALSE
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
      trace.it = FALSE
    )
    
    # Report AUC out of curiosity
    pred <- predict(propensity_model, newx = predictor_matrix, type = "response",s = "lambda.min")[,1]
    auroc <- (roc(c(sa_to_predict),c(pred),ci=TRUE))$auc
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

#--------------------------------------#
#--------------------------------------#
#### Open the dataset and impute it ####
#--------------------------------------#
#--------------------------------------#

# We want to do this analysis for both Diabetes and Lung cancer
for (WHICH_DISEASE in c("diabetes","lungcancer")) {
  
print(paste("Working on disease",WHICH_DISEASE))


# Open dataset
df <- read.csv(file = paste0("Data/",WHICH_DISEASE,"_df",which_output_path,".csv"))


# X is rownames
rownames(df) <- df$X
df <- df %>% select(-X)

# Turn factor values into factors
df$Ethnicity <- as.factor(df$Ethnicity)
df$Education <- as.factor(df$Education)
df$Smoker <- as.factor(df$Smoker)

# Make backup just in case
backup_df <- df
df <- backup_df

# Remove values not used in prediction or y or sa
df[c("Censored","Started","Death","Centre")] <- NULL

# Change to factors strings and binary variables
df[, sapply(df, is.character)] <- lapply(df[, sapply(df, is.character)], as.factor)

# Get quintiles of Townsend
quintiles_townsend <- quantile(df$Townsend,probs = c(0.2,0.4,0.6,0.8),na.rm=TRUE)

# Choose which ones are train and which ones are val
df_val <- df[df$is_val==1,]
df <- df[df$is_val!=1,]
df_val$is_val <- NULL
df$is_val <- NULL

# Impute non_val using mice
df_mice <- mice(df,maxit = MICE_MAXIT, m = MICE_M)
df <- complete(df_mice, action = "stacked")
df <- na.omit(df)

# For the validation data, we assume that we will have
# all present data at the time of use, so we do independent
# multiple imputation
# We want only true SA so we drop all for which df_val$Ethnicity or others is missing
df_val <- df_val[!is.na(df_val$Age),]
df_val <- df_val[!is.na(df_val$Female),]
df_val <- df_val[!is.na(df_val$Townsend),]
df_val <- df_val[!is.na(df_val$Ethnicity),]
df_val_mice <- mice(df_val,maxit = MICE_MAXIT, m = MICE_M)
df_val <- complete(df_val_mice, action = "stacked")
df_val <- na.omit(df_val)

# Define SA differently depending on disease
if (WHICH_DISEASE=="diabetes") {
  
  # Define sa according to ethnicity
  df$sa <- define_diabetes_sa(df$Ethnicity)
  df_val$sa <- define_diabetes_sa(df_val$Ethnicity)
  
} else if (WHICH_DISEASE=="lungcancer") {
  
  # Define sa according to ethnicity
  df$sa <- define_lungcancer_sa(df$Townsend,quintiles_townsend)
  df_val$sa <- define_lungcancer_sa(df_val$Townsend,quintiles_townsend)
  
}

# Save data for further work later
saveRDS(df, file = paste0(
  "Data/mice_imputed_",WHICH_DISEASE,"_training",which_output_path,".rds"
))
saveRDS(df_val, file = paste0(
  "Data/mice_imputed_",WHICH_DISEASE,"_validation",which_output_path,".rds"
))

# Load data for further work
df <- readRDS(paste0(
  "Data/mice_imputed_",WHICH_DISEASE,"_training",which_output_path,".rds"
))
df_val <- readRDS(paste0(
  "Data/mice_imputed_",WHICH_DISEASE,"_validation",which_output_path,".rds"
))

# Check for the entirety of the dataset
print(paste("Sample size calculations for the disease",WHICH_DISEASE))
expected_cstatistic <- if (WHICH_DISEASE=="diabetes") {0.75} else if (WHICH_DISEASE=="lungcancer") {0.8}
sample_size_needed <- pmsampsize(type = "b",
                                 parameters = ncol(df)-2,
                                 prevalence = mean(c(df$y,df_val$y)),
                                 cstatistic = expected_cstatistic,
                                 seed = 42)$sample_size
print(paste0(
  "For the overall population, the sample size needed is ",
  sample_size_needed,
  " while our population is ",
  nrow(df),
  " (",
  round((nrow(df)/sample_size_needed)*100,0),
  "% of Needed)"
))
# Check for each SA: Is the Sample Size enough
for (each_sa in levels(df$sa)) {
  sample_size_needed <- pmsampsize(type = "b",
                                   parameters = ncol(df)-2,
                                   prevalence = mean(c(df[df$sa==each_sa,"y"],df_val[df_val$sa==each_sa,"y"])),
                                   cstatistic = expected_cstatistic,
                                   seed = RANDOM_SEED)$sample_size
  print(paste0(
    "For the ",
    each_sa,
    " population, the sample size needed is ",
    sample_size_needed,
    " while our population is ",
    nrow(df[df$sa==each_sa,]),
    " (",
    round((nrow(df[df$sa==each_sa,])/sample_size_needed)*100,0),
    "% of Needed)"
  ))
}

#--------------------#
#--------------------#
#### Prepare data ####
#--------------------#
#--------------------#

# Turn validation and training data into right format
x_and_y <- prepareForGLMNET(df)
x <- x_and_y[["x"]]
y <- x_and_y[["y"]]
x_and_y_val <- prepareForGLMNET(df_val)
x_val <- x_and_y_val[["x"]]
y_val <- x_and_y_val[["y"]]

# What if we don't want SA information in our model?
if (WHICH_DISEASE=="diabetes") {
  x_and_y_nosa <- df %>% select(-Ethnicity) %>% prepareForGLMNET
  x_and_y_val_nosa <- df_val %>% select(-Ethnicity) %>% prepareForGLMNET
} else if (WHICH_DISEASE=="lungcancer") {
  x_and_y_nosa <- df %>% select(-Townsend) %>% prepareForGLMNET
  x_and_y_val_nosa <- df_val %>% select(-Townsend) %>% prepareForGLMNET
}
x_nosa <- x_and_y_nosa[["x"]]
y_nosa <- x_and_y_nosa[["y"]]
x_val_nosa <- x_and_y_val_nosa[["x"]]
y_val_nosa <- x_and_y_val_nosa[["y"]]

#--------------------------------------#
#--------------------------------------#
#### Calculate the propensity score ####
#--------------------------------------#
#--------------------------------------#

# Assign folds to train data
foldid <- sample(1:5,size=length(y),replace=TRUE)

# For propensity scores, clean df so that it doesn't have sas.
if (WHICH_DISEASE=="diabetes") {
  df_for_propensity <- df %>% select(-sa,-Ethnicity)
} else if (WHICH_DISEASE=="lungcancer") {
  df_for_propensity <- df %>% select(-sa,-Townsend)
}

# What are the sa levels?
sa_levels <- levels(df$sa)

# We calculate two sets of scores: the five-fold ones, and the full ones
# First, five folds
# (The five-folds ones are used for Cross-validation during hyperparameter
# selection...)
# Get propensity scores for weighting each group, in each inner fold
propensity_list_per_folds <- foreach(ii=1:5) %dopar% {
  propensity_score_result <- calculate_sensitive_propensity_score(
    df_for_propensity[foldid != ii,],df$sa[foldid != ii],sa_levels, s = "lambda.min"
  )
  return(propensity_score_result)
}
propensity_score_list_per_fold <- list()
propensity_models_list_per_fold <- list()
for (ii in 1:5) {
  propensity_score_list_per_fold[[ii]] <- propensity_list_per_folds[[ii]]$propensity
  propensity_models_list_per_fold[[ii]] <- propensity_list_per_folds[[ii]]$model
}
# Combine both items for saving
propensity_list_per_fold <- list()
propensity_list_per_fold$propensity <- propensity_score_list_per_fold
propensity_list_per_fold$model <- propensity_models_list_per_fold
propensity_list_per_fold$foldid <- foldid
propensity_list_per_fold$sa_levels <- sa_levels
# Save
saveRDS(propensity_list_per_fold,
        file=paste0("Output/propensity_",WHICH_DISEASE,"_folds",which_output_path,".rds")
)

# Calculate same scores and model for full data (not CV folds)
propensity_list_total <- calculate_sensitive_propensity_score(
  df_for_propensity,df$sa,sa_levels, s = "lambda.min"
)
training_weights <- propensity_list_total$propensity
propensity_training_model <-propensity_list_total$model
# Save
saveRDS(propensity_list_total,
        file=paste0("Output/propensity_",WHICH_DISEASE,"_total",which_output_path,".rds")
)

# Load scores
propensity_list_per_fold <- readRDS(paste0("Output/propensity_",WHICH_DISEASE,"_folds",which_output_path,".rds"))
propensity_list_total <- readRDS(paste0("Output/propensity_",WHICH_DISEASE,"_total",which_output_path,".rds"))

# Get propensity scores used for models
propensity_score_list_per_fold <- propensity_list_per_fold$propensity
propensity_score_list_per_total <- propensity_list_total$propensity
foldid <- propensity_list_per_fold$foldid
sa_levels <- propensity_list_per_fold$sa_levels

#--------------------------#
#--------------------------#
#### Running the models ####
#--------------------------#
#--------------------------#

# For each type of model, we will run the analysis
# And save the model into a list that we will keep loading
if (!file.exists(
  paste0("Output/Models/",WHICH_DISEASE,"_models",which_output_path,".rds")
)) {
  # If not, create the file with an empty R list
  saveRDS(list(), file = paste0("Output/Models/",WHICH_DISEASE,"_models",which_output_path,".rds"))
}

# Which file are we saving to?
file_to_save_models <- paste0("Output/Models/",WHICH_DISEASE,"_models",which_output_path,".rds")

print(paste(
  "We are saving the models to",file_to_save_models
))

#==================================#
#==================================#
#### XGBSingleConv ####
#==================================#
#==================================#

# Measure time to complete
start_time <- Sys.time()

# Scoring function calculates mean benefit for a particular lambda across all folds
scoring_function_xgbsingleconv <- function(eta,max_depth,subsample,colsample_bytree,nrounds) {
  
  all_logloss <- c()
  for (fold in 1:5) {
    
    # Train model
    xgb_model <- xgboost(data = x[foldid!=fold,], 
                         label = y[foldid!=fold], 
                         objective = "binary:logistic",
                         eta = eta,
                         max_depth = max_depth,
                         subsample = subsample,
                         colsample_bytree = colsample_bytree,
                         nrounds = nrounds,
                         verbose = 0)
    
    # Predictions of the fold
    outer_preds <- predict(xgb_model, x[foldid==fold,])
    
    # Calculate logloss
    negll <- -1*logLoss(y_true = y[foldid==fold], y_prob = outer_preds)
    
    # Append to folds
    all_logloss <- c(all_logloss,negll)
  }
  all_logloss <- as.numeric(all_logloss)
  
  # Calculate mean benefit and return as objective
  mean_logloss <- mean(all_logloss)
  return(list(Score = mean_logloss))
}

# What are the bounds of the parameters
random_hpam_set <- data.frame(
  eta = runif(5*NRUNS_PER_PARAM,min=0,max=1),
  max_depth = sample(1:6,5*NRUNS_PER_PARAM, replace=T),
  subsample = runif(5*NRUNS_PER_PARAM,min=0.1,max=1),
  colsample_bytree = runif(5*NRUNS_PER_PARAM,min=0.1,max=1),
  nrounds = sample(5:80,5*NRUNS_PER_PARAM, replace=T)
)

# For each row of the hparameter set, calculate metric
metric_per_hpam <- c()
for (ii in 1:nrow(random_hpam_set)) {
  print(paste("Exploring HPAM Set",ii))
  metric <- scoring_function_xgbsingleconv(random_hpam_set[ii,"eta"],
                                           random_hpam_set[ii,"max_depth"],
                                           random_hpam_set[ii,"subsample"],
                                           random_hpam_set[ii,"colsample_bytree"],
                                           random_hpam_set[ii,"nrounds"])
  metric_per_hpam <- c(metric_per_hpam,metric)
}

# What ii is best?
best_ii <- which.max(metric_per_hpam)

# Fetch best result
best_eta <- random_hpam_set[best_ii,"eta"]
best_max_depth <- random_hpam_set[best_ii,"max_depth"]
best_subsample <- random_hpam_set[best_ii,"subsample"]
best_colsample_bytree <- random_hpam_set[best_ii,"colsample_bytree"]
best_nrounds <- random_hpam_set[best_ii,"nrounds"]

# Train best model
xgb_model <- xgboost(data = x, 
                     label = y, 
                     objective = "binary:logistic",
                     eta = best_eta,
                     max_depth = best_max_depth,
                     subsample = best_subsample,
                     colsample_bytree = best_colsample_bytree,
                     nrounds = best_nrounds,
                     verbose = 0)

# Print AUROC just to check that things are good
pred <- predict(xgb_model, x)
pred_val <- predict(xgb_model, x_val)
auroc <- (roc(c(y_val),c(pred_val),ci=TRUE))$auc
print(
  paste(
    "The AUC of XGBSingleConv is",
    auroc
  )
)

# Measure time to complete
total_time <- Sys.time()-start_time

# Save model
model <- list()
model$Model <- xgb_model
model$Eta <- best_eta
model$Best_max_depth <- best_max_depth
model$Subsample <- best_subsample
model$Best_colsample_bytree <- best_colsample_bytree
model$Best_nrounds <- best_nrounds
model$Search$Hpams <- random_hpam_set
model$Search$Metric <- metric_per_hpam
model$Time <- total_time

# Save to path calculated before
models <- readRDS(file_to_save_models)
models$XGBSingleConv <- model
saveRDS(models, file = file_to_save_models)
rm(models)

#============================#
#============================#
#### XGBNoSA ####
#============================#
#============================#

# Measure time to complete
start_time <- Sys.time()

# Scoring function calculates mean benefit for a particular lambda across all folds
scoring_function_xgbsingleconv <- function(eta,max_depth,subsample,colsample_bytree,nrounds) {
  
  all_logloss <- c()
  for (fold in 1:5) {
    
    # Train model
    xgb_model <- xgboost(data = x_nosa[foldid!=fold,], 
                         label = y_nosa[foldid!=fold], 
                         objective = "binary:logistic",
                         eta = eta,
                         max_depth = max_depth,
                         subsample = subsample,
                         colsample_bytree = colsample_bytree,
                         nrounds = nrounds,
                         verbose = 0)
    
    # Predictions of the fold
    outer_preds <- predict(xgb_model, x_nosa[foldid==fold,])
    
    
    # Calculate logloss
    negll <- -1*logLoss(y_true = y_nosa[foldid==fold], y_prob = outer_preds)
    
    # Append to folds
    all_logloss <- c(all_logloss,negll)
  }
  all_logloss <- as.numeric(all_logloss)
  
  # Calculate mean benefit and return as objective
  mean_logloss <- mean(all_logloss)
  return(list(Score = mean_logloss))
}

# We get the random hpam_set from XGBSingleConv
models <- readRDS(file_to_save_models)
random_hpam_set <- models$XGBSingleConv$Search$Hpams
rm(models)

# For each row of the hparameter set, calculate metric
metric_per_hpam <- c()
for (ii in 1:nrow(random_hpam_set)) {
  print(paste("Exploring HPAM Set",ii))
  metric <- scoring_function_xgbsingleconv(random_hpam_set[ii,"eta"],
                                           random_hpam_set[ii,"max_depth"],
                                           random_hpam_set[ii,"subsample"],
                                           random_hpam_set[ii,"colsample_bytree"],
                                           random_hpam_set[ii,"nrounds"])
  metric_per_hpam <- c(metric_per_hpam,metric)
}

# What ii is best?
best_ii <- which.max(metric_per_hpam)

# Fetch best result
best_eta <- random_hpam_set[best_ii,"eta"]
best_max_depth <- random_hpam_set[best_ii,"max_depth"]
best_subsample <- random_hpam_set[best_ii,"subsample"]
best_colsample_bytree <- random_hpam_set[best_ii,"colsample_bytree"]
best_nrounds <- random_hpam_set[best_ii,"nrounds"]

# Train best model
xgb_model <- xgboost(data = x_nosa, 
                     label = y_nosa, 
                     objective = "binary:logistic",
                     eta = best_eta,
                     max_depth = best_max_depth,
                     subsample = best_subsample,
                     colsample_bytree = best_colsample_bytree,
                     nrounds = best_nrounds,
                     verbose = 0)

# Print AUROC just to check that things are good
pred <- predict(xgb_model, x_nosa)
pred_val <- predict(xgb_model, x_val_nosa)
auroc <- (roc(c(y_val_nosa),c(pred_val),ci=TRUE))$auc
print(
  paste(
    "The AUC of XGBSingleConv is",
    auroc
  )
)

# Measure time to complete
total_time <- Sys.time()-start_time

# Save model
model <- list()
model$Model <- xgb_model
model$Eta <- best_eta
model$Best_max_depth <- best_max_depth
model$Subsample <- best_subsample
model$Best_colsample_bytree <- best_colsample_bytree
model$Best_nrounds <- best_nrounds
model$Search$Hpams <- random_hpam_set
model$Search$Metric <- metric_per_hpam
model$Time <- total_time

# Save to path calculated before
models <- readRDS(file_to_save_models)
models$XGBNoSA <- model
saveRDS(models, file = file_to_save_models)
rm(models)

#===================================#
#===================================#
#### XGBMultiConv ####
#===================================#
#===================================#

# Measure time to complete
start_time <- Sys.time()

# for each group, we create a different model, propensity-weighted
for (each_sa in sa_levels) {
  
  
  # Scoring function calculates log loss for a particular lambda,forgetting factor across all folds
  scoring_function_xgbmulticonv <- function(eta,max_depth,subsample,colsample_bytree,nrounds,forget) {
    
    all_neglogloss <- c()
    for (fold in 1:5) {
      
      # To train model, create weights with prop score and forgetting factor
      training_weights <- propensity_score_list_per_fold[[fold]][[each_sa]]
      # Add forgetting to weights
      training_weights[df[foldid!=fold,"sa"] != each_sa] <- pmin(1,forget*training_weights[df[foldid!=fold,"sa"] != each_sa])
      
      # With weights and data, create matrix
      xgbMatrix <- xgb.DMatrix(x[foldid!=fold,], 
                               label = y[foldid!=fold], 
                               weight = training_weights)
      
      # Train model
      xgb_model <- xgboost(data = xgbMatrix, 
                           objective = "binary:logistic",
                           eta = eta,
                           max_depth = max_depth,
                           subsample = subsample,
                           colsample_bytree = colsample_bytree,
                           nrounds = nrounds,
                           verbose = 0)
      
      # Prediction of the fold
      inner_preds <- predict(xgb_model, x[(foldid!=fold)&(df$sa==each_sa),])
      outer_preds <- predict(xgb_model, x[(foldid==fold)&(df$sa==each_sa),])
      
      # No recalibration for log loss optimisation
      # outer_preds <- recalibratePredPerSA(y[foldid!=fold,],inner_preds, df[foldid!=fold,"sa"], outer_preds, df[foldid==fold,"sa"])
      
      # Calculate negative log loss
      negll <- -1*logLoss(y_true = y[(foldid==fold)&(df$sa==each_sa),],y_prob = outer_preds)
      
      # Append to folds
      all_neglogloss <- c(all_neglogloss,negll)
    }
    all_neglogloss <- as.numeric(all_neglogloss)
    
    # Calculate mean benefit and return as objective
    mean_negll <- mean(all_neglogloss)
    return(list(Score = mean_negll))
  }
  
  # We get the random hpam_set from XGBSingleConv
  models <- readRDS(file_to_save_models)
  random_hpam_set <- models$XGBSingleConv$Search$Hpams
  rm(models)
  # We add forgetting factor: all combinations from 0 to 1
  grid_forget <- expand.grid(1:nrow(random_hpam_set), forget_factor_range)
  random_hpam_set <- data.frame(random_hpam_set[grid_forget$Var1, ], forget = grid_forget$Var2)
  
  # For each row of the hparameter set, calculate metric
  metric_per_hpam <- c()
  for (ii in 1:nrow(random_hpam_set)) {
    print(paste("Exploring HPAM Set",ii,"with forget",random_hpam_set[ii,"forget"]))
    metric <- scoring_function_xgbmulticonv(random_hpam_set[ii,"eta"],
                                            random_hpam_set[ii,"max_depth"],
                                            random_hpam_set[ii,"subsample"],
                                            random_hpam_set[ii,"colsample_bytree"],
                                            random_hpam_set[ii,"nrounds"],
                                            random_hpam_set[ii,"forget"]
    )
    metric_per_hpam <- c(metric_per_hpam,metric)
  }
  
  # What ii is best?
  best_ii <- which.max(metric_per_hpam)
  
  # Fetch best result
  best_eta <- random_hpam_set[best_ii,"eta"]
  best_max_depth <- random_hpam_set[best_ii,"max_depth"]
  best_subsample <- random_hpam_set[best_ii,"subsample"]
  best_colsample_bytree <- random_hpam_set[best_ii,"colsample_bytree"]
  best_nrounds <- random_hpam_set[best_ii,"nrounds"]
  best_forget <- random_hpam_set[best_ii,"forget"]
  
  # Get best weights for SA
  training_weights <- propensity_score_list_per_total[[each_sa]]
  # Add forgetting to weights
  training_weights[df$sa != each_sa] <- pmin(1,best_forget*training_weights[df$sa != each_sa])
  
  # Create training matrix and train model
  xgbMatrix <- xgb.DMatrix(x, 
                           label = y, 
                           weight = training_weights)
  xgb_model <- xgboost(data = xgbMatrix, 
                       objective = "binary:logistic",
                       eta = best_eta,
                       max_depth = best_max_depth,
                       subsample = best_subsample,
                       colsample_bytree = best_colsample_bytree,
                       nrounds = best_nrounds,
                       verbose = 0)
  
  # Print AUROC just to check that things are good
  pred <- predict(xgb_model, x[df$sa==each_sa,])
  pred_val <- predict(xgb_model, x_val[df_val$sa==each_sa,])
  auroc <- (roc(c(y_val[df_val$sa==each_sa,]),c(pred_val),ci=TRUE))$auc
  print(
    paste(
      "The AUC of XGBMultiConv is",
      auroc,
      "for",
      each_sa
    )
  )
  
  # Save model
  model <- list()
  model <- list()
  model$Model <- xgb_model
  model$Eta <- best_eta
  model$Best_max_depth <- best_max_depth
  model$Subsample <- best_subsample
  model$Best_colsample_bytree <- best_colsample_bytree
  model$Best_nrounds <- best_nrounds
  model$Forget <- best_forget
  model$Search$Hpams <- random_hpam_set
  model$Search$Metric <- metric_per_hpam
  
  # Save to path calculated before
  models <- readRDS(file_to_save_models)
  models$XGBMultiConv[[each_sa]] <- model
  saveRDS(models, file = file_to_save_models)
  
}

# Calculate final time
total_time <- Sys.time()-start_time
models <- readRDS(file_to_save_models)
models$XGBMultiConv$Time <- total_time
saveRDS(models, file = file_to_save_models)


# Finish loop that goes over WHICH_DISEASE in c("diabetes","lungcancer")
}