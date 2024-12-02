# Define Source Directory
setwd("/mnt/bmh01-rds/Sperrin_UKBB_Fairness/Code/study_1_benefit/Manuscript_Ready_Code/FairSubgroupBenefit/")

# Set random seed
RANDOM_SEED <- 42
set.seed(RANDOM_SEED)

# Plotting colors
clrs <- c(
  "#E02C83","#1660F5","#04831B","#E02A00","#07F784"
)



#### Import libraries ####
library(tidyverse)  # For data managing
library(pROC)       # For calculating performance
library(glmnet)     # First model: Elastic-Net
library(xgboost)    # Second model: XGBoost
library(rBayesianOptimization)  # For HP tuning
library(caret)      # For confusion matrix
library(pracma)     # For trapezoidal integration
library(probably)  # For calibration plot
library(ggplot2) # For sample size calculations
library(KraljicMatrix) # For identifying Pareto
library(gtsummary) # For nice tables

# Get string of file to open
which_output_path <- "_randomsplit"

#### Misc. Functions ####



# Define lung cancer Sensitive attribute from Townsend
define_lungcancer_sa <- function(twnsd,quintiles_twnsd) {
  
  # get the intervals of each individual according to their Townsend
  twnsd_q <- as.factor(findInterval(twnsd,quintiles_twnsd))
  
  return(twnsd_q)
}

# Calculates the net benefit curve (Corrected as per our work)
calcNetBenefit <- function(clinical_thresholds,model_thresholds,y_true,s_pred,new_version=TRUE) {
  
  # First we create a matrix so that each column of the matrix is the binary predictions at given threshold
  binary_predictions <- outer(s_pred,model_thresholds,FUN=">")
  true_positives <- colMeans(binary_predictions&y_true)
  false_positives <- colMeans(binary_predictions&!(y_true))
  if (new_version==TRUE) {
    nb <- true_positives/clinical_thresholds - 1/(1-clinical_thresholds)*false_positives
  } else {
    nb <- true_positives - clinical_thresholds/(1-clinical_thresholds)*false_positives
  }
  
  return(nb)
}

# Calculates net benefit under a resource constraint
NB_under_constraints <- function(clinical_threshold,y_true,s_pred,resource_constraint,new_version=TRUE) {
  
  # If all s_pred are the same, we take a random sample of the patients and screen them, according to capacity
  # Unless s_pred is all 0, in that case we screen no one
  if (length(unique(s_pred))==1) {
    if (unique(s_pred)==0) {
      # Everyone is predicted to be zero
      constrained_benefit <- calcNetBenefit(clinical_threshold,clinical_threshold,y_true,s_pred*0,new_version)
    } else {
      # Random sample
      sample_to_screen <- sample(c(0,1),size=length(y_true),prob=c(1-capacity,capacity),replace=TRUE)
      constrained_benefit <- calcNetBenefit(clinical_threshold,clinical_threshold,y_true,sample_to_screen,new_version)
    }
  } else {
    
    # We assume the model is calibrated so the best threshold should be clinical_threshold
    # We then determine whether the resource_constraint is more restrictive than that
    resource_threshold <- quantile(s_pred,prob=(1-resource_constraint),na.rm=TRUE)
    
    # We go for the threshold that is the highest amongst both
    chosen_threshold <- pmax(clinical_threshold,resource_threshold)
    
    # The net benefit has the og clinical threshold but the new chosen threshold
    constrained_benefit <- calcNetBenefit(clinical_threshold,chosen_threshold,y_true,s_pred,new_version)
    
  }
  
  return(constrained_benefit)
  
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





#### NET BENEFIT PARAMETERS ####





# For diabetes, we need thresholds over which to integrate, and their weight
diabetes_threshold <- 0.15
diabetes_effect <- 0.58

# For lung cancer, we need threshold used and capacity
lung_cancer_threshold <- 0.015
capacity <- 1
lung_cancer_effect <- 0.20

# Calculates final expression of subgroup NB for diabetes
# In pp true negatives per patients
sNB_diabetes <- function(yy_true,yy_pred, pp = 1) {
  
  # Get prevalence
  prevalence <- mean(yy_true)
  
  # Get benefit
  bnft <- calcNetBenefit(diabetes_threshold,diabetes_threshold,yy_true,yy_pred,new_version=FALSE)
  
  
  # Get vinal sNB
  sNB <- (1 - prevalence + diabetes_effect*bnft)*pp
  
  return(sNB)
}
# Do the same for lung cancer
balance_lungcancer <- 1/lung_cancer_threshold
sNB_lungcancer <- function(yy_true,yy_pred, pp = 1, rc = capacity) {
  
  # Get prevalence
  prevalence <- mean(yy_true)
  
  # Get benefit
  bnft <- NB_under_constraints(lung_cancer_threshold,yy_true,yy_pred,rc,new_version=TRUE)
  
  # Get vinal sNB
  sNB <- (1 - prevalence + lung_cancer_effect*bnft/balance_lungcancer)*pp
  
  return(sNB)
}





#### Data Preparation: Load files/Predictions ####




## Diabetes

# Open training/validation data
df_train_diabetes <- readRDS(paste0(
  "Data/mice_imputed_diabetes_training",which_output_path,".rds"
))
df_val_diabetes <- readRDS(paste0(
  "Data/mice_imputed_diabetes_validation",which_output_path,".rds"
))

# Turn validation and training data into right format
x_and_y_train_diabetes <- prepareForGLMNET(df_train_diabetes)
x_train_diabetes <- x_and_y_train_diabetes[["x"]]
x_and_y_val_diabetes <- prepareForGLMNET(df_val_diabetes)
x_val_diabetes <- x_and_y_val_diabetes[["x"]]
# Without SA info
x_and_y_train_nosa_diabetes <- df_train_diabetes %>% select(-Ethnicity) %>% prepareForGLMNET
x_train_nosa_diabetes <- x_and_y_train_nosa_diabetes[["x"]]
x_and_y_val_nosa_diabetes <- df_val_diabetes %>% select(-Ethnicity) %>% prepareForGLMNET
x_val_nosa_diabetes <- x_and_y_val_nosa_diabetes[["x"]]
# SA info
sa_train_diabetes <- as.character(df_train_diabetes$Ethnicity) 
sa_val_diabetes <- as.character(df_val_diabetes$Ethnicity)
# Outcome info
y_train_diabetes <- as.numeric(x_and_y_train_diabetes[["y"]])
y_val_diabetes <- as.numeric(x_and_y_val_diabetes[["y"]])
# Open models
models_diabetes <- readRDS(paste0("Output/Models/diabetes_models",which_output_path,".rds"))


# Create predictions
# This will contain predictions for each model
predictions_diabetes <- list()

# 1) XGBNoSA
pred_val <- predict(models_diabetes$XGBNoSA$Model,x_val_nosa_diabetes)
predictions_diabetes$XGBNoSA <- as.numeric(pred_val)

# 2) XGBSingleSA
pred_val <- predict(models_diabetes$XGBSingleConv$Model,x_val_diabetes)
predictions_diabetes$XGBSingleSA <- as.numeric(pred_val)

# 3) XGBMultiSA
pred_val <- 0*y_val_diabetes
for (each_sa in c("Asian","Other","White","Black")) {
  
  # For each group, calculate the predictions for that group in train and val
  pred_val_per_sa <- predict(models_diabetes$XGBMultiConv[[each_sa]]$Model,x_val_diabetes[sa_val_diabetes==each_sa,], type = "response")
  
  # Save onto pred_val
  pred_val[sa_val_diabetes==each_sa] <- as.numeric(pred_val_per_sa)
}
predictions_diabetes$XGBMultiSA <- as.numeric(pred_val)


## Lung Cancer

# Open training/validation data
df_train_lungcancer <- readRDS(paste0(
  "Data/mice_imputed_lungcancer_training",which_output_path,".rds"
))
df_val_lungcancer <- readRDS(paste0(
  "Data/mice_imputed_lungcancer_validation",which_output_path,".rds"
))

# Turn validation and training data into right format
x_and_y_train_lungcancer <- prepareForGLMNET(df_train_lungcancer)
x_train_lungcancer <- x_and_y_train_lungcancer[["x"]]
x_and_y_val_lungcancer <- prepareForGLMNET(df_val_lungcancer)
x_val_lungcancer <- x_and_y_val_lungcancer[["x"]]
# Without SA info
x_and_y_train_nosa_lungcancer <- df_train_lungcancer %>% select(-Townsend) %>% prepareForGLMNET
x_train_nosa_lungcancer <- x_and_y_train_nosa_lungcancer[["x"]]
x_and_y_val_nosa_lungcancer <- df_val_lungcancer %>% select(-Townsend) %>% prepareForGLMNET
x_val_nosa_lungcancer <- x_and_y_val_nosa_lungcancer[["x"]]
# SA info
sa_train_lungcancer <- as.character(df_train_lungcancer$sa) 
sa_val_lungcancer <- as.character(df_val_lungcancer$sa)
# Outcome info
y_train_lungcancer <- as.numeric(x_and_y_train_lungcancer[["y"]])
y_val_lungcancer <- as.numeric(x_and_y_val_lungcancer[["y"]])
# Open models
models_lungcancer <- readRDS(paste0("Output/Models/lungcancer_models",which_output_path,".rds"))


# Create predictions
# This will contain predictions for each model
predictions_lungcancer <- list()

# 1) XGBNoSA
pred_val <- predict(models_lungcancer$XGBNoSA$Model,x_val_nosa_lungcancer)
predictions_lungcancer$XGBNoSA <- as.numeric(pred_val)

# 2) XGBSingleSA
pred_val <- predict(models_lungcancer$XGBSingleConv$Model,x_val_lungcancer)
predictions_lungcancer$XGBSingleSA <- as.numeric(pred_val)

# 3) XGBMultiSA
pred_val <- 0*y_val_lungcancer
for (each_sa in unique(sa_train_lungcancer)) {
  
  # For each group, calculate the predictions for that group in train and val
  pred_val_per_sa <- predict(models_lungcancer$XGBMultiConv[[each_sa]]$Model,x_val_lungcancer[sa_val_lungcancer==each_sa,], type = "response")
  
  # Save onto pred_val
  pred_val[sa_val_lungcancer==each_sa] <- as.numeric(pred_val_per_sa)
}
predictions_lungcancer$XGBMultiSA <- as.numeric(pred_val)



#### First section of resutls: Diabetes ####

#### Sup. Tab. 1: Diabetes table ####

# Open dataset
diabetes_table_df <- read.csv(
  file = paste0("Data/diabetes_df",which_output_path,".csv")
  )
rownames(diabetes_table_df) <- diabetes_table_df$X
diabetes_table_df <- diabetes_table_df %>% select(-X)
diabetes_table_df$Ethnicity <- as.factor(diabetes_table_df$Ethnicity)
diabetes_table_df$Education <- as.factor(diabetes_table_df$Education)
diabetes_table_df$Smoker <- as.factor(diabetes_table_df$Smoker)

# Only choose columns we care about
diabetes_table_df <- diabetes_table_df %>%
  select(c(Age,Female,Townsend,Ethnicity,Education,FamilyHxHipFracture,FamilyHxDiabetes,FamilyHxCVD,BMI,
           MedicationCholesterol,MedicationBloodPressure,MedicationHRT,FirstIncidenceO24,
           FirstIncidenceF20to29,FirstIncidenceF30to31,FirstIncidenceI10to11,y))

# Create a table of data for supplementary
table_summary <- diabetes_table_df %>%
  tbl_summary(by = Ethnicity)
table_summary

#### Sup. Fig. 2a: AUROC XGB table ####

metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Metric.Avg=double(),
                       Metric.Low=double(),
                       Metric.High=double(),
                       stringsAsFactors=FALSE)

# For each model, get the AUC (with 95% CI) of each subgroup
for (model_name in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
  
  # What model are we using
  model_predictions <- predictions_diabetes[[model_name]]
  
  # For each subgroup, calculate subgroup AUC
  for (each_sa in c("Asian","Black","Other","White")) {
    
    # Calculate subgroup AUC
    auroc <- roc(y_val_diabetes[sa_val_diabetes==each_sa],
                 model_predictions[sa_val_diabetes==each_sa],
                 ci=TRUE, plot=FALSE)
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,Model=model_name,
                                          Metric.Avg = auroc$auc,
                                          Metric.Low = auroc$ci[1],
                                          Metric.High = auroc$ci[3])
  }
  
  # Calculate overall AUC
  auroc <- roc(y_val_diabetes,model_predictions, ci=TRUE, plot=FALSE)
  
  # Attach overall model to df
  metric_df[nrow(metric_df)+1,] <- list(Subgroup="Overall",Model=model_name,
                                        Metric.Avg = auroc$auc,
                                        Metric.Low = auroc$ci[1],
                                        Metric.High = auroc$ci[3])
  
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
colnames(metric_df) <- c("Ethnicity","Model", "AUROC", "AUROC.Low", "AUROC.High")

# Plot
pd <- position_dodge(0.6) # move them .05 to the left and right
p <- ggplot(metric_df, aes(x=Ethnicity, y=AUROC, colour=Model,shape = Model)) + 
  geom_errorbar(aes(ymin=AUROC.Low, ymax=AUROC.High), width=.1, position=pd) +
  geom_point(position=pd) +
  xlab("Ethnicity") +
  ylab("C-Statistic") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "XGBNoSA" = 0,
                                "XGBSingleSA" = 5,
                                "XGBMultiSA" = 7),
                     name = "Model") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  )  +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  )

# Save
saveRDS(metric_df,"Output/Plots/Data_SupFig2a.rds")
ggsave("Output/Plots/SupFig2a.png", plot = p, width = 8, height = 4, units = "in", dpi = 400)



#### Sup. Fig. 2b: Slope XGB plot ####



metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Slope.Avg=double(),
                       Slope.Low=double(),
                       Slope.High=double(),
                       Intercept.Avg=double(),
                       Intercept.Low=double(),
                       Intercept.High=double(),
                       stringsAsFactors=FALSE)

# For each model, get the AUC (with 95% CI) of each subgroup
for (model_name in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
  
  # What model are we using
  model_predictions <- predictions_diabetes[[model_name]]
  linear_pred <- -1*log(
    (1-model_predictions)/
      model_predictions
  )
  
  # For each subgroup, calculate subgroup AUC
  for (each_sa in c("Asian","Black","Other","White")) {
    
    # Which group are we interested in
    sel <- sa_val_diabetes==each_sa
    
    # Calculate subgroup slope + intercept
    slope_model <- glm(y_val_diabetes[sel]~linear_pred[sel],family="binomial")
    intercept_model <- glm(y_val_diabetes[sel]~offset(linear_pred[sel])+1,family="binomial")
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,Model=model_name,
                                          Slope.Avg=slope_model$coefficients[2],
                                          Slope.Low=confint(slope_model)[2,1],
                                          Slope.High=confint(slope_model)[2,2],
                                          Intercept.Avg=intercept_model$coefficients[1],
                                          Intercept.Low=confint(intercept_model)[1],
                                          Intercept.High=confint(intercept_model)[2])
  }
  
  # Calculate overall slope + intercept
  slope_model <- glm(y_val_diabetes~linear_pred,family="binomial")
  intercept_model <- glm(y_val_diabetes~offset(linear_pred)+1,family="binomial")
  
  # Attach overall model to df
  metric_df[nrow(metric_df)+1,] <- list(Subgroup="Overall",Model=model_name,
                                        Slope.Avg=slope_model$coefficients[2],
                                        Slope.Low=confint(slope_model)[2,1],
                                        Slope.High=confint(slope_model)[2,2],
                                        Intercept.Avg=intercept_model$coefficients[1],
                                        Intercept.Low=confint(intercept_model)[1],
                                        Intercept.High=confint(intercept_model)[2])
  
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
colnames(metric_df) <- c("Ethnicity","Model", "Slope", "Slope.Low", "Slope.High",
                         "Intercept","Intercept.Low","Intercept.High")

# Plot Slope
pd <- position_dodge(0.6) # move them .05 to the left and right
p <- ggplot(metric_df, aes(x=Ethnicity, y=Slope, colour=Model,shape = Model)) + 
  geom_errorbar(aes(ymin=Slope.Low, ymax=Slope.High), width=.1, position=pd) +
  geom_point(position=pd) +
  xlab("Ethnicity") +
  ylab("Calibration Slope") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "XGBNoSA" = 0,
                                "XGBSingleSA" = 5,
                                "XGBMultiSA" = 7),
                     name = "Model") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.5) + # Add horizontal line
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  )  +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) +
  guides(color = "none",shape = "none")

# Save
saveRDS(metric_df,"Output/Plots/Data_SupFig2b.rds")
ggsave("Output/Plots/SupFig2b.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)



#### Sup. Fig. 2c: CITL XGB Plot ####



# Plot Intercept
pd <- position_dodge(0.6) # move them .05 to the left and right
p <- ggplot(metric_df, aes(x=Ethnicity, y=Intercept, colour=Model,shape = Model)) + 
  geom_errorbar(aes(ymin=Intercept.Low, ymax=Intercept.High), width=.1, position=pd) +
  geom_point(position=pd) +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "XGBNoSA" = 0,
                                "XGBSingleSA" = 5,
                                "XGBMultiSA" = 7),
                     name = "Model") +
  xlab("Ethnicity") +
  ylab("Calibration Intercept") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) + # Add horizontal line
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  )  +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) +
  guides(color = "none",shape = "none")

# Save
saveRDS(metric_df,"Output/Plots/Data_SupFig2c.rds")
ggsave("Output/Plots/SupFig2c.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)



#### Sup. Fig. 3: sNB XGB Comparison ####



metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Metric.Avg=double(),
                       Metric.Low=double(),
                       Metric.High=double(),
                       stringsAsFactors=FALSE)

# For the Treat.No.One, we add manually with estimates of proportion
for (each_sa in c("Asian","Black","Other","White","Overall")) {
  
  # What sum and n?
  if (each_sa=="Overall") {
    x <- sum(y_val_diabetes)
    n <- length(y_val_diabetes)
  } else {
    x <- sum(y_val_diabetes[sa_val_diabetes==each_sa])
    n <- length(y_val_diabetes[sa_val_diabetes==each_sa])
  }
  
  # Get test
  proportion_test <- prop.test(x,n)
  
  # Add to dataframe
  metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,Model="Treat.No.One",
                                        Metric.Avg = 10000*(1-x/n),
                                        Metric.Low = 10000*(1-proportion_test$conf.int[2]),
                                        Metric.High = 10000*(1-proportion_test$conf.int[1]))
}

# For random lottery, we have combination of treat none and treat everyone
for (each_sa in c("Asian","Black","Other","White","Overall")) {
  
  # What sum and n?
  if (each_sa=="Overall") {
    x <- sum(y_val_diabetes)
    n <- length(y_val_diabetes)
  } else {
    x <- sum(y_val_diabetes[sa_val_diabetes==each_sa])
    n <- length(y_val_diabetes[sa_val_diabetes==each_sa])
  }
  
  # Get test
  proportion_test <- prop.test(x,n)
  
  # Get values for mean and CI
  avg <- 1 - 0.05*diabetes_threshold/(1-diabetes_threshold) + (
    0.05/(1-diabetes_threshold) - 1
  )*x/n
  lowci <- 1 - 0.05*diabetes_threshold/(1-diabetes_threshold) + (
    0.05/(1-diabetes_threshold) - 1
  )*proportion_test$conf.int[2]
  highci <- 1 - 0.05*diabetes_threshold/(1-diabetes_threshold) + (
    0.05/(1-diabetes_threshold) - 1
  )*proportion_test$conf.int[1]
  
  # Add to dataframe
  metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,Model="Random 5%",
                                        Metric.Avg = 10000*avg,
                                        Metric.Low = 10000*lowci,
                                        Metric.High = 10000*highci
  )
}

# For each model...
for (model_name in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
  
  print(model_name)
  
  # What model are we using
  if (model_name == "Treat.No.One") {
    model_predictions <- 0*y_val_diabetes
  } else if (model_name == "Treat.All") {
    model_predictions <- 0*y_val_diabetes + 1
  } else {
    model_predictions <- predictions_diabetes[[model_name]]
  }
  
  # For each subgroup...
  for (each_sa in c("Asian","Black","Other","White","Overall")) {
    
    print(each_sa)
    
    # Which group are we interested in
    if (each_sa=="Overall") {
      sel <- rep(TRUE,length(sa_val_diabetes))
    } else {
      sel <- sa_val_diabetes==each_sa
    }
    y_pred_sel <- model_predictions[sel]
    y_true_sel <- y_val_diabetes[sel]
    
    # Create a vector of all the sNB bootstrapped
    bootstrapped_sNBs <- c()
    
    for (b in 1:500) {
      
      if (b %% 50 == 0) {print(b)}
      
      # Get bootstraps
      bootstrapped_selection <- sample(1:length(y_pred_sel),length(y_pred_sel),replace=TRUE)
      y_pred_boot <- y_pred_sel[bootstrapped_selection]
      y_true_boot <- y_true_sel[bootstrapped_selection]
      
      # Append sNB
      bootstrapped_sNBs <- c(bootstrapped_sNBs,sNB_diabetes(yy_true = y_true_boot,yy_pred = y_pred_boot,pp=10000))
    }
    
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,Model=model_name,
                                          Metric.Avg = mean(bootstrapped_sNBs),
                                          Metric.Low = quantile(x = bootstrapped_sNBs,probs = 0.05),
                                          Metric.High = quantile(x = bootstrapped_sNBs,probs = 0.95))
  }
  
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
colnames(metric_df) <- c("Ethnicity","Model", "sNB", "sNB.Low", "sNB.High")

# Plot sNB
pd <- position_dodge(0.3) # move them .05 to the left and right
p <- ggplot(metric_df, aes(x=Ethnicity, y=sNB, colour=Model,shape = Model)) + 
  geom_errorbar(aes(ymin=sNB.Low, ymax=sNB.High), width=.1, position=pd) +
  geom_point(position=pd) +
  xlab("Ethnicity") +
  ylab("Subgroup Net Benefit (TNs per 10,000 Patients)") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3],
                                 "Random 5%" = clrs[4]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "XGBNoSA" = 0,
                                "XGBSingleSA" = 5,
                                "XGBMultiSA" = 7,
                                "Random 5%" = 9),
                     name = "Model") +
  geom_hline(yintercept = 10000, linetype = "dashed", color = "black", size = 0.5) + # Add horizontal line
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  )  +
#  scale_y_continuous(breaks=seq(0, 10000, 50)) + # Set y-axis breaks every 100
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  )

# Save
saveRDS(metric_df,"Output/Plots/Data_SupFig3.rds")
ggsave("Output/Plots/SupFig3.png", plot = p, width = 8, height = 4, units = "in", dpi = 400)


#### Sup. Fig. 4: DCA for XGB ####



# Which thresholds do we plot?
thresholds <- seq(0,0.30,0.005)

# Net benefit metric
metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Clinical.Threshold=double(),
                       NB=double(),
                       stringsAsFactors=FALSE)

# For each model...
for (model_name in c("Treat.No.One","Treat.All","XGBNoSA","XGBSingleSA","XGBMultiSA")) {
  
  # What model are we using
  if (model_name == "Treat.No.One") {
    model_predictions <- 0*y_val_diabetes
  } else if (model_name == "Treat.All") {
    model_predictions <- 0*y_val_diabetes + 1
  } else {
    model_predictions <- predictions_diabetes[[model_name]]
  }
  
  # For each subgroup...
  for (each_sa in c("Asian","Black","Other","White","Overall")) {
    
    # Which group are we interested in
    if (each_sa=="Overall") {
      sel <- rep(TRUE,length(sa_val_diabetes))
    } else {
      sel <- sa_val_diabetes==each_sa
    }
    
    # For each threshold...
    nb <- 10000*calcNetBenefit(clinical_thresholds = thresholds,model_thresholds = thresholds,
                               y_true = y_val_diabetes[sel],s_pred = model_predictions[sel],new_version = FALSE)
    
    # Append to dataframe
    metric_df <- rbind(metric_df,
                       data.frame(Subgroup=each_sa,
                                  Model=model_name, 
                                  Clinical.Threshold=thresholds,
                                  NB=nb,
                                  stringsAsFactors=FALSE))
  }
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
saveRDS(metric_df,"Output/Plots/Data_SupFig4.rds")

# Figure 2a
df <- metric_df %>% filter(Subgroup == "Overall") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "XGBNoSA" = "solid",
                                   "XGBSingleSA" = "solid",
                                   "XGBMultiSA" = "solid"),
                        name = "Model") +
  coord_cartesian(ylim= c(0,NA)) +
  xlab("Threshold") +
  ylab("Net Benefit (TPs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) +
  geom_vline(xintercept = diabetes_threshold, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("Overall")

# Save
ggsave("Output/Plots/SupFig4a.png", plot = p, width = 8, height = 4, units = "in", dpi = 400)



# Figure 2b
df <- metric_df %>% filter(Subgroup == "Asian") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "XGBNoSA" = "solid",
                                   "XGBSingleSA" = "solid",
                                   "XGBMultiSA" = "solid"),
                        name = "Model") +
  coord_cartesian(ylim= c(0,NA)) +
  xlab("Threshold") +
  ylab("Net Benefit (TPs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) + guides(color = "none",shape = "none") + theme(legend.position="none") +
  geom_vline(xintercept = diabetes_threshold, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("Asian")

# Save
ggsave("Output/Plots/SupFig4b.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)


# Figure 2c
df <- metric_df %>% filter(Subgroup == "Black") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "XGBNoSA" = "solid",
                                   "XGBSingleSA" = "solid",
                                   "XGBMultiSA" = "solid"),
                        name = "Model") +
  coord_cartesian(ylim= c(0,NA)) +
  xlab("Threshold") +
  ylab("Net Benefit (TPs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) + guides(color = "none",shape = "none") + theme(legend.position="none") +
  geom_vline(xintercept = diabetes_threshold, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("Black")

# Save
ggsave("Output/Plots/SupFig4c.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)


# Figure 2d
df <- metric_df %>% filter(Subgroup == "Other") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "XGBNoSA" = "solid",
                                   "XGBSingleSA" = "solid",
                                   "XGBMultiSA" = "solid"),
                        name = "Model") +
  coord_cartesian(ylim= c(0,NA)) +
  xlab("Threshold") +
  ylab("Net Benefit (TPs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) + guides(color = "none",shape = "none") + theme(legend.position="none") +
  geom_vline(xintercept = diabetes_threshold, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("Other")

# Save
ggsave("Output/Plots/SupFig4d.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)



# Figure 2e
df <- metric_df %>% filter(Subgroup == "White") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "XGBNoSA" = "solid",
                                   "XGBSingleSA" = "solid",
                                   "XGBMultiSA" = "solid"),
                        name = "Model") +
  coord_cartesian(ylim= c(0,NA)) +
  xlab("Threshold") +
  ylab("Net Benefit (TPs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) + guides(color = "none",shape = "none") + theme(legend.position="none") +
  geom_vline(xintercept = diabetes_threshold, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("White")

# Save
ggsave("Output/Plots/SupFig4e.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)



### Second section: Lung cancer ####



#### Sup. Tab. 2: LungCancer table ####



# Open dataset
lungcancer_table_df <- read.csv(
  file = paste0("Data/lungcancer_df",which_output_path,".csv")
)
rownames(lungcancer_table_df) <- lungcancer_table_df$X
lungcancer_table_df <- lungcancer_table_df %>% select(-X)
lungcancer_table_df$Ethnicity <- as.factor(lungcancer_table_df$Ethnicity)
lungcancer_table_df$Education <- as.factor(lungcancer_table_df$Education)
lungcancer_table_df$Smoker <- as.factor(lungcancer_table_df$Smoker)
quintiles_townsend <- quantile(lungcancer_table_df$Townsend,probs = c(0.2,0.4,0.6,0.8),na.rm=TRUE)
lungcancer_table_df$Townsend.Quintile <- define_lungcancer_sa(lungcancer_table_df$Townsend,quintiles_townsend)
lungcancer_table_df$Townsend.Quintile <- as.character(lungcancer_table_df$Townsend.Quintile)
lungcancer_table_df$Townsend.Quintile[lungcancer_table_df$Townsend.Quintile=="0"] <- "1st (Lowest Deprivation)"
lungcancer_table_df$Townsend.Quintile[lungcancer_table_df$Townsend.Quintile=="1"] <- "2nd"
lungcancer_table_df$Townsend.Quintile[lungcancer_table_df$Townsend.Quintile=="2"] <- "3rd"
lungcancer_table_df$Townsend.Quintile[lungcancer_table_df$Townsend.Quintile=="3"] <- "4th"
lungcancer_table_df$Townsend.Quintile[lungcancer_table_df$Townsend.Quintile=="4"] <- "5th (Highest Deprivation)"
lungcancer_table_df$Townsend.Quintile <- as.factor(lungcancer_table_df$Townsend.Quintile)

# Only choose columns we care about
#lungcancer_table_df <- lungcancer_table_df %>%
#  select(c(Age,Female,Townsend,Ethnicity,Education,FamilyHxOtherCancer,FamilyHxLungCancer,
#           FamilyHxChronicBronchitis, Smoker, MedicationCholesterol,MedicationBloodPressure,
#           FirstIncidenceJ43to44,AnyCancer,y,Townsend.Quintile))
# If using synthetic data example, we don't have all of these variables
# so we plot simplified table
lungcancer_table_df <- lungcancer_table_df %>%
    select(c(Age,Female,Townsend,Ethnicity,Education,
             Smoker, MedicationCholesterol,MedicationBloodPressure,
             y,Townsend.Quintile))

# Create a table of data for supplementary
table_summary <- lungcancer_table_df %>%
  tbl_summary(by = Townsend.Quintile)
table_summary



#### Fig
#### Sup. Fig. 6a: XGB AUROC table ####

metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Metric.Avg=double(),
                       Metric.Low=double(),
                       Metric.High=double(),
                       stringsAsFactors=FALSE)

# For each model, get the AUC (with 95% CI) of each subgroup
for (model_name in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
  
  # What model are we using
  model_predictions <- predictions_lungcancer[[model_name]]
  
  # For each subgroup, calculate subgroup AUC
  for (each_sa in c("0","1","2","3","4")) {
    
    # Calculate subgroup AUC
    auroc <- roc(y_val_lungcancer[sa_val_lungcancer==each_sa],
                 model_predictions[sa_val_lungcancer==each_sa],
                 ci=TRUE, plot=FALSE)
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,Model=model_name,
                                          Metric.Avg = auroc$auc,
                                          Metric.Low = auroc$ci[1],
                                          Metric.High = auroc$ci[3])
  }
  
  # Calculate overall AUC
  auroc <- roc(y_val_lungcancer,model_predictions, ci=TRUE, plot=FALSE)
  
  # Attach overall model to df
  metric_df[nrow(metric_df)+1,] <- list(Subgroup="Overall",Model=model_name,
                                        Metric.Avg = auroc$auc,
                                        Metric.Low = auroc$ci[1],
                                        Metric.High = auroc$ci[3])
  
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
levels(metric_df$Subgroup) <- c("Q1","Q2","Q3","Q4","Q5","Overall")
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
colnames(metric_df) <- c("Townsend.Quintile","Model", "AUROC", "AUROC.Low", "AUROC.High")

# Plot
pd <- position_dodge(0.6) # move them .05 to the left and right
p <- ggplot(metric_df, aes(x=Townsend.Quintile, y=AUROC, colour=Model,shape = Model)) + 
  geom_errorbar(aes(ymin=AUROC.Low, ymax=AUROC.High), width=.1, position=pd) +
  geom_point(position=pd) +
  xlab("Townsend.Quintile") +
  ylab("C-Statistic") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "XGBNoSA" = 0,
                                "XGBSingleSA" = 5,
                                "XGBMultiSA" = 7),
                     name = "Model") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  )  +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  )

# Save
saveRDS(metric_df,"Output/Plots/Data_SupFig6a.rds")
ggsave("Output/Plots/SupFig6a.png", plot = p, width = 8, height = 4, units = "in", dpi = 400)



#### Sup. Fig. 6b: XGB Slope plot ####



metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Slope.Avg=double(),
                       Slope.Low=double(),
                       Slope.High=double(),
                       Intercept.Avg=double(),
                       Intercept.Low=double(),
                       Intercept.High=double(),
                       stringsAsFactors=FALSE)

# For each model, get the AUC (with 95% CI) of each subgroup
for (model_name in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
  
  # What model are we using
  model_predictions <- predictions_lungcancer[[model_name]]
  linear_pred <- -1*log(
    (1-model_predictions)/
      model_predictions
  )
  
  # For each subgroup, calculate subgroup AUC
  for (each_sa in c("0","1","2","3","4")) {
    
    # Which group are we interested in
    sel <- sa_val_lungcancer==each_sa
    
    # Calculate subgroup slope + intercept
    slope_model <- glm(y_val_lungcancer[sel]~linear_pred[sel],family="binomial")
    intercept_model <- glm(y_val_lungcancer[sel]~offset(linear_pred[sel])+1,family="binomial")
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,Model=model_name,
                                          Slope.Avg=slope_model$coefficients[2],
                                          Slope.Low=confint(slope_model)[2,1],
                                          Slope.High=confint(slope_model)[2,2],
                                          Intercept.Avg=intercept_model$coefficients[1],
                                          Intercept.Low=confint(intercept_model)[1],
                                          Intercept.High=confint(intercept_model)[2])
  }
  
  # Calculate overall slope + intercept
  slope_model <- glm(y_val_lungcancer~linear_pred,family="binomial")
  intercept_model <- glm(y_val_lungcancer~offset(linear_pred)+1,family="binomial")
  
  # Attach overall model to df
  metric_df[nrow(metric_df)+1,] <- list(Subgroup="Overall",Model=model_name,
                                        Slope.Avg=slope_model$coefficients[2],
                                        Slope.Low=confint(slope_model)[2,1],
                                        Slope.High=confint(slope_model)[2,2],
                                        Intercept.Avg=intercept_model$coefficients[1],
                                        Intercept.Low=confint(intercept_model)[1],
                                        Intercept.High=confint(intercept_model)[2])
  
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
levels(metric_df$Subgroup) <- c("Q1","Q2","Q3","Q4","Q5","Overall")
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
colnames(metric_df) <- c("Townsend.Quintile","Model", "Slope", "Slope.Low", "Slope.High",
                         "Intercept","Intercept.Low","Intercept.High")

# Plot Slope
pd <- position_dodge(0.6) # move them .05 to the left and right
p <- ggplot(metric_df, aes(x=Townsend.Quintile, y=Slope, colour=Model,shape = Model)) + 
  geom_errorbar(aes(ymin=Slope.Low, ymax=Slope.High), width=.1, position=pd) +
  geom_point(position=pd) +
  xlab("Townsend.Quintile") +
  ylab("Calibration Slope") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "XGBNoSA" = 0,
                                "XGBSingleSA" = 5,
                                "XGBMultiSA" = 7),
                     name = "Model") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.5) + # Add horizontal line
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  )  +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) +
  guides(color = "none",shape = "none")

# Save
saveRDS(metric_df,"Output/Plots/Data_SupFig6b.rds")
ggsave("Output/Plots/SupFig6b.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)



#### Sup. Fig. 6c: XGB CITL Plot ####



# Plot Intercept
pd <- position_dodge(0.6) # move them .05 to the left and right
p <- ggplot(metric_df, aes(x=Townsend.Quintile, y=Intercept, colour=Model,shape = Model)) + 
  geom_errorbar(aes(ymin=Intercept.Low, ymax=Intercept.High), width=.1, position=pd) +
  geom_point(position=pd) +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "XGBNoSA" = 0,
                                "XGBSingleSA" = 5,
                                "XGBMultiSA" = 7),
                     name = "Model") +
  xlab("Townsend.Quintile") +
  ylab("Calibration Intercept") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) + # Add horizontal line
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  )  +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) +
  guides(color = "none",shape = "none")

# Save
saveRDS(metric_df,"Output/Plots/Data_SupFig6c.rds")
ggsave("Output/Plots/SupFig6c.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)



#### Sup. Fig 7: XGB sNB Comparison Lung ####



metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Metric.Avg=double(),
                       Metric.Low=double(),
                       Metric.High=double(),
                       stringsAsFactors=FALSE)

# For the Treat.No.One, we add manually with estimates of proportion
for (each_sa in c("0","1","2","3","4","Overall")) {
  
  # What sum and n?
  if (each_sa=="Overall") {
    x <- sum(y_val_lungcancer)
    n <- length(y_val_lungcancer)
  } else {
    x <- sum(y_val_lungcancer[sa_val_lungcancer==each_sa])
    n <- length(y_val_lungcancer[sa_val_lungcancer==each_sa])
  }
  
  # Get test
  proportion_test <- prop.test(x,n)
  
  # Add to dataframe
  metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,Model="Treat.No.One",
                                        Metric.Avg = 10000*(1-x/n),
                                        Metric.Low = 10000*(1-proportion_test$conf.int[2]),
                                        Metric.High = 10000*(1-proportion_test$conf.int[1]))
}

# For random lottery, we have combination of treat none and treat everyone
for (each_sa in c("0","1","2","3","4","Overall")) {
  
  # What sum and n?
  if (each_sa=="Overall") {
    x <- sum(y_val_lungcancer)
    n <- length(y_val_lungcancer)
  } else {
    x <- sum(y_val_lungcancer[sa_val_lungcancer==each_sa])
    n <- length(y_val_lungcancer[sa_val_lungcancer==each_sa])
  }
  
  # Get test
  proportion_test <- prop.test(x,n)
  
  # Get values for mean and CI
  avg <- 1 - 0.05*lung_cancer_threshold/(1-lung_cancer_threshold) + (
    0.05/(1-lung_cancer_threshold) - 1
  )*x/n
  lowci <- 1 - 0.05*lung_cancer_threshold/(1-lung_cancer_threshold) + (
    0.05/(1-lung_cancer_threshold) - 1
  )*proportion_test$conf.int[2]
  highci <- 1 - 0.05*lung_cancer_threshold/(1-lung_cancer_threshold) + (
    0.05/(1-lung_cancer_threshold) - 1
  )*proportion_test$conf.int[1]
  
  # Add to dataframe
  metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,Model="Random 5%",
                                        Metric.Avg = 10000*avg,
                                        Metric.Low = 10000*lowci,
                                        Metric.High = 10000*highci
  )
}


# For each model...
for (model_name in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
  
  print(model_name)
  
  # What model are we using
  if (model_name == "Treat.No.One") {
    model_predictions <- 0*y_val_lungcancer
  } else if (model_name == "Treat.All") {
    model_predictions <- 0*y_val_lungcancer + 1
  } else {
    model_predictions <- predictions_lungcancer[[model_name]]
  }
  
  # For each subgroup...
  for (each_sa in c("0","1","2","3","4","Overall")) {
    
    print(each_sa)
    
    # Which group are we interested in
    if (each_sa=="Overall") {
      sel <- rep(TRUE,length(sa_val_lungcancer))
    } else {
      sel <- sa_val_lungcancer==each_sa
    }
    y_pred_sel <- model_predictions[sel]
    y_true_sel <- y_val_lungcancer[sel]
    
    # Create a vector of all the sNB bootstrapped
    bootstrapped_sNBs <- c()
    
    for (b in 1:500) {
      
      if (b %% 50 == 0) {print(b)}
      
      # Get bootstraps
      bootstrapped_selection <- sample(1:length(y_pred_sel),length(y_pred_sel),replace=TRUE)
      y_pred_boot <- y_pred_sel[bootstrapped_selection]
      y_true_boot <- y_true_sel[bootstrapped_selection]
      
      # Append sNB
      bootstrapped_sNBs <- c(bootstrapped_sNBs,sNB_lungcancer(yy_true = y_true_boot,yy_pred = y_pred_boot,pp=10000))
    }
    
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,Model=model_name,
                                          Metric.Avg = mean(bootstrapped_sNBs),
                                          Metric.Low = quantile(x = bootstrapped_sNBs,probs = 0.05),
                                          Metric.High = quantile(x = bootstrapped_sNBs,probs = 0.95))
  }
  
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
levels(metric_df$Subgroup) <- c("Q1","Q2","Q3","Q4","Q5","Overall")
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
colnames(metric_df) <- c("Townsend","Model", "sNB", "sNB.Low", "sNB.High")

# Plot sNB
pd <- position_dodge(0.3) # move them .05 to the left and right
p <- ggplot(metric_df, aes(x=Townsend, y=sNB, colour=Model,shape = Model)) + 
  geom_errorbar(aes(ymin=sNB.Low, ymax=sNB.High), width=.1, position=pd) +
  geom_point(position=pd) +
  xlab("Townsend") +
  ylab("Subgroup Net Benefit (TNs per 10,000 Patients)") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3],
                                 "Random 5%" = clrs[4]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "XGBNoSA" = 0,
                                "XGBSingleSA" = 5,
                                "XGBMultiSA" = 7,
                                "Random 5%" = 9),
                     name = "Model") +
  geom_hline(yintercept = 10000, linetype = "dashed", color = "black", size = 0.5) + # Add horizontal line
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  )  +
#  scale_y_continuous(breaks=seq(0, 10000, 50)) + # Set y-axis breaks every 100
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  )

# Save
saveRDS(metric_df,"Output/Plots/Data_SupFig7.rds")
ggsave("Output/Plots/SupFig7.png", plot = p, width = 8, height = 4, units = "in", dpi = 400)


#### Sup. Fig. 8: XGB DCA ####



# Which thresholds do we plot?
thresholds <- seq(0,0.05,0.0005)

# Net benefit metric
metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Clinical.Threshold=double(),
                       NB=double(),
                       stringsAsFactors=FALSE)

# For each model...
for (model_name in c("Treat.No.One","Treat.All","XGBNoSA","XGBSingleSA","XGBMultiSA","XGBNoSA","XGBSingleSA","XGBMultiSA")) {
  
  # What model are we using
  if (model_name == "Treat.No.One") {
    model_predictions <- 0*y_val_lungcancer
  } else if (model_name == "Treat.All") {
    model_predictions <- 0*y_val_lungcancer + 1
  } else {
    model_predictions <- predictions_lungcancer[[model_name]]
  }
  
  # For each subgroup...
  for (each_sa in c("0","1","2","3","4","Overall")) {
    
    # Which group are we interested in
    if (each_sa=="Overall") {
      sel <- rep(TRUE,length(sa_val_lungcancer))
    } else {
      sel <- sa_val_lungcancer==each_sa
    }
    
    # For each threshold...
    nb <- 10000*calcNetBenefit(clinical_thresholds = thresholds,model_thresholds = thresholds,
                               y_true = y_val_lungcancer[sel],s_pred = model_predictions[sel],new_version = FALSE)
    
    # Append to dataframe
    metric_df <- rbind(metric_df,
                       data.frame(Subgroup=each_sa,
                                  Model=model_name, 
                                  Clinical.Threshold=thresholds,
                                  NB=nb,
                                  stringsAsFactors=FALSE))
  }
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
levels(metric_df$Subgroup) <- c("Q1","Q2","Q3","Q4","Q5","Overall")
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
saveRDS(metric_df,"Output/Plots/Data_SupFig8.rds")

# Figure 2a
df <- metric_df %>% filter(Subgroup == "Overall") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "XGBNoSA" = "solid",
                                   "XGBSingleSA" = "solid",
                                   "XGBMultiSA" = "solid"),
                        name = "Model") +
  coord_cartesian(ylim= c(0,NA)) +
  xlab("Threshold") +
  ylab("Net Benefit (TPs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) +
  geom_vline(xintercept = 0.015, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("Overall")

# Save
ggsave("Output/Plots/SupFig8a.png", plot = p, width = 8, height = 4, units = "in", dpi = 400)



# Figure 2b
df <- metric_df %>% filter(Subgroup == "Q1") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "XGBNoSA" = "solid",
                                   "XGBSingleSA" = "solid",
                                   "XGBMultiSA" = "solid"),
                        name = "Model") +
  coord_cartesian(ylim= c(0,NA)) +
  xlab("Threshold") +
  ylab("Net Benefit (TPs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) + guides(color = "none",shape = "none") + theme(legend.position="none") +
  geom_vline(xintercept = 0.015, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("1st Quintile (Lowest)")

# Save
ggsave("Output/Plots/SupFig8b.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)


# Figure 2c
df <- metric_df %>% filter(Subgroup == "Q2") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "XGBNoSA" = "solid",
                                   "XGBSingleSA" = "solid",
                                   "XGBMultiSA" = "solid"),
                        name = "Model") +
  coord_cartesian(ylim= c(0,NA)) +
  xlab("Threshold") +
  ylab("Net Benefit (TPs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) + guides(color = "none",shape = "none") + theme(legend.position="none") +
  geom_vline(xintercept = 0.015, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("2nd Quintile")

# Save
ggsave("Output/Plots/SupFig8c.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)


# Figure 2d
df <- metric_df %>% filter(Subgroup == "Q3") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "XGBNoSA" = "solid",
                                   "XGBSingleSA" = "solid",
                                   "XGBMultiSA" = "solid"),
                        name = "Model") +
  coord_cartesian(ylim= c(0,NA)) +
  xlab("Threshold") +
  ylab("Net Benefit (TPs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) + guides(color = "none",shape = "none") + theme(legend.position="none") +
  geom_vline(xintercept = 0.015, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("3rd Quintile")

# Save
ggsave("Output/Plots/SupFig8d.png", plot = p, width = 4*2/3, height = 4, units = "in", dpi = 400)



# Figure 2e
df <- metric_df %>% filter(Subgroup == "Q4") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "XGBNoSA" = "solid",
                                   "XGBSingleSA" = "solid",
                                   "XGBMultiSA" = "solid"),
                        name = "Model") +
  coord_cartesian(ylim= c(0,NA)) +
  xlab("Threshold") +
  ylab("Net Benefit (TPs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) + guides(color = "none",shape = "none") + theme(legend.position="none") +
  geom_vline(xintercept = 0.015, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("4th Quintile")

# Save
ggsave("Output/Plots/SupFig8e.png", plot = p, width = 4*2/3, height = 4, units = "in", dpi = 400)

# Figure 2f

df <- metric_df %>% filter(Subgroup == "Q5") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "XGBNoSA" = "solid",
                                   "XGBSingleSA" = "solid",
                                   "XGBMultiSA" = "solid"),
                        name = "Model") +
  coord_cartesian(ylim= c(0,NA)) +
  xlab("Threshold") +
  ylab("Net Benefit (TPs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 10), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) + guides(color = "none",shape = "none") + theme(legend.position="none") +
  geom_vline(xintercept = 0.015, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("5th Quintile (Highest)")

# Save
ggsave("Output/Plots/SupFig8f.png", plot = p, width = 4*2/3, height = 4, units = "in", dpi = 400)




#### Sup. Fig 9: XGB sNB Resource constraints ####

# Currently this only 'works' on real data.
# As it is specific to the incidence and predictive nature of the
# real lung cancer example.

# First we need to create training predictions
predictions_lungcancer_training <- list()

# 1) XGBNoSA
pred_train <- predict(models_lungcancer$XGBNoSA$Model,x_train_nosa_lungcancer)
predictions_lungcancer_training$XGBNoSA <- as.numeric(pred_train)

# 2) XGBSingleSA
pred_train <- predict(models_lungcancer$XGBSingleConv$Model,x_train_lungcancer)
predictions_lungcancer_training$XGBSingleSA <- as.numeric(pred_train)

# 3) XGBMultiSA
pred_train <- 0*y_train_lungcancer
for (each_sa in unique(sa_train_lungcancer)) {
  
  # For each group, calculate the predictions for that group in train and train
  pred_train_per_sa <- predict(models_lungcancer$XGBMultiConv[[each_sa]]$Model,x_train_lungcancer[sa_train_lungcancer==each_sa,], type = "response")
  
  # Save onto pred_train
  pred_train[sa_train_lungcancer==each_sa] <- as.numeric(pred_train_per_sa)
}
predictions_lungcancer_training$XGBMultiSA <- as.numeric(pred_train)

# Get new grouping: lowest quintile?
grouping_lc_train <- as.factor(sa_train_lungcancer=="4")
grouping_lc_val <- as.factor(sa_val_lungcancer=="4")
two_groups <- c("Q1234","Q5")
levels(grouping_lc_train) <- two_groups
levels(grouping_lc_val) <- two_groups


### !!! HEAVY COMPUTATION

# How many points to explore
resolution <- 1000

# Save plots per capacity and per model
threshold_grid_per_model_per_capacity <- list()
bnft_grid_per_model_per_capacity <- list()
bnft_grid_per_model_per_capacity_validation <- list()

ii <- 1

for (capacity in c(0.10,0.03,0.01)) {
  
  print(capacity)
  
  # When do we reach that capacity in the LogSingleModel?  
  top_quantile <- 2*(quantile(predictions_lungcancer_training[["XGBSingleSA"]],capacity))
  # Which thresholds do we explore
  threshold_options <- seq(lung_cancer_threshold,0.2, length.out = resolution)
  
  # Total number of grid options?
  threshold_grid <- expand.grid(threshold_options,threshold_options)
  threshold_grid_per_model <- list()
  for (model in c("XGBSingleSA")) {threshold_grid_per_model[[model]] <- threshold_grid}
  
  # Group weights: What proportion is each weight in the population?
  sa_numbers <- list()
  sa_weights <- list()
  for (each_sa in two_groups) {
    sa_numbers[[each_sa]] <- sum(grouping_lc_train==each_sa)
    sa_weights[[each_sa]] <- mean(grouping_lc_train==each_sa)
  }
  
  # What is the benefit in each group for each threshold option?
  bnft_for_each_group <- list()
  bnft_grid_per_model <- list()
  capacity_for_each_group <- list()
  capacity_grid_per_model <- list()
  for (model in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
    
    print(model)
    
    bnft_for_each_group[[model]] <- list()
    capacity_for_each_group[[model]] <- list()
    
    for (each_sa in two_groups) {
      
      bnft_for_each_group[[model]][[each_sa]] <- c()
      capacity_for_each_group[[model]][[each_sa]] <- c()
      
      # Which group are we interested in
      sel <- (grouping_lc_train==each_sa)
      # who to select?
      y_pred_sel <- predictions_lungcancer_training[[model]][sel]
      y_true_sel <- y_train_lungcancer[sel]
      
      for (tt in threshold_options) {
        
        # What is NB (raw)?
        raw_nb <- calcNetBenefit(lung_cancer_threshold,tt,y_true_sel,y_pred_sel,new_version=TRUE)
        
        # Add to list
        bnft_for_each_group[[model]][[each_sa]] <- c(
          bnft_for_each_group[[model]][[each_sa]],
          (1 - mean(y_true_sel) + lung_cancer_effect*raw_nb/balance_lungcancer)*10000
        )
        
        # What is the number of positives at that threshold?
        capacity_for_each_group[[model]][[each_sa]] <- c(
          capacity_for_each_group[[model]][[each_sa]],
          sum(y_pred_sel>tt)
        )
        
      }
    }
    
    # What benefits correspond to each group?
    grid_df <- expand.grid(bnft_for_each_group[[model]][["Q1234"]],
                           bnft_for_each_group[[model]][["Q5"]])
    colnames(grid_df) <- c("Q1234","Q5")
    bnft_grid_per_model[[model]] <- grid_df
    
    
    # What capacity corresponds to each group?
    capacity_df <- expand.grid(capacity_for_each_group[[model]][["Q1234"]],
                               capacity_for_each_group[[model]][["Q5"]])
    colnames(capacity_df) <- c("Q1234","Q5")
    capacity_grid_per_model[[model]] <- capacity_df
    
  }
  
  # For each model, change df to have min and overall
  for (model in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
    
    grid_df <- bnft_grid_per_model[[model]]
    grid_df <- grid_df %>% 
      mutate(Overall = (sa_weights[["Q1234"]]*Q1234 +
                          sa_weights[["Q5"]]*Q5),
             Minimum = pmin(Q1234,Q5))
    bnft_grid_per_model[[model]] <- grid_df
    
    # For capacity, get maximum capacity
    capacity_df <- capacity_grid_per_model[[model]]
    capacity_df <- capacity_df %>% 
      mutate(Total.Positives = (Q1234+Q5))
    capacity_grid_per_model[[model]] <- capacity_df
  }
  
  # We only want models for which tot.positives is within capacity
  for (model in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {possible_choices <- capacity_grid_per_model[[model]]$Total.Positives < capacity*length(y_train_lungcancer)
  threshold_grid_per_model[[model]] <- threshold_grid_per_model[[model]][possible_choices,]
  bnft_grid_per_model[[model]] <- bnft_grid_per_model[[model]][possible_choices,]
  capacity_grid_per_model[[model]] <- capacity_grid_per_model[[model]][possible_choices,]}
  
  # Get Pareto frontier for each model
  pareto_frontier_per_model <- list()
  pareto_indices_per_model <- list()
  for (model in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
    pareto_frontier <- get_frontier(bnft_grid_per_model[[model]],
                                    x=Minimum,y=Overall,quadrant = "top.right")
    pareto_frontier_per_model[[model]] <- pareto_frontier
    pareto_indices <- rownames(pareto_frontier)
    pareto_indices_per_model[[model]] <- pareto_indices
  }
  
  
  # What is the benefit in each group IN VALIDATION for each threshold option?
  sa_weights_validation <- list()
  for (each_sa in two_groups) {
    sa_weights_validation[[each_sa]] <- mean(grouping_lc_val==each_sa)
  }
  bnft_for_each_group_validation <- list()
  bnft_grid_per_model_validation <- list()
  capacity_for_each_group_validation <- list()
  capacity_grid_per_model_validation <- list()
  for (model in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
    
    bnft_for_each_group_validation[[model]] <- list()
    capacity_for_each_group_validation[[model]] <- list()
    
    for (each_sa in two_groups) {
      
      bnft_for_each_group_validation[[model]][[each_sa]] <- c()
      capacity_for_each_group_validation[[model]][[each_sa]] <- c()
      
      # Which group are we interested in
      sel <- (grouping_lc_val==each_sa)
      # who to select?
      y_pred_sel <- predictions_lungcancer[[model]][sel]
      y_true_sel <- y_val_lungcancer[sel]
      
      for (tt in threshold_options) {
        
        # What is NB (raw)?
        raw_nb <- calcNetBenefit(lung_cancer_threshold,tt,y_true_sel,y_pred_sel,new_version=TRUE)
        
        # Add to list
        bnft_for_each_group_validation[[model]][[each_sa]] <- c(
          bnft_for_each_group_validation[[model]][[each_sa]],
          (1 - mean(y_true_sel) + lung_cancer_effect*raw_nb/balance_lungcancer)*10000
        )
        
        # What is the number of positives at that threshold?
        capacity_for_each_group_validation[[model]][[each_sa]] <- c(
          capacity_for_each_group_validation[[model]][[each_sa]],
          sum(y_pred_sel>tt)
        )
        
      }
    }
    
    # What benefits correspond to each group?
    grid_df <- expand.grid(bnft_for_each_group_validation[[model]][["Q1234"]],
                           bnft_for_each_group_validation[[model]][["Q5"]])
    colnames(grid_df) <- c("Q1234","Q5")
    bnft_grid_per_model_validation[[model]] <- grid_df
    
    
    # What capacity corresponds to each group?
    capacity_df <- expand.grid(capacity_for_each_group_validation[[model]][["Q1234"]],
                               capacity_for_each_group_validation[[model]][["Q5"]])
    colnames(capacity_df) <- c("Q1234","Q5")
    capacity_grid_per_model_validation[[model]] <- capacity_df
    
  }
  # For each model, change df to have min and overall
  for (model in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
    
    grid_df <- bnft_grid_per_model_validation[[model]]
    grid_df <- grid_df %>% 
      mutate(Overall = (sa_weights_validation[["Q1234"]]*Q1234 +
                          sa_weights_validation[["Q5"]]*Q5),
             Minimum = pmin(Q1234,Q5))
    bnft_grid_per_model_validation[[model]] <- grid_df
    
    # For capacity, get maximum capacity
    capacity_df <- capacity_grid_per_model_validation[[model]]
    capacity_df <- capacity_df %>% 
      mutate(Total.Positives = (Q1234+Q5))
    capacity_grid_per_model_validation[[model]] <- capacity_df
  }
  
  # Initialise empty lists
  bnft_grid_per_model_per_capacity[[ii]] <- list()
  bnft_grid_per_model_per_capacity_validation[[ii]] <- list()
  threshold_grid_per_model_per_capacity[[ii]] <- list()
  
  for (model in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
    #What is the plot?
    bnft_grid_per_model_per_capacity[[ii]][[model]] <- bnft_grid_per_model[[model]][pareto_indices_per_model[[model]],]
    
    # In validation?
    bnft_grid_per_model_per_capacity_validation[[ii]][[model]] <- bnft_grid_per_model_validation[[model]][pareto_indices_per_model[[model]],]
    
    # What thresholds are these?
    threshold_grid_per_model_per_capacity[[ii]][[model]] <- threshold_grid[pareto_indices_per_model[[model]],]
  }
  
  ii <- ii +1
}

# Save data file
fig_5_list <- list()
fig_5_list[["Training"]] <- bnft_grid_per_model_per_capacity
fig_5_list[["Validation"]] <- bnft_grid_per_model_per_capacity_validation
fig_5_list[["Thresholds"]] <- threshold_grid_per_model_per_capacity
saveRDS(fig_5_list,"Output/Plots/Data_Fig5.rds")

# Reload
fig_5_list <- readRDS("Output/Plots/Data_Fig5.rds")

# To plot Figure 5a, we take the training data and add it to a df file.
fig5_df <- list()
for (capacity_ii in c(1,2,3)) {
  for (model in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
    df_part <- fig_5_list$Training[[capacity_ii]][[model]][c("Overall","Minimum")]
    df_part$Capacity <- c("No Constraints","3%","1%")[capacity_ii]
    df_part$Model <- model
    fig5_df[[length(fig5_df)+1]] <- df_part
  }
}
fig5_df <- Reduce(full_join,fig5_df)
fig5_df$Subgroup <- factor(fig5_df$Model, levels = c("XGBNoSA","XGBSingleSA","XGBMultiSA"))

# What is the benefit when no patient is treated?
overall_treat_no_oone <- (1-mean(y_train_lungcancer==1))*10000
minimum_treat_no_one <- (1-mean(y_train_lungcancer[grouping_lc_train=="Q5"]))*10000

# Plot
p <- ggplot(fig5_df, aes(x = Minimum, y = Overall, colour = Model, linetype = Model, shape = Capacity)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c("XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Max. Capacity") +
  scale_linetype_manual(values = c("XGBNoSA" = 2,
                                   "XGBSingleSA" = 3,
                                   "XGBMultiSA" = 4),
                        name = "Model")  +
  geom_point(x = minimum_treat_no_one, y = overall_treat_no_oone, color = "#333a52", size = 1.5) +
  geom_text(x = minimum_treat_no_one, y = overall_treat_no_oone, label = "Treat No One",
            vjust = -0.6, hjust = 0, size = 3, color = "#333a52") +
  xlab("Net Benefit in Top Quintile of Deprivation (TNs per 10,000 patients)") +
  ylab("Overall Net Benefit (TNs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 8), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 8), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "transparent") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) +
  ggtitle("Training") +
  coord_cartesian(xlim = c(9933, 9940),
                  ylim = c(9960, 9964)) +
  theme(legend.position="none")

# Save
ggsave("Output/Plots/SupFig9a.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)



# For figure 5b, use validation data



fig5_df <- list()
for (capacity_ii in c(1,2,3)) {
  for (model in c("XGBNoSA","XGBSingleSA","XGBMultiSA")) {
    df_part <- fig_5_list$Validation[[capacity_ii]][[model]][c("Overall","Minimum")]
    df_part$Capacity <- c("No Constraints","3%","1%")[capacity_ii]
    df_part$Model <- model
    fig5_df[[length(fig5_df)+1]] <- df_part
  }
}
fig5_df <- Reduce(full_join,fig5_df)
fig5_df$Model <- factor(fig5_df$Model, levels = c("XGBNoSA","XGBSingleSA","XGBMultiSA"))

# What is the benefit when no patient is treated?
overall_treat_no_oone <- (1-mean(y_train_lungcancer==1))*10000
minimum_treat_no_one <- (1-mean(y_train_lungcancer[grouping_lc_train=="Q5"]))*10000

# Plot
p <- ggplot(fig5_df, aes(x = Minimum, y = Overall, colour = Model, linetype = Model, shape = Capacity)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c("XGBNoSA" = clrs[1],
                                 "XGBSingleSA" = clrs[2],
                                 "XGBMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("XGBNoSA" = 2,
                                   "XGBSingleSA" = 3,
                                   "XGBMultiSA" = 4),
                        name = "Model")  +
  coord_cartesian(xlim = c(9931, 9936),
                  ylim = c(9960, 9965)) +
  geom_point(x = minimum_treat_no_one, y = overall_treat_no_oone, color = "#333a52", size = 1.5) +
  geom_text(x = minimum_treat_no_one, y = overall_treat_no_oone, label = "Treat No One",
            vjust = -0.6, hjust = 0, size = 3, color = "#333a52") +
  xlab("Net Benefit in Top Quintile of Deprivation (TNs per 10,000 patients)") +
  ylab("Overall Net Benefit (TNs per 10,000 patients)") +
  theme_minimal(base_family = "Arial") + # Set font to Arial
  theme(
    text = element_text(size = 8), # Adjust text size
    axis.title = element_text(face = "bold"), # Make axis titles bold
    legend.position = "right", # Move legend to the right
    legend.title = element_text(face = "bold"), # Make legend title bold
    legend.text = element_text(size = 6), # Adjust legend text size
    panel.grid.major = element_line(colour = "gray90"), # Adjust grid line colour
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white"), # Set panel background colour
    legend.background = element_rect(fill = "white") # Set legend background to transparent
  ) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), # Set plot background to transparent
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Adjust plot margins
  ) +
  ggtitle("Validation")+
  theme(legend.position = c(0.2, 0.3)) +
  guides(
    shape = "none"
  )


# Save
ggsave("Output/Plots/SupFig9b.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)