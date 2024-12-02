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





#### Random_Split Outputs ####
which_output_path <- "_randomsplit"




#### Misc. Functions ####


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





#### NET BENEFIT PARAMETERS ####

# For diabetes, we need thresholds over which to integrate, and their weight
diabetes_threshold <- 0.15

# For lung cancer, we need threshold used and capacity
lung_cancer_threshold <- 0.015

#### First section of resutls: Diabetes ####

# Open training/validation data
df <- readRDS("Data/diabetes_df_randomsplit_imputed.rds")

# Turn validation and training data into right format
x_and_y <- prepareForGLMNET(df)
x <- x_and_y[["x"]]
y <- x_and_y[["y"]]
sa <- df$sa

# Open the diabetes performance
performance_dict <- readRDS("Output/performance_diabetes.rds")

#### Fig. 1a: AUROC table ####

# Open the diabetes performance
performance_dict <- readRDS("Output/performance_diabetes.rds")

metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Metric.Avg=double(),
                       Metric.Low=double(),
                       Metric.High=double(),
                       stringsAsFactors=FALSE)

# For each model, get the AUC (with 95% CI) of each subgroup
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  # For each subgroup, calculate subgroup AUC
  for (each_sa in c("Asian","Black","Other","White","Overall")) {
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,
                                          Model=model_name,
                                          Metric.Avg = performance_dict[[model_name]][[each_sa]][["AUROC"]][2],
                                          Metric.Low = performance_dict[[model_name]][[each_sa]][["AUROC"]][1],
                                          Metric.High = performance_dict[[model_name]][[each_sa]][["AUROC"]][3])
  }
  
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
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "LogNoSA" = 0,
                                "LogSingleSA" = 5,
                                "LogMultiSA" = 7),
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
  ) + 
  ggtitle("(a)")

# Save
saveRDS(metric_df,"Output/Plots/Data_Fig1a.rds")
ggsave("Output/Plots/Fig1a.png", plot = p, width = 8, height = 4, units = "in", dpi = 400)


#### Fig. 1b: Slope Plot ####

# Open the diabetes performance
performance_dict <- readRDS("Output/performance_diabetes.rds")

metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Metric.Avg=double(),
                       Metric.Low=double(),
                       Metric.High=double(),
                       stringsAsFactors=FALSE)

# For each model, get the AUC (with 95% CI) of each subgroup
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  # For each subgroup, calculate subgroup AUC
  for (each_sa in c("Asian","Black","Other","White","Overall")) {
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,
                                          Model=model_name,
                                          Metric.Avg = performance_dict[[model_name]][[each_sa]][["Slope"]][2],
                                          Metric.Low = performance_dict[[model_name]][[each_sa]][["Slope"]][1],
                                          Metric.High = performance_dict[[model_name]][[each_sa]][["Slope"]][3])
  }
  
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
colnames(metric_df) <- c("Ethnicity","Model", "Slope", "Slope.Low", "Slope.High")

# Plot Slope
pd <- position_dodge(0.6) # move them .05 to the left and right
p <- ggplot(metric_df, aes(x=Ethnicity, y=Slope, colour=Model,shape = Model)) + 
  geom_errorbar(aes(ymin=Slope.Low, ymax=Slope.High), width=.1, position=pd) +
  geom_point(position=pd) +
  xlab("Ethnicity") +
  ylab("Calibration Slope") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "LogNoSA" = 0,
                                "LogSingleSA" = 5,
                                "LogMultiSA" = 7),
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
  guides(color = "none",shape = "none")  + 
  ggtitle("(b)")

# Save
saveRDS(metric_df,"Output/Plots/Data_Fig1b.rds")
ggsave("Output/Plots/Fig1b.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)


#### Fig. 1c: Intercept Plot ####

# Open the diabetes performance
performance_dict <- readRDS("Output/performance_diabetes.rds")

metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Metric.Avg=double(),
                       Metric.Low=double(),
                       Metric.High=double(),
                       stringsAsFactors=FALSE)

# For each model, get the AUC (with 95% CI) of each subgroup
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  # For each subgroup, calculate subgroup AUC
  for (each_sa in c("Asian","Black","Other","White","Overall")) {
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,
                                          Model=model_name,
                                          Metric.Avg = performance_dict[[model_name]][[each_sa]][["Intercept"]][2],
                                          Metric.Low = performance_dict[[model_name]][[each_sa]][["Intercept"]][1],
                                          Metric.High = performance_dict[[model_name]][[each_sa]][["Intercept"]][3])
  }
  
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
colnames(metric_df) <- c("Ethnicity","Model", "Intercept", "Intercept.Low", "Intercept.High")

# Plot Intercept
pd <- position_dodge(0.6) # move them .05 to the left and right
p <- ggplot(metric_df, aes(x=Ethnicity, y=Intercept, colour=Model,shape = Model)) + 
  geom_errorbar(aes(ymin=Intercept.Low, ymax=Intercept.High), width=.1, position=pd) +
  geom_point(position=pd) +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "LogNoSA" = 0,
                                "LogSingleSA" = 5,
                                "LogMultiSA" = 7),
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
  guides(color = "none",shape = "none")  + 
  ggtitle("(c)")

# Save
saveRDS(metric_df,"Output/Plots/Data_Fig1c.rds")
ggsave("Output/Plots/Fig1c.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)




#### Fig. 2: sNB Comparison ####

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
    x <- sum(y)
    n <- length(y)
  } else {
    x <- sum(y[sa==each_sa])
    n <- length(y[sa==each_sa])
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
    x <- sum(y)
    n <- length(y)
  } else {
    x <- sum(y[sa==each_sa])
    n <- length(y[sa==each_sa])
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
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  print(model_name)
  
  # For each subgroup...
  for (each_sa in c("Asian","Black","Other","White","Overall")) {
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,
                                          Model=model_name,
                                          Metric.Avg = 10000*performance_dict[[model_name]][[each_sa]][["sNB"]][2],
                                          Metric.Low = 10000*performance_dict[[model_name]][[each_sa]][["sNB"]][1],
                                          Metric.High = 10000*performance_dict[[model_name]][[each_sa]][["sNB"]][3])
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
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3],
                                 "Random 5%" = clrs[4]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "LogNoSA" = 0,
                                "LogSingleSA" = 5,
                                "LogMultiSA" = 7,
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
saveRDS(metric_df,"Output/Plots/Data_Fig2.rds")
ggsave("Output/Plots/Fig2.png", plot = p, width = 8, height = 4, units = "in", dpi = 400)

# What's the difference between White and Other groups?
# For Treat No-One
mean_nw <- metric_df %>% filter((Ethnicity != "White")&(Ethnicity != "Overall")&(Model=="Treat.No.One")) %>%
  summarize(sNB = mean(sNB, na.rm = TRUE))
mean_w <- metric_df %>% filter((Ethnicity == "White")&(Model=="Treat.No.One")) %>%
  summarize(sNB = mean(sNB, na.rm = TRUE))
print(paste("For treating no-one, the White-Rest difference is",mean_w-mean_nw))
# For LogSingleSA
mean_nw <- metric_df %>% filter((Ethnicity != "White")&(Ethnicity != "Overall")&(Model=="LogSingleSA")) %>%
  summarize(sNB = mean(sNB, na.rm = TRUE))
mean_w <- metric_df %>% filter((Ethnicity == "White")&(Model=="LogSingleSA")) %>%
  summarize(sNB = mean(sNB, na.rm = TRUE))
print(paste("For LogSingleSA, the White-Rest difference is",mean_w-mean_nw))




#### Sup. Fig. 1: DCA ####

# Open data for Treatall/none
df <- readRDS("Data/diabetes_df_randomsplit_imputed.rds")
x_and_y <- prepareForGLMNET(df)
x <- x_and_y[["x"]]
y <- x_and_y[["y"]]
sa <- df$sa

# Net benefit metric
metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Clinical.Threshold=double(),
                       NB=double(),
                       stringsAsFactors=FALSE)

# For each model...
for (model_name in c("Treat.No.One","Treat.All","LogNoSA","LogSingleSA","LogMultiSA")) {
  
  # For each subgroup...
  for (each_sa in c("Asian","Black","Other","White","Overall")) {
    
    if (each_sa == "Overall") {
      sel <- rep(TRUE,length(sa))
    } else {
      sel <- sa==each_sa
    }
    
    # What model are we using
    if (model_name == "Treat.No.One") {
      nb <- calcNetBenefit(clinical_thresholds = performance_dict[["LogNoSA"]][[each_sa]][["DCA_x"]],
                           model_thresholds = performance_dict[["LogNoSA"]][[each_sa]][["DCA_x"]],
                           y_true = y[sel],
                           s_pred = 0*y[sel],
                           new_version = FALSE)
      metric_df <- rbind(metric_df,
                         data.frame(Subgroup=each_sa,
                                    Model=model_name, 
                                    Clinical.Threshold=performance_dict[["LogNoSA"]][[each_sa]][["DCA_x"]],
                                    NB=nb,
                                    NB.Low = nb,
                                    NB.High = nb,
                                    stringsAsFactors=FALSE))
    } else if (model_name == "Treat.All") {
      nb <- calcNetBenefit(clinical_thresholds = performance_dict[["LogNoSA"]][[each_sa]][["DCA_x"]],
                           model_thresholds = performance_dict[["LogNoSA"]][[each_sa]][["DCA_x"]],
                           y_true = y[sel],
                           s_pred = 0*y[sel]+1,
                           new_version = FALSE)
      metric_df <- rbind(metric_df,
                         data.frame(Subgroup=each_sa,
                                    Model=model_name, 
                                    Clinical.Threshold=performance_dict[["LogNoSA"]][[each_sa]][["DCA_x"]],
                                    NB=nb,
                                    NB.Low = nb,
                                    NB.High = nb,
                                    stringsAsFactors=FALSE))
    } else {
      # Append to dataframe
      metric_df <- rbind(metric_df,
                         data.frame(Subgroup=each_sa,
                                    Model=model_name, 
                                    Clinical.Threshold=performance_dict[[model_name]][[each_sa]][["DCA_x"]],
                                    NB=performance_dict[[model_name]][[each_sa]][["DCA"]][,2],
                                    NB.Low = performance_dict[[model_name]][[each_sa]][["DCA"]][,1],
                                    NB.High = performance_dict[[model_name]][[each_sa]][["DCA"]][,3],
                                    stringsAsFactors=FALSE))
    }
  }
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
# Remove  0 and 1 Clinical threshold
metric_df <- metric_df %>%
  filter(Clinical.Threshold != 0 & Clinical.Threshold != 1)
# Correct for /10,000
metric_df$NB <- 10000*metric_df$NB
metric_df$NB.Low <- 10000*metric_df$NB.Low
metric_df$NB.High <- 10000*metric_df$NB.High
# Save
saveRDS(metric_df,"Output/Plots/Data_SupFig1.rds")

# Figure 2a
df <- metric_df %>% filter(Subgroup == "Overall") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +  
  geom_ribbon(aes(ymin = NB.Low, ymax = NB.High, fill = Model), alpha = 0.1, colour = NA) +  
  scale_fill_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "LogNoSA" = "solid",
                                   "LogSingleSA" = "solid",
                                   "LogMultiSA" = "solid"),
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
  ggtitle("(a) Overall")

# Save
ggsave("Output/Plots/SupFig1a.png", plot = p, width = 8, height = 4, units = "in", dpi = 400)



# Figure 2b
df <- metric_df %>% filter(Subgroup == "Asian") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = NB.Low, ymax = NB.High, fill = Model), alpha = 0.1, colour = NA) +  
  scale_fill_manual(values = c("Treat.No.One" = "#333a52",
                               "Treat.All" = "#333a32",
                               "LogNoSA" = clrs[1],
                               "LogSingleSA" = clrs[2],
                               "LogMultiSA" = clrs[3]),
                    name = "Model") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "LogNoSA" = "solid",
                                   "LogSingleSA" = "solid",
                                   "LogMultiSA" = "solid"),
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
  ggtitle("(b) Asian")

# Save
ggsave("Output/Plots/SupFig1b.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)


# Figure 2c
df <- metric_df %>% filter(Subgroup == "Black") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = NB.Low, ymax = NB.High, fill = Model), alpha = 0.1, colour = NA) +  
  scale_fill_manual(values = c("Treat.No.One" = "#333a52",
                               "Treat.All" = "#333a32",
                               "LogNoSA" = clrs[1],
                               "LogSingleSA" = clrs[2],
                               "LogMultiSA" = clrs[3]),
                    name = "Model") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "LogNoSA" = "solid",
                                   "LogSingleSA" = "solid",
                                   "LogMultiSA" = "solid"),
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
  ggtitle("(c) Black")

# Save
ggsave("Output/Plots/SupFig1c.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)


# Figure 2d
df <- metric_df %>% filter(Subgroup == "Other") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = NB.Low, ymax = NB.High, fill = Model), alpha = 0.1, colour = NA) +  
  scale_fill_manual(values = c("Treat.No.One" = "#333a52",
                               "Treat.All" = "#333a32",
                               "LogNoSA" = clrs[1],
                               "LogSingleSA" = clrs[2],
                               "LogMultiSA" = clrs[3]),
                    name = "Model") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "LogNoSA" = "solid",
                                   "LogSingleSA" = "solid",
                                   "LogMultiSA" = "solid"),
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
  ggtitle("(d) Other")

# Save
ggsave("Output/Plots/SupFig1d.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)



# Figure 2e
df <- metric_df %>% filter(Subgroup == "White") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = NB.Low, ymax = NB.High, fill = Model), alpha = 0.1, colour = NA) +  
  scale_fill_manual(values = c("Treat.No.One" = "#333a52",
                               "Treat.All" = "#333a32",
                               "LogNoSA" = clrs[1],
                               "LogSingleSA" = clrs[2],
                               "LogMultiSA" = clrs[3]),
                    name = "Model") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "LogNoSA" = "solid",
                                   "LogSingleSA" = "solid",
                                   "LogMultiSA" = "solid"),
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
  ggtitle("(e) White")

# Save
ggsave("Output/Plots/SupFig1e.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)



### Second section: Lung cancer ####


# Open training/validation data
df <- readRDS("Data/lungcancer_df_randomsplit_imputed.rds")

# Turn validation and training data into right format
x_and_y <- prepareForGLMNET(df)
x <- x_and_y[["x"]]
y <- x_and_y[["y"]]
sa <- df$sa


# Open the lungcancer performance
performance_dict <- readRDS("Output/performance_lung_cancer.rds")


#### Fig. 3a: AUROC table ####

# Open the lungcancer performance
performance_dict <- readRDS("Output/performance_lung_cancer.rds")

metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Metric.Avg=double(),
                       Metric.Low=double(),
                       Metric.High=double(),
                       stringsAsFactors=FALSE)

# For each model, get the AUC (with 95% CI) of each subgroup
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  # For each subgroup, calculate subgroup AUC
  for (each_sa in c("0","1","2","3","4","Overall")) {
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,
                                          Model=model_name,
                                          Metric.Avg = performance_dict[[model_name]][[each_sa]][["AUROC"]][2],
                                          Metric.Low = performance_dict[[model_name]][[each_sa]][["AUROC"]][1],
                                          Metric.High = performance_dict[[model_name]][[each_sa]][["AUROC"]][3])
  }
  
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
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "LogNoSA" = 0,
                                "LogSingleSA" = 5,
                                "LogMultiSA" = 7),
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
  )  + 
  ggtitle("(a)")

# Save
saveRDS(metric_df,"Output/Plots/Data_Fig3a.rds")
ggsave("Output/Plots/Fig3a.png", plot = p, width = 8, height = 4, units = "in", dpi = 400)

#### Fig. 3b: Slope Plot ####

# Open the lungcancer performance
performance_dict <- readRDS("Output/performance_lung_cancer.rds")
metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Metric.Avg=double(),
                       Metric.Low=double(),
                       Metric.High=double(),
                       stringsAsFactors=FALSE)

# For each model, get the AUC (with 95% CI) of each subgroup
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  # For each subgroup, calculate subgroup AUC
  for (each_sa in c("0","1","2","3","4","Overall")) {
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,
                                          Model=model_name,
                                          Metric.Avg = performance_dict[[model_name]][[each_sa]][["Slope"]][2],
                                          Metric.Low = performance_dict[[model_name]][[each_sa]][["Slope"]][1],
                                          Metric.High = performance_dict[[model_name]][[each_sa]][["Slope"]][3])
  }
  
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
levels(metric_df$Subgroup) <- c("Q1","Q2","Q3","Q4","Q5","Overall")
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
colnames(metric_df) <- c("Townsend.Quintile","Model", "Slope", "Slope.Low", "Slope.High")

# Plot Slope
pd <- position_dodge(0.6) # move them .05 to the left and right
p <- ggplot(metric_df, aes(x=Townsend.Quintile, y=Slope, colour=Model,shape = Model)) + 
  geom_errorbar(aes(ymin=Slope.Low, ymax=Slope.High), width=.1, position=pd) +
  geom_point(position=pd) +
  xlab("Townsend.Quintile") +
  ylab("Calibration Slope") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "LogNoSA" = 0,
                                "LogSingleSA" = 5,
                                "LogMultiSA" = 7),
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
  guides(color = "none",shape = "none")  + 
  ggtitle("(b)")

# Save
saveRDS(metric_df,"Output/Plots/Data_Fig3b.rds")
ggsave("Output/Plots/Fig3b.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)


#### Fig. 3c: Intercept Plot ####

# Open the lungcancer performance
performance_dict <- readRDS("Output/performance_lung_cancer.rds")
metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Metric.Avg=double(),
                       Metric.Low=double(),
                       Metric.High=double(),
                       stringsAsFactors=FALSE)

# For each model, get the AUC (with 95% CI) of each subgroup
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  # For each subgroup, calculate subgroup AUC
  for (each_sa in c("0","1","2","3","4","Overall")) {
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,
                                          Model=model_name,
                                          Metric.Avg = performance_dict[[model_name]][[each_sa]][["Intercept"]][2],
                                          Metric.Low = performance_dict[[model_name]][[each_sa]][["Intercept"]][1],
                                          Metric.High = performance_dict[[model_name]][[each_sa]][["Intercept"]][3])
  }
  
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
levels(metric_df$Subgroup) <- c("Q1","Q2","Q3","Q4","Q5","Overall")
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
colnames(metric_df) <- c("Townsend.Quintile","Model", "Intercept", "Intercept.Low", "Intercept.High")

# Plot Intercept
pd <- position_dodge(0.6) # move them .05 to the left and right
p <- ggplot(metric_df, aes(x=Townsend.Quintile, y=Intercept, colour=Model,shape = Model)) + 
  geom_errorbar(aes(ymin=Intercept.Low, ymax=Intercept.High), width=.1, position=pd) +
  geom_point(position=pd) +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "LogNoSA" = 0,
                                "LogSingleSA" = 5,
                                "LogMultiSA" = 7),
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
  guides(color = "none",shape = "none")  + 
  ggtitle("(c)")

# Save
saveRDS(metric_df,"Output/Plots/Data_Fig3c.rds")
ggsave("Output/Plots/Fig3c.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)



#### Fig. 4: sNB Comparison ####

performance_dict <- readRDS("Output/performance_lung_cancer.rds")

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
    x <- sum(y)
    n <- length(y)
  } else {
    x <- sum(y[sa==each_sa])
    n <- length(y[sa==each_sa])
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
    x <- sum(y)
    n <- length(y)
  } else {
    x <- sum(y[sa==each_sa])
    n <- length(y[sa==each_sa])
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
for (model_name in c("LogNoSA","LogSingleSA","LogMultiSA")) {
  
  print(model_name)
  
  # For each subgroup...
  for (each_sa in c("0","1","2","3","4","Overall")) {
    
    # Attach to df
    metric_df[nrow(metric_df)+1,] <- list(Subgroup=each_sa,
                                          Model=model_name,
                                          Metric.Avg = 10000*performance_dict[[model_name]][[each_sa]][["sNB"]][2],
                                          Metric.Low = 10000*performance_dict[[model_name]][[each_sa]][["sNB"]][1],
                                          Metric.High = 10000*performance_dict[[model_name]][[each_sa]][["sNB"]][3])
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
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3],
                                 "Random 5%" = clrs[4]),
                      name = "Model") +
  scale_shape_manual(values = c("Treat.No.One" = 3,
                                "LogNoSA" = 0,
                                "LogSingleSA" = 5,
                                "LogMultiSA" = 7,
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
saveRDS(metric_df,"Output/Plots/Data_Fig4.rds")
metric_df <- readRDS("Output/Plots/Data_Fig4.rds")
ggsave("Output/Plots/Fig4.png", plot = p, width = 8, height = 4, units = "in", dpi = 400)



#### Sup. Fig. 5: DCA ####

# Open data for Treatall/none
df <- readRDS("Data/lungcancer_df_randomsplit_imputed.rds")
x_and_y <- prepareForGLMNET(df)
x <- x_and_y[["x"]]
y <- x_and_y[["y"]]
sa <- df$sa

# Net benefit metric
metric_df <-data.frame(Subgroup=character(),
                       Model=character(), 
                       Clinical.Threshold=double(),
                       NB=double(),
                       stringsAsFactors=FALSE)

# For each model...
for (model_name in c("Treat.No.One","Treat.All","LogNoSA","LogSingleSA","LogMultiSA")) {
  
  # For each subgroup...
  for (each_sa in c("0","1","2","3","4","Overall")) {
    
    if (each_sa == "Overall") {
      sel <- rep(TRUE,length(sa))
    } else {
      sel <- sa==each_sa
    }
    
    # What model are we using
    if (model_name == "Treat.No.One") {
      nb <- calcNetBenefit(clinical_thresholds = performance_dict[["LogNoSA"]][[each_sa]][["DCA_x"]],
                           model_thresholds = performance_dict[["LogNoSA"]][[each_sa]][["DCA_x"]],
                           y_true = y[sel],
                           s_pred = 0*y[sel],
                           new_version = FALSE)
      metric_df <- rbind(metric_df,
                         data.frame(Subgroup=each_sa,
                                    Model=model_name, 
                                    Clinical.Threshold=performance_dict[["LogNoSA"]][[each_sa]][["DCA_x"]],
                                    NB=nb,
                                    NB.Low = nb,
                                    NB.High = nb,
                                    stringsAsFactors=FALSE))
    } else if (model_name == "Treat.All") {
      nb <- calcNetBenefit(clinical_thresholds = performance_dict[["LogNoSA"]][[each_sa]][["DCA_x"]],
                           model_thresholds = performance_dict[["LogNoSA"]][[each_sa]][["DCA_x"]],
                           y_true = y[sel],
                           s_pred = 0*y[sel]+1,
                           new_version = FALSE)
      metric_df <- rbind(metric_df,
                         data.frame(Subgroup=each_sa,
                                    Model=model_name, 
                                    Clinical.Threshold=performance_dict[["LogNoSA"]][[each_sa]][["DCA_x"]],
                                    NB=nb,
                                    NB.Low = nb,
                                    NB.High = nb,
                                    stringsAsFactors=FALSE))
    } else {
      # Append to dataframe
      metric_df <- rbind(metric_df,
                         data.frame(Subgroup=each_sa,
                                    Model=model_name, 
                                    Clinical.Threshold=performance_dict[[model_name]][[each_sa]][["DCA_x"]],
                                    NB=performance_dict[[model_name]][[each_sa]][["DCA"]][,2],
                                    NB.Low = performance_dict[[model_name]][[each_sa]][["DCA"]][,1],
                                    NB.High = performance_dict[[model_name]][[each_sa]][["DCA"]][,3],
                                    stringsAsFactors=FALSE))
    }
  }
}

# Rename and factorise
metric_df$Subgroup <- factor(metric_df$Subgroup, levels = unique(metric_df$Subgroup))
levels(metric_df$Subgroup) <- c("Q1","Q2","Q3","Q4","Q5","Overall")
metric_df$Model <- factor(metric_df$Model, levels = unique(metric_df$Model))
# Remove  0 and 1 Clinical threshold
metric_df <- metric_df %>%
  filter(Clinical.Threshold != 0 & Clinical.Threshold != 1)
# Correct for /10,000
metric_df$NB <- 10000*metric_df$NB
metric_df$NB.Low <- 10000*metric_df$NB.Low
metric_df$NB.High <- 10000*metric_df$NB.High
# Save
saveRDS(metric_df,"Output/Plots/Data_SupFig1.rds")

# Figure 2a
df <- metric_df %>% filter(Subgroup == "Overall") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +  
  geom_ribbon(aes(ymin = NB.Low, ymax = NB.High, fill = Model), alpha = 0.1, colour = NA) +  
  scale_fill_manual(values = c("Treat.No.One" = "#333a52",
                               "Treat.All" = "#333a32",
                               "LogNoSA" = clrs[1],
                               "LogSingleSA" = clrs[2],
                               "LogMultiSA" = clrs[3]),
                    name = "Model") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "LogNoSA" = "solid",
                                   "LogSingleSA" = "solid",
                                   "LogMultiSA" = "solid"),
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
  geom_vline(xintercept = lung_cancer_threshold, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("(a) Overall")

# Save
ggsave("Output/Plots/SupFig5a.png", plot = p, width = 8, height = 4, units = "in", dpi = 400)



# Figure 2b
df <- metric_df %>% filter(Subgroup == "Q1") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = NB.Low, ymax = NB.High, fill = Model), alpha = 0.1, colour = NA) +  
  scale_fill_manual(values = c("Treat.No.One" = "#333a52",
                               "Treat.All" = "#333a32",
                               "LogNoSA" = clrs[1],
                               "LogSingleSA" = clrs[2],
                               "LogMultiSA" = clrs[3]),
                    name = "Model") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "LogNoSA" = "solid",
                                   "LogSingleSA" = "solid",
                                   "LogMultiSA" = "solid"),
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
  geom_vline(xintercept = lung_cancer_threshold, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("(b) Q1")

# Save
ggsave("Output/Plots/SupFig5b.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)


# Figure 2c
df <- metric_df %>% filter(Subgroup == "Q2") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = NB.Low, ymax = NB.High, fill = Model), alpha = 0.1, colour = NA) +  
  scale_fill_manual(values = c("Treat.No.One" = "#333a52",
                               "Treat.All" = "#333a32",
                               "LogNoSA" = clrs[1],
                               "LogSingleSA" = clrs[2],
                               "LogMultiSA" = clrs[3]),
                    name = "Model") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "LogNoSA" = "solid",
                                   "LogSingleSA" = "solid",
                                   "LogMultiSA" = "solid"),
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
  geom_vline(xintercept = lung_cancer_threshold, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("(c) Q2")

# Save
ggsave("Output/Plots/SupFig5c.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)


# Figure 2d
df <- metric_df %>% filter(Subgroup == "Q3") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = NB.Low, ymax = NB.High, fill = Model), alpha = 0.1, colour = NA) +  
  scale_fill_manual(values = c("Treat.No.One" = "#333a52",
                               "Treat.All" = "#333a32",
                               "LogNoSA" = clrs[1],
                               "LogSingleSA" = clrs[2],
                               "LogMultiSA" = clrs[3]),
                    name = "Model") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "LogNoSA" = "solid",
                                   "LogSingleSA" = "solid",
                                   "LogMultiSA" = "solid"),
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
  geom_vline(xintercept = lung_cancer_threshold, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("(d) Q3")

# Save
ggsave("Output/Plots/SupFig5d.png", plot = p, width = 4*2/3, height = 4, units = "in", dpi = 400)



# Figure 2e
df <- metric_df %>% filter(Subgroup == "Q4") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = NB.Low, ymax = NB.High, fill = Model), alpha = 0.1, colour = NA) +  
  scale_fill_manual(values = c("Treat.No.One" = "#333a52",
                               "Treat.All" = "#333a32",
                               "LogNoSA" = clrs[1],
                               "LogSingleSA" = clrs[2],
                               "LogMultiSA" = clrs[3]),
                    name = "Model") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "LogNoSA" = "solid",
                                   "LogSingleSA" = "solid",
                                   "LogMultiSA" = "solid"),
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
  geom_vline(xintercept = lung_cancer_threshold, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("(e) Q4")

# Save
ggsave("Output/Plots/SupFig5e.png", plot = p, width = 4*2/3, height = 4, units = "in", dpi = 400)



# Figure 2e
df <- metric_df %>% filter(Subgroup == "Q5") %>%
  select(-Subgroup)

# Plot
p <- ggplot(df, aes(x = Clinical.Threshold, y = NB, colour = Model, linetype = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = NB.Low, ymax = NB.High, fill = Model), alpha = 0.1, colour = NA) +  
  scale_fill_manual(values = c("Treat.No.One" = "#333a52",
                               "Treat.All" = "#333a32",
                               "LogNoSA" = clrs[1],
                               "LogSingleSA" = clrs[2],
                               "LogMultiSA" = clrs[3]),
                    name = "Model") +
  scale_colour_manual(values = c("Treat.No.One" = "#333a52",
                                 "Treat.All" = "#333a32",
                                 "LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Model") +
  scale_linetype_manual(values = c("Treat.No.One" = "twodash",
                                   "Treat.All" = "dashed",
                                   "LogNoSA" = "solid",
                                   "LogSingleSA" = "solid",
                                   "LogMultiSA" = "solid"),
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
  geom_vline(xintercept = lung_cancer_threshold, linetype = "dotted", colour = "#faaa1b", alpha = 0.8) +
  ggtitle("(f) Q5")

# Save
ggsave("Output/Plots/SupFig5f.png", plot = p, width = 4*2/3, height = 4, units = "in", dpi = 400)


#### Fig 5: sNB Resource constraints ####

# Currently this only 'works' on real data.
# As it is specific to the incidence and predictive nature of the
# real lung cancer example.

# Load Data on Pareto Optimality
pareto_dict <- readRDS("Output/pareto_with_ci_dict.rds")

# To plot Figure 5, we take the values and add it to a df file.
fig5_df <- list()
for (constraint in c("NoCon","3Perc","1Perc")) {
  for (model in c("LogNoSA","LogSingleSA","LogMultiSA")) {
    df_part <- pareto_dict[[model]][[constraint]]
    df_part$Capacity <- constraint
    df_part$Model <- model
    fig5_df[[length(fig5_df)+1]] <- df_part
  }
}
fig5_df <- Reduce(full_join,fig5_df)
fig5_df$Model <- factor(fig5_df$Model, levels = c("LogNoSA","LogSingleSA","LogMultiSA"))
fig5_df <-  fig5_df %>%
  mutate(column_name = recode(Capacity, 
                              "NoCon" = "No Constraints", 
                              "3Perc" = "3%", 
                              "1Perc" = "1%"))
fig5_df$TopDeprivOCNB <- 10000*fig5_df$TopDeprivOCNB
fig5_df$TopDeprivOCNB.Low <- 10000*fig5_df$TopDeprivOCNB.Low
fig5_df$TopDeprivOCNB.High <- 10000*fig5_df$TopDeprivOCNB.High
fig5_df$RestOCNB <- 10000*fig5_df$RestOCNB
fig5_df$RestOCNB.Low <- 10000*fig5_df$RestOCNB.Low
fig5_df$RestOCNB.High <- 10000*fig5_df$RestOCNB.High

# Calculate Treat No One
df <- readRDS("Data/lungcancer_df_randomsplit_imputed.rds")
overall_treat_no_oone <- (1-mean(df$y))*10000
minimum_treat_no_one <- (1-mean(df$y[df$sa=="4"]))*10000


# Plot
p <- ggplot(fig5_df, aes(x = TopDeprivOCNB, y = RestOCNB, colour = Model, linetype = Model, shape = Capacity)) +
  geom_line() +
  geom_point() +
#  geom_errorbar(aes(ymin = RestOCNB.Low, ymax = RestOCNB.High), width = 0.2, linetype = 1, alpha = 0.1) +
#  geom_errorbarh(aes(xmin = TopDeprivOCNB.Low, xmax = TopDeprivOCNB.High), height = 0.2, linetype = 1, alpha = 0.1) +
  scale_colour_manual(values = c("LogNoSA" = clrs[1],
                                 "LogSingleSA" = clrs[2],
                                 "LogMultiSA" = clrs[3]),
                      name = "Max. Capacity") +
  scale_linetype_manual(values = c("LogNoSA" = 2,
                                   "LogSingleSA" = 3,
                                   "LogMultiSA" = 4),
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
  ggtitle("") +
  coord_cartesian(xlim = c(9932.5, 9938),
                  ylim = c(9961, 9963.5)) +
#  scale_x_continuous(breaks = seq(9932.5, 9938, 1)) +
#  scale_y_continuous(breaks = seq(9961, 9963.5, 0.5)) +
  theme(legend.position="none")



# Save
ggsave("Output/Plots/Fig5.png", plot = p, width = 4, height = 4, units = "in", dpi = 400)

