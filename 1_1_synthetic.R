# Set working directory
setwd("/mnt/bmh01-rds/Sperrin_UKBB_Fairness/Code/study_1_benefit/Manuscript_Ready_Code/FairSubgroupBenefit/")

# Import libraries
library(dplyr)

# Set random seed
set.seed(42)

# For saving data
which_output_path = "_randomsplit"

#### Create synthetic diabetes data ####

# We choose 10000 rows for reproducibility
n <- 10000

# We populate a dataframe with reasonable looking random numbers.
df <- data.frame(
  X = 1:n,
  Censored = sample(2018:2024, n, replace = TRUE),
  Started = sample(2000:2020, n, replace = TRUE),
  Death = sample(2018:2024, n, replace = TRUE),
  Centre = sample(c("Reading", "Croydon", "Oxford", "Liverpool"), n, replace = TRUE),
  Age = sample(20:90, n, replace = TRUE),
  Female = sample(c(0, 1), n, replace = TRUE),
  Urban = sample(c(0, 1), n, replace = TRUE),
  Townsend = runif(n, -6, 6),
  Ethnicity = sample(c("White", "Black", "Asian", "Other"), n, replace = TRUE),
  FamilyHxHipFracture = sample(c(0, 1), n, replace = TRUE),
  FamilyHxDepression = sample(c(0, 1), n, replace = TRUE),
  FamilyHxDiabetes = sample(c(0, 1), n, replace = TRUE),
  FamilyHxHBP = sample(c(0, 1), n, replace = TRUE),
  FamilyHxCVD = sample(c(0, 1), n, replace = TRUE),
  DBP = runif(n, 60, 100),
  SBP = runif(n, 100, 180),
  BMI = runif(n, 18, 40),
  BodyFat = runif(n, 10, 50),
  WaistCircumference = runif(n, 60, 120),
  Education = sample(c("Level0", "Level1", "Level2", "Level3"), n, replace = TRUE),
  Smoker = sample(c("Previous", "Never", "Current"), n, replace = TRUE),
  PackYears = runif(n, 0, 60),
  AlcoholStatus = sample(c(0, 1), n, replace = TRUE),
  AlcoholHigh = sample(c(0, 1), n, replace = TRUE),
  HRT = sample(c(0, 1), n, replace = TRUE),
  Hysterectomy = sample(c(0, 1), n, replace = TRUE),
  Menopause = sample(c(0, 1), n, replace = TRUE),
  MedicationCholesterol = sample(c(0, 1), n, replace = TRUE),
  MedicationBloodPressure = sample(c(0, 1), n, replace = TRUE),
  MedicationHRT = sample(c(0, 1), n, replace = TRUE),
  MedicationContraceptive = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceM07 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceM08 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceL40 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceO24 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceO60 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceM15to19 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceE28 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceF20to29 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceM30to36 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceF32to39 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceF30to31 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceG47 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceE66 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceI25to73 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceI10to11 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceE78 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceE88 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceK05 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceF80to89 = sample(c(0, 1), n, replace = TRUE),
  y = sample(c(0, 1), n, replace = TRUE),
  is_val = rbinom(n=n, size=1, prob=0.15)
)

# We want to introduce a covariate shift for different Ethnicity
# and Townsend score in age, to make the propensity score catch
# some heterogeneity between sensitive attributes
df$Age <- ifelse(df$Ethnicity == "White", as.integer(1.1 * df$Age), df$Age)
df$Age <- ifelse(df$Ethnicity == "Asian", as.integer(0.7 * df$Age), df$Age)
df$Age <- ifelse(df$Ethnicity == "White", as.integer(1.3 * df$Age), df$Age)
df$Age <- ifelse(df$Townsend > 3, as.integer(1.5 * df$Age), df$Age)
df$Age <- ifelse(df$Townsend < -3, as.integer(0.7 * df$Age), df$Age)


# Turn strings into factors
df <- df %>%
  mutate_if(is.character, as.factor)

# Create a linear predictor of risk
# Add some non-linear/interaction effects
linear_predictor <- 0.1 * (df$Age - 50) +
  0.01 * (df$Age - 50)^2 +
  2 * df$Female +
  3 * df$FamilyHxDiabetes +
  1 * (df$Smoker == "Never") -
  2 * df$FirstIncidenceI10to11 +
  0.2 * df$FamilyHxDiabetes*df$Female +
  1.5 * (df$Ethnicity == "Asian") +
  1 * (df$Ethnicity == "Black") -
  10
probabilities <- 1 / (1 + exp(-linear_predictor))

# Use these probabilities to create an outcome
df$y <- rbinom(n=n, size=1, prob=probabilities)

# Add some reasonable looking missing data
df$Ethnicity[sample(n, size = as.integer(0.05 * n))] <- NA
df$Townsend[sample(n, size = as.integer(0.03 * n))] <- NA
df$Smoker[sample(n, size = as.integer(0.07 * n))] <- NA
df$FirstIncidenceI10to11[sample(n, size = as.integer(0.05 * n))] <- NA
# Add some extra missingness in non-white patients
df$FirstIncidenceI10to11[df$Ethnicity != "White" & runif(n) < 0.1] <- NA

# Save as data
write.csv(df,paste0("Data/diabetes_df",which_output_path,".csv"))


#### Create synthetic lung cancer data ####

# This is going to be almost the same as the diabetes synthetic example
# with same covariates for simplicity, even though the real data does not
# use the same columns exactly.

# We choose 10000 rows for reproducibility
n <- 10000

# We populate a dataframe with reasonable looking random numbers.
df <- data.frame(
  X = 1:n,
  Censored = sample(2018:2024, n, replace = TRUE),
  Started = sample(2000:2020, n, replace = TRUE),
  Death = sample(2018:2024, n, replace = TRUE),
  Centre = sample(c("Reading", "Croydon", "Oxford", "Liverpool"), n, replace = TRUE),
  Age = sample(20:90, n, replace = TRUE),
  Female = sample(c(0, 1), n, replace = TRUE),
  Urban = sample(c(0, 1), n, replace = TRUE),
  Townsend = runif(n, -6, 6),
  Ethnicity = sample(c("White", "Black", "Asian", "Other"), n, replace = TRUE),
  FamilyHxHipFracture = sample(c(0, 1), n, replace = TRUE),
  FamilyHxDepression = sample(c(0, 1), n, replace = TRUE),
  FamilyHxDiabetes = sample(c(0, 1), n, replace = TRUE),
  FamilyHxHBP = sample(c(0, 1), n, replace = TRUE),
  FamilyHxCVD = sample(c(0, 1), n, replace = TRUE),
  DBP = runif(n, 60, 100),
  SBP = runif(n, 100, 180),
  BMI = runif(n, 18, 40),
  BodyFat = runif(n, 10, 50),
  WaistCircumference = runif(n, 60, 120),
  Education = sample(c("Level0", "Level1", "Level2", "Level3"), n, replace = TRUE),
  Smoker = sample(c("Previous", "Never", "Current"), n, replace = TRUE),
  PackYears = runif(n, 0, 60),
  AlcoholStatus = sample(c(0, 1), n, replace = TRUE),
  AlcoholHigh = sample(c(0, 1), n, replace = TRUE),
  HRT = sample(c(0, 1), n, replace = TRUE),
  Hysterectomy = sample(c(0, 1), n, replace = TRUE),
  Menopause = sample(c(0, 1), n, replace = TRUE),
  MedicationCholesterol = sample(c(0, 1), n, replace = TRUE),
  MedicationBloodPressure = sample(c(0, 1), n, replace = TRUE),
  MedicationHRT = sample(c(0, 1), n, replace = TRUE),
  MedicationContraceptive = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceM07 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceM08 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceL40 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceO24 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceO60 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceM15to19 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceE28 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceF20to29 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceM30to36 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceF32to39 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceF30to31 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceG47 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceE66 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceI25to73 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceI10to11 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceE78 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceE88 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceK05 = sample(c(0, 1), n, replace = TRUE),
  FirstIncidenceF80to89 = sample(c(0, 1), n, replace = TRUE),
  y = sample(c(0, 1), n, replace = TRUE),
  is_val = rbinom(n=n, size=1, prob=0.15)
)

# We want to introduce a covariate shift for different Ethnicity
# and Townsend score in age, to make the propensity score catch
# some heterogeneity between sensitive attributes
df$Age <- ifelse(df$Ethnicity == "White", as.integer(1.1 * df$Age), df$Age)
df$Age <- ifelse(df$Ethnicity == "Asian", as.integer(0.7 * df$Age), df$Age)
df$Age <- ifelse(df$Ethnicity == "White", as.integer(1.3 * df$Age), df$Age)
df$Age <- ifelse(df$Townsend > 3, as.integer(1.5 * df$Age), df$Age)
df$Age <- ifelse(df$Townsend < -3, as.integer(0.7 * df$Age), df$Age)

# Turn strings into factors
df <- df %>%
  mutate_if(is.character, as.factor)

# Create a linear predictor of risk
# Add some non-linear/interaction effects
linear_predictor <- 0.1 * (df$Age - 50) +
  0.01 * (df$Age - 50)^2 -
  1 * df$Female -
  5 * (df$Smoker == "Never") -
  2 * df$FirstIncidenceI10to11 +
  1.5 * (df$Smoker == "Current")*(df$Female==0) +
  1.5 * (df$Ethnicity == "Asian") +
  1 * (df$Ethnicity == "Black") -
  10
probabilities <- 1 / (1 + exp(-linear_predictor))

# Use these probabilities to create an outcome
df$y <- rbinom(n=n, size=1, prob=probabilities)

# Add some reasonable looking missing data
df$Ethnicity[sample(n, size = as.integer(0.05 * n))] <- NA
df$Townsend[sample(n, size = as.integer(0.03 * n))] <- NA
df$Smoker[sample(n, size = as.integer(0.07 * n))] <- NA
df$FirstIncidenceI10to11[sample(n, size = as.integer(0.05 * n))] <- NA
# Add some extra missingness in non-white patients
df$FirstIncidenceI10to11[df$Ethnicity != "White" & runif(n) < 0.1] <- NA

# Save as data
write.csv(df,paste0("Data/lungcancer_df",which_output_path,".csv"))


