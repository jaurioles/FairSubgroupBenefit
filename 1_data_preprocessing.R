# Define Source Directory
# Change if needed/trying to reproduce
setwd("/mnt/bmh01-rds/Sperrin_UKBB_Fairness/Code/study_1_benefit/Manuscript_Ready_Code/FairSubgroupBenefit/")

# The purpose of this file is to pre-process the UK Biobank raw data into
# a dataset to fit the clinical prediction models. 

#------------------------#
#------------------------#
#### Import Libraries ####
#------------------------#
#------------------------#

library(tidyverse)
library(dplyr)

# This function is provided by UKBiobank to pre-process data
source("ukb676287.r")


#------------------#
#------------------#
#### Parameters ####
#------------------#
#------------------#

# Get output path
which_output_path <- "_randomsplit"

#-----------------#
#-----------------#
#### Functions ####
#-----------------#
#-----------------#

get_first_dates <- function(df) {
  
  #For the ICD-10 first incidence data, we map each code of the UKB columns to
  #their corresponding ICD-10 code. For each of these codes, we return the first
  #date this happened for a patient.
  
  ocurrences_guide <- read.csv(file="Data/all_first_ocurrences.csv", header = FALSE)
  ocurrences_guide <- ocurrences_guide[grepl("Date", ocurrences_guide$V2),]
  get_ICD <- function(x) {
    code <- strsplit(x," ")[[1]][2]
    return(code)
  }
  ocurrences_guide$V3 <- sapply(ocurrences_guide$V2,get_ICD)
  ocurrences_guide$V4 <- sapply(ocurrences_guide$V1,function(x) {
    paste("f.",x,".0.0",sep="")
  })
  ocurrences_guide$V5 <- sapply(ocurrences_guide$V3,function(x) {
    paste("FirstIncidence",x,sep="")
  })
  
  all_dates <- df[,c()]
  for (ii in 1:dim(ocurrences_guide)[1]) {
    if (ocurrences_guide$V4[ii] %in% colnames(df)) {
      all_dates[,ocurrences_guide$V5[ii]] <- df[,ocurrences_guide$V4[ii]]
    }
  }
  
  return(all_dates)
  
}

#-----------------------------------#
#-----------------------------------#
#### Open the UK Biobank Dataset ####
#-----------------------------------#
#-----------------------------------#

# Open real data
df <- read.table("Data/ukb676287.tab", header = TRUE, sep="\t")

#---------------------------#
#---------------------------#
#### 1) DIABETES DATASET ####
#---------------------------#
#---------------------------#

# Define matrix we will base data on
output <- df[c()]

# Follow-Up and Date at which assesment was made #

output$Censored <- df$f.191.0.0
output$Started  <- df$f.53.0.0
output$Death    <- coalesce(df$f.40000.0.0)

# Assesment centre #

assesment_centres <- list(
"10003"="Stockport (pilot)",
"11001"="Manchester",
"11002"="Oxford",
"11003"="Cardiff",
"11004"="Glasgow",
"11005"="Edinburgh",
"11006"="Stoke",
"11007"="Reading",
"11008"="Bury",
"11009"="Newcastle",
"11010"="Leeds",
"11011"="Bristol",
"11012"="Barts",
"11013"="Nottingham",
"11014"="Sheffield",
"11016"="Liverpool",
"11017"="Middlesborough",
"11018"="Hounslow",
"11020"="Croydon",
"11021"="Birmingham",
"11022"="Swansea",
"11023"="Wrexham",
"11024"="Cheadle (revisit)",
"11025"="Cheadle (imaging)",
"11026"="Reading (imaging)",
"11027"="Newcastle (imaging)",
"11028"="Bristol (imaging)"
)
assesment_column <- df["f.54.0.0"]
assesment_column <- assesment_centres[sapply(df["f.54.0.0"],as.character)]
output$Centre <- unlist(unname(assesment_column))
output$Centre <- as.factor(output$Centre)

# Demographic Characteristics #

output$Age <- df$f.21022.0
output$Female <- as.numeric(df$f.31.0.0=="Female")

# Does the individual live in rural area? #

output$Urban <- as.numeric(grepl("urban",df$f.20118,ignore.case=TRUE))
output[as.character(output$f.20118) == "Postcode not linkable","Urban"] <- NA

# Townsend score to signify deprivation (IMD not consistent between nations) #
output$Townsend <- df$f.22189.0.0

# For Ethnicity, keep as factor #
output$Ethnicity <- NA
output[as.character(df$f.21000.0.0) %in% c("White","British","Irish","Any other white background"),"Ethnicity"] <- "White"
output[as.character(df$f.21000.0.0) %in% c("Asian or Asian British","Indian","Pakistani","Bangladeshi","Any other Asian background", "Chinese"),"Ethnicity"] <- "Asian"
output[as.character(df$f.21000.0.0) %in% c("Black or Black British","Caribbean","African","Any other Black background"),"Ethnicity"] <- "Black"
output[as.character(df$f.21000.0.0) %in% c("Other ethnic group","Do not know","Prefer not to answer","Mixed","White and Black Caribbean","White and Asian","Any other mixed background"),"Ethnicity"] <- "Other"
output$Ethnicity <- as.factor(output$Ethnicity)

# Predictors #

# Get if parents had history of diseases
family_history_columns <- colnames(df)[grepl("f.20107.0.",colnames(df))|
                                         grepl("f.20110.0.",colnames(df))|
                                         grepl("f.20111.0.",colnames(df))]
family_history <- df[family_history_columns]
family_history <- sapply(family_history,as.numeric)
family_history[is.na(family_history)] <- 0
# Get matrix of all family histories
family_history_matrix <- output[c()]
family_history_matrix$FamilyHxHipFracture <- as.numeric(apply(family_history==14,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxProstateCancer <- as.numeric(apply(family_history==13,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxDepression <- as.numeric(apply(family_history==12,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxParkinson <- as.numeric(apply(family_history==11,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxAlzheimer <- as.numeric(apply(family_history==10,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxDiabetes <- as.numeric(apply(family_history==9,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxHBP <- as.numeric(apply(family_history==8,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxChronicBronchitis <- as.numeric(apply(family_history==6,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxBreastCancer <- as.numeric(apply(family_history==5,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxBowelCancer <- as.numeric(apply(family_history==4,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxLungCancer <- as.numeric(apply(family_history==3,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxStroke <- as.numeric(apply(family_history==2,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxHeartDisease <- as.numeric(apply(family_history==1,FUN=any,MARGIN=1,na.rm=FALSE))
# For diabetes prediction, we care about family history in: Hip Fracture, Depression, Diabetes, HBP, CVD (Stroke & Heart Disease)
output$FamilyHxHipFracture <- family_history_matrix$FamilyHxHipFracture
output$FamilyHxDepression <- family_history_matrix$FamilyHxDepression
output$FamilyHxDiabetes <- family_history_matrix$FamilyHxDiabetes
output$FamilyHxHBP <- family_history_matrix$FamilyHxHBP
output$FamilyHxCVD <- as.numeric(family_history_matrix$FamilyHxStroke|family_history_matrix$FamilyHxHeartDisease)

# Blood pressure
dbp_automatic <- rowMeans(df[,c("f.4079.0.0","f.4079.0.1")],na.rm = FALSE)
dbp_manual <- rowMeans(df[,c("f.94.0.0","f.94.0.1")],na.rm = FALSE)
output$DBP <- coalesce(dbp_automatic,dbp_manual)
sbp_automatic <- rowMeans(df[,c("f.4080.0.0","f.4080.0.1")],na.rm = FALSE)
sbp_manual <- rowMeans(df[,c("f.93.0.0","f.93.0.1")],na.rm = FALSE)
output$SBP <- coalesce(sbp_automatic,sbp_manual)

# Body measures
output$BMI <- coalesce(df$f.21001.0.0,df$f.23099.0.0)
output$BodyFat <- df$f.23099.0.0
output$WaistCircumference <- df$f.48.0.0

# Education
all_education_columns_1 <- colnames(df)[grepl("f.6138.0.",colnames(df))]
all_education_columns_2 <- colnames(df)[grepl("f.10722.0.",colnames(df))]
all_education_columns <- c(all_education_columns_1,all_education_columns_2)
output$Education <- NA
output[apply(
  sapply(df[all_education_columns],FUN=function(x) {grepl("None",x)}),
  MARGIN = 1,
  FUN = any
),"Education"] <- "Level0"
output[apply(
  sapply(df[all_education_columns],FUN=function(x) {
    grepl("NQV",x)|grepl("CSE",x)|grepl("GCSE",x)
    }),
  MARGIN = 1,
  FUN = any
),"Education"] <- "Level1"
output[apply(
  sapply(df[all_education_columns],FUN=function(x) {
    grepl("AS levels",x)
  }),
  MARGIN = 1,
  FUN = any
),"Education"] <- "Level2"
output[apply(
  sapply(df[all_education_columns],FUN=function(x) {
    grepl("College",x)
  }),
  MARGIN = 1,
  FUN = any
),"Education"] <- "Level3"
output[apply(
  sapply(df[all_education_columns],FUN=function(x) {
    grepl("Prefer",x)
  }),
  MARGIN = 1,
  FUN = any
),"Education"] <- NA
output$Education <- as.factor(output$Education)

# Smoking status + Packyears
output$Smoker <- df$f.20116.0.0
output$Smoker[output$Smoker=="Prefer not to answer"] <- NA
output$PackYears <- df$f.20161.0.0
output$PackYears[output$Smoker == "Never"] <- 0

# Alcohol intake
output$AlcoholStatus <- as.character(df$f.20117.0.0)
output$AlcoholStatus[output$AlcoholStatus=="Current"] <- 1
output$AlcoholStatus[output$AlcoholStatus!=1] <- 0
output$AlcoholStatus <- as.numeric(output$AlcoholStatus)

# Drinks More than Three times a Week
output$AlcoholHigh[df$f.1558.0.0 == "Daily or almost daily"] <- 1
output$AlcoholHigh[grepl("Three",df$f.1558.0.0)] <- 1
output$AlcoholHigh[output$AlcoholHigh != 1] <- 0
output$AlcoholHigh[is.na(output$AlcoholHigh)] <- 0

# Menopause/HRT. If male just 0. If missing, Hysteroctomy and HRT set to 0.
output$HRT <- as.numeric(df$f.2814.0.0=="Yes")
output$HRT[is.na(output$HRT)] <- 0
output$Hysterectomy <- as.numeric(df$f.3591.0.0=="Yes")
output$Hysterectomy[is.na(output$Hysterectomy)] <- 0
output$Menopause <- as.numeric(df$f.2724.0.0=="Yes")
# If missing for menopause, set to 0 if under 45 and to 1 if over 55
# As given by https://www.nhs.uk/conditions/menopause/
output$Menopause[is.na(output$Menopause)&output$Female==0] <- 0
output$Menopause[is.na(output$Menopause)&output$Age<45] <- 0
output$Menopause[is.na(output$Menopause)&output$Age>55] <- 1

# Current Medications
medication_columns_1 <- colnames(df)[grepl("f.6177.0.",colnames(df))]
medication_columns_2 <- colnames(df)[grepl("f.6153.0.",colnames(df))]
medication_columns_3 <- colnames(df)[grepl("f.6154.0.",colnames(df))]
medication_columns_4 <- colnames(df)[grepl("f.10004.0.",colnames(df))]
medication_columns_5 <- colnames(df)[grepl("f.10005.0.",colnames(df))]
medication_columns <- c(
  medication_columns_1, medication_columns_2, medication_columns_3,
  medication_columns_4, medication_columns_5
)
output$MedicationCholesterol <- 0
output$MedicationBloodPressure <- 0
output$MedicationInsulin <- 0
output$MedicationHRT <- 0
output$MedicationContraceptive <- 0
output[apply(
  sapply(
  df[medication_columns],
  FUN=function(x) {
    grepl("cholesterol",as.character(x),ignore.case = TRUE)
  }
  ),
  MARGIN = 1,
  FUN = any
  ),"MedicationCholesterol"] <- 1
output[apply(
  sapply(
    df[medication_columns],
    FUN=function(x) {
      grepl("blood pressure",as.character(x),ignore.case = TRUE)
    }
  ),
  MARGIN = 1,
  FUN = any
),"MedicationBloodPressure"] <- 1
output[apply(
  sapply(
    df[medication_columns],
    FUN=function(x) {
      grepl("insulin",as.character(x),ignore.case = TRUE)
    }
  ),
  MARGIN = 1,
  FUN = any
),"MedicationInsulin"] <- 1
output[apply(
  sapply(
    df[medication_columns],
    FUN=function(x) {
      grepl("hormone replacement",as.character(x),ignore.case = TRUE)
    }
  ),
  MARGIN = 1,
  FUN = any
),"MedicationHRT"] <- 1
output[apply(
  sapply(
    df[medication_columns],
    FUN=function(x) {
      grepl("contraceptive",as.character(x),ignore.case = TRUE)
    }
  ),
  MARGIN = 1,
  FUN = any
),"MedicationContraceptive"] <- 1
# We don't actually want Insulin takers in our cohort, so we remove variable
# but keep it separate for exclusion
taking_insulin_to_exclude <- output$MedicationInsulin
output$MedicationInsulin <- NULL

# All first diagnosis
all_dates <- get_first_dates(df)
dateIndeces <- !is.na(all_dates)
all_dates[dateIndeces] <- as.Date(all_dates[dateIndeces])
ever_had_icd10 <- as.Date(df$f.53.0.0)
before_recruitment_diagnosis <- all_dates
for (ii in 1:dim(all_dates)[2]) {
  before_recruitment_diagnosis[,ii] <- all_dates[,ii]<=df$f.53.0.0
}
before_recruitment_diagnosis[is.na(before_recruitment_diagnosis)] <- 0
# For first diagnosis, only keep those relevant to diabetes
# Taken from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5860745/pdf/pone.0194127.pdf
# And https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8050730/pdf/main.pdf
# https://bmcmedicine.biomedcentral.com/counter/pdf/10.1186/1741-7015-9-103.pdf
# https://qdiabetes.org/
# Learning disabilities
# Psoriatic Arthritis
output$FirstIncidenceM07 <- before_recruitment_diagnosis$FirstIncidenceM07
# Juvenile Arthritis
output$FirstIncidenceM08 <- before_recruitment_diagnosis$FirstIncidenceM08
# Psoriasis
output$FirstIncidenceL40 <- before_recruitment_diagnosis$FirstIncidenceL40
# Diabetes during pregnancy
output$FirstIncidenceO24 <- before_recruitment_diagnosis$FirstIncidenceO24
# Preterm Labour
output$FirstIncidenceO60 <- before_recruitment_diagnosis$FirstIncidenceO60
# Arthrosis
output["FirstIncidenceM15to19"] <- before_recruitment_diagnosis %>%
  mutate(FirstIncidenceM15to19 = as.numeric(FirstIncidenceM15|FirstIncidenceM16|FirstIncidenceM17|
    FirstIncidenceM18|FirstIncidenceM19)) %>%
  select(FirstIncidenceM15to19) # Arthrosis
# Ovarian dysfunction (For PCOS)
output$FirstIncidenceE28 <- before_recruitment_diagnosis$FirstIncidenceE28 
# Schizophrenia, schizotypal and delusional disorders
output["FirstIncidenceF20to29"] <- before_recruitment_diagnosis %>%
  mutate(FirstIncidenceF20to29 = as.numeric(FirstIncidenceF20|FirstIncidenceF21|FirstIncidenceF22|
    FirstIncidenceF23|FirstIncidenceF24|FirstIncidenceF25|
    FirstIncidenceF28|FirstIncidenceF29)) %>%
  select(FirstIncidenceF20to29)
# Systemic connective tissue disorders (For Major Cell Arteritis)
output["FirstIncidenceM30to36"] <- before_recruitment_diagnosis %>%
  mutate(FirstIncidenceM30to36 = as.numeric(FirstIncidenceM30|FirstIncidenceM31|FirstIncidenceM32|
    FirstIncidenceM33|FirstIncidenceM34|FirstIncidenceM35|
    FirstIncidenceM36)) %>%
  select(FirstIncidenceM30to36)
# Mood disorders (For Depression)
output["FirstIncidenceF32to39"] <- before_recruitment_diagnosis %>%
  mutate(FirstIncidenceF32to39 = as.numeric(FirstIncidenceF32|FirstIncidenceF33|FirstIncidenceF34|
    FirstIncidenceF38|FirstIncidenceF39)) %>%
  select(FirstIncidenceF32to39)
# Bipolar disorder or Mania
output["FirstIncidenceF30to31"] <- before_recruitment_diagnosis %>%
  mutate(FirstIncidenceF30to31 = as.numeric(FirstIncidenceF30|FirstIncidenceF31)) %>%
  select(FirstIncidenceF30to31)
# Sleep disorder
output$FirstIncidenceG47 <- before_recruitment_diagnosis$FirstIncidenceG47
# Obesity
output$FirstIncidenceE66 <- before_recruitment_diagnosis$FirstIncidenceE66
# These CVD Codes are gotten from https://www.sonoraquest.com/media/1305/icd-10cardio_0417-2.pdf
# Cardiovascular and Ischaemic Disease, Circulatory System Diseases
output["FirstIncidenceI25to73"] <- before_recruitment_diagnosis %>%
  mutate(FirstIncidenceI25to73 = as.numeric(FirstIncidenceI20|FirstIncidenceI21|FirstIncidenceI25|
  FirstIncidenceI48|FirstIncidenceI50|FirstIncidenceI63|FirstIncidenceI64|FirstIncidenceI65|FirstIncidenceI67|
  FirstIncidenceI73)) %>%
  select(FirstIncidenceI25to73)
# Hypertensive Disease
output["FirstIncidenceI10to11"] <- before_recruitment_diagnosis %>%
  mutate(FirstIncidenceI10to11 = as.numeric(FirstIncidenceI10|FirstIncidenceI11)) %>%
  select(FirstIncidenceI10to11)
# Hypercholesterolemia, Hypertriglyceridemia, Hyperlipidemia
output["FirstIncidenceE78"] <- before_recruitment_diagnosis$FirstIncidenceE78
# Other metabolic disorders (Metabolic disorder, unspecified)
output$FirstIncidenceE88 <- before_recruitment_diagnosis$FirstIncidenceE88
# Gingivitis and periodontal diseases
output$FirstIncidenceK05 <- before_recruitment_diagnosis$FirstIncidenceK05
#  Disorders of psychological development (Learning disability)
output["FirstIncidenceF80to89"] <- before_recruitment_diagnosis %>%
  mutate(FirstIncidenceF80to89 = as.numeric(FirstIncidenceF80|FirstIncidenceF81|
    FirstIncidenceF83|FirstIncidenceF84|FirstIncidenceF88|FirstIncidenceF89)) %>%
  select(FirstIncidenceF80to89)


# For Cancer, identify if there are any diagnosis
# At the time of diagnosis
all_cancer_diagnoses_columns <- colnames(df)[grepl("f.40006.",colnames(df))]
all_cancer_dates_columns <- colnames(df)[grepl("f.40005.",colnames(df))]
all_cancer_diagnoses <- df[,all_cancer_diagnoses_columns]
all_cancer_dates <- df[,all_cancer_dates_columns]
dateIndeces <- !is.na(all_cancer_dates)
all_cancer_dates[dateIndeces] <- as.Date(all_cancer_dates[dateIndeces])
enrollment_dates <- matrix(rep((df$f.53.0.0),each=dim(all_cancer_dates)[2]),
                           ncol=dim(all_cancer_dates)[2], byrow=TRUE)
cancer_before_attendance <- sapply(data.frame(all_cancer_dates <= enrollment_dates),as.numeric)
cancer_before_attendance[is.na(cancer_before_attendance)] <- 0
# For each positive cancer_before_attendance, fetch code
previous_cancers <- df[c()]
where_is_cancer <- which(cancer_before_attendance == 1,arr.ind=TRUE)
for (ii in 1:dim(where_is_cancer)[1]) {
  icd_cancer_code <- substr(all_cancer_diagnoses[where_is_cancer[ii,1],where_is_cancer[ii,2]], start = 1, stop = 3)
  if (!is.na(icd_cancer_code)) {
    previous_cancers[where_is_cancer[ii,1],paste("Cancer",icd_cancer_code,sep="")] <- 1
  }
}

# Final refactoring
output <- droplevels(output)

#-----------------------------#
#-----------------------------#
#### Definition of outcome ####
#-----------------------------#
#-----------------------------#

# Which are the outputs ICD10 codes?
columns_of_output <- c("FirstIncidenceE11")

# Identify all ICD-10 diagnosis
all_dates <- get_first_dates(df)[columns_of_output]
within_5_years <- all_dates
for (ii in 1:dim(all_dates)[2]) {
  within_5_years[,ii] <- ((all_dates[,ii]>=df$f.53.0.0)&(all_dates[,ii]<=(df$f.53.0.0+1825)))
}
within_5_years[is.na(within_5_years)] <- 0

# Sum over all columns (if multiple)
within_5_years <- sapply(within_5_years,as.logical)
outcome_of_interest <- apply(within_5_years,MARGIN=1,FUN=any)
outcome_of_interest <- sapply(outcome_of_interest,as.numeric)

# Define it within dataframe
output$y <- outcome_of_interest

#-----------------------------------#
#-----------------------------------#
#### Cohort definition/Exclusion ####
#-----------------------------------#
#-----------------------------------#

# Define exclusion criteria in exclusion, which we will filter out at the end

# Only 40 to 70 year olds
exclusion <- (output$Age<40)|(output$Age>70)

# Without T2 or T1 diabetes before entrance
exclusion <- exclusion|(before_recruitment_diagnosis$FirstIncidenceE10==1)|(before_recruitment_diagnosis$FirstIncidenceE11==1)

# Exclude those taking insulin
exclusion <- exclusion|(taking_insulin_to_exclude==1)

# Exclude if death before 5 year follow-up
start_of_follow_up <- data.frame(as.Date(output$Started))
death_dates <- data.frame(sapply(output$Death, FUN = function(x) {
  return(as.Date(x, optional = TRUE))
}))
censor_dates <- data.frame(sapply(output$Censored, FUN = function(x) {
  return(as.Date(x, optional = TRUE))
}))
exclude_1 <- (death_dates<(start_of_follow_up+1825))&(output$y==0)
exclude_1[is.na(exclude_1)] <- FALSE
exclude_2 <- (censor_dates<(start_of_follow_up+1825))&(output$y==0)
exclude_2[is.na(exclude_2)] <- FALSE
exclude <- exclude_1|exclude_2
exclusion <- exclusion|exclude

# Save "Detailed" Output for Tables later on
output_detailed <- output
output_detailed$DetailedEthnicity <- as.character(df$f.21000.0.0)
output_detailed <- output_detailed[exclusion==0,]
write.csv(output_detailed,
          paste0("Data/detailed_diabetes_df",which_output_path,".csv"))


# APPLY EXCLUSION: If exclusion==1, then filter out
output <- output[exclusion==0,]


# Reshuffle
output <- output[sample(nrow(output)), ]

#------------------------#
#------------------------#
#### Define Groupings ####
#------------------------#
#------------------------#

# Set random seed
set.seed(42)

# Choose a random 15% set of population
output$is_val <- rbinom(n=nrow(output), size=1, prob=0.15)

#-----------------#
#-----------------#
#### Save Data ####
#-----------------#
#-----------------#

# Save data
write.csv(output,paste0("Data/diabetes_df",which_output_path,".csv"))

















#------------------------------#
#------------------------------#
#### 2) LUNG CANCER DATASET ####
#------------------------------#
#------------------------------#

# Define matrix we will base data on
output <- df[c()]

# Follow-Up and Date at which assesment was made #

output$Censored <- df$f.191.0.0
output$Started  <- df$f.53.0.0
output$Death    <- coalesce(df$f.40000.0.0)

# Assesment centre #

assesment_centres <- list(
  "10003"="Stockport (pilot)",
  "11001"="Manchester",
  "11002"="Oxford",
  "11003"="Cardiff",
  "11004"="Glasgow",
  "11005"="Edinburgh",
  "11006"="Stoke",
  "11007"="Reading",
  "11008"="Bury",
  "11009"="Newcastle",
  "11010"="Leeds",
  "11011"="Bristol",
  "11012"="Barts",
  "11013"="Nottingham",
  "11014"="Sheffield",
  "11016"="Liverpool",
  "11017"="Middlesborough",
  "11018"="Hounslow",
  "11020"="Croydon",
  "11021"="Birmingham",
  "11022"="Swansea",
  "11023"="Wrexham",
  "11024"="Cheadle (revisit)",
  "11025"="Cheadle (imaging)",
  "11026"="Reading (imaging)",
  "11027"="Newcastle (imaging)",
  "11028"="Bristol (imaging)"
)
assesment_column <- df["f.54.0.0"]
assesment_column <- assesment_centres[sapply(df["f.54.0.0"],as.character)]
output$Centre <- unlist(unname(assesment_column))
output$Centre <- as.factor(output$Centre)


# Demographic Characteristics #

output$Age <- df$f.21022.0
output$Female <- as.numeric(df$f.31.0.0=="Female")
# Does the individual live in rural area?
output$Urban <- as.numeric(grepl("urban",df$f.20118,ignore.case=TRUE))
output[as.character(output$f.20118) == "Postcode not linkable","Urban"] <- NA

# Townsend score to signify deprivation (IMD not consistent between nations)
output$Townsend <- df$f.22189.0.0

# For Ethnicity, keep as factor
output$Ethnicity <- NA
output[as.character(df$f.21000.0.0) %in% c("White","British","Irish","Any other white background"),"Ethnicity"] <- "White"
output[as.character(df$f.21000.0.0) %in% c("Asian or Asian British","Indian","Pakistani","Bangladeshi","Any other Asian background", "Chinese"),"Ethnicity"] <- "Asian"
output[as.character(df$f.21000.0.0) %in% c("Black or Black British","Caribbean","African","Any other Black background"),"Ethnicity"] <- "Black"
output[as.character(df$f.21000.0.0) %in% c("Other ethnic group","Do not know","Prefer not to answer","Mixed","White and Black Caribbean","White and Asian","Any other mixed background"),"Ethnicity"] <- "Other"
output$Ethnicity <- as.factor(output$Ethnicity)

# Predictors #

# Get if parents had history of diseases
family_history_columns <- colnames(df)[grepl("f.20107.0.",colnames(df))|
                                         grepl("f.20110.0.",colnames(df))|
                                         grepl("f.20111.0.",colnames(df))]
family_history <- df[family_history_columns]
family_history <- sapply(family_history,as.numeric)
family_history[is.na(family_history)] <- 0
# Get matrix of all family histories
family_history_matrix <- output[c()]
family_history_matrix$FamilyHxHipFracture <- as.numeric(apply(family_history==14,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxProstateCancer <- as.numeric(apply(family_history==13,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxDepression <- as.numeric(apply(family_history==12,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxParkinson <- as.numeric(apply(family_history==11,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxAlzheimer <- as.numeric(apply(family_history==10,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxDiabetes <- as.numeric(apply(family_history==9,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxHBP <- as.numeric(apply(family_history==8,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxChronicBronchitis <- as.numeric(apply(family_history==6,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxBreastCancer <- as.numeric(apply(family_history==5,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxBowelCancer <- as.numeric(apply(family_history==4,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxLungCancer <- as.numeric(apply(family_history==3,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxStroke <- as.numeric(apply(family_history==2,FUN=any,MARGIN=1,na.rm=FALSE))
family_history_matrix$FamilyHxHeartDisease <- as.numeric(apply(family_history==1,FUN=any,MARGIN=1,na.rm=FALSE))
# For diabetes prediction, we care about family history in: Lung Cancer, Other Cancers, HBP, Chronic Bronchitis, CVD (Stroke & Heart Disease)
output$FamilyHxOtherCancer <- family_history_matrix$FamilyHxProstateCancer|family_history_matrix$FamilyHxBreastCancer|family_history_matrix$FamilyHxBowelCancer
output$FamilyHxLungCancer <- family_history_matrix$FamilyHxLungCancer
output$FamilyHxHBP <- family_history_matrix$FamilyHxHBP
output$FamilyHxChronicBronchitis <- family_history_matrix$FamilyHxChronicBronchitis
output$FamilyHxCVD <- as.numeric(family_history_matrix$FamilyHxStroke|family_history_matrix$FamilyHxHeartDisease)

# Blood pressure
dbp_automatic <- rowMeans(df[,c("f.4079.0.0","f.4079.0.1")],na.rm = FALSE)
dbp_manual <- rowMeans(df[,c("f.94.0.0","f.94.0.1")],na.rm = FALSE)
output$DBP <- coalesce(dbp_automatic,dbp_manual)
sbp_automatic <- rowMeans(df[,c("f.4080.0.0","f.4080.0.1")],na.rm = FALSE)
sbp_manual <- rowMeans(df[,c("f.93.0.0","f.93.0.1")],na.rm = FALSE)
output$SBP <- coalesce(sbp_automatic,sbp_manual)

# Body measures
output$BMI <- coalesce(df$f.21001.0.0,df$f.23099.0.0)
output$BodyFat <- df$f.23099.0.0
output$WaistCircumference <- df$f.48.0.0

# Education
all_education_columns_1 <- colnames(df)[grepl("f.6138.0.",colnames(df))]
all_education_columns_2 <- colnames(df)[grepl("f.10722.0.",colnames(df))]
all_education_columns <- c(all_education_columns_1,all_education_columns_2)
output$Education <- NA
output[apply(
  sapply(df[all_education_columns],FUN=function(x) {grepl("None",x)}),
  MARGIN = 1,
  FUN = any
),"Education"] <- "Level0"
output[apply(
  sapply(df[all_education_columns],FUN=function(x) {
    grepl("NQV",x)|grepl("CSE",x)|grepl("GCSE",x)
  }),
  MARGIN = 1,
  FUN = any
),"Education"] <- "Level1"
output[apply(
  sapply(df[all_education_columns],FUN=function(x) {
    grepl("AS levels",x)
  }),
  MARGIN = 1,
  FUN = any
),"Education"] <- "Level2"
output[apply(
  sapply(df[all_education_columns],FUN=function(x) {
    grepl("College",x)
  }),
  MARGIN = 1,
  FUN = any
),"Education"] <- "Level3"
output[apply(
  sapply(df[all_education_columns],FUN=function(x) {
    grepl("Prefer",x)
  }),
  MARGIN = 1,
  FUN = any
),"Education"] <- NA
output$Education <- as.factor(output$Education)

# Smoking status + Packyears
output$Smoker <- df$f.20116.0.0
output$Smoker[output$Smoker=="Prefer not to answer"] <- NA
output$PackYears <- df$f.20161.0.0
output$PackYears[output$Smoker == "Never"] <- 0

# Alcohol intake
output$AlcoholStatus <- as.character(df$f.20117.0.0)
output$AlcoholStatus[output$AlcoholStatus=="Current"] <- 1
output$AlcoholStatus[output$AlcoholStatus!=1] <- 0
output$AlcoholStatus <- as.numeric(output$AlcoholStatus)

# Drinks More than Three times a Week
output$AlcoholHigh[df$f.1558.0.0 == "Daily or almost daily"] <- 1
output$AlcoholHigh[grepl("Three",df$f.1558.0.0)] <- 1
output$AlcoholHigh[output$AlcoholHigh != 1] <- 0
output$AlcoholHigh[is.na(output$AlcoholHigh)] <- 0

# Menopause/HRT. If male just 0. If missing, Hysteroctomy and HRT set to 0.
output$HRT <- as.numeric(df$f.2814.0.0=="Yes")
output$HRT[is.na(output$HRT)] <- 0
output$Hysterectomy <- as.numeric(df$f.3591.0.0=="Yes")
output$Hysterectomy[is.na(output$Hysterectomy)] <- 0
output$Menopause <- as.numeric(df$f.2724.0.0=="Yes")
# If missing for menopause, set to 0 if under 45 and to 1 if over 55
# As given by https://www.nhs.uk/conditions/menopause/
output$Menopause[is.na(output$Menopause)&output$Female==0] <- 0
output$Menopause[is.na(output$Menopause)&output$Age<45] <- 0
output$Menopause[is.na(output$Menopause)&output$Age>55] <- 1

# Current Medications
medication_columns_1 <- colnames(df)[grepl("f.6177.0.",colnames(df))]
medication_columns_2 <- colnames(df)[grepl("f.6153.0.",colnames(df))]
medication_columns_3 <- colnames(df)[grepl("f.6154.0.",colnames(df))]
medication_columns_4 <- colnames(df)[grepl("f.10004.0.",colnames(df))]
medication_columns_5 <- colnames(df)[grepl("f.10005.0.",colnames(df))]
medication_columns <- c(
  medication_columns_1, medication_columns_2, medication_columns_3,
  medication_columns_4, medication_columns_5
)
output$MedicationCholesterol <- 0
output$MedicationBloodPressure <- 0
output$MedicationInsulin <- 0
output$MedicationHRT <- 0
output$MedicationContraceptive <- 0
output[apply(
  sapply(
    df[medication_columns],
    FUN=function(x) {
      grepl("cholesterol",as.character(x),ignore.case = TRUE)
    }
  ),
  MARGIN = 1,
  FUN = any
),"MedicationCholesterol"] <- 1
output[apply(
  sapply(
    df[medication_columns],
    FUN=function(x) {
      grepl("blood pressure",as.character(x),ignore.case = TRUE)
    }
  ),
  MARGIN = 1,
  FUN = any
),"MedicationBloodPressure"] <- 1
output[apply(
  sapply(
    df[medication_columns],
    FUN=function(x) {
      grepl("insulin",as.character(x),ignore.case = TRUE)
    }
  ),
  MARGIN = 1,
  FUN = any
),"MedicationInsulin"] <- 1
output[apply(
  sapply(
    df[medication_columns],
    FUN=function(x) {
      grepl("hormone replacement",as.character(x),ignore.case = TRUE)
    }
  ),
  MARGIN = 1,
  FUN = any
),"MedicationHRT"] <- 1
output[apply(
  sapply(
    df[medication_columns],
    FUN=function(x) {
      grepl("contraceptive",as.character(x),ignore.case = TRUE)
    }
  ),
  MARGIN = 1,
  FUN = any
),"MedicationContraceptive"] <- 1

# All first diagnosis
all_dates <- get_first_dates(df)
dateIndeces <- !is.na(all_dates)
all_dates[dateIndeces] <- as.Date(all_dates[dateIndeces])
ever_had_icd10 <- as.Date(df$f.53.0.0)
before_recruitment_diagnosis <- all_dates
for (ii in 1:dim(all_dates)[2]) {
  before_recruitment_diagnosis[,ii] <- all_dates[,ii]<=df$f.53.0.0
}
before_recruitment_diagnosis[is.na(before_recruitment_diagnosis)] <- 0
# For first diagnosis, only keep those relevant to lung cancer
# Taken from https://www.nejm.org/doi/pdf/10.1056/NEJMoa1211776?articleTools=true (PLCOm2012)
# And https://www.sciencedirect.com/science/article/pii/S1525730415002715?via%3Dihub (Review)
# https://www.thelancet.com/action/showPdf?pii=S2213-2600%2823%2900050-4
# Asthma
output$FirstIncidenceJ45 <- before_recruitment_diagnosis$FirstIncidenceJ45
# COPD
output["FirstIncidenceJ43to44"] <- before_recruitment_diagnosis %>%
  mutate(FirstIncidenceJ43to44 = as.numeric(FirstIncidenceJ43|FirstIncidenceJ44)) %>%
  select(FirstIncidenceJ43to44)
# Pneumonia
output["FirstIncidenceJ12to18"] <- before_recruitment_diagnosis %>%
  mutate(FirstIncidenceJ12to18 = as.numeric(FirstIncidenceJ12|FirstIncidenceJ13|
                                              FirstIncidenceJ14|FirstIncidenceJ15|
                                              FirstIncidenceJ16|FirstIncidenceJ17|
                                              FirstIncidenceJ18)) %>%
  select(FirstIncidenceJ12to18)
# Respiratory Tuberculosis
output["FirstIncidenceA15to16"] <- before_recruitment_diagnosis %>%
  mutate(FirstIncidenceA15to16 = as.numeric(FirstIncidenceA15|FirstIncidenceA16)) %>%
  select(FirstIncidenceA15to16)

# For Cancer, identify if there are any diagnosis
# At the time of diagnosis
all_cancer_diagnoses_columns <- colnames(df)[grepl("f.40006.",colnames(df))]
all_cancer_dates_columns <- colnames(df)[grepl("f.40005.",colnames(df))]
all_cancer_diagnoses <- df[,all_cancer_diagnoses_columns]
all_cancer_dates <- df[,all_cancer_dates_columns]
dateIndeces <- !is.na(all_cancer_dates)
all_cancer_dates[dateIndeces] <- as.Date(all_cancer_dates[dateIndeces])
enrollment_dates <- matrix(rep((df$f.53.0.0),each=dim(all_cancer_dates)[2]),
                           ncol=dim(all_cancer_dates)[2], byrow=TRUE)
cancer_before_attendance <- sapply(data.frame(all_cancer_dates <= enrollment_dates),as.numeric)
cancer_before_attendance[is.na(cancer_before_attendance)] <- 0
# For each positive cancer_before_attendance, fetch code
previous_cancers <- df[c()]
where_is_cancer <- which(cancer_before_attendance == 1,arr.ind=TRUE)
for (ii in 1:dim(where_is_cancer)[1]) {
  icd_cancer_code <- substr(all_cancer_diagnoses[where_is_cancer[ii,1],where_is_cancer[ii,2]], start = 1, stop = 3)
  if (!is.na(icd_cancer_code)) {
    previous_cancers[where_is_cancer[ii,1],paste("Cancer",icd_cancer_code,sep="")] <- 1
  }
}
previous_cancers[is.na(previous_cancers)] <- 0
# We only include if there is ANY malignant cancer
any_cancer <- apply(previous_cancers[ , grepl( "CancerC" , names( previous_cancers ) ) ],1,max)
output$AnyCancer <- any_cancer

# Final refactoring
output <- droplevels(output)

#-----------------------------#
#-----------------------------#
#### Definition of outcome ####
#-----------------------------#
#-----------------------------#

# Which are the outputs ICD10 codes?
columns_of_output <- c("CancerC34")

# Identify all Cancer diagnoses within 5 years
all_cancer_diagnoses_columns <- colnames(df)[grepl("f.40006.",colnames(df))]
all_cancer_dates_columns <- colnames(df)[grepl("f.40005.",colnames(df))]
all_cancer_diagnoses <- df[,all_cancer_diagnoses_columns]
all_cancer_dates <- df[,all_cancer_dates_columns]
dateIndeces <- !is.na(all_cancer_dates)
all_cancer_dates[dateIndeces] <- as.Date(all_cancer_dates[dateIndeces])
enrollment_dates <- matrix(rep((df$f.53.0.0),each=dim(all_cancer_dates)[2]),
                           ncol=dim(all_cancer_dates)[2], byrow=TRUE)
cancer_before_5y <- sapply(data.frame((all_cancer_dates >= enrollment_dates)&(all_cancer_dates <= enrollment_dates+2190)),as.numeric)
cancer_before_5y[is.na(cancer_before_5y)] <- 0
# For each positive cancer_before_5y, fetch code
fiveyear_cancers <- df[c()]
where_is_cancer <- which(cancer_before_5y == 1,arr.ind=TRUE)
for (ii in 1:dim(where_is_cancer)[1]) {
  icd_cancer_code <- substr(all_cancer_diagnoses[where_is_cancer[ii,1],where_is_cancer[ii,2]], start = 1, stop = 3)
  if (!is.na(icd_cancer_code)) {
    fiveyear_cancers[where_is_cancer[ii,1],paste("Cancer",icd_cancer_code,sep="")] <- 1
  }
}
fiveyear_cancers[is.na(fiveyear_cancers)] <- 0

# Assign as positive only if lung cancer
outcome_of_interest <- fiveyear_cancers[ ,columns_of_output]

# Define it within dataframe
output$y <- outcome_of_interest

#-----------------------------------#
#-----------------------------------#
#### Cohort definition/Exclusion ####
#-----------------------------------#
#-----------------------------------#

# Define exclusion criteria, which we will filter out at the end

# Only 40 to 70 year olds
exclusion <- (output$Age<40)|(output$Age>70)

# If there is a cancer in tha past 5 years, excluding melanoma, exclude
# At the time of diagnosis
all_cancer_diagnoses_columns <- colnames(df)[grepl("f.40006.",colnames(df))]
all_cancer_dates_columns <- colnames(df)[grepl("f.40005.",colnames(df))]
all_cancer_diagnoses <- df[,all_cancer_diagnoses_columns]
all_cancer_dates <- df[,all_cancer_dates_columns]
dateIndeces <- !is.na(all_cancer_dates)
all_cancer_dates[dateIndeces] <- as.Date(all_cancer_dates[dateIndeces])
enrollment_dates <- matrix(rep((df$f.53.0.0),each=dim(all_cancer_dates)[2]),
                           ncol=dim(all_cancer_dates)[2], byrow=TRUE)
cancer_before_attendance <- sapply(data.frame((all_cancer_dates >= enrollment_dates - 1825)&(all_cancer_dates <= enrollment_dates)),as.numeric)
cancer_before_attendance[is.na(cancer_before_attendance)] <- 0
# For each positive cancer_before_attendance, fetch code
cancer_last5years <- df[c()]
where_is_cancer <- which(cancer_before_attendance == 1,arr.ind=TRUE)
for (ii in 1:dim(where_is_cancer)[1]) {
  icd_cancer_code <- substr(all_cancer_diagnoses[where_is_cancer[ii,1],where_is_cancer[ii,2]], start = 1, stop = 3)
  if (!is.na(icd_cancer_code)) {
    cancer_last5years[where_is_cancer[ii,1],paste("Cancer",icd_cancer_code,sep="")] <- 1
  }
}
cancer_last5years[is.na(cancer_last5years)] <- 0
# Exclude melanomas
cancer_last5years$CancerC44 <- NULL
cancer_last5years$CancerC45 <- NULL
# We only exclude if there is ANY malignant cancer in last 5 years, excluding melanoma
any_cancer_last5years <- apply(cancer_last5years[ , grepl( "CancerC" , names( cancer_last5years ) ) ],1,max)
# Add cancers in the last 5 years as exclusion
exclusion <- exclusion|(any_cancer_last5years==1)

# Exclude if death before 6 year follow-up
start_of_follow_up <- data.frame(as.Date(output$Started))
death_dates <- data.frame(sapply(output$Death, FUN = function(x) {
  return(as.Date(x, optional = TRUE))
}))
censor_dates <- data.frame(sapply(output$Censored, FUN = function(x) {
  return(as.Date(x, optional = TRUE))
}))
exclude_1 <- (death_dates<(start_of_follow_up+2190))&(output$y==0)
exclude_1[is.na(exclude_1)] <- FALSE
exclude_2 <- (censor_dates<(start_of_follow_up+2190))&(output$y==0)
exclude_2[is.na(exclude_2)] <- FALSE
exclude <- exclude_1|exclude_2
exclusion <- exclusion|exclude

mean(output$y)
mean(output$y[exclusion==0])

# APPLY EXCLUSION: If exclusion==1, then filter out
output <- output[exclusion==0,]

# Reshuffle
output <- output[sample(nrow(output)), ]

#------------------------#
#------------------------#
#### Define Groupings ####
#------------------------#
#------------------------#

# Set random seed
set.seed(42)

# Choose a random 15% set of population
output$is_val <- rbinom(n=nrow(output), size=1, prob=0.15)

#-----------------#
#-----------------#
#### Save Data ####
#-----------------#
#-----------------#

write.csv(output,paste0("Data/Synth_UKBB/Clean/lungcancer_df",which_output_path,".csv"))

