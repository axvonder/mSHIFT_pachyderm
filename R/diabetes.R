#################### functions for diabetes calculations ###################

# Sigmoid function
Sigmoid <- function(x) {
  z <- 1 / (1 + exp(-x))
  return(z)
}

# Load Diabetes_dict
Diabetes_dict <- list(
  CARDIA = list(
    `Age Group` = 0.295,
    Black = -0.055,
    Male = -0.958,
    BMI = 0.083,
    `Parental History` = 0.507,
    Smoker = -0.13,
    `High SBP` = 1.347,
    `High Cholesterol` = 0.431,
    Time = 10,
    Constant = -5.171
  ),
  `CARDIA-10` = list(
    `Age Group` = 0.217,
    Black = 0.342,
    Male = 0.322,
    BMI = 0.134,
    `Parental History` = 0.857,
    Smoker = -0.013,
    `High SBP` = 0.09,
    `High Cholesterol` = 0.328,
    Time = 10,
    Constant = -7.516
  ),
  ARIC = list(
    `Age Group` = 0.076,
    Black = 0.280,
    Male = 0.454,
    BMI = 0.130,
    `Parental History` = 0.626,
    Smoker = 0.305,
    `High SBP` = 0.386,
    `High Cholesterol` = 0.002,
    Time = 9,
    Constant = -6.662
  ),
  CHS = list(
    `Age Group` = -0.164,
    Black = 0.235,
    Male = 0.414,
    BMI = 0.135,
    `Parental History` = 0.281,
    Smoker = 0.181,
    `High SBP` = 0.635,
    `High Cholesterol` = -0.054,
    Time = 7,
    Constant = -7.405
  )
)



# Diabetes_risk_Alva function (from diabetes.R)
Diabetes_risk_Alva <- function(row, Diabetes_dict) {
  Age <- as.numeric(row['Age'])
  Ethnicity <- as.numeric(row['Ethnicity'])
  Sex <- as.numeric(row['Sex'])
  BMI <- as.numeric(row['BMI'])
  Parental_History <- as.numeric(row['Parental Diabetes History'])
  Smoker <- as.numeric(row['Current Smoker'])
  SBP <- as.numeric(row['Systolic Blood Pressure'])
  Cholesterol <- as.numeric(row['Total Cholesterol'])
  HR_CM <- as.numeric(row['PM_HR_Diabetes'])
  HR_RM <- as.numeric(row['RM_HR_Diabetes'])
  
  if (Age < 35) {
    dict <- Diabetes_dict$CARDIA
  } else if (Age < 55) {
    dict <- Diabetes_dict$`CARDIA-10`
  } else if (Age < 75) {
    dict <- Diabetes_dict$ARIC
  } else {
    dict <- Diabetes_dict$CHS
  }
  
  S <- 0
  S <- S + dict$`Age Group`
  
  if (Ethnicity == 4) {
    S <- S + dict$Black
  }
  if (Sex == 1) {
    S <- S + dict$Male
  }
  
  S <- S + BMI * dict$BMI
  S <- S + Parental_History * dict$`Parental History`
  S <- S + Smoker * dict$Smoker
  
  if (SBP > 140) {
    S <- S + dict$`High SBP`
  }
  if (Cholesterol > 240) {
    S <- S + dict$`High Cholesterol`
  }
  
  S <- S + dict$Constant
  P <- Sigmoid(S)
  Diabetes_rate <- -log(1 - P) / dict$Time
  Diabetes_risk <- 1 - exp(-Diabetes_rate)
  
  Diabetes_risk <- Diabetes_risk * HR_CM * HR_RM
  
  # Calibration
  if (Age < 45) {
    Diabetes_risk <- Diabetes_risk / 2.37
  } else if (Age < 65) {
    Diabetes_risk <- Diabetes_risk / 1.73
  } else {
    Diabetes_risk <- Diabetes_risk / 1.62
  }
  
  return(Diabetes_risk)
}

# Hazard_ratio_CM_Diabetes function
Hazard_ratio_CM_Diabetes <- function(row, x_PM) {
  row['Processed meat intake'] <- x_PM * as.numeric(row['Processed meat intake'])
  
  PM_intake <- as.numeric(row['Processed meat intake'])
  if (PM_intake > 0) {
    if (PM_intake < 15) {
      HR_best_fit <- 0.02 * PM_intake + 1
    } else {
      HR_best_fit <- 0.0022 * PM_intake + 1.27
    }
    
    if (PM_intake < 15) {
      HR_2sigma <- 0.03 * PM_intake + 1
    } else {
      HR_2sigma <- 0.00366 * PM_intake + 1.28
    }
    
    sigma <- (HR_2sigma - HR_best_fit) / 2
    Hazard_ratio <- max(rnorm(1, mean = HR_best_fit, sd = sigma), 1)
  } else {
    Hazard_ratio <- 1
  }
  
  return(Hazard_ratio)
}

# Hazard_ratio_RM_Diabetes function
Hazard_ratio_RM_Diabetes <- function(row, x_RM) {
  row['Red meat intake'] <- x_RM * as.numeric(row['Red meat intake'])
  
  RM_intake <- as.numeric(row['Red meat intake'])
  if (RM_intake > 0) {
    if (RM_intake < 40) {
      HR_best_fit <- 0.004 * RM_intake + 1
      HR_2sigma <- 0.0078 * RM_intake + 1
    } else {
      HR_best_fit <- 0.004 * RM_intake + 1.0
      HR_2sigma <- 0.004 * RM_intake + 1.19
    }
    
    sigma <- (HR_2sigma - HR_best_fit) / 2
    Hazard_ratio <- max(rnorm(1, mean = HR_best_fit, sd = sigma), 1)
  } else {
    Hazard_ratio <- 1
  }
  
  return(Hazard_ratio)
}
