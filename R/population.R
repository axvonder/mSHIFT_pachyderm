

#####################       population.py       #####################



# Function to perform demographic binning (if needed)
demographic_bins <- function(df) {
  sex_mapping <- c('1' = 'Male', '2' = 'Female')
  ethnicity_mapping <- c(
    '1' = 'Mexican American',
    '2' = 'Other Hispanic',
    '3' = 'Non-hispanic white',
    '4' = 'Non-hispanic black',
    '6' = 'Non-hispanic Asian',
    '7' = 'Other, including multi-racial'
  )
  income_mapping <- c(
    '1' = "AHI<25k", '2' = "AHI<25k", '3' = "AHI<25k", '4' = "AHI<25k", '5' = "AHI<25k",
    '6' = "25k<AHI<55k", '7' = "25k<AHI<55k", '8' = "25k<AHI<55k",
    '9' = "55k<AHI<100k", '10' = "55k<AHI<100k", '11' = "55k<AHI<100k",
    '12' = "AHI>100k"
  )
  df$Sex <- sex_mapping[as.character(df$Sex)]
  df$Ethnicity <- ethnicity_mapping[as.character(df$Ethnicity)]
  df$`Annual Household Income` <- income_mapping[as.character(df$`Annual Household Income`)]
  df$Age <- cut(
    df$Age,
    breaks = c(-Inf, 49, 64, 79, Inf),
    labels = c("18-49", "50-64", "65-79", ">80")
  )
  return(df)
}



# Include necessary functions from other scripts
# Since we are converting all code, we'll include all dependencies here.

# CRC_family_history function (from colorectal_cancer.R)
CRC_family_history <- function(row) {
  family_history_prob <- 0.034
  sample <- runif(1)
  if (sample < family_history_prob) {
    return(1)
  } else {
    return(0)
  }
}

# Hazard_ratio_CM_Diabetes function (from diabetes.R)
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

# Hazard_ratio_RM_Diabetes function (from diabetes.R)
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

# Hazard_ratio_CM_CVD function (from CVD.R)
Hazard_ratio_CM_CVD <- function(row, x_PM) {
  row['Processed meat intake'] <- x_PM * as.numeric(row['Processed meat intake'])
  serving_size <- as.numeric(row['Processed meat intake']) / 30  # Serving size in Zhong et al.
  
  if (serving_size == 0) {
    Hazard_ratio <- 1
  } else if (serving_size < 0.1) {
    Hazard_ratio <- rnorm(1, mean = 1.0267, sd = 0.005)
  } else if (serving_size < 0.2) {
    Hazard_ratio <- rnorm(1, mean = 1.052, sd = 0.015)
  } else if (serving_size < 0.3) {
    Hazard_ratio <- rnorm(1, mean = 1.0759, sd = 0.02)
  } else if (serving_size < 0.4) {
    Hazard_ratio <- rnorm(1, mean = 1.0981, sd = 0.025)
  } else if (serving_size < 0.5) {
    Hazard_ratio <- rnorm(1, mean = 1.1186, sd = 0.035)
  } else if (serving_size < 0.6) {
    Hazard_ratio <- rnorm(1, mean = 1.1372, sd = 0.035)
  } else if (serving_size < 0.7) {
    Hazard_ratio <- rnorm(1, mean = 1.1539, sd = 0.04)
  } else if (serving_size < 0.8) {
    Hazard_ratio <- rnorm(1, mean = 1.1685, sd = 0.045)
  } else if (serving_size < 0.9) {
    Hazard_ratio <- rnorm(1, mean = 1.181, sd = 0.045)
  } else if (serving_size < 1.0) {
    Hazard_ratio <- rnorm(1, mean = 1.1912, sd = 0.05)
  } else if (serving_size < 1.1) {
    Hazard_ratio <- rnorm(1, mean = 1.1992, sd = 0.05)
  } else if (serving_size < 1.2) {
    Hazard_ratio <- rnorm(1, mean = 1.2049, sd = 0.05)
  } else if (serving_size < 1.3) {
    Hazard_ratio <- rnorm(1, mean = 1.2082, sd = 0.055)
  } else if (serving_size < 1.4) {
    Hazard_ratio <- rnorm(1, mean = 1.2092, sd = 0.055)
  } else if (serving_size < 1.5) {
    Hazard_ratio <- rnorm(1, mean = 1.2077, sd = 0.055)
  } else if (serving_size < 1.6) {
    Hazard_ratio <- rnorm(1, mean = 1.2039, sd = 0.056)
  } else if (serving_size < 1.7) {
    Hazard_ratio <- rnorm(1, mean = 1.1978, sd = 0.056)
  } else if (serving_size < 1.8) {
    Hazard_ratio <- rnorm(1, mean = 1.1893, sd = 0.06)
  } else if (serving_size < 1.9) {
    Hazard_ratio <- rnorm(1, mean = 1.1786, sd = 0.07)
  } else if (serving_size < 2.0) {
    Hazard_ratio <- rnorm(1, mean = 1.1657, sd = 0.075)
  } else {
    Hazard_ratio <- rnorm(1, mean = 1.1507, sd = 0.075)
  }
  
  return(Hazard_ratio)
}

# Hazard_ratio_RM_CVD function (from CVD.R)
Hazard_ratio_RM_CVD <- function(row, x_RM) {
  row['Red meat intake'] <- x_RM * as.numeric(row['Red meat intake'])
  serving_size <- (as.numeric(row['Red meat intake']) * 0.035274) / 4  # Convert grams to oz and normalize
  
  if (serving_size == 0) {
    Hazard_ratio <- 1
  } else if (serving_size < 0.1) {
    Hazard_ratio <- rnorm(1, mean = 1.0176, sd = 0.006)
  } else if (serving_size < 0.2) {
    Hazard_ratio <- rnorm(1, mean = 1.0347, sd = 0.015)
  } else if (serving_size < 0.3) {
    Hazard_ratio <- rnorm(1, mean = 1.0512, sd = 0.022)
  } else if (serving_size < 0.4) {
    Hazard_ratio <- rnorm(1, mean = 1.0672, sd = 0.025)
  } else if (serving_size < 0.5) {
    Hazard_ratio <- rnorm(1, mean = 1.0825, sd = 0.035)
  } else if (serving_size < 0.6) {
    Hazard_ratio <- rnorm(1, mean = 1.0972, sd = 0.035)
  } else if (serving_size < 0.7) {
    Hazard_ratio <- rnorm(1, mean = 1.1112, sd = 0.045)
  } else if (serving_size < 0.8) {
    Hazard_ratio <- rnorm(1, mean = 1.1246, sd = 0.045)
  } else if (serving_size < 0.9) {
    Hazard_ratio <- rnorm(1, mean = 1.1371, sd = 0.05)
  } else if (serving_size < 1.0) {
    Hazard_ratio <- rnorm(1, mean = 1.1489, sd = 0.05)
  } else if (serving_size < 1.1) {
    Hazard_ratio <- rnorm(1, mean = 1.1599, sd = 0.05)
  } else if (serving_size < 1.2) {
    Hazard_ratio <- rnorm(1, mean = 1.1701, sd = 0.06)
  } else if (serving_size < 1.3) {
    Hazard_ratio <- rnorm(1, mean = 1.1794, sd = 0.065)
  } else if (serving_size < 1.4) {
    Hazard_ratio <- rnorm(1, mean = 1.1879, sd = 0.075)
  } else if (serving_size < 1.5) {
    Hazard_ratio <- rnorm(1, mean = 1.1955, sd = 0.075)
  } else if (serving_size < 1.6) {
    Hazard_ratio <- rnorm(1, mean = 1.2022, sd = 0.075)
  } else if (serving_size < 1.7) {
    Hazard_ratio <- rnorm(1, mean = 1.208, sd = 0.075)
  } else if (serving_size < 1.8) {
    Hazard_ratio <- rnorm(1, mean = 1.2128, sd = 0.085)
  } else if (serving_size < 1.9) {
    Hazard_ratio <- rnorm(1, mean = 1.2167, sd = 0.085)
  } else if (serving_size < 2.0) {
    Hazard_ratio <- rnorm(1, mean = 1.2196, sd = 0.085)
  } else {
    Hazard_ratio <- rnorm(1, mean = 1.2216, sd = 0.1)
  }
  
  return(Hazard_ratio)
}

# Relative_risk_Zhao_CM function (from colorectal_cancer.R)
Relative_risk_Zhao_CM <- function(row, x_PM) {
  row['Processed meat intake'] <- x_PM * as.numeric(row['Processed meat intake'])
  PM_intake <- as.numeric(row['Processed meat intake'])
  
  if (PM_intake > 0) {
    gradient_processed_meat <- rnorm(1, mean = 0.0044, sd = 0.002)
    RR_processed_meat <- max(PM_intake * gradient_processed_meat + 1, 1)
  } else {
    RR_processed_meat <- 1
  }
  
  return(RR_processed_meat)
}

# Relative_risk_Zhao_RM function (from colorectal_cancer.R)
Relative_risk_Zhao_RM <- function(row, x_RM) {
  row['Red meat intake'] <- x_RM * as.numeric(row['Red meat intake'])
  RM_intake <- as.numeric(row['Red meat intake'])
  
  if (RM_intake > 0) {
    gradient_red_meat <- rnorm(1, mean = 0.0016, sd = 0.0005)
    RR_red_meat <- max(RM_intake * gradient_red_meat + 1, 1)
  } else {
    RR_red_meat <- 1
  }
  
  return(RR_red_meat)
}

# Run sampling function
run_sampling <- function(df, x_PM, x_RM, seed) {
  set.seed(seed)
  
  df$`CRC family history` <- apply(df, 1, CRC_family_history)
  df$PM_HR_Diabetes <- apply(df, 1, function(row) Hazard_ratio_CM_Diabetes(row, x_PM))
  df$RM_HR_Diabetes <- apply(df, 1, function(row) Hazard_ratio_RM_Diabetes(row, x_RM))
  df$PM_HR_CVD <- apply(df, 1, function(row) Hazard_ratio_CM_CVD(row, x_PM))
  df$RM_HR_CVD <- apply(df, 1, function(row) Hazard_ratio_RM_CVD(row, x_RM))
  df$PM_HR_CRC <- apply(df, 1, function(row) Relative_risk_Zhao_CM(row, x_PM))
  df$RM_HR_CRC <- apply(df, 1, function(row) Relative_risk_Zhao_RM(row, x_RM))
  
  return(df)
}

# update_smoking_pack_years function
update_smoking_pack_years <- function(row) {
  smoking_pack_years <- 0
  
  if (as.numeric(row['Smoking pack years']) != 0) {
    smoking_pack_years <- as.numeric(row['Smoking pack years']) + as.numeric(row['Daily Cigarettes']) / 20
  }
  
  return(smoking_pack_years)
}

# update_SBP function
update_SBP <- function(row) {
  SBP <- as.numeric(row['Systolic Blood Pressure'])
  
  if (as.numeric(row['Sex']) == 0) {
    if (as.numeric(row['Age']) < 50) {
      SBP <- SBP + 0.288
    } else {
      SBP <- SBP + 0.308
    }
  } else if (as.numeric(row['Sex']) == 1) {
    if (as.numeric(row['Age']) < 50) {
      SBP <- SBP + 0.430
    } else {
      SBP <- SBP + 0.613
    }
  }
  
  return(SBP)
}

# update_df function
update_df <- function(df) {
  df$Age <- df$Age + 1
  df$`Smoking pack years` <- apply(df, 1, update_smoking_pack_years)
  df$`Systolic Blood Pressure` <- apply(df, 1, update_SBP)
  return(df)
}

# Sigmoid function (from diabetes.R)
Sigmoid <- function(x) {
  z <- 1 / (1 + exp(-x))
  return(z)
}

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

# Load Diabetes_dict (from diabetes.R)
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

# CVD_Framingham_coeff_2008_update (from CVD.R)
CVD_Framingham_coeff_2008_update <- list(
  `1` = list(
    Age = 3.06117,
    TC = 1.12370,
    HDL = -0.93263,
    SBP_treated = 1.99881,
    SBP_untreated = 1.93303,
    Diabetes = 0.57367,
    Smoker = 0.65451,
    Baseline = 0.88936
  ),
  `2` = list(
    Age = 2.32888,
    TC = 1.20904,
    HDL = -0.70833,
    SBP_treated = 2.82263,
    SBP_untreated = 2.76157,
    Diabetes = 0.69154,
    Smoker = 0.52873,
    Baseline = 0.95012
  )
)

# updated_Framingham_CVD_risk function (from CVD.R)
updated_Framingham_CVD_risk <- function(row, with_diabetes) {
  age <- as.numeric(row['Age'])
  sex <- as.numeric(row['Sex'])
  tot_cholesterol <- as.numeric(row['Total Cholesterol'])
  HDL_cholesterol <- as.numeric(row['HDL Cholesterol'])
  taking_blood_pressure_medication <- as.numeric(row['Taking BPM'])
  systolic_blood_pressure <- as.numeric(row['Systolic Blood Pressure'])
  smoker <- as.numeric(row['Current Smoker'])
  
  if (with_diabetes) {
    diabetes <- 1
  } else {
    diabetes <- as.numeric(row['Diabetes'])
  }
  
  CVD_dict <- CVD_Framingham_coeff_2008_update[[as.character(sex)]]
  
  S <- 0
  S <- S + log(age) * CVD_dict$Age
  S <- S + log(tot_cholesterol) * CVD_dict$TC
  S <- S + log(HDL_cholesterol) * CVD_dict$HDL
  if (taking_blood_pressure_medication == 1) {
    S <- S + log(systolic_blood_pressure) * CVD_dict$SBP_treated
  } else {
    S <- S + log(systolic_blood_pressure) * CVD_dict$SBP_untreated
  }
  if (diabetes == 1) {
    S <- S + CVD_dict$Diabetes
  }
  if (smoker == 1) {
    S <- S + CVD_dict$Smoker
  }
  
  if (sex == 2) {
    CVD_risk_10yr <- 1 - 0.95012 ^ (exp(S - 26.1931))
  } else if (sex == 1) {
    CVD_risk_10yr <- 1 - 0.88936 ^ (exp(S - 23.9802))
  }
  
  CVD_rate <- -log(1 - CVD_risk_10yr) / 10
  CVD_risk <- 1 - exp(-CVD_rate)
  
  HR_RM <- as.numeric(row['RM_HR_CVD'])
  HR_CM <- as.numeric(row['PM_HR_CVD'])
  
  CVD_risk <- CVD_risk * HR_CM * HR_RM
  
  CVD_risk <- CVD_risk / 1.2  # Calibration
  
  return(CVD_risk)
}

# CRC_risk function (from colorectal_cancer.R)
CRC_risk <- function(row, with_diabetes) {
  age <- as.numeric(row['Age'])
  sex <- as.numeric(row['Sex'])
  ethnicity <- as.numeric(row['Ethnicity'])
  pack_years <- as.numeric(row['Smoking pack years'])
  alcohol_drinks <- as.numeric(row['Daily Drinks'])
  BMI <- as.numeric(row['BMI'])
  years_education <- as.numeric(row['Years of Education'])
  previously_taken_aspirin <- as.numeric(row['Taking Aspirin'])
  taking_aspirin <- as.numeric(row['Taking Aspirin'])
  CRC_family_history <- as.numeric(row['CRC family history'])
  taking_multivitamins <- as.numeric(row['Taking Multivitamins'])
  hours_moderate_PA <- as.numeric(row['Hours moderate physical activity'])
  taking_estrogen <- as.numeric(row['Taking Estrogen'])
  previously_taken_estrogen <- as.numeric(row['Taking Estrogen'])
  
  if (with_diabetes) {
    diabetes <- 1
  } else {
    diabetes <- as.numeric(row['Diabetes'])
  }
  
  if (sex == 1) {
    intercept <- -6.6419738
    risk_age <- 0.091669179 * age - 3.7411814e-05 * max(age - 47, 0) ^ 3 + 7.794128e-05 * max(age - 60, 0) ^ 3 - 4.0529466e-05 * max(age - 72, 0) ^ 3
    if (ethnicity == 1 || ethnicity == 2) {
      risk_ethnicity <- -0.13659953
    } else if (ethnicity == 3) {
      risk_ethnicity <- -0.16728044
    } else if (ethnicity == 6) {
      risk_ethnicity <- 0.25353936
    } else {
      risk_ethnicity <- 0
    }
    risk_pack_years <- 0.00022581331 * pack_years + 1.1341047e-05 * max(pack_years, 0) ^ 3 - 1.3522018e-05 * max(pack_years - 6.375, 0) ^ 3 + 2.1809706e-06 * max(pack_years - 39.525, 0) ^ 3
    risk_alcohol <- 0.28379769 * alcohol_drinks - 0.21424251 * max(alcohol_drinks, 0) ^ 3 + 0.22570057 * max(alcohol_drinks - 0.14457189, 0) ^ 3 - 0.011458065 * max(alcohol_drinks - 2.8477722, 0) ^ 3
    risk_BMI <- 0.018020786 * BMI + 9.4715899e-05 * max(BMI - 22.047175, 0) ^ 3 - 0.00015791645 * max(BMI - 25.941735, 0) ^ 3 + 6.3200548e-05 * max(BMI - 31.778341, 0) ^ 3
    risk_years_education <- 0.072052428 * years_education - 0.00060634342 * max(years_education - 7, 0) ^ 3 + 0.0016674444 * max(years_education - 14, 0) ^ 3 - 0.001061101 * max(years_education - 18, 0) ^ 3
    risk_aspirin <- -0.032284161 * previously_taken_aspirin - 0.20960315 * taking_aspirin
    risk_family_history <- 0.24250922 * CRC_family_history
    risk_multivitamins <- -0.19175375 * taking_multivitamins
    if (diabetes == 1) {
      risk_diabetes <- 0.11020556
    } else {
      risk_diabetes <- 0
    }
    risk_physical_activity <- -0.090669913 * hours_moderate_PA + 0.0093816671 * max(hours_moderate_PA - 0.10714286, 0) ^ 3 - 0.011850527 * max(hours_moderate_PA - 0.82142857, 0) ^ 3 + 0.0024688598 * max(hours_moderate_PA - 3.5357143, 0) ^ 3
    Sum <- intercept + risk_age + risk_ethnicity + risk_pack_years + risk_alcohol + risk_BMI + risk_years_education + risk_aspirin + risk_family_history + risk_multivitamins + risk_diabetes + risk_physical_activity
    P <- 1 - 0.9846654 ^ exp(Sum)
  } else if (sex == 2) {
    intercept <- -5.9026635
    risk_age <- 0.090012542 * age - 4.4217156e-05 * max(age - 47, 0) ^ 3 + 9.2119076e-05 * max(age - 60, 0) ^ 3 - 4.7901919e-05 * max(age - 72, 0) ^ 3
    if (ethnicity == 1 || ethnicity == 2) {
      risk_ethnicity <- -0.39669678
    } else if (ethnicity == 3) {
      risk_ethnicity <- -0.34056094
    } else if (ethnicity == 6) {
      risk_ethnicity <- 0.014010715
    } else {
      risk_ethnicity <- 0
    }
    risk_pack_years <- 0.062703176 * pack_years - 0.002446026 * max(pack_years, 0) ^ 3 + 0.003038396 * max(pack_years - 1.25, 0) ^ 3 - 0.00059134632 * max(pack_years - 6.375, 0) ^ 3 - 1.023614e-06 * max(pack_years - 27.5125, 0) ^ 3
    risk_alcohol <- -0.08856241 * alcohol_drinks + 0.62375456 * max(alcohol_drinks, 0) ^ 3 - 0.7191129 * max(alcohol_drinks - 0.10740682, 0) ^ 3 + 0.095358349 * max(alcohol_drinks - 0.8099724, 0) ^ 3
    risk_BMI <- 0.0075233925 * BMI + 6.7918662e-05 * max(BMI - 20.371336, 0) ^ 3 - 0.00011039091 * max(BMI - 25.508027, 0) ^ 3 + 4.2472244e-05 * max(BMI - 33.722266, 0) ^ 3
    risk_years_education <- 0.07443905 * years_education - 0.00062546554 * max(years_education - 7, 0) ^ 3 + 0.0017200302 * max(years_education - 14, 0) ^ 3 - 0.0010945647 * max(years_education - 18, 0) ^ 3
    risk_estrogen <- -0.24500762 * taking_estrogen - 0.044320489 * previously_taken_estrogen
    risk_pain_med <- -0.046383323 * previously_taken_aspirin - 0.236997 * taking_aspirin
    risk_family_history <- 0.31589053 * CRC_family_history
    risk_multivitamins <- -0.1665365 * taking_multivitamins
    if (diabetes == 1) {
      risk_diabetes <- 0.23328937
    } else {
      risk_diabetes <- 0
    }
    Sum <- intercept + risk_age + risk_ethnicity + risk_pack_years + risk_alcohol + risk_BMI + risk_years_education + risk_pain_med + risk_family_history + risk_multivitamins + risk_diabetes
    P <- 1 - 0.9901043 ^ exp(Sum)
  }
  
  CRC_rate <- -log(1 - P) / 10
  CRC_risk <- 1 - exp(-CRC_rate)
  
  RR_red_meat <- as.numeric(row['RM_HR_CRC'])
  RR_processed_meat <- as.numeric(row['PM_HR_CRC'])
  
  CRC_risk <- CRC_risk * RR_processed_meat * RR_red_meat
  
  # Calibration
  if (sex == 1) {
    if (age < 50) {
      CRC_risk <- CRC_risk / 1.56
    } else if (age < 65) {
      CRC_risk <- CRC_risk / 1.5
    } else if (age < 80) {
      CRC_risk <- CRC_risk / 1.79
    } else {
      CRC_risk <- CRC_risk / 1.27
    }
  } else {
    if (age < 50) {
      CRC_risk <- CRC_risk / 1.67
    } else if (age < 65) {
      CRC_risk <- CRC_risk / 1.78
    } else if (age < 80) {
      CRC_risk <- CRC_risk / 1.92
    } else {
      CRC_risk <- CRC_risk / 0.93
    }
  }
  
  return(CRC_risk)
}

# calculate_risks function
calculate_risks <- function(df, mortality_table, seed) {
  set.seed(seed)
  
  df$`Diabetes risk` <- apply(df, 1, function(row) Diabetes_risk_Alva(row, Diabetes_dict))
  df$`CVD risk no diabetes` <- apply(df, 1, function(row) updated_Framingham_CVD_risk(row, with_diabetes = FALSE))
  df$`CVD risk with diabetes` <- apply(df, 1, function(row) updated_Framingham_CVD_risk(row, with_diabetes = TRUE))
  
  df$`CRC risk no diabetes` <- apply(df, 1, function(row) CRC_risk(row, with_diabetes = FALSE))
  df$`CRC risk with diabetes` <- apply(df, 1, function(row) CRC_risk(row, with_diabetes = TRUE))
  
  df$`Diabetes and CVD risk` <- df$`Diabetes risk` * df$`CVD risk no diabetes`
  df$`Diabetes and CRC risk` <- df$`Diabetes risk` * df$`CRC risk no diabetes`
  df$`CVD and CRC risk with diabetes` <- df$`CVD risk with diabetes` * df$`CRC risk with diabetes`
  df$`CVD and CRC risk no diabetes` <- df$`CVD risk no diabetes` * df$`CRC risk no diabetes`
  df$`Diabetes and CVD and CRC risk` <- df$`Diabetes risk` * df$`CVD risk no diabetes` * df$`CRC risk no diabetes`
  
  df$`diabetes mortality risk` <- apply(df, 1, function(row) mortality_prob(row, mortality_table, with_diabetes = TRUE, with_CVD = FALSE))
  df$`CVD mortality risk` <- apply(df, 1, function(row) mortality_prob(row, mortality_table, with_diabetes = FALSE, with_CVD = TRUE))
  df$`diabetes and CVD mortality risk` <- apply(df, 1, function(row) mortality_prob(row, mortality_table, with_diabetes = TRUE, with_CVD = TRUE))
  df$`CRC mortality risk` <- apply(df, 1, CRC_mortality_prob)
  df$`Healthy mortality risk` <- apply(df, 1, function(row) mortality_prob(row, mortality_table, with_diabetes = FALSE, with_CVD = FALSE))
  
  return(df)
}

# mortality_prob function (from mortalities.R)
mortality_prob <- function(row, mortality_table, with_diabetes, with_CVD) {
  age <- as.numeric(row['Age'])
  sex <- as.numeric(row['Sex'])
  
  if (age < 100) {
    mortality_prob <- mortality_table[[as.character(sex)]][[as.character(floor(age))]]
  } else {
    mortality_prob <- 0.99
  }
  
  if (with_CVD) {
    RR_CVD <- rnorm(1, mean = 2.0, sd = 0.05)
    mortality_prob <- mortality_prob * RR_CVD
  }
  
  if (with_diabetes) {
    if (age < 55) {
      RR_Diabetes <- rnorm(1, mean = 2.35, sd = 0.085)
    } else if (age < 65) {
      RR_Diabetes <- rnorm(1, mean = 1.79, sd = 0.03)
    } else if (age < 75) {
      RR_Diabetes <- rnorm(1, mean = 1.46, sd = 0.015)
    } else {
      RR_Diabetes <- rnorm(1, mean = 1.19, sd = 0.005)
    }
    mortality_prob <- mortality_prob * RR_Diabetes
  }
  
  return(mortality_prob)
}

# CRC_mortality_prob function (from mortalities.R)
CRC_mortality_prob <- function(row) {
  age <- as.numeric(row['Age'])
  sex <- as.numeric(row['Sex'])
  
  if (sex == 1) {
    if (age < 50) {
      mortality_probability <- 0.0236
    } else if (age < 65) {
      mortality_probability <- 0.247
    } else if (age < 80) {
      mortality_probability <- 0.0226
    } else {
      mortality_probability <- 0.0796
    }
  } else if (sex == 2) {
    if (age < 50) {
      mortality_probability <- 0.0426
    } else if (age < 65) {
      mortality_probability <- 0.0118
    } else if (age < 80) {
      mortality_probability <- 0.0347
    } else {
      mortality_probability <- 0.0308
    }
  }
  
  return(mortality_probability)
}

# mortality_table (from mortalities.R)
# This is a nested list in R
mortality_table <- list(
  `1` = list(
    `18` = 0.000811, `19` = 0.000945, `20` = 0.001082, `21` = 0.001214, `22` = 0.001327, `23` = 0.001413,
    `24` = 0.001476, `25` = 0.001531, `26` = 0.001584, `27` = 0.001633, `28` = 0.001681, `29` = 0.001739,
    `30` = 0.001779, `31` = 0.001829, `32` = 0.001888, `33` = 0.001957, `34` = 0.002032, `35` = 0.002119,
    `36` = 0.002209, `37` = 0.002286, `38` = 0.002346, `39` = 0.002401, `40` = 0.002468, `41` = 0.002565,
    `42` = 0.0027, `43` = 0.002876, `44` = 0.003084, `45` = 0.003318, `46` = 0.003572, `47` = 0.00385,
    `48` = 0.004161, `49` = 0.004515, `50` = 0.004895, `51` = 0.005362, `52` = 0.005835, `53` = 0.006438,
    `54` = 0.007098, `55` = 0.007765, `56` = 0.008432, `57` = 0.009126, `58` = 0.00987, `59` = 0.01067,
    `60` = 0.011534, `61` = 0.012431, `62` = 0.013332, `63` = 0.014219, `64` = 0.015117, `65` = 0.016078,
    `66` = 0.017216, `67` = 0.018401, `68` = 0.019666, `69` = 0.021099, `70` = 0.022544, `71` = 0.024099,
    `72` = 0.026447, `73` = 0.028617, `74` = 0.03139, `75` = 0.034322, `76` = 0.03797, `77` = 0.041944,
    `78` = 0.045881, `79` = 0.050573, `80` = 0.055675, `81` = 0.061704, `82` = 0.068389, `83` = 0.075732,
    `84` = 0.085241, `85` = 0.094489, `86` = 0.104787, `87` = 0.117465, `88` = 0.131319, `89` = 0.146372,
    `90` = 0.162625, `91` = 0.18006, `92` = 0.198626, `93` = 0.218248, `94` = 0.23882, `95` = 0.260206,
    `96` = 0.282243, `97` = 0.304747, `98` = 0.327517, `99` = 0.350342
  ),
  `2` = list(
    `18` = 0.000329, `19` = 0.000365, `20` = 0.000402, `21` = 0.000441, `22` = 0.000481, `23` = 0.000521,
    `24` = 0.00056, `25` = 0.000598, `26` = 0.000635, `27` = 0.000675, `28` = 0.000718, `29` = 0.000765,
    `30` = 0.000818, `31` = 0.000872, `32` = 0.000928, `33` = 0.000983, `34` = 0.001037, `35` = 0.001097,
    `36` = 0.00116, `37` = 0.001274, `38` = 0.001274, `39` = 0.00133, `40` = 0.001396, `41` = 0.00148,
    `42` = 0.001584, `43` = 0.001709, `44` = 0.001849, `45` = 0.002002, `46` = 0.002167, `47` = 0.002348,
    `48` = 0.002551, `49` = 0.002781, `50` = 0.003028, `51` = 0.003299, `52` = 0.003615, `53` = 0.003974,
    `54` = 0.004357, `55` = 0.004746, `56` = 0.005135, `57` = 0.00553, `58` = 0.00594, `59` = 0.006376,
    `60` = 0.006849, `61` = 0.007354, `62` = 0.007879, `63` = 0.008425, `64` = 0.009008, `65` = 0.009638,
    `66` = 0.010386, `67` = 0.011235, `68` = 0.012237, `69` = 0.013393, `70` = 0.014731, `71` = 0.01608,
    `72` = 0.017952, `73` = 0.019637, `74` = 0.021744, `75` = 0.023929, `76` = 0.026629, `77` = 0.029547,
    `78` = 0.032857, `79` = 0.036519, `80` = 0.040589, `81` = 0.045639, `82` = 0.051261, `83` = 0.057423,
    `84` = 0.065035, `85` = 0.072862, `86` = 0.080121, `87` = 0.09065, `88` = 0.102322, `89` = 0.115196,
    `90` = 0.129319, `91` = 0.144719, `92` = 0.161402, `93` = 0.179347, `94` = 0.198504, `95` = 0.218788,
    `96` = 0.240082, `97` = 0.262235, `98` = 0.285066, `99` = 0.308369
  )
)

# Assume all previous code and functions are already defined in population.R

# Function to update mortalities
update_mortalities <- function(df, year, pre_new_cases) {
  if (pre_new_cases) {
    df[[paste0('Healthy mortalities year ', year, ' pre')]] <- apply(df, 1, expected_healthy_mortalities)
    df[[paste0('diabetes mortalities year ', year, ' pre')]] <- apply(df, 1, expected_diabetes_mortalities)
    df[[paste0('CVD mortalities year ', year, ' pre')]] <- apply(df, 1, expected_CVD_mortalities)
    df[[paste0('CRC mortalities year ', year, ' pre')]] <- apply(df, 1, expected_CRC_mortality)
    df[[paste0('diabetes and CVD mortalities year ', year, ' pre')]] <- apply(df, 1, expected_diabetes_CVD_mortalities)
    df[[paste0('diabetes and CRC mortalities year ', year, ' pre')]] <- apply(df, 1, expected_diabetes_CRC_mortalities)
    df[[paste0('CVD and CRC mortalities year ', year, ' pre')]] <- apply(df, 1, expected_CVD_CRC_mortalities)
    df[[paste0('diabetes and CVD and CRC mortalities year ', year, ' pre')]] <- apply(df, 1, expected_diabetes_CVD_CRC_mortalities)
    
    df[[paste0('Total mortalities year ', year, ' pre')]] <- df[[paste0('Healthy mortalities year ', year, ' pre')]] +
      df[[paste0('diabetes mortalities year ', year, ' pre')]] +
      df[[paste0('CVD mortalities year ', year, ' pre')]] +
      df[[paste0('CRC mortalities year ', year, ' pre')]] -
      df[[paste0('diabetes and CVD mortalities year ', year, ' pre')]] -
      df[[paste0('diabetes and CRC mortalities year ', year, ' pre')]] -
      df[[paste0('CVD and CRC mortalities year ', year, ' pre')]] +
      df[[paste0('diabetes and CVD and CRC mortalities year ', year, ' pre')]]
  } else {
    df[[paste0('Healthy mortalities year ', year, ' post')]] <- apply(df, 1, expected_healthy_mortalities)
    df[[paste0('diabetes mortalities year ', year, ' post')]] <- apply(df, 1, expected_diabetes_mortalities)
    df[[paste0('CVD mortalities year ', year, ' post')]] <- apply(df, 1, expected_CVD_mortalities)
    df[[paste0('CRC mortalities year ', year, ' post')]] <- apply(df, 1, expected_CRC_mortality)
    df[[paste0('diabetes and CVD mortalities year ', year, ' post')]] <- apply(df, 1, expected_diabetes_CVD_mortalities)
    df[[paste0('diabetes and CRC mortalities year ', year, ' post')]] <- apply(df, 1, expected_diabetes_CRC_mortalities)
    df[[paste0('CVD and CRC mortalities year ', year, ' post')]] <- apply(df, 1, expected_CVD_CRC_mortalities)
    df[[paste0('diabetes and CVD and CRC mortalities year ', year, ' post')]] <- apply(df, 1, expected_diabetes_CVD_CRC_mortalities)
    
    df[[paste0('Total mortalities year ', year, ' post')]] <- df[[paste0('Healthy mortalities year ', year, ' post')]] +
      df[[paste0('diabetes mortalities year ', year, ' post')]] +
      df[[paste0('CVD mortalities year ', year, ' post')]] +
      df[[paste0('CRC mortalities year ', year, ' post')]] -
      df[[paste0('diabetes and CVD mortalities year ', year, ' post')]] -
      df[[paste0('diabetes and CRC mortalities year ', year, ' post')]] -
      df[[paste0('CVD and CRC mortalities year ', year, ' post')]] +
      df[[paste0('diabetes and CVD and CRC mortalities year ', year, ' post')]]
  }
  return(df)
}

# Function to update cases
update_cases <- function(df, year) {
  df[[paste0('New diabetes cases year ', year)]] <- apply(df, 1, expected_new_diabetes_cases)
  df[[paste0('New CVD cases year ', year)]] <- apply(df, 1, expected_new_CVD_cases)
  df[[paste0('New CRC cases year ', year)]] <- apply(df, 1, expected_new_CRC_cases)
  df[[paste0('New diabetes and CVD cases year ', year)]] <- apply(df, 1, expected_new_diabetes_CVD_cases)
  df[[paste0('New diabetes and CRC cases year ', year)]] <- apply(df, 1, expected_new_diabetes_CRC_cases)
  df[[paste0('New CVD and CRC cases year ', year)]] <- apply(df, 1, expected_new_CVD_CRC_cases)
  df[[paste0('New diabetes and CVD and CRC cases year ', year)]] <- apply(df, 1, expected_new_diabetes_CVD_CRC_cases)
  
  # Update cumulative cases
  df[['New diabetes cases']] <- df[['New diabetes cases']] + df[[paste0('New diabetes cases year ', year)]]
  df[['New CVD cases']] <- df[['New CVD cases']] + df[[paste0('New CVD cases year ', year)]]
  df[['New CRC cases']] <- df[['New CRC cases']] + df[[paste0('New CRC cases year ', year)]]
  df[['New diabetes and CVD cases']] <- df[['New diabetes and CVD cases']] + df[[paste0('New diabetes and CVD cases year ', year)]]
  df[['New diabetes and CRC cases']] <- df[['New diabetes and CRC cases']] + df[[paste0('New diabetes and CRC cases year ', year)]]
  df[['New CVD and CRC cases']] <- df[['New CVD and CRC cases']] + df[[paste0('New CVD and CRC cases year ', year)]]
  df[['New diabetes and CVD and CRC cases']] <- df[['New diabetes and CVD and CRC cases']] + df[[paste0('New diabetes and CVD and CRC cases year ', year)]]
  
  # Update healthy population
  df[['healthy']] <- apply(df, 1, healthy_pop)
  
  return(df)
}

# Function to apply new mortalities
apply_new_mortalities <- function(df, year) {
  df[['Sample Weight']] <- apply(df, 1, function(row) update_sample_weight(row, year))
  df[['New diabetes cases']] <- apply(df, 1, function(row) update_new_diabetes_cases(row, year))
  df[['New CVD cases']] <- apply(df, 1, function(row) update_new_CVD_cases(row, year))
  df[['New CRC cases']] <- apply(df, 1, function(row) update_new_CRC_cases(row, year))
  df[['New diabetes and CVD cases']] <- apply(df, 1, function(row) update_new_diabetes_CVD_cases(row, year))
  df[['New diabetes and CRC cases']] <- apply(df, 1, function(row) update_new_diabetes_CRC_cases(row, year))
  df[['New CVD and CRC cases']] <- apply(df, 1, function(row) update_new_CVD_CRC_cases(row, year))
  df[['New diabetes and CVD and CRC cases']] <- apply(df, 1, function(row) update_new_diabetes_CVD_CRC_cases(row, year))
  
  return(df)
}

# Function to update sample weight
update_sample_weight <- function(row, year) {
  SW <- as.numeric(row[['Sample Weight']])
  total_mortalities <- as.numeric(row[[paste0('Total mortalities year ', year, ' post')]])
  SW <- SW - total_mortalities
  return(SW)
}

# Function to update new diabetes cases
update_new_diabetes_cases <- function(row, year) {
  new_cases <- as.numeric(row[['New diabetes cases']])
  
  if (as.numeric(row[['Diabetes']]) != 1) {
    mortalities <- as.numeric(row[[paste0('diabetes mortalities year ', year, ' post')]])
    new_cases <- new_cases - mortalities
    if (new_cases < 0) {
      new_cases <- 0
    }
  }
  return(new_cases)
}

# Function to update new CVD cases
update_new_CVD_cases <- function(row, year) {
  new_cases <- as.numeric(row[['New CVD cases']])
  
  if (as.numeric(row[['CVD']]) != 1) {
    mortalities <- as.numeric(row[[paste0('CVD mortalities year ', year, ' post')]])
    new_cases <- new_cases - mortalities
    if (new_cases < 0) {
      new_cases <- 0
    }
  }
  return(new_cases)
}

# Function to update new CRC cases
update_new_CRC_cases <- function(row, year) {
  new_cases <- as.numeric(row[['New CRC cases']])
  
  if (as.numeric(row[['CRC']]) != 1) {
    mortalities <- as.numeric(row[[paste0('CRC mortalities year ', year, ' post')]])
    new_cases <- new_cases - mortalities
    if (new_cases < 0) {
      new_cases <- 0
    }
  }
  return(new_cases)
}

# Function to update new diabetes and CVD cases
update_new_diabetes_CVD_cases <- function(row, year) {
  new_cases <- as.numeric(row[['New diabetes and CVD cases']])
  
  if (as.numeric(row[['Diabetes']]) != 1 && as.numeric(row[['CVD']]) != 1) {
    mortalities <- as.numeric(row[[paste0('diabetes and CVD mortalities year ', year, ' post')]])
    new_cases <- new_cases - mortalities
    if (new_cases < 0) {
      new_cases <- 0
    }
  }
  return(new_cases)
}

# Function to update new diabetes and CRC cases
update_new_diabetes_CRC_cases <- function(row, year) {
  new_cases <- as.numeric(row[['New diabetes and CRC cases']])
  
  if (as.numeric(row[['Diabetes']]) != 1 && as.numeric(row[['CRC']]) != 1) {
    mortalities <- as.numeric(row[[paste0('diabetes and CRC mortalities year ', year, ' post')]])
    new_cases <- new_cases - mortalities
    if (new_cases < 0) {
      new_cases <- 0
    }
  }
  return(new_cases)
}

# Function to update new CVD and CRC cases
update_new_CVD_CRC_cases <- function(row, year) {
  new_cases <- as.numeric(row[['New CVD and CRC cases']])
  
  if (as.numeric(row[['CVD']]) != 1 && as.numeric(row[['CRC']]) != 1) {
    mortalities <- as.numeric(row[[paste0('CVD and CRC mortalities year ', year, ' post')]])
    new_cases <- new_cases - mortalities
    if (new_cases < 0) {
      new_cases <- 0
    }
  }
  return(new_cases)
}

# Function to update new diabetes, CVD, and CRC cases
update_new_diabetes_CVD_CRC_cases <- function(row, year) {
  new_cases <- as.numeric(row[['New diabetes and CVD and CRC cases']])
  
  if (as.numeric(row[['Diabetes']]) != 1 && as.numeric(row[['CVD']]) != 1 && as.numeric(row[['CRC']]) != 1) {
    mortalities <- as.numeric(row[[paste0('diabetes and CVD and CRC mortalities year ', year, ' post')]])
    new_cases <- new_cases - mortalities
    if (new_cases < 0) {
      new_cases <- 0
    }
  }
  return(new_cases)
}

# Function to calculate expected healthy mortalities
expected_healthy_mortalities <- function(row) {
  if (as.numeric(row[['Diabetes']]) == 1 || as.numeric(row[['CVD']]) == 1 || as.numeric(row[['CRC']]) == 1) {
    return(0)
  } else {
    return(as.numeric(row[['healthy']]) * as.numeric(row[['Healthy mortality risk']]))
  }
}

# Function to calculate expected diabetes mortalities
expected_diabetes_mortalities <- function(row) {
  exp_mortalities <- 0
  SW <- as.numeric(row[['Sample Weight']])
  Diabetes <- as.numeric(row[['Diabetes']])
  CVD <- as.numeric(row[['CVD']])
  CRC <- as.numeric(row[['CRC']])
  
  # Extract necessary values
  diabetes_mortality_risk <- as.numeric(row[['diabetes mortality risk']])
  diabetes_and_CVD_mortality_risk <- as.numeric(row[['diabetes and CVD mortality risk']])
  CRC_mortality_risk <- as.numeric(row[['CRC mortality risk']])
  
  New_CVD_cases <- as.numeric(row[['New CVD cases']])
  New_CRC_cases <- as.numeric(row[['New CRC cases']])
  New_diabetes_cases <- as.numeric(row[['New diabetes cases']])
  New_diabetes_and_CVD_cases <- as.numeric(row[['New diabetes and CVD cases']])
  New_diabetes_and_CRC_cases <- as.numeric(row[['New diabetes and CRC cases']])
  New_CVD_and_CRC_cases <- as.numeric(row[['New CVD and CRC cases']])
  New_diabetes_and_CVD_and_CRC_cases <- as.numeric(row[['New diabetes and CVD and CRC cases']])
  
  if (CRC == 1 && Diabetes == 1) {
    exp_mortalities <- SW * CRC_mortality_risk
  } else if (CRC == 1 && Diabetes != 1) {
    exp_mortalities <- as.numeric(row[['New diabetes and CRC cases']]) * CRC_mortality_risk
  } else if (Diabetes == 1 && CVD != 1) {
    SW_diabetes_only <- SW - New_CRC_cases - New_CVD_cases + New_CVD_and_CRC_cases
    exp_mortalities <- SW_diabetes_only * diabetes_mortality_risk
    exp_mortalities <- exp_mortalities + (New_CVD_cases - New_CVD_and_CRC_cases) * diabetes_and_CVD_mortality_risk
    exp_mortalities <- exp_mortalities + New_CRC_cases * CRC_mortality_risk
  } else if (CVD == 1 && Diabetes != 1) {
    exp_mortalities <- (New_diabetes_cases - New_diabetes_and_CRC_cases) * diabetes_and_CVD_mortality_risk
    exp_mortalities <- exp_mortalities + New_diabetes_and_CRC_cases * CRC_mortality_risk
  } else if (Diabetes == 1 && CVD == 1) {
    exp_mortalities <- (SW - New_CRC_cases) * diabetes_and_CVD_mortality_risk
    exp_mortalities <- exp_mortalities + New_CRC_cases * CRC_mortality_risk
  } else if (Diabetes != 1 && CVD != 1) {
    SW_diabetes_only <- New_diabetes_cases - New_diabetes_and_CRC_cases - New_diabetes_and_CVD_cases + New_diabetes_and_CVD_and_CRC_cases
    exp_mortalities <- SW_diabetes_only * diabetes_mortality_risk
    exp_mortalities <- exp_mortalities + (New_diabetes_and_CVD_cases - New_diabetes_and_CVD_and_CRC_cases) * diabetes_and_CVD_mortality_risk
    exp_mortalities <- exp_mortalities + New_diabetes_and_CRC_cases * CRC_mortality_risk
  }
  return(exp_mortalities)
}

# Function to calculate expected CVD mortalities
expected_CVD_mortalities <- function(row) {
  exp_mortalities <- 0
  SW <- as.numeric(row[['Sample Weight']])
  Diabetes <- as.numeric(row[['Diabetes']])
  CVD <- as.numeric(row[['CVD']])
  CRC <- as.numeric(row[['CRC']])
  
  # Extract necessary values
  CVD_mortality_risk <- as.numeric(row[['CVD mortality risk']])
  diabetes_and_CVD_mortality_risk <- as.numeric(row[['diabetes and CVD mortality risk']])
  CRC_mortality_risk <- as.numeric(row[['CRC mortality risk']])
  
  New_diabetes_cases <- as.numeric(row[['New diabetes cases']])
  New_CRC_cases <- as.numeric(row[['New CRC cases']])
  New_CVD_cases <- as.numeric(row[['New CVD cases']])
  New_diabetes_and_CRC_cases <- as.numeric(row[['New diabetes and CRC cases']])
  New_CVD_and_CRC_cases <- as.numeric(row[['New CVD and CRC cases']])
  New_diabetes_and_CVD_cases <- as.numeric(row[['New diabetes and CVD cases']])
  New_diabetes_and_CVD_and_CRC_cases <- as.numeric(row[['New diabetes and CVD and CRC cases']])
  
  if (CRC == 1 && CVD == 1) {
    exp_mortalities <- SW * CRC_mortality_risk
  } else if (CRC == 1 && CVD != 1) {
    exp_mortalities <- New_CVD_and_CRC_cases * CRC_mortality_risk
  } else if (CVD == 1 && Diabetes != 1) {
    SW_CVD_only <- SW - New_CRC_cases - New_diabetes_cases + New_diabetes_and_CRC_cases
    exp_mortalities <- SW_CVD_only * CVD_mortality_risk
    exp_mortalities <- exp_mortalities + (New_diabetes_cases - New_diabetes_and_CRC_cases) * diabetes_and_CVD_mortality_risk
    exp_mortalities <- exp_mortalities + New_CRC_cases * CRC_mortality_risk
  } else if (Diabetes == 1 && CVD != 1) {
    exp_mortalities <- (New_CVD_cases - New_CVD_and_CRC_cases) * diabetes_and_CVD_mortality_risk
    exp_mortalities <- exp_mortalities + New_CVD_and_CRC_cases * CRC_mortality_risk
  } else if (Diabetes == 1 && CVD == 1) {
    exp_mortalities <- (SW - New_CRC_cases) * diabetes_and_CVD_mortality_risk
    exp_mortalities <- exp_mortalities + New_CRC_cases * CRC_mortality_risk
  } else if (Diabetes != 1 && CVD != 1 && CRC != 1) {
    SW_CVD_only <- New_CVD_cases - New_CVD_and_CRC_cases - New_diabetes_and_CVD_cases + New_diabetes_and_CVD_and_CRC_cases
    exp_mortalities <- SW_CVD_only * CVD_mortality_risk
    exp_mortalities <- exp_mortalities + (New_diabetes_and_CVD_cases - New_diabetes_and_CVD_and_CRC_cases) * diabetes_and_CVD_mortality_risk
    exp_mortalities <- exp_mortalities + New_CVD_and_CRC_cases * CRC_mortality_risk
  }
  return(exp_mortalities)
}

# Function to calculate expected CRC mortalities
expected_CRC_mortality <- function(row) {
  CRC <- as.numeric(row[['CRC']])
  SW <- as.numeric(row[['Sample Weight']])
  CRC_mortality_risk <- as.numeric(row[['CRC mortality risk']])
  
  if (CRC == 1) {
    exp_mortalities <- SW * CRC_mortality_risk
  } else {
    exp_mortalities <- as.numeric(row[['New CRC cases']]) * CRC_mortality_risk
  }
  return(exp_mortalities)
}

# Function to calculate expected diabetes and CVD mortalities
expected_diabetes_CVD_mortalities <- function(row) {
  exp_mortalities <- 0
  SW <- as.numeric(row[['Sample Weight']])
  Diabetes <- as.numeric(row[['Diabetes']])
  CVD <- as.numeric(row[['CVD']])
  CRC <- as.numeric(row[['CRC']])
  
  # Extract necessary values
  diabetes_and_CVD_mortality_risk <- as.numeric(row[['diabetes and CVD mortality risk']])
  CRC_mortality_risk <- as.numeric(row[['CRC mortality risk']])
  
  New_diabetes_and_CVD_cases <- as.numeric(row[['New diabetes and CVD cases']])
  New_diabetes_and_CVD_and_CRC_cases <- as.numeric(row[['New diabetes and CVD and CRC cases']])
  New_CRC_cases <- as.numeric(row[['New CRC cases']])
  
  if (CRC == 1 && Diabetes == 1 && CVD == 1) {
    exp_mortalities <- SW * CRC_mortality_risk
  } else if (CRC == 1) {
    exp_mortalities <- New_diabetes_and_CVD_cases * CRC_mortality_risk
  } else if (Diabetes == 1 && CVD == 1) {
    exp_mortalities <- (SW - New_CRC_cases) * diabetes_and_CVD_mortality_risk
    exp_mortalities <- exp_mortalities + New_CRC_cases * CRC_mortality_risk
  } else {
    exp_mortalities <- (New_diabetes_and_CVD_cases - New_diabetes_and_CVD_and_CRC_cases) * diabetes_and_CVD_mortality_risk
    exp_mortalities <- exp_mortalities + New_diabetes_and_CVD_and_CRC_cases * CRC_mortality_risk
  }
  return(exp_mortalities)
}

# Function to calculate expected diabetes and CRC mortalities
expected_diabetes_CRC_mortalities <- function(row) {
  Diabetes <- as.numeric(row[['Diabetes']])
  CRC <- as.numeric(row[['CRC']])
  SW <- as.numeric(row[['Sample Weight']])
  CRC_mortality_risk <- as.numeric(row[['CRC mortality risk']])
  
  if (CRC == 1 && Diabetes == 1) {
    exp_mortalities <- SW * CRC_mortality_risk
  } else {
    exp_mortalities <- as.numeric(row[['New diabetes and CRC cases']]) * CRC_mortality_risk
  }
  return(exp_mortalities)
}

# Function to calculate expected CVD and CRC mortalities
expected_CVD_CRC_mortalities <- function(row) {
  CVD <- as.numeric(row[['CVD']])
  CRC <- as.numeric(row[['CRC']])
  SW <- as.numeric(row[['Sample Weight']])
  CRC_mortality_risk <- as.numeric(row[['CRC mortality risk']])
  
  if (CRC == 1 && CVD == 1) {
    exp_mortalities <- SW * CRC_mortality_risk
  } else {
    exp_mortalities <- as.numeric(row[['New CVD and CRC cases']]) * CRC_mortality_risk
  }
  return(exp_mortalities)
}

# Function to calculate expected diabetes, CVD, and CRC mortalities
expected_diabetes_CVD_CRC_mortalities <- function(row) {
  Diabetes <- as.numeric(row[['Diabetes']])
  CVD <- as.numeric(row[['CVD']])
  CRC <- as.numeric(row[['CRC']])
  SW <- as.numeric(row[['Sample Weight']])
  CRC_mortality_risk <- as.numeric(row[['CRC mortality risk']])
  
  if (Diabetes == 1 && CVD == 1 && CRC == 1) {
    exp_mortalities <- SW * CRC_mortality_risk
  } else {
    exp_mortalities <- as.numeric(row[['New diabetes and CVD and CRC cases']]) * CRC_mortality_risk
  }
  return(exp_mortalities)
}

# Function to calculate expected new diabetes cases
expected_new_diabetes_cases <- function(row) {
  if (as.numeric(row[['Diabetes']]) != 1) {
    exp_cases <- as.numeric(row[['Diabetes risk']]) * (as.numeric(row[['Sample Weight']]) - as.numeric(row[['New diabetes cases']]))
  } else {
    exp_cases <- 0
  }
  return(exp_cases)
}

# Function to calculate expected new CVD cases
expected_new_CVD_cases <- function(row) {
  exp_cases <- 0
  Diabetes <- as.numeric(row[['Diabetes']])
  CVD <- as.numeric(row[['CVD']])
  SW <- as.numeric(row[['Sample Weight']])
  New_diabetes_cases <- as.numeric(row[['New diabetes cases']])
  New_CVD_cases <- as.numeric(row[['New CVD cases']])
  New_diabetes_and_CVD_cases <- as.numeric(row[['New diabetes and CVD cases']])
  
  if (CVD != 1) {
    if (Diabetes != 1) {
      exp_cases <- exp_cases + as.numeric(row[['CVD risk no diabetes']]) * (SW - New_diabetes_cases - New_CVD_cases + New_diabetes_and_CVD_cases)
      exp_cases <- exp_cases + as.numeric(row[['CVD risk with diabetes']]) * (New_diabetes_cases - New_diabetes_and_CVD_cases)
    } else {
      exp_cases <- exp_cases + as.numeric(row[['CVD risk with diabetes']]) * (SW - New_CVD_cases)
    }
  }
  return(exp_cases)
}

# Function to calculate expected new CRC cases
expected_new_CRC_cases <- function(row) {
  exp_cases <- 0
  Diabetes <- as.numeric(row[['Diabetes']])
  CRC <- as.numeric(row[['CRC']])
  SW <- as.numeric(row[['Sample Weight']])
  New_diabetes_cases <- as.numeric(row[['New diabetes cases']])
  New_CRC_cases <- as.numeric(row[['New CRC cases']])
  New_diabetes_and_CRC_cases <- as.numeric(row[['New diabetes and CRC cases']])
  
  if (CRC != 1) {
    if (Diabetes != 1) {
      exp_cases <- exp_cases + as.numeric(row[['CRC risk no diabetes']]) * (SW - New_diabetes_cases - New_CRC_cases + New_diabetes_and_CRC_cases)
      exp_cases <- exp_cases + as.numeric(row[['CRC risk with diabetes']]) * (New_diabetes_cases - New_diabetes_and_CRC_cases)
    } else {
      exp_cases <- exp_cases + as.numeric(row[['CRC risk with diabetes']]) * (SW - New_CRC_cases)
    }
  }
  return(exp_cases)
}



# Function to calculate expected new diabetes and CVD cases
expected_new_diabetes_CVD_cases <- function(row) {
  exp_cases <- 0
  
  Diabetes <- as.numeric(row[['Diabetes']])
  CVD <- as.numeric(row[['CVD']])
  SW <- as.numeric(row[['Sample Weight']])
  
  New_diabetes_cases <- as.numeric(row[['New diabetes cases']])
  New_CVD_cases <- as.numeric(row[['New CVD cases']])
  New_diabetes_and_CVD_cases <- as.numeric(row[['New diabetes and CVD cases']])
  
  Diabetes_risk <- as.numeric(row[['Diabetes risk']])
  CVD_risk_with_diabetes <- as.numeric(row[['CVD risk with diabetes']])
  Diabetes_and_CVD_risk <- as.numeric(row[['Diabetes and CVD risk']])
  
  if (CVD != 1 && Diabetes != 1) {
    exp_cases <- exp_cases +
      Diabetes_and_CVD_risk * 
      (SW - New_diabetes_cases - New_CVD_cases + New_diabetes_and_CVD_cases)
    
    exp_cases <- exp_cases +
      Diabetes_risk * 
      (New_CVD_cases - New_diabetes_and_CVD_cases)
    
    exp_cases <- exp_cases +
      CVD_risk_with_diabetes * 
      (New_diabetes_cases - New_diabetes_and_CVD_cases)
    
  } else if (CVD == 1 && Diabetes != 1) {
    exp_cases <- exp_cases +
      Diabetes_risk * 
      (SW - New_diabetes_cases)
    
  } else if (CVD != 1 && Diabetes == 1) {
    exp_cases <- exp_cases +
      CVD_risk_with_diabetes * 
      (SW - New_CVD_cases)
    
  } else {
    # Both CVD and Diabetes == 1; no new cases expected
  }
  
  return(exp_cases)
}

# Function to calculate expected new diabetes and CRC cases
expected_new_diabetes_CRC_cases <- function(row) {
  exp_cases <- 0
  
  Diabetes <- as.numeric(row[['Diabetes']])
  CRC <- as.numeric(row[['CRC']])
  SW <- as.numeric(row[['Sample Weight']])
  
  New_diabetes_cases <- as.numeric(row[['New diabetes cases']])
  New_CRC_cases <- as.numeric(row[['New CRC cases']])
  New_diabetes_and_CRC_cases <- as.numeric(row[['New diabetes and CRC cases']])
  
  Diabetes_risk <- as.numeric(row[['Diabetes risk']])
  CRC_risk_with_diabetes <- as.numeric(row[['CRC risk with diabetes']])
  Diabetes_and_CRC_risk <- as.numeric(row[['Diabetes and CRC risk']])
  
  if (CRC != 1 && Diabetes != 1) {
    exp_cases <- exp_cases +
      Diabetes_and_CRC_risk * 
      (SW - New_diabetes_cases - New_CRC_cases + New_diabetes_and_CRC_cases)
    
    exp_cases <- exp_cases +
      Diabetes_risk * 
      (New_CRC_cases - New_diabetes_and_CRC_cases)
    
    exp_cases <- exp_cases +
      CRC_risk_with_diabetes * 
      (New_diabetes_cases - New_diabetes_and_CRC_cases)
    
  } else if (CRC == 1 && Diabetes != 1) {
    exp_cases <- exp_cases +
      Diabetes_risk * 
      (SW - New_diabetes_cases)
    
  } else if (CRC != 1 && Diabetes == 1) {
    exp_cases <- exp_cases +
      CRC_risk_with_diabetes * 
      (SW - New_CRC_cases)
    
  } else {
    # Both CRC and Diabetes == 1; no new cases expected
  }
  
  return(exp_cases)
}

# Function to calculate expected new CVD and CRC cases
expected_new_CVD_CRC_cases <- function(row) {
  exp_cases <- 0
  
  Diabetes <- as.numeric(row[['Diabetes']])
  CVD <- as.numeric(row[['CVD']])
  CRC <- as.numeric(row[['CRC']])
  SW <- as.numeric(row[['Sample Weight']])
  
  New_diabetes_cases <- as.numeric(row[['New diabetes cases']])
  New_CVD_cases <- as.numeric(row[['New CVD cases']])
  New_CRC_cases <- as.numeric(row[['New CRC cases']])
  New_CVD_and_CRC_cases <- as.numeric(row[['New CVD and CRC cases']])
  New_diabetes_and_CVD_cases <- as.numeric(row[['New diabetes and CVD cases']])
  New_diabetes_and_CRC_cases <- as.numeric(row[['New diabetes and CRC cases']])
  New_diabetes_and_CVD_and_CRC_cases <- as.numeric(row[['New diabetes and CVD and CRC cases']])
  
  CVD_risk_no_diabetes <- as.numeric(row[['CVD risk no diabetes']])
  CRC_risk_no_diabetes <- as.numeric(row[['CRC risk no diabetes']])
  CVD_risk_with_diabetes <- as.numeric(row[['CVD risk with diabetes']])
  CRC_risk_with_diabetes <- as.numeric(row[['CRC risk with diabetes']])
  CVD_and_CRC_risk_no_diabetes <- as.numeric(row[['CVD and CRC risk no diabetes']])
  CVD_and_CRC_risk_with_diabetes <- as.numeric(row[['CVD and CRC risk with diabetes']])
  
  if (Diabetes == 1) {
    if (CRC != 1 && CVD != 1) {
      exp_cases <- exp_cases +
        CVD_and_CRC_risk_with_diabetes * 
        (SW - New_CVD_cases - New_CRC_cases + New_CVD_and_CRC_cases)
      
      exp_cases <- exp_cases +
        CVD_risk_with_diabetes * 
        (New_CRC_cases - New_CVD_and_CRC_cases)
      
      exp_cases <- exp_cases +
        CRC_risk_with_diabetes * 
        (New_CVD_cases - New_CVD_and_CRC_cases)
      
    } else if (CRC == 1 && CVD != 1) {
      exp_cases <- exp_cases +
        CVD_risk_with_diabetes * 
        (SW - New_CVD_cases)
      
    } else if (CRC != 1 && CVD == 1) {
      exp_cases <- exp_cases +
        CRC_risk_with_diabetes * 
        (SW - New_CRC_cases)
      
    } else {
      # Both CRC and CVD == 1; no new cases expected
    }
    
  } else {
    if (CRC != 1 && CVD != 1) {
      exp_cases <- exp_cases +
        CVD_and_CRC_risk_no_diabetes * 
        as.numeric(row[['healthy']])
      
      exp_cases <- exp_cases +
        CVD_and_CRC_risk_with_diabetes * 
        (New_diabetes_cases - New_diabetes_and_CVD_cases - New_diabetes_and_CRC_cases + New_diabetes_and_CVD_and_CRC_cases)
      
      exp_cases <- exp_cases +
        CRC_risk_no_diabetes * 
        (New_CVD_cases - New_CVD_and_CRC_cases - New_diabetes_and_CVD_cases + New_diabetes_and_CVD_and_CRC_cases)
      
      exp_cases <- exp_cases +
        CRC_risk_with_diabetes * 
        (New_diabetes_and_CVD_cases - New_diabetes_and_CVD_and_CRC_cases)
      
      exp_cases <- exp_cases +
        CVD_risk_no_diabetes * 
        (New_CRC_cases - New_CVD_and_CRC_cases - New_diabetes_and_CRC_cases + New_diabetes_and_CVD_and_CRC_cases)
      
      exp_cases <- exp_cases +
        CVD_risk_with_diabetes * 
        (New_diabetes_and_CRC_cases - New_diabetes_and_CVD_and_CRC_cases)
      
    } else if (CRC == 1 && CVD != 1) {
      exp_cases <- exp_cases +
        CVD_risk_no_diabetes * 
        (SW - New_diabetes_cases - New_CVD_cases + New_diabetes_and_CVD_cases)
      
      exp_cases <- exp_cases +
        CVD_risk_with_diabetes * 
        (New_diabetes_cases - New_diabetes_and_CVD_cases)
      
    } else if (CRC != 1 && CVD == 1) {
      exp_cases <- exp_cases +
        CRC_risk_no_diabetes * 
        (SW - New_diabetes_cases - New_CRC_cases + New_diabetes_and_CRC_cases)
      
      exp_cases <- exp_cases +
        CRC_risk_with_diabetes * 
        (New_diabetes_cases - New_diabetes_and_CRC_cases)
      
    } else {
      # Both CRC and CVD == 1; no new cases expected
    }
  }
  
  return(exp_cases)
}

# Function to calculate expected new diabetes, CVD, and CRC cases
expected_new_diabetes_CVD_CRC_cases <- function(row) {
  exp_cases <- 0
  
  Diabetes <- as.numeric(row[['Diabetes']])
  CVD <- as.numeric(row[['CVD']])
  CRC <- as.numeric(row[['CRC']])
  SW <- as.numeric(row[['Sample Weight']])
  
  New_diabetes_cases <- as.numeric(row[['New diabetes cases']])
  New_CVD_cases <- as.numeric(row[['New CVD cases']])
  New_CRC_cases <- as.numeric(row[['New CRC cases']])
  New_diabetes_and_CVD_cases <- as.numeric(row[['New diabetes and CVD cases']])
  New_diabetes_and_CRC_cases <- as.numeric(row[['New diabetes and CRC cases']])
  New_CVD_and_CRC_cases <- as.numeric(row[['New CVD and CRC cases']])
  New_diabetes_and_CVD_and_CRC_cases <- as.numeric(row[['New diabetes and CVD and CRC cases']])
  
  Diabetes_risk <- as.numeric(row[['Diabetes risk']])
  CVD_risk_with_diabetes <- as.numeric(row[['CVD risk with diabetes']])
  CRC_risk_with_diabetes <- as.numeric(row[['CRC risk with diabetes']])
  Diabetes_and_CVD_and_CRC_risk <- as.numeric(row[['Diabetes and CVD and CRC risk']])
  Diabetes_and_CVD_risk <- as.numeric(row[['Diabetes and CVD risk']])
  Diabetes_and_CRC_risk <- as.numeric(row[['Diabetes and CRC risk']])
  CVD_and_CRC_risk_with_diabetes <- as.numeric(row[['CVD and CRC risk with diabetes']])
  
  if (Diabetes == 1) {
    if (CRC != 1 && CVD != 1) {
      exp_cases <- exp_cases +
        CVD_and_CRC_risk_with_diabetes * 
        (SW - New_CVD_cases - New_CRC_cases + New_CVD_and_CRC_cases)
      
    } else if (CRC == 1 && CVD != 1) {
      exp_cases <- exp_cases +
        CVD_risk_with_diabetes * 
        (SW - New_CVD_cases)
      
    } else if (CRC != 1 && CVD == 1) {
      exp_cases <- exp_cases +
        CRC_risk_with_diabetes * 
        (SW - New_CRC_cases)
      
    } else {
      # All diseases present; no new cases expected
    }
    
  } else {
    if (CRC != 1 && CVD != 1) {
      exp_cases <- exp_cases +
        Diabetes_and_CVD_and_CRC_risk * 
        as.numeric(row[['healthy']])
      
      exp_cases <- exp_cases +
        Diabetes_risk * 
        (New_CVD_and_CRC_cases - New_diabetes_and_CVD_and_CRC_cases)
      
      exp_cases <- exp_cases +
        CRC_risk_with_diabetes * 
        (New_diabetes_and_CVD_cases - New_diabetes_and_CVD_and_CRC_cases)
      
      exp_cases <- exp_cases +
        CVD_risk_with_diabetes * 
        (New_diabetes_and_CRC_cases - New_diabetes_and_CVD_and_CRC_cases)
      
      exp_cases <- exp_cases +
        Diabetes_and_CRC_risk * 
        (New_CVD_cases - New_CVD_and_CRC_cases - New_diabetes_and_CVD_cases + New_diabetes_and_CVD_and_CRC_cases)
      
      exp_cases <- exp_cases +
        CVD_and_CRC_risk_with_diabetes * 
        (New_diabetes_cases - New_diabetes_and_CRC_cases - New_diabetes_and_CVD_cases + New_diabetes_and_CVD_and_CRC_cases)
      
      exp_cases <- exp_cases +
        Diabetes_and_CVD_risk * 
        (New_CRC_cases - New_CVD_and_CRC_cases - New_diabetes_and_CRC_cases + New_diabetes_and_CVD_and_CRC_cases)
      
    } else if (CRC == 1 && CVD != 1) {
      exp_cases <- exp_cases +
        Diabetes_and_CVD_risk * 
        (SW - New_diabetes_cases - New_CVD_cases + New_diabetes_and_CVD_cases)
      
      exp_cases <- exp_cases +
        Diabetes_risk * 
        (New_CVD_cases - New_diabetes_and_CVD_cases)
      
      exp_cases <- exp_cases +
        CVD_risk_with_diabetes * 
        (New_diabetes_cases - New_diabetes_and_CVD_cases)
      
    } else if (CRC != 1 && CVD == 1) {
      exp_cases <- exp_cases +
        Diabetes_and_CRC_risk * 
        (SW - New_diabetes_cases - New_CRC_cases + New_diabetes_and_CRC_cases)
      
      exp_cases <- exp_cases +
        Diabetes_risk * 
        (New_CRC_cases - New_diabetes_and_CRC_cases)
      
      exp_cases <- exp_cases +
        CRC_risk_with_diabetes * 
        (New_diabetes_cases - New_diabetes_and_CRC_cases)
      
    } else {
      # Both CRC and CVD == 1; no new cases expected
    }
  }
  
  return(exp_cases)
}





# Function to calculate healthy population
healthy_pop <- function(row) {
  if (as.numeric(row[['Diabetes']]) == 1 || as.numeric(row[['CVD']]) == 1 || as.numeric(row[['CRC']]) == 1) {
    return(0)
  } else {
    with_disease <- as.numeric(row[['New diabetes cases']]) +
      as.numeric(row[['New CVD cases']]) +
      as.numeric(row[['New CRC cases']]) -
      as.numeric(row[['New diabetes and CVD cases']]) -
      as.numeric(row[['New diabetes and CRC cases']]) -
      as.numeric(row[['New CVD and CRC cases']]) +
      as.numeric(row[['New diabetes and CVD and CRC cases']])
    return(as.numeric(row[['Sample Weight']]) - with_disease)
  }
}



