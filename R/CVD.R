#################### functions for CVD calculations ###################

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

# updated_Framingham_CVD_risk function
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

# Hazard_ratio_CM_CVD function
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

# Hazard_ratio_RM_CVD function
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
