#################### functions for moralities calculations ###################

# mortality_table (this is a nested list)
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

# mortality_prob function
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

# CRC_mortality_prob function
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