

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



