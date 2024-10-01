
# Define the run_simulation function
run_simulation <- function(dat, reductions, start_seed, end_seed, red_meat, processed_meat, years, test_mode) {
  for (reduction in reductions) {
    for (seed in start_seed:end_seed) {
      if (test_mode) {
        output_directory <- 'Output/Tests/'
      } else {
        output_directory <- 'Output/Final/'
      }
      
      if (!dir.exists(output_directory)) {
        dir.create(output_directory, recursive = TRUE)
      }
      
      if (red_meat && !processed_meat) {
        mdir <- paste0(reduction, '_reduction_RM_alone')
        file_name <- paste0('df_seed_', seed, '_RM.csv')
      } else if (processed_meat && !red_meat) {
        mdir <- paste0(reduction, '_reduction_PM_alone')
        file_name <- paste0('df_seed_', seed, '_PM.csv')
      } else if (red_meat && processed_meat) {
        mdir <- paste0(reduction, '_reduction')
        file_name <- paste0('df_seed_', seed, '.csv')
      } else if (reduction == 0) {
        mdir <- paste0(reduction, '_reduction')
        file_name <- paste0('df_seed_', seed, '.csv')
      } else {
        stop('Invalid combination of red_meat and processed_meat')
      }
      
      if (test_mode) {
        mdir <- paste0(mdir, '_test/')
      } else {
        mdir <- paste0(mdir, '/')
      }
      
      path <- file.path(output_directory, mdir)
      if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
      }
      
      cat('\n')
      cat(rep('#',20), '\n')
      cat(sprintf('REDUCTION: %s%%, Random Seed: %s, Processed meat: %s, Red meat: %s\n', reduction, seed, processed_meat, red_meat))
      cat(rep('#',20), '\n')
      cat('\n')
      
      if (red_meat) {
        x_RM <- (100 - reduction) / 100
      } else {
        x_RM <- 1
      }
      
      if (processed_meat) {
        x_PM <- (100 - reduction) / 100
      } else {
        x_PM <- 1
      }
      
      df_sim <- dat
      
      if (test_mode) {
        df_sim <- head(df_sim, n = 50)
      }
      
      # Run the simulation functions
      df_sim <- run_sampling(df_sim, x_PM = x_PM, x_RM = x_RM, seed = seed)
      df_sim <- calculate_risks(df_sim, mortality_table = mortality_table, seed = seed)
      
      # Initialize cumulative case columns if they don't exist
      df_sim$`New diabetes cases` <- 0
      df_sim$`New CVD cases` <- 0
      df_sim$`New CRC cases` <- 0
      df_sim$`New diabetes and CVD cases` <- 0
      df_sim$`New diabetes and CRC cases` <- 0
      df_sim$`New CVD and CRC cases` <- 0
      df_sim$`New diabetes and CVD and CRC cases` <- 0
      df_sim$healthy <- 0
      
      for (year in 1:years) {
        df_sim <- update_mortalities(df_sim, year = year, pre_new_cases = TRUE)
        df_sim <- update_cases(df_sim, year = year)
        df_sim <- update_mortalities(df_sim, year = year, pre_new_cases = FALSE)
        df_sim <- apply_new_mortalities(df_sim, year = year)
        df_sim <- update_df(df_sim)
        df_sim <- calculate_risks(df_sim, mortality_table = mortality_table, seed = seed)
      }
      
      # Write the data frame to CSV
      write.csv(df_sim, file = file.path(path, file_name), row.names = FALSE)
    }
  }
}

# Main analysis function
run_analysis <- function(
    diseases,
    reductions,
    red_meat,
    processed_meat,
    mortalities,
    years,
    data_dir,
    with_filter,
    start_seed,
    end_seed,
    test_mode
) {
  for (reduction in reductions) {
    for (NCD in diseases) {
      # Determine the intervention path based on meat reduction
      if (red_meat && processed_meat) {
        int_dir <- paste0(reduction, '_reduction')
      } else if (red_meat && !processed_meat) {
        int_dir <- paste0(reduction, '_reduction_RM_alone')
      } else if (!red_meat && processed_meat) {
        int_dir <- paste0(reduction, '_reduction_PM_alone')
      } else if (reduction == 0) {
        int_dir <- paste0(reduction, '_reduction')
      } else {
        stop('Must specify a reduction in either red or processed meat')
      }
      
      if (test_mode) {
        baseline_dir <- paste0('0_reduction_test')
        int_dir <- paste0(int_dir, '_test')
      } else {
        baseline_dir <- paste0('0_reduction')
      }
      
      baseline_path <- file.path(data_dir, baseline_dir)
      intervention_path <- file.path(data_dir, int_dir)
      
      cases_prevented_total <- c()
      
      for (seed in start_seed:end_seed) {
        # Read baseline and intervention data
        baseline_file <- file.path(baseline_path, paste0('df_seed_', seed, '.csv'))
        if (red_meat && processed_meat) {
          int_file <- file.path(intervention_path, paste0('df_seed_', seed, '.csv'))
        } else if (red_meat && !processed_meat) {
          int_file <- file.path(intervention_path, paste0('df_seed_', seed, '_RM.csv'))
        } else if (!red_meat && processed_meat) {
          int_file <- file.path(intervention_path, paste0('df_seed_', seed, '_PM.csv'))
        }
        
        # Check if files exist
        if (!file.exists(baseline_file)) {
          stop(sprintf('Baseline file not found: %s', baseline_file))
        }
        if (!file.exists(int_file)) {
          stop(sprintf('Intervention file not found: %s', int_file))
        }
        
        baseline_data <- read.csv(baseline_file)
        int_data <- read.csv(int_file)
        
        # Adjust age if necessary
        baseline_data$Age <- baseline_data$Age - years
        int_data$Age <- int_data$Age - years
        
        # Apply demographic filter if required
        if (with_filter) {
          # Example filter: Only males
          baseline_data <- subset(baseline_data, Sex == 'Male')
          int_data <- subset(int_data, Sex == 'Male')
          # Add additional filters as needed
        }
        
        cases_prevented_total_seed <- 0
        
        for (year in 1:years) {
          # Adjust column names for CSV files (replace spaces and special characters with dots)
          if (mortalities) {
            baseline_col <- paste0(NCD, '.mortalities.year.', year, '.post')
            int_col <- paste0(NCD, '.mortalities.year.', year, '.post')
          } else {
            baseline_col <- paste0('New.', NCD, '.cases.year.', year)
            int_col <- paste0('New.', NCD, '.cases.year.', year)
          }
          
          # Check if columns exist
          if (!(baseline_col %in% names(baseline_data))) {
            stop(sprintf('Column not found in baseline data: %s', baseline_col))
          }
          if (!(int_col %in% names(int_data))) {
            stop(sprintf('Column not found in intervention data: %s', int_col))
          }
          
          baseline_cases <- sum(baseline_data[[baseline_col]], na.rm = TRUE)
          int_cases <- sum(int_data[[int_col]], na.rm = TRUE)
          
          cases_prevented_year <- baseline_cases - int_cases
          cases_prevented_total_seed <- cases_prevented_total_seed + cases_prevented_year
        }
        cases_prevented_total <- c(cases_prevented_total, cases_prevented_total_seed)
      }
      
      mean_cases <- mean(cases_prevented_total)
      lower <- quantile(cases_prevented_total, probs = 0.025)
      upper <- quantile(cases_prevented_total, probs = 0.975)
      
      # Print the results
      cat(sprintf(
        '%d%% reduction in %s %s:\n',
        reduction,
        ifelse(red_meat, 'red meat', ''),
        ifelse(processed_meat, 'and processed meat', '')
      ))
      cat(sprintf(
        'Total %s %s prevented over %d years: %.2f (95%% CI: %.2f - %.2f)\n\n',
        NCD,
        ifelse(mortalities, 'mortalities', 'cases'),
        years,
        mean_cases,
        lower,
        upper
      ))
    }
  }
}
