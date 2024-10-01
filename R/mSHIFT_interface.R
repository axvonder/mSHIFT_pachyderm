# Inside R/mSHIFT_interface.R

#' Run mSHIFT Simulation
#'
#' This function runs the mSHIFT simulation with specified parameters.
#'
#' @param data Input data frame for the simulation.
#' @param reductions Numeric vector of reduction percentages.
#' @param start_seed Starting seed number.
#' @param end_seed Ending seed number.
#' @param red_meat Logical indicating whether to reduce red meat.
#' @param processed_meat Logical indicating whether to reduce processed meat.
#' @param years Number of years to run the simulation.
#' @param test_mode Logical indicating test mode.
#' @return None
#' @export
run_sim <- function(data, reductions = c(0, 25, 50, 75, 100), start_seed = 1, end_seed = 5,
                                  red_meat = TRUE, processed_meat = TRUE, years = 10, test_mode = TRUE) {

  # Call the internal function to run the simulation
  run_simulation(
    dat = data,
    reductions = reductions,
    start_seed = start_seed,
    end_seed = end_seed,
    red_meat = red_meat,
    processed_meat = processed_meat,
    years = years,
    test_mode = test_mode
  )
}


# Inside R/mSHIFT_interface.R

#' Run mSHIFT Analysis
#'
#' This function runs the mSHIFT analysis with specified parameters.
#'
#' @param diseases Character vector of diseases to analyze.
#' @param reductions Numeric vector of reduction percentages.
#' @param red_meat Logical indicating whether red meat was reduced.
#' @param processed_meat Logical indicating whether processed meat was reduced.
#' @param mortalities Logical indicating whether to analyze mortalities.
#' @param years Number of years to analyze.
#' @param data_dir Directory containing simulation outputs.
#' @param with_filter Logical indicating whether to apply demographic filters.
#' @param start_seed Starting seed number.
#' @param end_seed Ending seed number.
#' @param test_mode Logical indicating test mode.
#' @return None
#' @export
run_analysis <- function(diseases = c('diabetes', 'CVD', 'CRC'), reductions = c(0, 25, 50, 75, 100),
                                red_meat = TRUE, processed_meat = TRUE, mortalities = TRUE,
                                years = 10, data_dir = NULL, with_filter = FALSE,
                                start_seed = 1, end_seed = 5, test_mode = TRUE) {
  if (is.null(data_dir)) {
    data_dir <- ifelse(test_mode, 'Output/Tests/', 'Output/Final/')
  }
  # Call the internal function to run the analysis
  run_analysis(
    diseases = diseases,
    reductions = reductions,
    red_meat = red_meat,
    processed_meat = processed_meat,
    mortalities = mortalities,
    years = years,
    data_dir = data_dir,
    with_filter = with_filter,
    start_seed = start_seed,
    end_seed = end_seed,
    test_mode = test_mode
  )
}



#' Load Sample Data
#'
#' This function loads the sample data included in the package.
#'
#' @return Data frame containing sample data.
#' @examples
#' \dontrun{
#' sample_data <- load_sample_data()
#' }
#' @export
load_sample_data <- function() {
  data_path <- system.file("extdata", "mSHIFT_two_day_recall_data.plk", package = "mSHIFT")
  # Get the path to the Python script
  python_script_path <- system.file("python", "pickle_reader.py", package = "mSHIFT")
  # Check and install 'pandas' if necessary
  if (!reticulate::py_module_available("pandas")) {
    reticulate::py_install("pandas")
  }
  # Source the Python script
  reticulate::source_python(python_script_path)
  # Read the pickle file using the Python function
  data <- read_pickle_file(data_path)
  return(data)
}




