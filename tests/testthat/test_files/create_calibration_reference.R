# This script creates the reference results to test
# that radiometric calibration results remain the same
# despite changes to the interface

# Loading the package
# The data objects 'temperature' and 'panels' are loaded along with the package
# functions
library(thermal)

# Getting a vector of absolute paths to the test TIR images
thermal_files <- dir(system.file("extdata", package = "thermal"),
                     pattern = "thermal\\.tiff$", full.names = TRUE)

# Reading metadata from EXIF tags
tmeta_df <- read_metadata(thermal_files,
                          camera_tz = "Etc/GMT+5",
                          display_tz = "Etc/GMT+4",
                          tags = "minimal")

# Proviving IDs for each picture as row names
rownames(tmeta_df) <- substr(basename(tmeta_df$SourceFile), 1, 8)

# Adding the temperature data to the metadata
tmeta_df <- add_temp_metadata(tmeta_df, temperature, tolerance = as.difftime(5, units = "secs"))

# Joining the pixel and temperature data in a format that can easily be used for fitting models
pixel_data <- join_thermal(tmeta_df, panels)

# Computing the linear models of temperature as a function of digital numbers
thermal_models <- thermal_lm(metadata = tmeta_df,
                             pixel_values = pixel_data,
                             summary_functions = list(black = median, gray = median, white = median),
                             use_panels = c("black", "gray", "white"))

thermal_models$metadata$prediction <- thermal_predict(thermal_models$metadata, dn = 3500)

# These two objects should be sufficient to test the stability of the calibration output
saveRDS(pixel_data, "pixel_data.rds")
saveRDS(thermal_models, "thermal_models.rds")
