# Tests that verify whether the functions related to radiometric calibration
# work as expected
test_that("Radiometric calibration works as expected", {

		  # Proviving IDs for each picture as row names
		  rownames(tmeta) <- substr(basename(tmeta$SourceFile), 1, 8)

		  # Adding the temperature data to the metadata
		  tmeta <- add_temp_metadata(tmeta, temperature, tolerance = as.difftime(5, units = "secs"))

		  # Joining the pixel and temperature data in a format that can easily be used for fitting models
		  pixel_data <- join_thermal(tmeta, panels, columns = c("black", "gray", "white"))

		  # Computing the linear models of temperature as a function of digital numbers
		  thermal_models <- thermal_lm(metadata = tmeta,
					       pixel_values = pixel_data,
					       summary_functions = list(black = median, gray = median, white = median),
					       surfaces = c("black", "gray", "white"))

		  # Generating the prediction from the models
		  thermal_models$metadata$prediction <- thermal_predict(thermal_models$metadata, dn = 3500)

		  thermal_reference <- readRDS(test_path("test_files", "thermal_models.rds"))
		  reference_meta <- thermal_reference$metadata[, -1]

		  expect_identical(pixel_data, readRDS(test_path("test_files", "pixel_data.rds")))
		  expect_identical(thermal_models$metadata[, -1], reference_meta)
})
