# Tests that verify whether the functions related to radiometric calibration
# work as expected
test_that("Radiometric calibration works as expected", {

		  # Proviving IDs for each picture as row names
		  rownames(tmeta) <- substr(basename(tmeta$SourceFile), 1, 8)

		  # Adding the temperature data to the metadata
		  tmeta <- add_temp_metadata(tmeta, temperature, tolerance = as.difftime(5, units = "secs"))

		  # Computing the linear models of temperature as a function of digital numbers
		  thermal_models <- thermal_lm(metadata = tmeta,
					       polygons = panels,
					       summary_functions = list(black = median, gray = median, white = median),
					       surfaces = c("black", "gray", "white"),
					       ncores = 2)

		  thermal_reference <- readRDS(test_path("test_files", "thermal_models.rds"))

		  # Testing whether the pixels extracted from the panels are the same
		  expect_identical(thermal_models$pixels, thermal_reference$pixels)

		  # Testing whether the updated metadata (with linear model parameters) is the same
		  reference_meta <- thermal_reference$metadata[, -1]

		  # We allow for some numerical differences between model estimates by rounding to the 12th decimal
		  thermal_models$metadata$intercept <- round(thermal_models$metadata$intercept, 12)
		  reference_meta$intercept <- round(reference_meta$intercept, 12)

		  thermal_models$metadata$slope <- round(thermal_models$metadata$slope, 12)
		  reference_meta$slope <- round(reference_meta$slope, 12)

		  expect_identical(thermal_models$metadata[, -1], reference_meta)
		  
		  # We want to check that thermal_lm works even when supplying a single value for summary_functions
		  single_func_models <- thermal_lm(metadata = tmeta,
						   polygons = panels,
						   summary_functions = median,
						   surfaces = c("black", "gray", "white"),
						   ncores = 2)
		  
		  single_func_models$metadata$intercept <- round(single_func_models$metadata$intercept, 12)
		  single_func_models$metadata$slope <- round(single_func_models$metadata$slope, 12)

		  expect_identical(single_func_models$metadata[, -1], reference_meta)

		  # Also testing global_lm functionality
		  global_model <- global_lm(thermal_models)
		  global_model$metadata$intercept <- round(global_model$metadata$intercept, 12)
		  global_model$metadata$slope <- round(global_model$metadata$slope, 12)

		  global_reference <- readRDS(test_path("test_files", "global_model.rds"))
		  global_reference$metadata$intercept <- round(global_reference$metadata$intercept, 12)
		  global_reference$metadata$slope <- round(global_reference$metadata$slope, 12)

		  # We do not test the models themselves because they contain complex data
		  # However, we do test both the metadata and pixels
		  expect_identical(global_model$metadata[, -1], global_reference$metadata[, -1])
		  expect_identical(global_model$pixels, global_reference$pixels)

		  # We also test that the predictions made using this model remain the same
		  reference_temp <- terra::rast(test_path("test_files", "surface_temp.tiff"), noflip = TRUE)
		  test_rast <- terra::rast(tmeta[1, "SourceFile"], noflip = TRUE)

		  # We test both ways to pass the model parameters to thermal_predict
		  # Up to some numerical precision
		  expect_identical(round(terra::values(reference_temp), 6),
				   round(terra::values(thermal_predict(x = test_rast,
								       model = global_model$models$global)), 6))

		  expect_identical(round(terra::values(reference_temp), 6),
				   round(terra::values(thermal_predict(x = test_rast,
								       slope = global_model$metadata$slope[1],
								       intercept = global_model$metadata$intercept[1])), 6))

		  # We test that thermal_predict generates the expected warnings and error
		  expect_warning(thermal_predict(x = test_rast,
						 model = global_model$models$global,
						 slope = global_model$metadata$slope[1],
						 intercept = global_model$metadata$intercept[1]))

		  expect_warning(thermal_predict(x = test_rast,
						 model = global_model$models$global,
						 intercept = global_model$metadata$intercept[1]))

		  expect_warning(thermal_predict(x = test_rast,
						 model = global_model$models$global,
						 slope = global_model$metadata$slope[1]))

		  expect_error(thermal_predict(x = test_rast))
})
