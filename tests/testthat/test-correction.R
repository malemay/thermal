# Tests that verify that the output of the correction methods remains
# the same as the package is developed

# Testing the overall mean correction method
test_that("overall mean correction works as expected", {

		  overall_tmpdir <- withr::local_tempdir(pattern = "overall_test")

		  overall_tmeta <- correct_thermal(base_data = tmeta,
						   correction_type = "overall",
						   output_dir = overall_tmpdir,
						   camera_tz = "Etc/GMT+5",
						   display_tz = "Etc/GMT+4",
						   tags = "minimal",
						   overwrite_dst = TRUE,
						   verbose = FALSE)

		  overall_tmeta$mean <- thermal_mean(overall_tmeta)

		  # Comparing the results without the SourceFile column which is likely to vary
		  expect_identical(overall_tmeta[, -1], readRDS(test_path("test_files", "overall_tmeta.rds"))[, -1])
})

# Testing the linear model correction method
test_that("Linear model correction works as expected", {

		  lm_tmpdir <- withr::local_tempdir(pattern = "lm_test")

		  lm_tmeta <- correct_thermal(base_data = tmeta,
					      correction_type = "lm",
					      output_dir = lm_tmpdir,
					      camera_tz = "Etc/GMT+5",
					      display_tz = "Etc/GMT+4",
					      tags = "minimal",
					      nuc_threshold = 50,
					      overwrite_dst = TRUE,
					      verbose = FALSE)

		  lm_tmeta$mean <- thermal_mean(lm_tmeta)

		  # Comparing the results without the SourceFile column which is likely to vary
		  expect_identical(lm_tmeta[, -1], readRDS(test_path("test_files", "lm_tmeta.rds"))[, -1])
})

# Testing the spline correction method
test_that("Spline correction works as expected", {

		  spline_tmpdir <- withr::local_tempdir(pattern = "spline_test")

		  spline_tmeta <- correct_thermal(base_data = tmeta,
						  correction_type = "spline",
						  output_dir = spline_tmpdir,
						  camera_tz = "Etc/GMT+5",
						  display_tz = "Etc/GMT+4",
						  tags = "minimal",
						  nuc_threshold = 50,
						  overwrite_dst = TRUE,
						  verbose = FALSE)

		  spline_tmeta$mean <- thermal_mean(spline_tmeta)

		  expect_identical(spline_tmeta[, -1], readRDS(test_path("test_files", "spline_tmeta.rds"))[, -1])
})

# Testing the overlap correction method
test_that("Overlap correction works as expected", {

		  # This list will store the output of optimize transform for each image pair
		  tparams <- list()

		  # We only process the first 30 images for testing because otherwise this is very long
		  for(i in 2:30) {
			  # Computing the initial guess for theta
			  theta_guess <- -(tmeta[i, "GimbalYawDegree"] - tmeta[i - 1, "GimbalYawDegree"])

			  # Filling the tparams list with the parameter estimates
			  tparams[[i - 1]] <- optimize_transform(x = terra::rast(tmeta[i - 1, "SourceFile"]),
								 y = terra::rast(tmeta[i, "SourceFile"]),
								 theta1 = theta_guess,
								 min_cor = 0.6,
								 fact = 4,
								 cores = if(.Platform$OS.type == "unix") 2 else 1,
								 reltol = 10^-2)
		  }

		  # Performing the overlap correction
		  overlap_tmpdir <- withr::local_tempdir(pattern = "overlap_test")

		  overlap_tmeta <- correct_thermal(base_data = tmeta[1:30, ],
						   correction_type = "overlap",
						   output_dir = overlap_tmpdir,
						   camera_tz = "Etc/GMT+5",
						   display_tz = "Etc/GMT+4",
						   tags = "minimal",
						   overwrite_dst = TRUE,
						   tparams = tparams,
						   verbose = FALSE)

		  overlap_tmeta$mean <- thermal_mean(overlap_tmeta)

		  expect_identical(overlap_tmeta[, -1], readRDS(test_path("test_files", "overlap_tmeta.rds"))[, -1])

})

