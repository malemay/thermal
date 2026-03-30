# Tests that verify that the output of the correction methods remains
# the same as the package is developed

# Testing the overall mean correction method
test_that("overall mean correction works as expected", {

		  overall_tmpdir <- withr::local_tempdir(pattern = "overall_test")

		  overall_tmeta <- correct_thermal(metadata = tmeta,
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

		  lm_tmeta <- correct_thermal(metadata = tmeta,
					      correction_type = "lm",
					      output_dir = lm_tmpdir,
					      camera_tz = "Etc/GMT+5",
					      display_tz = "Etc/GMT+4",
					      tags = "minimal",
					      nuc_threshold = 80,
					      overwrite_dst = TRUE,
					      verbose = FALSE)

		  lm_tmeta$mean <- thermal_mean(lm_tmeta)

		  # Comparing the results without the SourceFile column which is likely to vary
		  expect_identical(lm_tmeta[, -1], readRDS(test_path("test_files", "lm_tmeta.rds"))[, -1])
})

# Testing the spline correction method
test_that("Spline correction works as expected", {

		  spline_tmpdir <- withr::local_tempdir(pattern = "spline_test")

		  spline_tmeta <- correct_thermal(metadata = tmeta,
						  correction_type = "spline",
						  output_dir = spline_tmpdir,
						  camera_tz = "Etc/GMT+5",
						  display_tz = "Etc/GMT+4",
						  tags = "minimal",
						  spline_spar = 0.5,
						  nuc_threshold = 80,
						  overwrite_dst = TRUE,
						  verbose = FALSE)

		  spline_tmeta$mean <- thermal_mean(spline_tmeta)

		  expect_identical(spline_tmeta[, -1], readRDS(test_path("test_files", "spline_tmeta.rds"))[, -1])
})

# Testing the overlap correction method
test_that("Overlap correction works as expected", {

		  # Computing the overlap transform parameters
		  tmeta <- compute_overlaps(tmeta,
					    theta_guess = TRUE,
					    min_cor = 0.6,
					    fact = 4,
					    cores = if(.Platform$OS.type == "unix") 2 else 1,
					    verbose = FALSE,
					    reltol = 10^-2)

		  # Performing the overlap correction
		  overlap_tmpdir <- withr::local_tempdir(pattern = "overlap_test")

		  overlap_tmeta <- correct_thermal(metadata = tmeta,
						   correction_type = "overlap",
						   output_dir = overlap_tmpdir,
						   camera_tz = "Etc/GMT+5",
						   display_tz = "Etc/GMT+4",
						   tags = "minimal",
						   overwrite_dst = TRUE,
						   verbose = FALSE)

		  overlap_tmeta$mean <- thermal_mean(overlap_tmeta)

		  expect_identical(overlap_tmeta[, -1], readRDS(test_path("test_files", "overlap_tmeta.rds"))[, -1])
})

