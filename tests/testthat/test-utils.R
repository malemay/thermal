# Test the functions in the file R/utils.R

# Testing the thermal_mean function
test_that("thermal_mean works properly", {
    # We test mean computation in a different way from the first 5 images of the test set
    means <- numeric(length(thermal_files))

    for(i in 1:length(thermal_files)) {
	    i_matrix <- as.matrix(terra::rast(thermal_files[i]))
	    means[i] <- sum(i_matrix) / length(i_matrix)
    }

    # EXPECTATIONS
    # Comparing the home-made solution to thermal_mean
    expect_identical(thermal_mean(tmeta), means)
    # Also testing with multiple cores
    expect_identical(thermal_mean(tmeta, ncores = if(.Platform$OS.type == "unix") 2 else 1), means)
})

# Testing the compute_vignetting function
test_that("compute_vignetting works properly", {
    # Computing vignetting differently from what is done in the package as a double-check
    rasters <- lapply(thermal_files, terra::rast)
    vignetting <- sum(do.call(c, rasters)) / length(rasters)

    # Testing that the computed values are the same
    expect_identical(as.numeric(terra::values(vignetting)), as.numeric(terra::values(compute_vignetting(tmeta))))
})

