# Test the functions in the file R/utils.R

# Testing the thermal_mean function
test_that("thermal_mean works properly", {
    # We test mean computation in a different way from the first 5 images of the test set
    means <- numeric(length(thermal_files))

    for(i in 1:length(thermal_files)) {
	    i_matrix <- as.matrix(terra::rast(thermal_files[i]))
	    means[i] <- sum(i_matrix) / length(i_matrix)
    }

    tmeta <- read_metadata(thermal_files, camera_tz = "Etc/GMT+5", tags = "minimal")

    # EXPECTATIONS
    # Comparing the home-made solution to thermal_mean
    expect_identical(thermal_mean(tmeta), means)
    # Also testing with multiple cores
    expect_identical(thermal_mean(tmeta, ncores = 2), means)
})
