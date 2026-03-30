# This script creates the files that will be used as reference
# for correction when comparing the output of the functions
# to the expected results

# Loading the required packages
library(thermal)

# Loading the metadata
thermal_files <- dir(system.file("extdata", package = "thermal"),
                     pattern = "thermal\\.tiff$", full.names = TRUE)

# Reading the metadata
tmeta <- read_metadata(thermal_files, camera_tz = "Etc/GMT+5", display_tz = "Etc/GMT+4", tags = "minimal")
 
### RESULTS OF THE OVERALL MEAN CORRECTION
overall_tmpdir <- withr::local_tempdir(pattern = "overall_test")

overall_tmeta <- correct_thermal(metadata = tmeta,
                                 correction_type = "overall",
                                 output_dir = overall_tmpdir,
                                 camera_tz = "Etc/GMT+5",
                                 display_tz = "Etc/GMT+4",
                                 tags = "minimal",
                                 overwrite_dst = TRUE)

overall_tmeta$mean <- thermal_mean(overall_tmeta)

### RESULTS OF THE LINEAR MODEL CORRECTION
lm_tmpdir <- withr::local_tempdir(pattern = "lm_test")

lm_tmeta <- correct_thermal(metadata = tmeta,
                            correction_type = "lm",
                            output_dir = lm_tmpdir,
                            camera_tz = "Etc/GMT+5",
                            display_tz = "Etc/GMT+4",
                            tags = "minimal",
			    nuc_threshold = 80,
                            overwrite_dst = TRUE)

lm_tmeta$mean <- thermal_mean(lm_tmeta)

### RESULTS OF THE SPLINE CORRECTION
spline_tmpdir <- withr::local_tempdir(pattern = "spline_test")

spline_tmeta <- correct_thermal(metadata = tmeta,
				correction_type = "spline",
				output_dir = spline_tmpdir,
				camera_tz = "Etc/GMT+5",
				display_tz = "Etc/GMT+4",
				tags = "minimal",
				spline_spar = 0.5,
				nuc_threshold = 80,
				overwrite_dst = TRUE)

spline_tmeta$mean <- thermal_mean(spline_tmeta)

### RESULTS OF THE OVERLAP CORRECTION

# This list will store the output of optimize transform for each image pair
tparams <- list()

for(i in 2:nrow(tmeta)) {
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

overlap_tmeta <- correct_thermal(metadata = tmeta,
                                 correction_type = "overlap",
                                 output_dir = overlap_tmpdir,
                                 camera_tz = "Etc/GMT+5",
                                 display_tz = "Etc/GMT+4",
                                 tags = "minimal",
                                 overwrite_dst = TRUE,
                                 tparams = tparams)

overlap_tmeta$mean <- thermal_mean(overlap_tmeta)

# Saving all the results to file
saveRDS(overall_tmeta, file = "overall_tmeta.rds")
saveRDS(lm_tmeta, file = "lm_tmeta.rds")
saveRDS(spline_tmeta, file = "spline_tmeta.rds")
saveRDS(overlap_tmeta, file = "overlap_tmeta.rds")

