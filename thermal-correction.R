## ----options, include = FALSE-------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6,
  fig.align = "center")

## ----thermal-files------------------------------------------------------------
# Loading the package
library(thermal)

# Getting a vector of absolute paths to the test TIR images
thermal_files <- dir(system.file("extdata", package = "thermal"),
                     pattern = "thermal\\.tiff$", full.names = TRUE)

head(thermal_files)

## ----raster-images, echo = FALSE, fig.align = "center"------------------------
par(mfrow = c(2, 2))

for(i in thermal_files[1:4]) {
	terra::plot(terra::rast(i), mar = rep(0.5, 4), legend = FALSE, axes = FALSE)
}

## ----read-metadata------------------------------------------------------------
# Reading metadata from EXIF tags
tmeta_df <- read_metadata(thermal_files,
                          camera_tz = "Etc/GMT+5",
                          display_tz = "Etc/GMT+4",
                          tags = "minimal")

str(tmeta_df)

## ----position-plot------------------------------------------------------------
# If base_crs is not set, it will be set to WGS84 (EPSG:4326) with a warning
plot_metadata(tmeta_df, base_crs = sf::st_crs("EPSG:4326"))

## ----position-plot-crs--------------------------------------------------------
plot_metadata(tmeta_df, base_crs = sf::st_crs("EPSG:4326"), new_crs = sf::st_crs("EPSG:32198"))

## ----position-plot-arrows-----------------------------------------------------
plot_metadata(tmeta_df, base_crs = sf::st_crs("EPSG:4326"),
              arrows = TRUE, arrow_angle = 20,
              arrow_length = grid::unit(0.15, "inches"))

## ----position-plot-factor-----------------------------------------------------
bearing <- ifelse(abs(tmeta_df$GimbalYawDegree) > 150, "south",
                  ifelse(abs(tmeta_df$GimbalYawDegree) < 50, "north", "west"))

plot_metadata(tmeta_df, base_crs = sf::st_crs("EPSG:4326"), color = as.factor(bearing))

## ----position-plot-numeric----------------------------------------------------
plot_metadata(tmeta_df, base_crs = sf::st_crs("EPSG:4326"), color = tmeta_df$GimbalYawDegree)

## ----thermal-mean-------------------------------------------------------------
# We store the mean of each image as a new column in the metadata data.frame
tmeta_df$mean <- thermal_mean(tmeta_df)

## ----raw-mean-plot------------------------------------------------------------
plot_drift(tmeta_df)

## ----vignetting---------------------------------------------------------------
vignetting <- compute_vignetting(tmeta_df)
terra::plot(vignetting)

## ----overall-model-plot-------------------------------------------------------
plot_drift(tmeta_df)

## ----overall-correction-------------------------------------------------------
# Creating a temporary directory to write the corrected images to
# Normally users will use a directory on persistent storage for retrieval later
tmpdir <- tempdir()

# Running correct_thermal on the test data
# Time zone arguments need to be provided again because EXIF tags are
# transferred unmodified from the source files to the corrected files
overall_tmeta <- correct_thermal(metadata = tmeta_df,
                                 correction_type = "overall",
                                 output_dir = tmpdir,
                                 camera_tz = "Etc/GMT+5",
                                 display_tz = "Etc/GMT+4",
                                 tags = "minimal",
                                 overwrite_dst = TRUE)

## ----overall-mean-plot--------------------------------------------------------
plot_drift(overall_tmeta)

## ----nuc-computation----------------------------------------------------------
plot_drift(tmeta_df, nuc_threshold = 80)

## ----lm-model-plot------------------------------------------------------------
plot_fit(tmeta_df, nuc_threshold = 80, method = "lm")

## ----lm-correction------------------------------------------------------------
lm_tmeta <- correct_thermal(metadata = tmeta_df,
                            correction_type = "lm",
                            output_dir = tmpdir,
                            camera_tz = "Etc/GMT+5",
                            display_tz = "Etc/GMT+4",
                            tags = "minimal",
                            nuc_threshold = 80,
                            overwrite_dst = TRUE,
                            verbose = FALSE)

## ----lm-mean-plot-------------------------------------------------------------
plot_drift(lm_tmeta)

## ----spline-model-plot--------------------------------------------------------
plot_fit(tmeta_df, nuc_threshold = 80, method = "spline")

## ----spline-model-plot-spar1--------------------------------------------------
plot_fit(tmeta_df, nuc_threshold = 80, method = "spline", spline_spar = 0.5)

## ----spline-correction--------------------------------------------------------
spline_tmeta <- correct_thermal(metadata = tmeta_df,
                                correction_type = "spline",
                                spline_spar = 0.5,
                                output_dir = tmpdir,
                                camera_tz = "Etc/GMT+5",
                                display_tz = "Etc/GMT+4",
                                tags = "minimal",
                                nuc_threshold = 80,
                                overwrite_dst = TRUE,
                                verbose = FALSE)

## ----spline-mean-plot---------------------------------------------------------
plot_drift(spline_tmeta)

## ----optimize-example---------------------------------------------------------
# We can provide a first guess of what the rotation will be based on the
# difference in gimbal yaw angle recorded by the camera at the two time points
theta_guess <- -(tmeta_df[2, "GimbalYawDegree"] - tmeta_df[1, "GimbalYawDegree"])

trans_params <- optimize_transform(x = terra::rast(tmeta_df[1, "SourceFile"]),
                                   y = terra::rast(tmeta_df[2, "SourceFile"]),
                                   theta1 = theta_guess,
                                   reltol = 10^-3)

## ----check-overlaps, fig.width = 8, fig.height = 6----------------------------
# Setting the seed for reproducibility of the vignette
set.seed(1)

# Verifying the transformation parameters non-interactively
check_transform(x = terra::rast(tmeta_df[1, "SourceFile"]),
                y = terra::rast(tmeta_df[2, "SourceFile"]),
                params = trans_params, n = 40)

## ----optimize-all-------------------------------------------------------------
# Adding the transformation parameters to the metadata
tmeta_df <- compute_overlaps(tmeta_df,
                             min_cor = 0.6,
                             fact = 4,
                             cores = if(.Platform$OS.type == "unix") 2 else 1,
                             verbose = FALSE,
                             reltol = 10^-2)

str(tmeta_df)

## ----corr-hist, fig.width = 6, fig.height = 6---------------------------------
hist(tmeta_df$corr, breaks = 10, main = "Histogram of image overlap correlation values",
     ylab = "Count", xlab = "Correlation")

## ----check-overlaps-mincorr, fig.width = 8, fig.height = 6--------------------
# Setting the seed for reproducibility of the vignette
set.seed(1)

# Finding the image index with the lowest correlation
min_image <- which.min(tmeta_df$corr)

# Verifying the transformation parameters non-interactively
check_transform(x = terra::rast(tmeta_df[min_image - 1, "SourceFile"]),
                y = terra::rast(tmeta_df[min_image, "SourceFile"]),
                params = trans_params, n = 40)

## ----overlap-correction-------------------------------------------------------
overlap_tmeta <- correct_thermal(metadata = tmeta_df,
                                 correction_type = "overlap",
                                 output_dir = tmpdir,
                                 camera_tz = "Etc/GMT+5",
                                 display_tz = "Etc/GMT+4",
                                 tags = "minimal",
                                 overwrite_dst = TRUE,
                                 verbose = FALSE)

## ----overlap-mean-plot--------------------------------------------------------
plot_drift(overlap_tmeta)

## ----session-info-------------------------------------------------------------
sessionInfo()

