# This code prepares the test images that will be distributed
# along with the package for the vignette and examples

# Loading the required packages
library(thermal)
library(terra)
library(exifr)

### Setting some parameters for the generation of images

# The directory containing the images to analyse
flight_dir <- "~/Documents/drone_data/drone_pictures/drone_test_01092023/12h"
# We will use the first n images from the flight
# The first 90 images are from the first of four subflights
n <- 90 
# The aggregation factor used for reducing the size of the thermal images
thermal_agg_fact <- 4
# The aggregation factor used for reducing the size of the visible images
visible_agg_fact <- 10
# The time zone that the images were originally taken in
camera_tz <- "Etc/GMT+5"
# The time zone to use for the analyses
display_tz <- "Etc/GMT+4"


# Loading the metadata of thermal and visible images
thermal_images <- dir(flight_dir, pattern = "R.JPG$", full.names = TRUE)
visible_images <- dir(flight_dir, pattern = ".jpg$", full.names = TRUE)
tmeta <- read_metadata(thermal_images, camera_tz = camera_tz, display_tz = display_tz, tags = "minimal")
vmeta <- read_metadata(visible_images, camera_tz = camera_tz, display_tz = display_tz, tags = "minimal")

# The first n images are from the first subflight - these are the images we will use as test images
tmeta <- tmeta[1:n, ]
vmeta <- vmeta[1:n, ]

# Checking that all thermal images are in the same order as their match
stopifnot(all(as.numeric(substr(basename(vmeta$SourceFile), 5, 8)) - 
	      as.numeric(substr(basename(tmeta$SourceFile), 5, 8)) == 1))

# We will create aggregated images from the original ones
for(i in tmeta$SourceFile) {
	message("Processing image ", i)

	# Aggregating the original raster and adjusting its extent so that it matches the resolution
	i_rast <- aggregate(rast(i), fact = thermal_agg_fact)
	ext(i_rast) <- c(0, ncol(i_rast), 0, nrow(i_rast))

	# Writing the aggregated raster to file and copying the original EXIF tags to it
	output_name <- paste0("../inst/extdata/", sub(pattern = "_R.JPG$", "_thermal.tiff", basename(i)))
	writeRaster(i_rast, filename = output_name, datatype = "INT2U", overwrite = TRUE)
	exiftool_call(args = paste0("-overwrite_original -tagsFromFile ", i), fnames = output_name)
}

# We also aggregate the visible images similarly
for(i in vmeta$SourceFile) {
	message("Processing image ", i)

	# Aggregating the original raster and adjusting its extent so that it matches the resolution
	i_rast <- aggregate(rast(i), fact = visible_agg_fact)
	ext(i_rast) <- c(0, ncol(i_rast), 0, nrow(i_rast))

	# We also need to explicitly set the RGB channels
	RGB(i_rast) <- c(1, 2, 3)

	# Writing the aggregated raster to file and copying the original EXIF tags to it
	output_name <- paste0("../inst/extdata/", sub(pattern = ".jpg$", "_visible.tiff", basename(i)))
	writeRaster(i_rast, filename = output_name, datatype = "INT1U", overwrite = TRUE)
	exiftool_call(args = paste0("-overwrite_original -tagsFromFile ", i), fnames = output_name)

	# Removing the auxiliary file created by terra
	unlink(paste0(output_name, ".aux.xml"))
}

# Checking that the images make sense and they they match correctly

# Reading the metadata of the test files
tmeta_test <- read_metadata(dir("../inst/extdata", pattern = "_thermal.tiff$", full.names = TRUE),
			    camera_tz = camera_tz, display_tz = display_tz)
vmeta_test <- read_metadata(dir("../inst/extdata", pattern = "_visible.tiff$", full.names = TRUE),
			    camera_tz = camera_tz, display_tz = display_tz)

# Matching the visible and thermal images
# We use a simplified matching function because we expect the images to be sorted
matching_function <- function(visible, thermal) 1:nrow(visible)

# Still we check that they matched correctly
file_matches <- match_images(vmeta_test, tmeta_test, match_func = matching_function,
			     max_difftime = as.difftime(1.5, units = "secs"))

stopifnot(all(as.numeric(substr(basename(vmeta_test$SourceFile), 5, 8)) - 
	      as.numeric(substr(basename(tmeta_test$SourceFile), 5, 8)) == 1))

par(mfrow = c(1, 2))

for(i in 1:n) {
	plot(rast(file_matches[i, "visible_file"]), mar = c(0.5, 0.5, 0.5, 0.5), axes = FALSE)
	plot(rast(file_matches[i, "thermal_file"]), mar = c(1.3, 1.3, 1.3, 1.3), axes = FALSE)
	Sys.sleep(0.8)
}

