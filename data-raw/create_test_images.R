# This code prepares the test images that will be distributed
# along with the package for the vignette and examples

# Loading the required packages
library(thermal)
library(terra)
library(exifr)

### Setting some parameters for the generation of images

# The directory containing the images to analyse
flight_dir <- "~/Documents/large_data/drone_data/drone_pictures/drone_test_01092023/12h"
# We will use the first n images from the flight
# The first 90 images are from the first of four subflights
n <- 90 
# The aggregation factor used for reducing the size of the thermal images
thermal_agg_fact <- 8
# The time zone that the images were originally taken in
camera_tz <- "Etc/GMT+5"
# The time zone to use for the analyses
display_tz <- "Etc/GMT+4"

# Loading the metadata of thermal images
thermal_images <- dir(flight_dir, pattern = "R.JPG$", full.names = TRUE)
tmeta <- read_metadata(thermal_images, camera_tz = camera_tz, display_tz = display_tz, tags = "minimal")

# The first n images are from the first subflight - these are the images we will use as test images
tmeta <- tmeta[1:n, ]

# We will create aggregated images from the original ones
# We keep only one image out of three to save disk space and processing time
for(i in tmeta$SourceFile[c(TRUE, FALSE, FALSE)]) {
	message("Processing image ", i)

	# Aggregating the original raster and adjusting its extent so that it matches the resolution
	i_rast <- aggregate(rast(i), fact = thermal_agg_fact)
	ext(i_rast) <- c(0, ncol(i_rast), 0, nrow(i_rast))

	# Writing the aggregated raster to file and copying the original EXIF tags to it
	output_name <- paste0("../inst/extdata/", sub(pattern = "_R.JPG$", "_thermal.tiff", basename(i)))
	writeRaster(i_rast, filename = output_name, datatype = "INT2U", overwrite = TRUE)
	exiftool_call(args = paste0("-overwrite_original -tagsFromFile ", i), fnames = output_name)
}

