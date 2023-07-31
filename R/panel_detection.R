#' Convert a visible image to a binary raster of uniform surfaces
#'
#' To complete
#'
#' @param input_file A character. The path of tfe visible image to read.
#' @param output_file A character. The path of the image to output
#' @param canny_param A character. A geometry string used as input to image_canny.
#' @param kernel A character. The kernel to use as input to image_morphology at the
#'               step used to identify homogeneous surfaces.
#' @param output_format A character. The image format to use for output. Defaults to png
#'                      because it is not lossy.
#'
#' @return NULL, invisibly. This function is called for its file writing effect
#'
#' @examples
#' NULL
#'
#' @export
find_surfaces <- function(input_file, output_file, canny_param = "5x1+5%+10%", kernel = "Square:5", output_format = "png") {
	# First we need to read in the image with colorspace=RGB (instead of sRGB)
	img <- magick::image_convert(magick::image_read(input_file), colorspace = "rgb")

	# Then we carry canny edge identification and threshold the image so that it is binary (edge/no edge)
	img <- magick::image_threshold(magick::image_canny(img, geometry = canny_param), threshold = "50%")

	# Then we negate the image and identify homogeneous regions of the image using a square kernel
	img <- magick::image_morphology(magick::image_negate(img), method = "Erode", kernel = kernel)

	# We write the resulting image to file
	magick::image_write(img, path = output_file, format = output_format)

	# Destroying the image because we no longer need it and have to free memory
	magick::image_destroy(img)

	# Also using garbage collection because destroying the image does not seem to free the memory in itself
	gc()

	invisible(NULL)
}

#' Identify panels from a set of polygons
#'
#' To complete
#'
#' @param filename A character. The path to a file containing candidate polygons.
#' @param min_area A numeric. The minimum area (in pixels) for a polygon to be considered.
#' @param max_area A numeric. The maximum area (in pixels) for a polygon to be considered.
#' @param max_ratio A numeric. The maximum ratio of the polygon perimeter to the perimeter
#'                  of a perfect square that would have the same area as the polygon. The
#'                  lower this value is, the more the shapes kept are close to a perfect square.
#' @param max_distance A numeric. The maximum distance (as the sum of all values in the distance
#'                     matrix) between the final set of three polygons for them to be considered.
#'
#' @return A list of three elements containing polygons at various filtering steps:
#'         all_polygons: All polygons input to the function
#'	   filtered: The polygons that passed the area and perimeter ratio filters
#'         selected: The three filtered polygons closest to each other, or NULL if no solution could be found
#'
#' @examples
#' NULL
#'
#' @export
find_panels <- function(filename, min_area = 800, max_area = 3000, max_ratio = 1.5, max_distance = 2000) {
	# Reading the polygons
	polygons <- sf::st_read(filename, quiet = TRUE)

	# Subsetting to only those with a value of 1
	polygons <- polygons[polygons$DN == 1, ]

	# Computing their area and perimeter
	polygons$area <- sf::st_area(polygons)
	polygons$perimeter <- lwgeom::st_perimeter(polygons)

	# Computing the ratio of the perimeter to the perimeter if the polygon were a square
	polygons$peri_ratio <- polygons$perimeter / (sqrt(polygons$area) * 4)

	# Creating a copy of all polygons for debugging purposes
	all_polygons <- polygons

	# Removing polygons that are too small or too big to be the panels
	polygons <- polygons[polygons$area >= min_area & polygons$area <= max_area, ]

	# Keeping only those with an acceptable ratio
	polygons <- polygons[polygons$peri_ratio <= max_ratio, ]

	filtered_polygons <- polygons

	# If there are fewer than 3 polygons remaining at this point, we did not find a solution
	n_polygons <- nrow(polygons)

	if(n_polygons < 3) return(list(all_polygons = all_polygons, filtered = filtered_polygons, selected = NULL))

	# Finding the three such polygons with the smallest distance to each other
	to_test <- t(combn(1:n_polygons, 3) )

	distance_matrix <- sf::st_distance(polygons)

	distance_sum <- apply(to_test, 1, function(x) sum(distance_matrix[x, x]))

	output <- polygons[to_test[which.min(distance_sum), ],]

	if(min(distance_sum) > max_distance) output <- NULL

	return(list(all_polygons = all_polygons, filtered = filtered_polygons, selected = output))
}

#' Generate a set of plots to debug panel identification
#'
#' To complete
#'
#' @param image A named character vector with the names corresponding to the names of the
#'              polygons and the values being the paths to their corresponding visible image
#'              files.
#' @param polygons A list of lists of polygons such as returned by \code{\link{find_panels}}.
#'                 Each element of the list corresponds to a different picture. The names of
#'                 the list elements must match the names in the image vector, as they will
#'                 be used to match polygons with their target image.
#' @param output_dir The directory in which to create the plots. Defaults to the current
#'                   directory.
#' @param all_poly A logical. Whether to output the plots showing all polygons
#' @param filtered_poly A logical. Whether to output the plots showing filtered polygons
#' @param selected_poly A logical. Whether to output the plots showing selected polygons
#' @param verbose A logical. Whether the user should be kept informed of progress.
#' @param res. A numeric. The resolution of the output images, passed to the png device.
#'
#' @return NULL, invisibly
#'
#' @examples
#' NULL
#'
#' @export
debug_panels <- function(images, polygons, output_dir = ".",
			 all_poly = TRUE, filtered_poly = TRUE, selected_poly = TRUE,
			 verbose = TRUE, res = 200) {

	# Sanity check
	stopifnot(all(names(polygons) %in% names(images)))

	# Looping over all the polygon lists
	for(prefix in names(polygons)) {

		if(verbose) message("Processing image ", prefix)
		i_polygons <- polygons[[prefix]]

		# Plotting the selected panels over the image itself
		png(paste0(output_dir, "/", prefix, "_1panels.png"), width = 4, height = 4, units = "in", res = res)
		terra::plot(rast(images[prefix]), main = if(is.null(i_polygons$selected)) "NO PANELS" else "PANELS FOUND")
		if(!is.null(i_polygons$selected)) {
			plot(i_polygons$selected, add = TRUE, col = "transparent",
			     border = if(is.null(i_polygons$selected$color)) "black" else i_polygons$selected$color)
		}
		dev.off()

		# Generating the area and perimeter ratio plots for all polygons if requested
		if(all_poly) {
			png(paste0(output_dir, "/", prefix, "_2all_area.png"), width = 4, height = 4, units = "in", res = res)
			plot(i_polygons$all_polygons[, "area"])
			dev.off()

			png(paste0(output_dir, "/", prefix, "_3all_ratio.png"), width = 4, height = 4, units = "in", res = res)
			plot(i_polygons$all_polygons[, "peri_ratio"])
			dev.off()
		}

		# Generating the area and perimeter ratio plots for filtered polygons if requested and if there are any
		if(filtered_poly && nrow(i_polygons$filtered)) {
			png(paste0(output_dir, "/", prefix, "_4filtered_area.png"), width = 4, height = 4, units = "in", res = res)
			plot(i_polygons$filtered[, "area"])
			dev.off()

			png(paste0(output_dir, "/", prefix, "_5filtered_ratio.png"), width = 4, height = 4, units = "in", res = res)
			plot(i_polygons$filtered[, "peri_ratio"])
			dev.off()
		}

		# Generating the area and perimeter ratio plots for selected polygons if requested and if there are any
		if(selected_poly && !is.null(i_polygons$selected)) {
			png(paste0(output_dir, "/", prefix, "_6selected_area.png"), width = 4, height = 4, units = "in", res = res)
			plot(i_polygons$selected[, "area"])
			dev.off()

			png(paste0(output_dir, "/", prefix, "_7selected_ratio.png"), width = 4, height = 4, units = "in", res = res)
			plot(i_polygons$selected[, "peri_ratio"])
			dev.off()
		}
	}

	return(invisible(NULL))
}

#' Classify polygons according to panel color
#'
#' To complete
#'
#' @param picture A raster representing a visible image.
#' @param polygons A set of three polygons representing the panels found
#'
#' @return A set of polygons similar to the input ones, but with added
#'         columns for the mean value of the pixel sums over the polygon
#'         and the color identified for the polygon (black, gray or white). 
#'
#' @examples
#' NULL
#'
#' @export
classify_panels <- function(picture, polygons) {

	# Extracting the pixels from the visible picture
	# at the location of the panels and computing the
	# sum of the three channels
	pixel_values <- terra::extract(picture, polygons)
	pixel_values$sum <- apply(pixel_values[, 2:4], 1, sum)

	# Creating a column in the polygons object with the sum
	polygons$pixel_mean <- 0
	polygon_means <- tapply(pixel_values$sum, pixel_values$ID, mean)
	polygons$pixel_mean[as.integer(names(polygon_means))] <- unname(polygon_means)

	# Reordering the polygons according to pixel means and attributing the color accordingly
	polygons <- polygons[order(polygons$pixel_mean), ]
	polygons$color <- c("black", "gray", "white")

	polygons
}

