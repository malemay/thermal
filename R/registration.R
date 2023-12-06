#' Identify matching points on visible and thermal images
#'
#' To be completed
#'
#' @param visible A terra raster with a visible image
#' @param thermal A terra raster with the corresponding thermal image
#' @param npoints Integer. The number of values to query.
#'
#' @return A data.frame with the pixel coordinates of matching points
#'         on both pictures
#'
#' @examples
#' NULL
#'
#' @export
matching_points <- function(visible, thermal, npoints = 10) {
	vlist <- list()
	tlist <- list()

	for(i in 1:npoints) {
		terra::plot(thermal)
		tlist[[i]] <- terra::click(thermal, n = 1, xy = TRUE)

		terra::plot(visible)
		vlist[[i]] <- terra::click(visible, n = 1, xy = TRUE)
	}

	vlist <- do.call("rbind", vlist)
	tlist <- do.call("rbind", tlist)

	list(vlist, tlist)
}

#' Align a thermal and a visible image
#'
#' To be completed
#'
#' @param visible a terra raster representing a visible image
#' @param thermal a terra raster representing a thermal image
#' @param start_values a numeric of length 4, the starting values for
#'          the optimization procedure. The algorithm optimizes a linear
#'          model that translates coordinates on the visible image to
#'          coordinates on the thermal image. The four parameters to
#'          optimize are the slope and intercept of the translation
#'          in x and the slope and intercept of the translation in y.
#' @param the method to use for the optim function.
#'
#' @return A list of optimization results, as returned by optim.
#'
#' @examples
#' NULL
#'
#' @export
align_images <- function(visible, thermal, start_values, distortion_center = c(2000, 1500),
			 aggregate_factor = 1, crop_values = c(0, 0), min_overlap = 10000,
			 method = c("Nelder-Mead", "BFGS"), reltol = 10^-8) {

	# Getting the sum of the values from the visible raster
	visible_sum <- sum(visible)

	# Aggregating the values if provided
	if(aggregate_factor != 1) visible_sum <- terra::aggregate(visible_sum, fact = aggregate_factor)

	# Setting the positions to query in the visible raster
	vcoords <- xyFromCell(visible_sum, 1:ncell(visible_sum))

	# Cropping some of the coordinates if provided
	if(! all(crop_values == 0)) {
		vcoords <- vcoords[vcoords[, 1] > crop_values[1] & vcoords[, 1] < (xmax(visible_sum) - crop_values[1]) &
				   vcoords[, 2] > crop_values[2] & vcoords[, 2] < (ymax(visible_sum) - crop_values[2]), ]
	}

	# And the corresponding values
	vvalues <- values(visible_sum, mat = FALSE)[cellFromXY(visible_sum, vcoords)]
	tvalues <- values(thermal, mat = FALSE)

	# Pre-compute r2 distances prior to input to convert_coords_optim
	r2 <- (vcoords[, 1] - distortion_center[1])^2 + (vcoords[, 2] - distortion_center[2])^2

	# Also pre-transform the visible coordinates
	vcoords[, 1] <- vcoords[, 1] - distortion_center[1]
	vcoords[, 2] <- vcoords[, 2] - distortion_center[2]

	# We return the results of the optimization
	optim(start_values, function(x, vcoords, r2, distortion_center, vvalues, tvalues, nrows, ncols, extent, min_overlap) {
		      -assess_registration_cpp(x, vcoords, r2, distortion_center, vvalues, tvalues, nrows, ncols, extent, min_overlap)
			 },
		      vcoords = vcoords, r2 = r2,
		      distortion_center = distortion_center,
		      vvalues = vvalues, tvalues = tvalues,
		      nrows = nrow(thermal), ncols = ncol(thermal),
		      extent = unlist(as.list(ext(thermal))),
		      min_overlap = min_overlap,
		      method = method,
		      control = list(trace = 3, reltol = reltol))
}

#' Convert visible image coordinates to thermal image coordinates
#'
#' To complete
#'
#' @param coords A two-column (x-y) matrix of coordinates to convert.
#' @param params A numeric vector of transformation parameters (slope, bx, by, k).
#'
#' @return A matrix similar to the input with converted coordinates.
#'
#' @examples
#' NULL
#'
#' @export
convert_coordinates <- function(coords, params, distortion_center = c(2000, 1500)) {

	# Extract the parameters of the optimized model
	slope <- params[1]
	bx <- params[2]
	by <- params[3]
	k <- params[4]

	# Removing the distortion in the visible image coordinates
	r2 <- (coords[, 1] - distortion_center[1])^2 + (coords[, 2] - distortion_center[2])^2
	denom <- 1 + k * r2
	coords[, 1] <- distortion_center[1] + (coords[, 1] - distortion_center[1]) / denom
	coords[, 2] <- distortion_center[2] + (coords[, 2] - distortion_center[2]) / denom

	# Adjusting for the thermal image coordinates
	coords[, 1] <- coords[, 1] * slope + bx
	coords[, 2] <- coords[, 2] * slope + by

	return(coords)
}

#' Compare aligned images interactively
#'
#' To complete
#'
#' @param visible a terra raster representing a visible image
#' @param thermal a terra raster representing a thermal image
#' @param params a list of optimization results as returned by
#'          align_images.
#' @param nclicks the number of clicks before the function returns.
#'
#' @return NULL, invisibly. This function is invoked for its plotting.
#'
#' @examples
#' NULL
#'
#' @export
thermclick <- function(visible, thermal, params, nclicks = 1, distortion_center = c(2000, 1500)) {

	dev.new()
	vdev <- dev.cur()

	terra::plot(visible)

	dev.new()
	tdev <- dev.cur()
	terra::plot(thermal)

	for(i in 1:nclicks) {
		dev.set(vdev)
		vpoint <- terra::click(visible, n = 1, xy = TRUE, col = "red")

		tcoords <- convert_coordinates(coords = as.matrix(vpoint[, c("x", "y")]),
					       params = params,
					       distortion_center = distortion_center)

		dev.set(tdev)
		terra::points(data.frame(x = tcoords[, 1], y = tcoords[, 2]), pch = 16, col = "blue")
	}
}

