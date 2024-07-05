#' Identify matching points on visible and thermal images
#'
#' To complete. Functionality should eventually be added to fit a coordinate
#' conversion model to manually identified points.
#'
#' @param visible A terra raster representing a visible image.
#' @param thermal A terra raster representing the corresponding thermal image
#' @param npoints an integer. The number of values to query.
#'
#' @return A list with the pixel coordinates of matching points
#'         on both pictures
#'
#' @export
#' @examples
#' NULL
#'
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
#' This function aims to identify transformation parameters that allow to
#' convert the coordinates on a visible raster to corresponding coordinates
#' on a thermal raster.
#'
#' The coordinate transformation model implemented is relatively simple
#' and involves four parameters:
#' - slope: this parameter adjusts the resolution of the visible image
#' to that of the thermal image, and is the same in the x- and y- directions.
#' - bx: this is the offset of the left margin of the thermal image relative
#' to that of the visible image
#' - by: similar to bx, but for the top margin
#' - k: a distortion parameter as implemented in the Brown-Conrady model of
#' image distortion
#'
#' In the coordinate transformation framework, the visible image is first
#' adjusted for distortion using the k parameter and the Brown-Conrady model.
#' x- and y- coordinates are then adjusted based on simple linear models
#' that relate the visible and thermal coordinates according to the following
#' equations:
#' thermal_x = visible_x * slope + bx
#' thermal_y = visible_y * slope + by
#'
#' Optimization is performed by finding the parameters that maximize the
#' correlation between the sum of all visible channels and the thermal values.
#'
#' The coordinate transformation of the optimization procedure described here
#' is implemented in C++ for maximum speed. For more general R code to compute
#' transformed coordinates based on a given set of parameters, see \code{\link{convert_coordinates}}.
#'
#' @param visible a terra raster representing a visible image.
#' @param thermal a terra raster representing a thermal image.
#' @param start_values a numeric vector of length 4 containing the starting
#' values for the optimization procedure. The elements are in the following
#' order: slope, bx, by, k (see details for explanations).
#' @param distortion_center A numeric vector of length two containing coordinates
#' of the center of the visible image for the purposes of computing distortion
#' parameters.
#' @param aggregate_factor A numeric of length one. If the value is different
#' from 1 (the default value), the visible image will be aggregated according
#' to this factor (passed to \code{\link[terra]{aggregate}}) to speed up computation
#' and allow consistency with the resolution of the thermal raster.
#' @param crop_values A numeric vector of length 2 specifying the extent of coordinates
#' to crop from the visible image in the horizontal (x) direction and vertical (y)
#' direction. Cropping can be useful to speed up computation in cases where the
#' thermal raster is expected to cover a reduced area compared to the visible image.
#' If both values are 0 (the default), then no cropping is done.
#' @param min_overlap A numeric of length one. The minimal number of visible and thermal
#' raster values that must be compared for the correlation to be valid. Otherwise,
#' correlation is set to 0 and the algorithm therefore will not converge towards parameters
#' that result in a low number of overlapping raster cells.
#' @param method the method to use for the \code{\link{optim}} function.
#' @param reltol A numeric of length one. The relative tolerance for the optimization
#' procedure implemented by \code{\link{optim}}. Lower values should yield a more
#' precise result at the cost of longer computation.
#'
#' @return A list of optimization results, as returned by optim.
#'
#' @export
#' @examples
#' NULL
#'
align_images <- function(visible, thermal, start_values, distortion_center = c(2000, 1500),
			 aggregate_factor = 1, crop_values = c(0, 0), min_overlap = 10000,
			 method = "Nelder-Mead", reltol = 10^-8) {

	# Getting the sum of the values from the visible raster
	visible_sum <- sum(visible)

	# Aggregating the values if provided
	if(aggregate_factor != 1) visible_sum <- terra::aggregate(visible_sum, fact = aggregate_factor)

	# Setting the positions to query in the visible raster
	vcoords <- terra::xyFromCell(visible_sum, 1:(terra::ncell(visible_sum)))

	# Cropping some of the coordinates if provided
	if(! all(crop_values == 0)) {
		vcoords <- vcoords[vcoords[, 1] > crop_values[1] & vcoords[, 1] < (terra::xmax(visible_sum) - crop_values[1]) &
				   vcoords[, 2] > crop_values[2] & vcoords[, 2] < (terra::ymax(visible_sum) - crop_values[2]), ]
	}

	# And the corresponding values
	vvalues <- terra::values(visible_sum, mat = FALSE)[terra::cellFromXY(visible_sum, vcoords)]
	tvalues <- terra::values(thermal, mat = FALSE)

	# Pre-compute r2 distances prior to input to convert_coords_optim
	r2 <- (vcoords[, 1] - distortion_center[1])^2 + (vcoords[, 2] - distortion_center[2])^2

	# Also pre-transform the visible coordinates
	vcoords[, 1] <- vcoords[, 1] - distortion_center[1]
	vcoords[, 2] <- vcoords[, 2] - distortion_center[2]

	# We return the results of the optimization
	stats::optim(start_values, function(x, vcoords, r2, distortion_center, vvalues, tvalues, nrows, ncols, extent, min_overlap) {
		      -assess_registration_cpp(x, vcoords, r2, distortion_center, vvalues, tvalues, nrows, ncols, extent, min_overlap)
			 },
		      vcoords = vcoords, r2 = r2,
		      distortion_center = distortion_center,
		      vvalues = vvalues, tvalues = tvalues,
		      nrows = terra::nrow(thermal), ncols = terra::ncol(thermal),
		      extent = unlist(as.list(terra::ext(thermal))),
		      min_overlap = min_overlap,
		      method = method,
		      control = list(trace = 3, reltol = reltol))
}

#' Convert visible image coordinates to thermal image coordinates
#'
#' This function converts visible image coordinates to their corresponding
#' position on a matching thermal image using a 4-parameter transformation
#' model described in detail in \code{\link{align_images}}.
#'
#' @param coords A two-column (x-y) matrix of coordinates to convert.
#' @param params A numeric vector of length 4 representing transformation
#' parameters (slope, bx, by, k).
#' @param distortion_center A numeric vector of length two containing coordinates
#' of the center of the visible image for the purposes of asjusting distortion.
#'
#' @return A matrix similar to the input with converted coordinates.
#'
#' @export
#' @examples
#' NULL
#'
convert_coordinates <- function(coords, params, distortion_center = c(2000, 1500)) {

	# Extract the parameters of the optimized model
	stopifnot(length(params) == 4)
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
#' Given a set of transformation parameters from visible image to thermal image
#' coordinates, this function allows to interactively and graphically check the
#' quality of transformation parameters. Both the visible and thermal image
#' will be displayed, and the user is then prompted to click on the visible
#' image. The corresponding position on the thermal raster will then be
#' displayed automatically.
#'
#' @param visible a terra raster representing a visible image.
#' @param thermal a terra raster representing a thermal image.
#' @param params A numeric vector of length 4 representing transformation
#' parameters (slope, bx, by, k).
#' @param distortion_center A numeric vector of length two containing coordinates
#' of the center of the visible image for the purposes of asjusting distortion.
#' @param nclicks the number of clicks that the user is prompted for before the
#' function returns.
#'
#' @return NULL, invisibly. This function is invoked for its plotting side-effect.
#'
#' @export
#' @examples
#' NULL
#'
thermclick <- function(visible, thermal, params, nclicks = 1, distortion_center = c(2000, 1500)) {

	grDevices::dev.new()
	vdev <- grDevices::dev.cur()

	terra::plot(visible)

	grDevices::dev.new()
	tdev <- grDevices::dev.cur()
	terra::plot(thermal)

	for(i in 1:nclicks) {
		grDevices::dev.set(vdev)

		message("Click on the visible image")
		vpoint <- terra::click(visible, n = 1, xy = TRUE, col = "red")

		tcoords <- convert_coordinates(coords = as.matrix(vpoint[, c("x", "y")]),
					       params = params,
					       distortion_center = distortion_center)

		grDevices::dev.set(tdev)
		terra::points(data.frame(x = tcoords[, 1], y = tcoords[, 2]), pch = 16, col = "blue")
	}
}

