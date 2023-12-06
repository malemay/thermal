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

#' Create linear models to convert from visible to thermal coordinates
#'
#' To complete
#'
#' @param coords A list of matching points in visible and thermal images, such as returned by
#'               the function matching_points.
#'
#' @return A list of two linear models, the first one allowing to convert x coordinates
#'         and the other allowing to convert y coordinates
#'
#' @examples
#' NULL
#'
#' @export
match_lm <- function(coords) {
	list(x = lm(coords[[2]]$x ~ coords[[1]]$x),
	     x = lm(coords[[2]]$y ~ coords[[1]]$y))
}

#' Plot the results of a linear model that converts from visible to thermal coordinates
#'
#' To be completed
#'
#' @param coord_lm A list of two linear models such as returned by match_lm.
#'
#' @return NULL. This function is only invoked for its plotting side-effect.
#'
#' @examples
#' NULL
#'
#' @export
plot_lm <- function(lmod) {
	par(mfrow = c(1, length(lmod)))

	for(i in 1:length(lmod)) {
		plot(lmod[[i]]$model[[2]], lmod[[i]]$model[[1]],
		     xlab = "Visible image coordinates",
		     ylab = "Thermal image coordinates",
		     main = "")
		     abline(lmod[[i]], lty = 2)

		     slope <- as.character(round(coef(lmod[[i]])[2], 3))
		     intercept <- as.character(round(coef(lmod[[i]])[1], 3))
		     rsquared <- as.character(round(summary(lmod[[i]])$r.squared, 2))
		     resid_sum <- as.character(round(mean(abs(residuals(lmod[[i]]))), 3))

		     mtext(text = paste0("y = ", slope, "x + ", intercept, " r2 = ", rsquared, " Mean of resid. = ", resid_sum),
			   side = 3)
	}

	invisible(NULL)
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
#' @param x A numeric vector of x pixel coordinates
#' @param y A numeric vector of y pixel coordinates
#' @param optimout The output of an alignment model optimization, as returned
#'                 by the function align_images.
#'
#' @return A list of two vectors, x and y, with the coordinates of the
#'         corresponding pixels in the thermal image.
#'
#' @examples
#' NULL
#'
#' @export
convert_coordinates <- function(x, y, optimout, distortion_center = c(2000, 1500)) {

	# Extract the parameters of the optimized model
	slope <- optimout[1]
	bx <- optimout[2]
	by <- optimout[3]
	k <- optimout[4]

	# Removing the distortion in the visible image coordinates
	r2 <- (x - distortion_center[1])^2 + (y - distortion_center[2])^2
	xout <- distortion_center[1] + (x - distortion_center[1]) / (1 + k * r2)
	yout <- distortion_center[2] + (y - distortion_center[2]) / (1 + k * r2)

	# Adjusting for the thermal image coordinates
	xout <- xout * slope + bx
	yout <- yout * slope + by

	return(cbind(xout, yout))
}

#' Faster visible to thermal coordinate transform for performance-critical code
convert_coords_optim <- function(coords, optimout, r2, distortion_center = c(2000, 1500)) {

	# Extract the parameters of the optimized model
	slope <- optimout[1]
	bx <- optimout[2]
	by <- optimout[3]
	k <- optimout[4]

	# Removing the distortion in the visible image coordinates
	# and adjusting for the thermal image coordinates
	denom <- 1 + k * r2
	coords[, 1] <- (distortion_center[1] + coords[, 1] / denom) * slope + bx
	coords[, 2] <- (distortion_center[2] + coords[, 2] / denom) * slope + by

	return(coords)
}

#' Plot aligned images with corresponding points
#'
#' To complete
#'
#' @param visible a terra raster representing a visible image
#' @param thermal a terra raster representing a thermal image
#' @param optimout a list of optimization results as returned by
#'          align_images.
#' @param points_df A data.frame with the location of points
#'          to plot, with a column named x and other named y.
#'
#' @return NULL, invisibly. This function is invoked for plotting.
#'
#' @examples
#' NULL
#' 
#' @export
plot_alignment <- function(visible, thermal, optimout, points_df, distortion_center = c(2000, 1500)) {

	thermal_coords <- convert_coordinates(x = points_df[, 1],
					      y = points_df[, 2],
					      optimout = optimout,
					      distortion_center = distortion_center)

	thermal_points <- points_df
	thermal_points[, 1] <- thermal_coords$x
	thermal_points[, 2] <- thermal_coords$y

	terra::plot(visible)
	plot_colors <- colors()[sample(length(colors()), nrow(points_df))]
	points(points_df, pch = 16, col = plot_colors)
	dev.new()
	terra::plot(thermal)
	points(thermal_points, pch = 16, col = plot_colors)

	invisible(NULL)
}

#' Compare aligned images interactively
#'
#' To complete
#'
#' @param visible a terra raster representing a visible image
#' @param thermal a terra raster representing a thermal image
#' @param optimout a list of optimization results as returned by
#'          align_images.
#' @param nclicks the number of clicks before the function returns.
#'
#' @return NULL, invisibly. This function is invoked for its plotting.
#'
#' @examples
#' NULL
#'
#' @export
thermclick <- function(visible, thermal, optimout, nclicks = 1, distortion_center = c(2000, 1500)) {

	dev.new()
	vdev <- dev.cur()

	terra::plot(visible)

	dev.new()
	tdev <- dev.cur()
	terra::plot(thermal)


	for(i in 1:nclicks) {
		dev.set(vdev)
		vpoint <- terra::click(visible, n = 1, xy = TRUE, col = "red")

		tcoords <- convert_coordinates(x = vpoint$x,
					       y = vpoint$y,
					       optimout = optimout$par,
					       distortion_center = distortion_center)

		dev.set(tdev)
		terra::points(data.frame(x = tcoords[, 1], y = tcoords[, 2]), pch = 16, col = "blue")
	}
}

