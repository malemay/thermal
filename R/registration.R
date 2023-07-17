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
align_images <- function(visible, thermal, start_values, method = c("Nelder-Mead", "BFGS")) {
	# Getting the sum of the values from the visible raster
	visible_sum <- sum(visible)

	# Setting the positions to query in the thermal raster
	tcoords <- xyFromCell(trast, 1:ncell(trast))
	# And the corresponding values
	tvalues <- extract(trast, tcoords)[[1]]

	ofunct <- function(x, vsum, tcoords, tvalues) {
		ax <- x[1]
		bx <- x[2]
		ay <- x[3]
		by <- x[4]

		# Getting the corresponding position in visible raster given params
		vcoords <- tcoords
		vcoords[, 1] <- (vcoords[, 1] - bx) / ax
		vcoords[, 2] <- (vcoords[, 2] - by) / ay

		# And extracting the corresponding values
		vvalues <- extract(vsum, vcoords)[[1]]

		# We return 0 if the correlation would be NA
		if(all(is.na(tvalues) | is.na(vvalues))) {
			return(0)
		} else {
			# Getting the negative correlation because we want to minimize
			return(-cor(tvalues, vvalues, use = "complete.obs"))
		}
	}

	# We return the results of the optimization
	optim(start_values, ofunct, method = method, control = list(trace = 3), vsum = visible_sum, tcoords = tcoords, tvalues = tvalues)
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
plot_alignment <- function(visible, thermal, optimout, points_df) {
	ax <- optimout$par[1]
	bx <- optimout$par[2]
	ay <- optimout$par[3]
	by <- optimout$par[4]

	thermal_points <- points_df
	thermal_points[, 1] <- thermal_points[, 1] * ax + bx
	thermal_points[, 2] <- thermal_points[, 2] * ay + by

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
thermclick <- function(visible, thermal, optimout, nclicks = 1) {

	# Extracting the parameters from the optimization
	ax <- optimout$par[1]
	bx <- optimout$par[2]
	ay <- optimout$par[3]
	by <- optimout$par[4]

	dev.new()
	vdev <- dev.cur()

	terra::plot(visible)

	# Find the corners of the thermal image on the visible and plot them
	xmin <- -bx / ax
	xmax <- (640 - bx) / ax
	ymin <- -by / ay
	ymax <- (512 - by) / ay
	polygon(c(xmin, xmin, xmax, xmax), c(ymin, ymax, ymax, ymin))

	dev.new()
	tdev <- dev.cur()
	terra::plot(thermal)


	for(i in 1:nclicks) {
		dev.set(vdev)
		vpoint <- terra::click(visible, n = 1, xy = TRUE, col = "red")
		tx <- vpoint$x * ax + bx
		ty <- vpoint$y * ay + by

		dev.set(tdev)
		terra::points(data.frame(x = tx, y = ty), pch = 16, col = "blue")
	}
}
