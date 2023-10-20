#' Create a rotation matrix from an angle theta
#'
#' @param theta
#'
#' @return A 2x2 matrix that can be used in matrix product to rotate
#' x-y coordinates in a counterclockwise direction.
#'
#' @export
#'
#' @examples
#' NULL
create_rmat <- function(theta) {
	theta_rad <- theta * pi / 180
	rot_matrix <- matrix(c(cos(theta_rad), sin(theta_rad), -sin(theta_rad), cos(theta_rad)), nrow = 2, ncol = 2)
	rot_matrix
}

#' Rotate x-y coordinates using matrix algebra
#'
#' @param coords A 2-column matrix with x-coordinates as the first column
#' and y-coordinates as the second column
#' @param theta A numeric. The angle (in degrees) to rotate counterclockwise.
#' @param center A numeric vector of length 2 indicating the x, y coordinates
#' of the rotation center
#'
#' @return A matrix of coordinates similar to the input one, but with
#' rotated coordinates.
#'
#' @export
#'
#' @examples
#' NULL
rotate_coords <- function(coords, theta = 0, center = c(0, 0)) {

	# Creating the rotation matrix
	rot_matrix <- create_rmat(theta)

	# We transform the coordinates to rotate around the center of the image
	coords[, 1] <- coords[, 1] - center[1]
	coords[, 2] <- coords[, 2] - center[2]

	# Rotating the coordinates through matrix multiplication and converting them back to the original coordinate system
	coords <- t(rot_matrix %*% t(coords))
	coords[, 1] <- coords[, 1] + center[1]
	coords[, 2] <- coords[, 2] + center[2]

	coords
}

#' A function that performs coordinate translation
#'
#' @param coords A 2-column matrix with x-coordinates as the first column
#' and y-coordinates as the second column
#' @param htrans A numeric. The translation along the x-axis.
#' @param vtrans A numeric. The translation along the y-axis.
#'
#' @return A matrix of coordinates similar to the input one, but with
#' rotated coordinates.
#'
#' @export
#'
#' @examples
#' NULL
translate_coords <- function(coords, htrans = 0, vtrans = 0) {

	# Translating the coordinates
	coords[, 1] <- coords[, 1] + htrans
	coords[, 2] <- coords[, 2] + vtrans

	coords
}

#' Apply a coordinate transformation to raster values
#'
#' @param x A template raster whose values will be changed.
#' @param new_coords A matrix with two columns containing the new coordinates
#' of the raster values, with x-coordinates in the first column and y-coordinates
#' in the second column.
#'
#' @return a raster with similar extent to the previous one, with values modified
#' according to the new coordinates.
#'
#' @export
#'
#' @examples
#' NULL
apply_transform <- function(x, new_coords) {

	# Filtering out coordinates that are outside the boundaries of the raster
	coords_ok <- new_coords[, 1] > xmin(x) & new_coords[, 1] < xmax(x) & new_coords[, 2] > ymin(x) & new_coords[, 2] < ymax(x)

	# Removing rows that are now outside the boundaries
	new_coords <- new_coords[coords_ok, ]
	new_values <- values(x)[coords_ok]

	# Assigning the new values, taking into account the fact that:
	# - x coordinates index into columns, and y-coordinates into rows
	# - the direction of y coordinates is inverted relative to the direction of indexing
	new_coords[, 2] <- ymax(x) - new_coords[, 2]

	# Replacing the values of the raster
	x[] <- NA
	x[new_coords[, c(2, 1)]] <- new_values
	x
}

#' Rotate the values in a raster
#'
#' @param x A raster whose values will be rotated
#' @param theta A numeric. The number of degrees to rotate the raster
#' values, counterclockwise.
#'
#' @return A raster similar to the input one, but whose values have been rotated.
#'
#' @export
#'
#' @examples
#' NULL
rotate_raster <- function(x, theta) {

	# We extract the x-y coordinates of the cells
	new_coords <- rotate_coords(xyFromCell(x, 1:ncell(x)),
				    theta = theta,
				    center = c((xmax(x) - xmin(x)) / 2, (ymax(x) - ymin(x)) / 2))

	apply_transform(x, new_coords)
}

#' Translates the values of a raster
#'
#' @param x A raster whose values will be rotated
#' @param htrans A numeric. The translation along the x-axis.
#' @param vtrans A numeric. The translation along the y-axis.
#'
#' @return A raster similar to the input one, but whose values have been rotated.
#'
#' @export
#'
#' @examples
#' NULL
translate_raster <- function(x, htrans, vtrans) {

	# We extract the coordinates of the cells
	new_coords <- translate_coords(xyFromCell(x, 1:ncell(x)),
				       htrans = htrans,
				       vtrans = vtrans)

	apply_transform(x, new_coords)
}

#' Transform a raster by rotation and translation
#'
#' @param x A raster whose values will be rotated
#' @param theta A numeric. The number of degrees to rotate the raster
#' values, counterclockwise.
#' @param htrans A numeric. The translation along the x-axis.
#' @param vtrans A numeric. The translation along the y-axis.
#'
#' @return A raster similar to the input one, but whose values have been rotated.
#'
#' @export
#'
#' @examples
#' NULL
transform_thermal <- function(x, theta, htrans, vtrans) {

	# We extract the x-y coordinates of the cells
	new_coords <- transform_coords(xyFromCell(x, 1:ncell(x)),
				       theta = theta,
				       center = c((xmax(x) - xmin(x)) / 2, (ymax(x) - ymin(x)) / 2),
				       htrans = htrans,
				       vtrans = vtrans)

	apply_transform(x, new_coords)
}

#' Transform coordinates according to a rotation/translation model
#'
#' @param coords A 2-column matrix with x-coordinates as the first column
#' and y-coordinates as the second column
#' @param theta A numeric. The number of degrees to rotate the raster
#' values, counterclockwise.
#' @param center A numeric vector of length 2 indicating the x, y coordinates
#' of the rotation center
#' @param htrans A numeric. The translation along the x-axis.
#' @param vtrans A numeric. The translation along the y-axis.
#'
#' @return A raster similar to the input one, but whose values have been rotated.
#'
#' @export
#'
#' @examples
#' NULL
transform_coords <- function(coords, theta, center, htrans, vtrans) {

	# Creating the rotation matrix
	rot_matrix <- create_rmat(theta)

	# Appending a translation to it
	transform_matrix <- cbind(rot_matrix, c(htrans, vtrans))
	transform_matrix <- rbind(transform_matrix, c(0, 0, 1))

	# We transform the coordinates to rotate around the center of the image
	coords[, 1] <- coords[, 1] - center[1]
	coords[, 2] <- coords[, 2] - center[2]

	# We embed the translation back to the original coordinates in the transformation matrix<
	transform_matrix[1, 3] <- transform_matrix[1, 3] + center[1]
	transform_matrix[2, 3] <- transform_matrix[2, 3] + center[2]

	# Rotating the coordinates through matrix multiplication and converting them back to the original coordinate system
	(cbind(coords, rep(1, nrow(coords))) %*% t(transform_matrix))[, 1:2, drop = FALSE]
}

#' Assess accuracy of raster transformation by comparing to target raster
#'
#' @param x The target raster.
#' @param y A raster whose coordinates are to be transformed to match the
#' target raster.
#' @param theta A numeric. The number of degrees to rotate the raster
#' values, counterclockwise.
#' @param htrans A numeric. The translation along the x-axis.
#' @param vtrans A numeric. The translation along the y-axis.
#'
#' @return The correlation between the values of the transformed raster
#' and the target raster.
#'
#' @export
#'
#' @examples
#' NULL
assess_transform <- function(x, y, theta, htrans, vtrans) {

	# We extract the x-y coordinates of the cells
	new_coords <- transform_coords(xyFromCell(y, 1:ncell(y)),
				       theta = theta,
				       center = c((xmax(y) - xmin(y)) / 2, (ymax(y) - ymin(y)) / 2),
				       htrans = htrans,
				       vtrans = vtrans)

	# Filtering out coordinates that are outside the boundaries of the target raster
	coords_ok <- new_coords[, 1] > xmin(x) & new_coords[, 1] < xmax(x) & new_coords[, 2] > ymin(x) & new_coords[, 2] < ymax(x)

	# Computing the correlation
	cor(terra::values(x)[terra::cellFromXY(x, new_coords[coords_ok, ])], terra::values(y)[coords_ok, ])
}

#' Optimize the transformation for a raster to match its target
#'
#' @param x The target raster.
#' @param y The raster to transform.
#' @param theta1 A numeric. The initial value for theta.
#' @param htrans1 A numeric. The initial value for htrans.
#' @param vtrans1 A numeric. The initial value for vtrans.
#' @param method A character. The method to use for optimization, passed to optim.
#' @param reltol A numeric. The relative tolerance for convergence in optim.
#'
#' @return The result of running the optim function on the transformation
#' of y to the target raster x.
#'
#' @export
#'
#' @examples
#' NULL
optimize_transform <- function(x, y, theta1 = 0, htrans1 = 0, vtrans1 = 0, method = "Nelder-Mead", reltol = 10^-8) {

	optim_coords <- terra::xyFromCell(y, 1:ncell(y))
	center <- dim(y)[2:1] / 2
	optim_coords <- cbind(optim_coords, rep(1, nrow(optim_coords)))
	optim_coords[, 1] <- optim_coords[, 1] - center[1]
	optim_coords[, 2] <- optim_coords[, 2] - center[2]

	optim(c(theta1, htrans1, vtrans1),
	      fn = function(param, xval, yval, coords, nrows, ncols, extent) {
		      theta <- param[1]
		      htrans <- param[2]
		      vtrans <- param[3]
		      -assess_transform_cpp(xval, yval, coords, theta * pi / 180, htrans + center[1], vtrans + center[2], nrows, ncols, extent)
	      },
	      xval = values(x),
	      yval = values(y),
	      coords = optim_coords,
	      nrows = nrow(x),
	      ncols = ncol(x),
	      extent = unlist(as.list(ext(x))),
	      method = method,
	      control = list(trace = 1, reltol = reltol))
}

#' Verify the accuracy of a raster transformation
#'
#' Interactively click on a raster whose coordinates are to be transformed
#' to match another raster and verify where the clicks land on the target
#' raster.
#'
#' @inheritParams assess_transform
#' @param n An integer. The number of clicks to query for.
#'
#' @return NULL, invisibly. This function is invoked for its plotting side-effect.
#' 
#' @export
#'
#' @examples
#' NULL
check_transform <- function(x, y, theta, htrans, vtrans, n = 10) {
	# Initalizing the plotting regions
	terra::plot(x)
	dev.new()
	terra::plot(y)
	dev.set(dev.prev())

	# Lopping for as many clicks as required
	for(i in 1:n) {
		# Asking for clicks on the untransformed raster
		raw_coords <- click(x, n = 1, xy = TRUE)[, c("x", "y")]
		# Plotting the transformed coordinates on the transformed raster
		dev.set(dev.next())
		transformed_coords <- transform_coords(as.matrix(raw_coords),
						       theta = theta,
						       center = dim(x)[2:1] / 2,
						       htrans = htrans,
						       vtrans = vtrans)

		points(x = transformed_coords[, 1], y = transformed_coords[, 2])
		dev.set(dev.prev())
	}

	invisible(NULL)
}

#' Find initial values for thermal image overlap optimization using a brute-force approach
#'
#' This function will test several combinations of theta, htrans and vtrans for
#' determining the optimal transformation between thermal images to find the best
#' overlap. All combinations of the theta, htrans and vtrans values provided to
#' the function will be tested, and the best combination returned.
#'
#' @param x The target raster.
#' @param y A raster that will be transformed to match the target raster x.
#' @param theta_vect A numeric vector of values to test for theta.
#' @param htrans_vect A numeric vector of values to test for htrans.
#' @param vtrans_vect A numeric vector of values to test for vtrans.
#' @param fact A factor of aggregation for the rasters, passed on to
#' \code{\link[terra]{aggregate}}.  Larger values will speed up computations
#' but will result in less accurate estimates. Default value is 1, i.e. no aggregation.
#' @param cores The number of cores to use when parallelizing the iterations,
#' based on \code{\link[parallel]{mclapply}}. Defaults to 1, i.e. no parallel computing.
#'
#' @return A numeric vector of length 3, with the first value corresponding to theta,
#' the second to htrans, and the third to vtrans.
#'
#' @export
#'
#' @examples
#' NULL
find_initial <- function(x, y, theta_vect, htrans_vect, vtrans_vect, fact = 1, cores = 1) {

	# Reducing the size of the rasters
	x <- terra::aggregate(x, fact = fact)
	terra::ext(x) <- c(0, terra::ncol(x), 0, terra::nrow(x))

	y <- terra::aggregate(y, fact = fact)
	terra::ext(y) <- c(0, terra::ncol(y), 0, terra::nrow(y))

	# Creating the data.frame of combinations to test and adjusting the values of htrans and vtrans
	to_test <- expand.grid(theta = theta_vect, htrans = htrans_vect / fact, vtrans = vtrans_vect / fact)

	# Pre-processing the data for assess_transform_cpp
	optim_coords <- terra::xyFromCell(y, 1:ncell(y))
	center <- dim(y)[2:1] / 2
	optim_coords <- cbind(optim_coords, rep(1, nrow(optim_coords)))
	optim_coords[, 1] <- optim_coords[, 1] - center[1]
	optim_coords[, 2] <- optim_coords[, 2] - center[2]

	# Initializing a vector for the correlations
	correlations <- parallel::mclapply(1:nrow(to_test), function(i, x, y, params) {
						   assess_transform_cpp(xval = values(x),
									yval = values(y),
									coords = optim_coords,
									theta = params[i, 1] * pi / 180,
									htrans = params[i, 2] + center[1],
									vtrans = params[i, 3] + center[2],
									nrows = nrow(x), ncols = ncol(x),
									extent = unlist(as.list(ext(x))))
						       }, x = x, y = y, params = to_test,
						   mc.cores = cores)

	# Extracting the combination of parameters that yielded the best results
	output <- to_test[which.max(unlist(correlations)), ]
	output$htrans <- output$htrans * fact
	output$vtrans <- output$vtrans * fact

	as.numeric(output)
}

