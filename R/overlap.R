#' Create a rotation matrix from an angle theta
#'
#' @param theta The angle (in degrees) to rotate the coordinates with.
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
	coords_ok <- new_coords[, 1] > terra::xmin(x) & new_coords[, 1] < terra::xmax(x) & new_coords[, 2] > terra::ymin(x) & new_coords[, 2] < terra::ymax(x)

	# Removing rows that are now outside the boundaries
	new_coords <- new_coords[coords_ok, ]
	new_values <- terra::values(x)[coords_ok]

	# Assigning the new values, taking into account the fact that:
	# - x coordinates index into columns, and y-coordinates into rows
	# - the direction of y coordinates is inverted relative to the direction of indexing
	new_coords[, 2] <- terra::ymax(x) - new_coords[, 2]

	# Replacing the values of the raster
	x[] <- NA
	x[new_coords[, c(2, 1)]] <- new_values
	x
}

#' Transform a raster by rotation and translation
#'
#' @param x A raster whose values will be rotated
#' @param theta A numeric. The number of degrees to rotate the raster
#' values, counterclockwise.
#' @param htrans A numeric. The translation along the x-axis.
#' @param vtrans A numeric. The translation along the y-axis.
#' @param reverse A logical. Whether the transformation (translation and rotation)
#' should be inversed relative to the input parameters. If TRUE, then the inverse
#' of the transformation matrix is used instead of the matrix itself.
#'
#' @return A raster similar to the input one, but whose values have been rotated.
#'
#' @export
#'
#' @examples
#' NULL
transform_raster <- function(x, theta, htrans, vtrans, reverse = FALSE) {

	# We extract the x-y coordinates of the cells
	new_coords <- transform_coords(terra::xyFromCell(x, 1:(terra::ncell(x))),
				       theta = theta,
				       center = c((terra::xmax(x) - terra::xmin(x)) / 2, (terra::ymax(x) - terra::ymin(x)) / 2),
				       htrans = htrans,
				       vtrans = vtrans,
				       reverse = reverse)

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
#' @param reverse A logical. Whether the transformation (translation and rotation)
#' should be inversed relative to the input parameters. If TRUE, then the inverse
#' of the transformation matrix is used instead of the matrix itself.
#'
#' @return A raster similar to the input one, but whose values have been rotated.
#'
#' @export
#'
#' @examples
#' NULL
transform_coords <- function(coords, theta, center, htrans, vtrans, reverse = FALSE) {

	# Creating the rotation matrix
	rot_matrix <- create_rmat(theta)

	# Appending a translation to it
	transform_matrix <- cbind(rot_matrix, c(htrans, vtrans))
	transform_matrix <- rbind(transform_matrix, c(0, 0, 1))
	if(reverse) transform_matrix <- solve(transform_matrix)

	# We transform the coordinates to rotate around the center of the image
	coords[, 1] <- coords[, 1] - center[1]
	coords[, 2] <- coords[, 2] - center[2]

	# We embed the translation back to the original coordinates in the transformation matrix<
	transform_matrix[1, 3] <- transform_matrix[1, 3] + center[1]
	transform_matrix[2, 3] <- transform_matrix[2, 3] + center[2]

	# Rotating the coordinates through matrix multiplication and converting them back to the original coordinate system
	(cbind(coords, rep(1, nrow(coords))) %*% t(transform_matrix))[, 1:2, drop = FALSE]
}

#' Optimize the transformation for a raster to match its target
#'
#' @param x The target raster.
#' @param y The raster to transform.
#' @param theta1 A numeric. The initial value for theta.
#' @param htrans1 A numeric. The initial value for htrans.
#' @param vtrans1 A numeric. The initial value for vtrans.
#' @param theta_range A numeric. The ± range of values to explore around
#' theta1 should a call to find_initial be necessary.
#' @param theta_length A numeric. The number of values of theta to test initially
#' if a call to find_initial is required.
#' @param htrans_range Similar to theta_range for the htrans parameter.
#' @param htrans_length Similar to theta_length for the htrans parameter.
#' @param vtrans_range Similar to theta_range for the vtrans parameter
#' @param vtrans_length Similar to theta_length for the vtrans parameter.
#' @param min_overlap A numeric value of length one. The minimum number of values
#' that must be compared between two rasters for the correlation to be considered
#' valid. Otherwise, a correlation of 0 is returned.
#' @param fact A numeric. The factor to use to reduce the resolution of the
#' rasters to speed up the estimation of the best parameters in find_initial.
#' Defaults to 1 (no aggregation is performed).
#' @param cores The number of cores to use for parallel computing in
#' find_initial. Defaults to 1.
#' @param min_cor A numeric. The minimum correlation that must be obtained
#' from initial guesses to go on with optim.
#' @param maxiter An integer. The maximum number of iterations of find_initial
#' to go through. The number of parameter combinations texted increases by
#' 8 times for each iteration, so this number should be kept low otherwise
#' computation time will become extremely long.
#' @param method A character. The method to use for optimization, passed to optim.
#' @param reltol A numeric. The relative tolerance for convergence in optim.
#' @param verbose A logical value. If TRUE, progress on the iterations
#' in the search for parameters will be output to the console.
#' @param trace A numeric value passed to optim to determine the level of
#' output verbosity. By default, it is set to the same value as verbose
#' (1 if TRUE and 0 if FALSE), but this parameter allows to control it
#' separately.
#'
#' @return The result of running the optim function on the transformation
#' of y to the target raster x.
#'
#' @export
#'
#' @examples
#' NULL
optimize_transform <- function(x, y, theta1 = 0, htrans1 = 0, vtrans1 = 0,
			       theta_range = 5, theta_length = 9,
			       htrans_range = 20, htrans_length = 5,
			       vtrans_range = 20, vtrans_length = 5,
			       min_overlap = 100,
			       fact = 1, cores = 1, min_cor = 0.8, maxiter = 2,
			       method = "Nelder-Mead", reltol = 10^-5,
			       verbose = TRUE, trace = verbose) {

	# Pre-processing the data that does not vary from one call to another
	optim_coords <- terra::xyFromCell(y, 1:(terra::ncell(y)))
	center <- c(terra::xmax(y) - terra::xmin(y), terra::ymax(y) - terra::ymin(y)) / 2
	optim_coords <- cbind(optim_coords, rep(1, nrow(optim_coords)))
	optim_coords[, 1] <- optim_coords[, 1] - center[1]
	optim_coords[, 2] <- optim_coords[, 2] - center[2]

	best_params <- c(theta1, htrans1, vtrans1)

	best_cor <- assess_transform_cpp(xval = terra::values(x), yval = terra::values(y),
					 coords = optim_coords,
					 theta = theta1 * pi / 180, htrans = htrans1 + center[1], vtrans = vtrans1 + center[2],
					 nrows = nrow(x), ncols = ncol(x),
					 extent = unlist(as.list(terra::ext(x))),
					 min_overlap = min_overlap)

	if(verbose) {
		message("Initial guess: cor = ", best_cor)
		message("Current best parameters: theta = ", best_params[1], ", htrans = ", best_params[2], ", vtrans = ", best_params[3])
	}

	# Keep track of the number of find_initial iterations that we have been through
	niter <- 0

	while(best_cor < min_cor && niter < maxiter) {
		# Increment the number of iterations
		niter <- niter + 1

		# Create the vector of values to test
		theta_vect <- seq(theta1 - theta_range, theta1 + theta_range, length.out = theta_length)
		htrans_vect <- seq(htrans1 - htrans_range, htrans1 + htrans_range, length.out = htrans_length)
		vtrans_vect <- seq(vtrans1 - vtrans_range, vtrans1 + vtrans_range, length.out = vtrans_length)

		if(verbose) message("Running iteration ", niter, " of find_initial on ", theta_length * htrans_length * vtrans_length, " values.")

		fi_output <- find_initial(x, y, theta_vect, htrans_vect, vtrans_vect, fact = fact, cores = cores)

		if(verbose) message("find_initial iteration ", niter, ": cor = ", fi_output$value)

		if(fi_output$value > best_cor) {
			best_cor <- fi_output$value
			best_params <- fi_output$param 
		}

		if(verbose) message("Current best parameters: theta = ", best_params[1], ", htrans = ", best_params[2], ", vtrans = ", best_params[3])

		# Update the range and length values
		theta_range <- theta_range * 2
		htrans_range <- htrans_range * 2
		vtrans_range <- vtrans_range * 2

		theta_length <- theta_length * 2
		htrans_length <- htrans_length * 2
		vtrans_length <- vtrans_length * 2
	}

	optim_output <- 
		stats::optim(best_params,
			     fn = function(param, xval, yval, coords, nrows, ncols, extent, min_overlap) {
				     theta <- param[1]
				     htrans <- param[2]
				     vtrans <- param[3]
				     -assess_transform_cpp(xval, yval, coords, theta * pi/180, htrans + center[1], vtrans + center[2], nrows, ncols, extent, min_overlap)
			     },
			     xval = terra::values(x),
			     yval = terra::values(y),
			     coords = optim_coords,
			     nrows = nrow(x),
			     ncols = ncol(x),
			     extent = unlist(as.list(terra::ext(x))),
			     min_overlap = min_overlap,
			     method = method,
			     control = list(trace = trace, reltol = reltol))

	if(verbose) {
		message("Final best parameters: theta = ", optim_output$par[1], ", htrans = ", optim_output$par[2], ", vtrans = ", optim_output$par[3])
		message("Final cor: ", -optim_output$value)
	}

	# We return comprehensive information to allow debugging
	list(optim = optim_output,
	     init = c(theta1, htrans1, vtrans1),
	     fi_params = best_params,
	     niter = niter)
}

#' Verify the accuracy of a raster transformation
#'
#' This function accepts two rasters and transformation parameters between the
#' two, and generates a plot showing the correspondence between features in both
#' rasters. It supports interactive analysis by clicking on a raster whose
#' coordinates are to be transformed to match another raster and plotting where
#' the clicks land on the target raster. Interactive analysis opens two distinct
#' display devices whereas non-interactive analysis splits a display device
#' into two parts.
#'
#' @param x,y Two rast objects whose transform is to be verified.
#' @param params A list of optimized parameters, as output by
#' \code{\link{optimize_transform}}.
#' @param reverse A logical. Whether the transformation (translation and
#' rotation) should be inversed relative to the input parameters. If TRUE, then
#' the inverse of the transformation matrix is used instead of the matrix
#' itself. This is TRUE by default because transformation parameters internally
#' describe the transformation from one image to the previous one, but the
#' other way around is more intuitive. When TRUE, x is the image that was taken
#' first, and y is the image that was taken after.
#' @param interactive A logical value: should the transform be checked
#' interactively?
#' @param n An integer. The number of points to plot or clicks to query for.
#' Defaults to 10.
#'
#' @return NULL, invisibly. This function is invoked for its plotting side-effect.
#' 
#' @export
#'
#' @examples
#' NULL
check_transform <- function(x, y, params, reverse = TRUE, interactive = FALSE, n = 10) {
	# Extracting the parameters from the list provided
	params <- params$optim$par
	theta <- params[1]
	htrans <- params[2]
	vtrans <- params[3]

	# Get the locations of the target raster corners on the source one
	x_corners <- get_corners(x = y,
				 theta = theta,
				 htrans = htrans,
				 vtrans = vtrans,
				 reverse = FALSE)

	# Get the locations of the source raster corners on the target one
	y_corners <- get_corners(x = x,
				 theta = theta,
				 htrans = htrans,
				 vtrans = vtrans,
				 reverse = TRUE)

	if(interactive) {
		# Initializing the plotting regions and plot the outline of the rasters
		terra::plot(x)
		plot(x_corners, border = "red", add = TRUE)

		grDevices::dev.new()
		terra::plot(y)
		plot(y_corners, border = "red", add = TRUE)

		grDevices::dev.set(grDevices::dev.prev())

		# Lopping for as many clicks as required
		for(i in 1:n) {
			# Asking for clicks on the untransformed raster
			raw_coords <- terra::click(x, n = 1, xy = TRUE, pch = 1)[, c("x", "y")]

			# Transforming the coordinates
			transformed_coords <- transform_coords(as.matrix(raw_coords),
							       theta = theta,
							       center = dim(x)[2:1] / 2,
							       htrans = htrans,
							       vtrans = vtrans,
							       reverse = reverse)

			# Plotting the transformed coordinates on the target raster
			grDevices::dev.set(grDevices::dev.next())
			graphics::points(x = transformed_coords[, 1], y = transformed_coords[, 2])
			grDevices::dev.set(grDevices::dev.prev())
		}
	} else {

		# Split the plotting region in two
		graphics::par(mfrow = c(1, 2))

		# Generate points randomly from the raster to transform
		rpoints <- terra::xyFromCell(x, 1:(terra::ncell(x)))[sample(terra::ncell(x), n), ]

		# Sampling some random colors as well
		rcolors <- sample(grDevices::colors(), n)

		# Plot the x raster with points on top of it
		terra::plot(x, main = "x")
		graphics::points(rpoints, col = rcolors)

		# Plot the corners of the target raster on this one
		plot(x_corners, border = "red", add = TRUE)

		# Plotting the raster to transform with points on top
		terra::plot(y, main = "y")

		# With transformed coordinates on top of it
		graphics::points(transform_coords(rpoints,
						  theta = theta, htrans = htrans, vtrans = vtrans,
						  center = c((terra::xmax(x) - terra::xmin(x)) / 2, (terra::ymax(x) - terra::ymin(x)) / 2),
						  reverse = TRUE),
				 col = rcolors)


		# Plot the corners of the source raster on this one
		plot(y_corners, border = "red", add = TRUE)
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
#' @param min_overlap A numeric value of length one. The minimum number of values
#' that must be compared between two rasters for the correlation to be considered
#' valid. Otherwise, a correlation of 0 is returned.
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
find_initial <- function(x, y, theta_vect, htrans_vect, vtrans_vect, min_overlap = 100, fact = 1, cores = 1) {

	# Reducing the size of the rasters
	if(fact != 1) x <- terra::aggregate(x, fact = fact)
	terra::ext(x) <- c(0, terra::ncol(x), 0, terra::nrow(x))

	if(fact != 1) y <- terra::aggregate(y, fact = fact)
	terra::ext(y) <- c(0, terra::ncol(y), 0, terra::nrow(y))

	# Creating the data.frame of combinations to test and adjusting the values of htrans and vtrans
	to_test <- expand.grid(theta = theta_vect, htrans = htrans_vect / fact, vtrans = vtrans_vect / fact)

	# Pre-processing the data for assess_transform_cpp
	optim_coords <- terra::xyFromCell(y, 1:(terra::ncell(y)))
	center <- c(terra::xmax(y) - terra::xmin(y), terra::ymax(y) - terra::ymin(y)) / 2
	optim_coords <- cbind(optim_coords, rep(1, nrow(optim_coords)))
	optim_coords[, 1] <- optim_coords[, 1] - center[1]
	optim_coords[, 2] <- optim_coords[, 2] - center[2]

	# Initializing a vector for the correlations
	correlations <- parallel::mclapply(1:nrow(to_test), function(i, x, y, params, min_overlap) {
						   assess_transform_cpp(xval = terra::values(x),
									yval = terra::values(y),
									coords = optim_coords,
									theta = params[i, 1] * pi / 180,
									htrans = params[i, 2] + center[1],
									vtrans = params[i, 3] + center[2],
									nrows = nrow(x), ncols = ncol(x),
									extent = unlist(as.list(terra::ext(x))),
									min_overlap = min_overlap)
						       },
						   x = x, y = y, params = to_test, min_overlap = min_overlap,
						   mc.cores = cores)

	# Extracting the combination of parameters that yielded the best results
	output <- to_test[which.max(unlist(correlations)), ]
	output$htrans <- output$htrans * fact
	output$vtrans <- output$vtrans * fact

	list(param = as.numeric(output), value = max(unlist(correlations), na.rm = TRUE))
}

#' Add overlap transform parameters to flight metadata
#'
#' This function computes transform parameters for all images in a dataset to
#' allow overlap correction in downstream analyses. These parameters, as well
#' as the resulting correlation between successive images, are added as columns
#' to the input metadata. This function is basically a convenient wrapper for
#' the workhorse function \code{\link{optimize_transform}} and the helper
#' \code{\link{add_tparams}}. Users who want to understand how these functions
#' work or gain finer control over the analysis should read their respective
#' documentation.
#'
#' @param metadata A data.frame of metadata on a flight, such as read by
#' \code{\link{read_metadata}}.
#' @param theta_guess A logical value indicating whether initial guesses for
#' theta (transform rotation angle) should be guessed from the metadata. This
#' requires the column 'GimbalYawDegree' (which specifies the yaw angle of the
#' gimbal when the picture was taken) to be found in the metadata. If
#' theta_guess is FALSE, then the initial guess is equal to the value of
#' theta1, if not NULL. If theta_guess is FALSE and theta1 is NULL, then theta1
#' is set to 0 with a warning.
#' @inheritParams optimize_transform
#'
#' @return A metadata data.frame similar to the input one, but with added
#' columns "theta", "htrans", "vtrans" and "corr" related to the coordinate
#' transformation.
#'
#' @export
#' @examples
#' NULL
compute_overlaps <- function(metadata, theta_guess = TRUE,
			     theta1 = NULL, htrans1 = 0, vtrans1 = 0,
			     theta_range = 5, theta_length = 9,
			     htrans_range = 20, htrans_length = 5,
			     vtrans_range = 20, vtrans_length = 5,
			     min_overlap = 100,
			     fact = 1, cores = 1, min_cor = 0.8, maxiter = 2,
			     method = "Nelder-Mead", reltol = 10^-5,
			     verbose = TRUE, trace = verbose) {

	# Checking if theta_guess is set to TRUE
	if(theta_guess) {
		# Some sanity checks
		if(!is.null(theta1)) warning("theta_guess is set to TRUE, theta1 value ignored")

		# Checking that the values needed to compute the guesses are found in the metadata
		if(! "GimbalYawDegree" %in% colnames(metadata)) {
			stop("theta values cannot be guessed from metadata if column 'GimbalYawDegree' is not provided")
		}
	} else {
		if(is.null(theta1)) {
			warning("theta_guess set to FALSE and theta1 = NULL: theta1 value set to 0")
			theta1 <- 0
		}
	}

	# This list will store the output of optimize_transform for each image pair
	tparams <- list()

	# Looping over all metadata rows but the first (there is no transform for the first picture)
	for(i in 2:nrow(metadata)) {
		# Computing the initial guess for theta
		if(theta_guess) theta1 <- -(metadata[i, "GimbalYawDegree"] - metadata[i - 1, "GimbalYawDegree"])

		# Filling the tparams list with the parameter estimates
		tparams[[i - 1]] <- optimize_transform(x = terra::rast(metadata[i - 1, "SourceFile"]),
						       y = terra::rast(metadata[i, "SourceFile"]),
						       theta1 = theta1, theta_range = theta_range, theta_length = theta_length,
						       htrans1 = htrans1, htrans_range = htrans_range, htrans_length = htrans_length,
						       vtrans1 = vtrans1, vtrans_range = vtrans_range, vtrans_length = vtrans_length,
						       min_overlap = min_overlap, min_cor = min_cor,
						       fact = fact, cores = cores, maxiter = maxiter,
						       method = method, reltol = reltol,
						       verbose = verbose, trace = trace)
	}

	# Add the parameters to the metadata
	metadata <- add_tparams(metadata, tparams)

	return(metadata)
}
