#' Compute the mean value of all thermal pictures in a dataset
#'
#' @param metadata A data.frame of metadata on a set of thermal pictures.
#' @param ncores An integer. The number of processes to launch (passed to \code{\link[parallel]{mclapply}}).
#'
#' @return a numeric vector of mean thermal values for all files in the dataset.
#' 
#' @export
#' @examples
#' NULL
thermal_mean <- function(metadata, ncores = 1) {
	# Sanity check
	stopifnot("SourceFile" %in% names(metadata))

	# Computing the means in parallel
	output <- parallel::mclapply(metadata$SourceFile, function(x) terra::global(suppressWarnings(terra::rast(x)), mean)$mean,
				     mc.cores = ncores)
	
	unlist(output)
}

#' Extract the vignetting pattern of a thermal dataset
#'
#' @param metadata A data.frame of metadata on a thermal picture dataset
#' 
#' @return A terra raster with the mean value of each pixel over the whole flight
#'
#' @export
#' @examples
#' NULL
compute_vignetting <- function(metadata) {

	# Extracting the file names from the metadata
	stopifnot("SourceFile" %in% names(metadata))
	filenames <- metadata$SourceFile

	output <- suppressWarnings(terra::rast(filenames[1]))

	if(length(filenames) == 1) return(output)

	for(i in 2:length(filenames)) {
		output <- output + suppressWarnings(terra::rast(filenames[i]))
	}

	output / length(filenames)
}

#' Add transformation parameters to a thermal picture dataset
#'
#' @param metadata A data.frame of metadata on a set of thermal pictures.
#' @param tparams A set of transformation parameters, as determined by
#' \code{\link{optimize_transform}}.
#'
#' @return A data.frame of metadata similar to the input one, with
#' added columns "theta", "htrans" and "vtrans" describing the transformation
#' parameters to go from one image to the previous one. Therefore, the parameters
#' for the very first image in the dataset are NA.
#'
#' @export
#' @examples
#' NULL
add_tparams <- function(metadata, tparams) {
	# A few sanity checks
	if(!(length(tparams) == nrow(metadata) - 1)) stop("Length of tparams object is not the number of metadata rows - 1")

	# Coercing the list of parameters to a data.frame and naming the columns
	params <- as.data.frame(t(sapply(tparams, function(x) x$optim$par)))
	colnames(params) <- c("theta", "htrans", "vtrans")

	# Appending the parameters to the input metadata
	metadata <- cbind(metadata, rbind(NA, params))
}

#' Plots the location of pictures from flight metadata
#'
#' This function can be used to get a first glance of the layout of a flight.
#' The color parameter can be used to specify different colors for points
#' depending on characteristics of the metadata.
#'
#' @param metadata a data.frame of metadata on a flight, such as returned by
#' \code{\link{read_metadata}}.  Must minimally contain columns GPSLongitude and
#' GPSLatitude for longitude and latitude coordinates in the WGS84 datum.
#' @param new_crs A coordinate reference system (such as returned by
#' \code{\link[sf]{st_crs}} to project the points to before display. If NULL
#' (the default), then data is displayed in the original coordinate system.
#' @param panel_pos An sf object representing the position of the panels (or
#' any point of interest), which will be transformed if new_crs is not NULL.
#' If NULL, (the default), then panels are not plotted.
#' @param panel_color The color to use for the panels.
#' @param color A vector used to specify the colors of the points for the
#' position of the images. If NULL (the default), then all points are colored
#' black.
#' @param title A single character string. The title of the plot.
#' @param cex.axis A character expansion value applied to the axis labels.
#' @param cex.points A character expansion value applied to the points.
#'
#' @return NULL, invisibly. This function is invoked for its plotting
#' side-effect.
#'
#' @export
#' @examples
#' NULL
position_plot <- function(metadata, new_crs = NULL, panel_pos = NULL, panel_color = "black", color = NULL, title = NULL, cex.axis = 1, cex.points = 1) {
	# Coercing to a spatial object
	metadata <- sf::st_as_sf(metadata,
				 coords = c("GPSLongitude", "GPSLatitude"),
				 crs = sf::st_crs("WGS84"),
				 dim = "XY")

	# Transforming the CRS if needed
	if(!is.null(new_crs)) {
		metadata <- sf::st_transform(metadata, crs = new_crs)

		if(!is.null(panel_pos)) panel_pos <- sf::st_transform(panel_pos, crs = new_crs)
	}

	# Setting up the viewports
	grid::pushViewport(grid::plotViewport())
	grid::pushViewport(grid::dataViewport(xData = sf::st_coordinates(metadata)[, "X"], yData = sf::st_coordinates(metadata)[, "Y"]))

	# Plotting the data and axes
	grid::grid.rect()
	grid::grid.points(sf::st_coordinates(metadata)[, "X"],
			  sf::st_coordinates(metadata)[, "Y"],
			  gp = grid::gpar(col = if(is.null(color)) "black" else color, cex = cex.points),
			  default.units = "native")

	# Adding the position of the panels (if provided)
	if(!is.null(panel_pos)) {
		grid::grid.draw(sf::st_as_grob(sf::st_geometry(panel_pos),
					       gp = grid::gpar(col = panel_color),
					       pch = 16))
	}

	grid::grid.xaxis(gp = grid::gpar(cex = cex.axis))
	grid::grid.yaxis(gp = grid::gpar(cex = cex.axis))

	# Adding a title (if requested)
	if(!is.null(title)) {
		grid::grid.text(label = title, y = grid::unit(1, "npc") + grid::unit(1, "lines"))
	}

	# Going back to the top viewport
	grid::upViewport(2)

	invisible(NULL)
}

#' Plot the fit of a drift model
#'
#' This function fits a statistical model of drift to flight metadata
#' under the specified parameters. It is to be used interactively so
#' as to visualize the impact of various methodological choices on
#' the model that will be used for thermal drift correction.
#'
#' @param metadata A data.frame of metadata on a given flight. Must minimally
#' contain the columns "mean" and "DateTimeOriginal"
#' @param nuc_threshold The threshold to use for declaring non-uniformity
#' correction events. Will be used to split the dataset into different segments
#' on which distinct models will be fitted
#' @param method A character. The type of statistical model to use for fitting.
#' One of "lm" or "spline".
#'
#' @return A data.frame similar to the input data with an added column
#' "fitted" for the fitted values. The value is returned invisibly as this
#' functon is mostly invoked for its plotting side-effect.
#'
#' @export
#' @examples
#' NULL
plot_fit <- function(metadata, nuc_threshold, method = c("lm", "spline")) {
	# The heavy lifting is done by fit_model
	models <- fit_model(x = metadata, nuc_threshold = nuc_threshold, model_type = method)

	# We compute the fitted values from those models
	metadata$fitted <- unlist(lapply(models, stats::fitted))

	# Base plotting is sufficient for this use case
	plot(x = metadata$DateTimeOriginal, y = metadata$mean)
	lines(x = metadata$DateTimeOriginal, y = metadata$fitted, col = "blue", lty = 2)

	invisible(metadata)
}
