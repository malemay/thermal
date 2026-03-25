#' Set reference surface coordinates through interactive clicking
#'
#' Panel corners are to be clicked on clockwise.
#'
#' @param x A SpatRaster object to query for coordinates.
#'
#' @return An polygon object of class sfc that represents the coordinates of a
#' surface in a reference coordinate system.
#'
#' @noRd
set_coords <- function(x, nclicks = 4) {

	# Interactively getting the coordinates of the panel
	coords <- terra::click(x, n = nclicks, xy = TRUE)

	# Formatting the coordinates and closing the shape by repeating the first row
	coords <- rbind(coords[, c("x", "y")], coords[1, c("x", "y")])

	# Creating an object of the sfc class from the supplied coordinates
	output <- sf::st_as_sf(sf::st_as_sfc(sf::st_as_binary(sf::st_polygon(list(as.matrix(coords))))))

	# Setting the coordinate reference system (same as the raster)
	sf::st_crs(output) <- sf::st_crs(x)

	output
}

#' Interactively identify reference surface coordinates from a raster
#'
#' This function is meant as a helper to interactively define the coordinates
#' of reference surfaces to be used in radiometric calibration with package
#' thermal. Calling this function will zoom onto an area of the input image
#' where reference surfaces are located and will query for as many polygons as
#' the number of values in the vector 'ids'. The user will be asked to click in
#' clockwise order on as many vertices as specified by argument 'nclicks' to
#' define the shape of each surface.  The resulting shapes can be written to
#' disk storage using \code{\link[sf]{st_write}}.
#'
#' @param image Either a SpatRaster object loaded using
#' \code{\link[terra]{rast}} or a character value which is interpreted as the
#' path to a file which can be read using \code{\link[terra]{rast}}.
#' @param ids A charactor vector of reference surface identifiers.
#' @param nclicks A numeric value corresponding to the number of clicks to
#' query for each surface. At the moment, the number of clicks is the same for
#' all surfaces in the vector 'ids'. Defaults to 4, assuming a square or
#' rectangular shape.
#'
#' @return An sf object with the coordinates of the reference surfaces,
#' including column 'id' for the identifier of the surface. The value returned
#' by this function can be combined in a named list with other sf objects to
#' use as input to function \code{\link{thermal_lm}}.
#'
#' @export
#' @examples
#' NULL
create_surfaces <- function(image, ids, nclicks = 4) {

	if(!interactive()) stop("create_surfaces can only be called in interactive mode")

	# Reading a picture if provided with a file path
	if(is.character(image)) image <- terra::rast(image, noflip = TRUE)

	# Querying for an area to zoom on the picture
	message("Click on the image twice to zoom on an area of interest (top left and botton right corners)")
	terra::plot(image)
	panel_region <- terra::zoom(image)

	# Initializing an object to return as output
	output_shapes <- list()

	# Looping over the ids
	for(i in ids) {
		# Getting the coordinates for each of the panels
		message("Click clockwise on ", nclicks, " vertices of surface ", i)
		output_shapes[[i]] <- set_coords(terra::crop(image, panel_region), nclicks = nclicks)
	}

	# Joining all the shapes together and assigning the id column
	output_shapes <- do.call("rbind", output_shapes)
	output_shapes$ID <- ids

	plot(output_shapes, add = TRUE)

	return(output_shapes)
}

#' Plot reference surface temperature data over a given time range
#'
#' This function plots the temperature data of reference surfaces over a given
#' time range. Since this function uses base R graphics functionality, further
#' lines or points can be drawn on top of the graph produced by this one after
#' the function returns.
#'
#' @param tempdata A list with a data.frame of temperature data for each of the
#' reference surfaces. Each data.frame must contain the columns "time" for the
#' time stamps in POSIXct format and "temp" for the surface temperature. The
#' names of the list elements must correspond to the names of the reference
#' surfaces.
#' @param xrange A pair of dates identifying the minimum and maximum values
#' used for the x-axis. If NULL, the whole input range is plotted.
#' @param at A date object indicating the location of the x-axis labels.  If
#' NULL, the default is used.
#' @param main A character. The title of the plot.
#' @param xlab A character string. The title of the x-axis.
#' @param ylab A character string. The title of the y-axis.
#' @param lcol A named character vector indicating the colors to use for
#' plotting the temperature values of each reference surface. If NULL (the
#' default), then R's default integer colors are used.
#' @param legend.pos A character keyword specifying the position of the legend
#' in the plot. If NULL (default), then no legend is drawn. See
#' \code{\link[graphics]{legend}} for accepted keyword values.
#' @param lty The line type to use for plotting. Defaults to a solid line.
#' @param ... Further graphical parameters passed to the drawing functions.
#'
#' @return NULL, invisibly. This function is invoked for its plotting
#' side-effect.
#'
#' @export
#' @examples
#' NULL
plot_temp <- function(tempdata, xrange = NULL, at = NULL, main = NULL,
		      xlab = "Time", ylab = "Temperature (\u00b0C)",
		      lcol = NULL, legend.pos = NULL, lty = 1, ...) {
	
	# Subsetting the input data to the range of interest
	if(!is.null(xrange)) {
		stopifnot(length(xrange) == 2)
		stopifnot(xrange[1] < xrange[2])

		tempdata <- lapply(tempdata,
				   function(x, xrange) x[x$time >= xrange[1] & x$time <= xrange[2], ],
				   xrange = xrange)
	}

	# Getting the limits of the x- and y-axes and setting the plotting region accordingly
	xlim <- range(do.call("rbind", tempdata)$time)
	ylim <- range(do.call("rbind", tempdata)$temp, na.rm = TRUE)

	# Creating a blank canvas for the plot
	plot(1, type = "n",
	     xlim = as.numeric(xlim), ylim = ylim,
	     xlab = xlab, ylab = ylab, main = main,
	     xaxt = "n", ...)

	# Setting the colors for each reference surface if they are not provided
	if(is.null(lcol)) {
		lcol <- 1:length(tempdata)
		names(lcol) <- names(tempdata)
	} else {
		stopifnot(length(tempdata) == length(lcol))
		stopifnot(all(names(lcol) %in% names(tempdata)) && all(names(tempdata) %in% names(lcol)))
	}
	
	# Adding lines for each of the surfaces
	for(i in names(tempdata)) graphics::lines(tempdata[[i]]$time, tempdata[[i]]$temp, col = lcol[i], lty = lty, ...)

	# Adding an axis for the time
	if(is.null(at)) {
		graphics::axis.POSIXct(1, tempdata[[1]]$time, ...)
	} else {
		graphics::axis.POSIXct(1, tempdata[[1]]$time, at = at, ...)
	}

	if(!is.null(legend.pos)) {
		graphics::legend(x = legend.pos, legend = names(tempdata), lty = lty, col = lcol[names(tempdata)])
	}

	invisible(NULL)
}

#' Extract the temperature at the nearest time point to a thermal image
#'
#' This function can be used to populate the columns in a flight metadata
#' data.frame with the temperatures associated to a given reference surface.
#'
#' @param metadata A data.frame containing metadata on a set of thermal images,
#' as returned by \code{\link{read_metadata}}. Must minimally contain a column
#' called "DateTimeOriginal" which indicates when the picture was taken.
#' @param temperature A data.frame of reference surface temperature data.
#' Should contain a column called "time" to allow matching the time stamps of
#' both datasets and a "temp" column for the temperature.
#' @param tolerance A difftime object of length 1 indicating the maximum time
#' difference acceptable between a picture and a temperature measurement to
#' allow both values to be matched.
#'
#' @return A numeric vector containing the temperature at the nearest time
#' point for every row in the metadata input.
#'
#' @noRd
extract_temp <- function(metadata, temperature, tolerance = as.difftime(10, units = "secs")) {

	# Creating a vector of matches in the temperature data.frame
	indices <- numeric(nrow(metadata))

	# Subsetting the temperature data to the range in the metadata
	time_range <- range(metadata$DateTimeOriginal) + c(-tolerance, tolerance)
	temperature <- temperature[temperature$time >= time_range[1] & temperature$time <= time_range[2], ]

	# Performing operations on numeric vectors because it is faster
	metadata_times <- as.numeric(metadata$DateTimeOriginal)
	temperature_times <- as.numeric(temperature$time)

	# Getting the minimum match for each row of metadata
	for(i in 1:length(indices)) indices[i] <- which.min(abs(metadata_times[i] - temperature_times))

	# Adding a sanity check to ensure that no time difference is greater than the tolerance
	if(any(abs(metadata$DateTimeOriginal - temperature[indices, "time"]) > tolerance)) stop("No temperature available given tolerance value ", tolerance)

	temperature[indices, "temp"]
}

#' Add reference surface temperature readings to a metadata dataset
#'
#' This is a convenience function that adds columns with the temperature of
#' reference surfaces at given time points to the flight metadata.
#'
#' @param metadata A data.frame containing metadata on a set of thermal images,
#' as returned by \code{\link{read_metadata}}. Must minimally contain a column
#' called "DateTimeOriginal" which indicates when the picture was taken.
#' @param temperature_list A named list of data.frames with reference surface
#' temperature data. The names of the list will be added as columns to the
#' output with their respective temperatures at each time. Each data.frame
#' should contain a column called "time" to allow matching the time stamps of
#' both datasets and a "temp" column for the temperature.
#' @param tolerance A difftime object of length 1 indicating the maximum time
#' difference acceptable between a picture and a temperature measurement to
#' allow both values to be matched.
#'
#' @return A metadata data.frame similar to the input one, but with added
#' columns containing the temperatures of each reference surface at the time
#' points when each picture was taken.
#'
#' @export
#' @examples
#' NULL
add_temp_metadata <- function(metadata, temperature_list, tolerance = as.difftime(10, units = "secs")) {

	# Adding the temperature of each of the surfaces to the metadata
	for(i in names(temperature_list)) {
		# Check that there is a "time" column in the temperature data
		if(! "time" %in% colnames(temperature_list[[i]]) || !inherits(temperature_list[[i]]$time, "POSIXct")) {
			stop("Temperature data must contain POSIX")
		}

		# Warn if any of the surface names already correspond to columns in metadata
		if(i %in% colnames(metadata)) warning("Column ", i, " already present in metadata. Overwriting its contents.")

		# Add the temperature from that surface to the metadata
		metadata[[i]] <- extract_temp(metadata, temperature_list[[i]], tolerance = tolerance)
	}

	# We return the modified data.frame
	metadata
}

#' Join thermal values to surface temperatures
#'
#' This generates data.frames relating thermal pixel values to temperature data
#' for the reference surfaces for each row in a metadata input data.frame.
#'
#' @param metadata A data.frame with metadata on thermal pictures taken during
#' a flight. The data.frames must have row names that correspond to names in
#' the polygons list to link the two. Must also have columns with the
#' temperatures of reference surfaces (see the argument \code{columns})
#' at the moment when the pictures were taken.
#' @param polygons A list of sf polygons representing the location of reference
#' surfaces. An "ID" column should be used to identify each of the reference
#' surfaces Each element of the list must be named according to the row names
#' of the metadata.
#' @param surfaces A character vector of columns in the metadata that represent
#' names of reference surfaces which are found in the column "ID" of the
#' \code{polygons} argument.
#' @param ncores An integer indicating the number of cores to use. Passed to
#' \code{\link[parallel]{mclapply}}. Defaults to 1 (no parallel processing).
#'
#' @return A list of data.frames suitable for plotting and modelling with one
#' row per pixel and the following columns:
#' \itemize{
#'   \item ID: the name of the reference surface
#'   \item thermal: the value of the pixel
#'   \item temp: the temperature of the pixel (fixed for a given surface ID)
#'   }
#' Each element of the list returned is named after the ID of the picture,
#' which also matches the row names of the metadata.
#' 
#' @noRd
join_thermal <- function(metadata, polygons, surfaces, ncores = 1) {

	# We verify that the geometry of all polygons is valid and therefore force loading the sf namespace
	stopifnot(all(sapply(polygons, sf::st_is_valid)))

	# Extracting the values of the thermal pixels based on the polygons
	output <- parallel::mclapply(names(polygons), function(i, polygons, mdata, surfaces){
					     # We subset the polygons to the requested surfaces
					     i_poly <- polygons[[i]][polygons[[i]]$ID %in% surfaces, ]

					     # Checking that the IDs of the polygons match the surfaces
					     if(!any(i_poly$ID %in% surfaces)) {
						     warning("None of the polygon IDs in polygon ", i, " correspond to the requested surfaces.")
						     return(data.frame(ID = character(), thermal = numeric(), temp = numeric()))
					     }

					     # Extracting the values from the thermal raster and naming the column
					     stopifnot(file.exists(mdata[i, "SourceFile"]))
					     irast <- suppressWarnings(terra::rast(mdata[i, "SourceFile"], noflip = TRUE))
					     pix_values <- terra::extract(irast, i_poly)
					     names(pix_values)[2] <- "thermal"

					     # Using the ID to extract temperature values
					     # by using surface IDs to index into temperature values
					     pix_values$ID <- i_poly$ID[pix_values$ID]
					     temperature_lookup <- as.numeric(mdata[i, surfaces])
					     names(temperature_lookup) <- surfaces
					     pix_values$temp <- temperature_lookup[pix_values$ID]

					     # A sanity check before returning
					     if(any(!stats::complete.cases(pix_values))) warning("NA values found in pixel data for polygon ", i)

					     # Returning the data.frame of thermal pixel/temperature values
					     return(pix_values)
	     }, polygons = polygons, mdata = metadata, surfaces = surfaces, mc.cores = ncores)

	# Naming the list elements after the polygon IDs
	names(output) <- names(polygons)

	return(output)
}

#' Create linear models of temperature as a function of thermal digital numbers
#'
#' This function creates thermal models for each dataset of temperature and
#' thermal digital numbers provided. This may be a subset of the provided
#' metadata (if for example some pictures did not include reference surfaces),
#' but all elements of the pixel_values list must be found within the input
#' metadata. Computation of thermal models from multiple pixel values at once
#' could eventually be supported, but at the moment every surface is summarized
#' into a single value which is used for fitting the linear models.
#'
#' @param metadata A data.frame of metadata on a given flight, such as returned
#' by \code{\link{read_metadata}}. Row names of this data.frame must correspond
#' to the names of the elements in pixel_values.
#' @param polygons A list of sf polygons representing the location of reference
#' surfaces. An "ID" column should be used to identify each of the reference
#' surfaces Each element of the list must be named according to the row names
#' of the metadata.
#' @param summary_functions A single function or a named list with functions
#' used to summarize the pixel values for each given reference surface. If a
#' list, the names of the list elements should correspond to the IDs of the
#' reference surfaces. The functions used should return a single value
#' summarizing the distribution of pixel values for each surface (for example,
#' the mean or median) or a numeric vector of values that have been filtered or
#' even left as is (if the function \code{\link[base]{identity}} is used, for
#' example). By default, the function \code{\link{max_density}} is used.
#' @param surfaces A character vector with the names of the reference surfaces
#' to use for computing linear models. Not all available surfaces may be used
#' if any should be excluded from computation.
#' @param ncores An integer indicating the number of cores to use. Passed to
#' \code{\link[parallel]{mclapply}}. Defaults to 1 (no parallel processing).
#'
#' @return A list of three elements: 
#' \itemize{
#' \item metadata: a data.frame of flight metadata updated with critical model
#' information ('slope' and 'intercept' for the linear models' slope and
#' intercept, respectively). If the summary functions used return a single
#' value summarizing the distribution for a given surface, columns following
#' the template 'surface name + pix' will be added with the summarized value
#' for that surface.
#' \item models: a list of the fitted linear models of radiometric calibration.
#' \item pixels: a list of data.frames (one for each element in polygons) with
#' the thermal pixels and temperature corresponding to each surface.
#' }
#'
#' @export
#' @examples
#' NULL
thermal_lm <- function(metadata, polygons, surfaces, summary_functions = max_density, ncores = 1) {

	# Ensuring that the metadata contains temperature data
	if(!all(surfaces %in% colnames(metadata))) {
		stop("metadata input does not contain the columns: ", paste(surfaces, collapse = ","))
	}

	# Ensuring that all polygons correspond to rows in the metadata
	if(!all(names(polygons) %in% rownames(metadata))) {
		stop("All input polygon names must correspond to row names in metadata.")
	}

	# Checking that summary_functions is accurately specified
	if(length(summary_functions) == 1 && is.function(summary_functions)) {
		# If this is a single function then we create a list in which every surface uses this function
		# This allows for a unified interface to computation no matter how summary_functions is specified
		temp_func <- summary_functions
		summary_functions <- list()
		for(i in surfaces) summary_functions[[i]] <- temp_func
	} else if(is.list(summary_functions) && all(sapply(summary_functions, is.function))) {
		if(!all(surfaces %in% names(summary_functions))) stop("You need to provide a summary function for each reference surface.")
	} else {
		stop("summary_functions must be either a single function or a named list with names corresponding to reference surfaces")
	}

	# Extracting data.frames relating thermal pixel values to temperature for each polygon
	pixel_values <- join_thermal(metadata, polygons, surfaces, ncores = ncores)

	# Computing model data and the model itself for all elements in pixel_values
	models <- parallel::mclapply(names(pixel_values), function(i, pixels, surfaces, functions) {

				 # Extracting the pixels for this picture
				 i_pixels <- pixels[[i]]
				 
				 # We only keep the surfaces that are explicitly asked for
				 i_pixels <- i_pixels[i_pixels$ID %in% surfaces, ]

				 # Formatting the data for fitting the model by summarizing the pixels from each surface
				 model_data <- lapply(split(i_pixels, i_pixels$ID), function(x, functions) {
							      # Extracting the ID and temperature for that image and surface
							      id <- unique(x$ID)
							      temp <- unique(x$temp)
							      stopifnot(length(id) == 1 && length(temp) == 1)

							      # Formatting the data.frame by summarizing the data
							      data.frame(ID = id,
									 pixel = functions[[id]](x$thermal),
									 temp = temp)
							  }, functions = functions)

				 # Combining the data from all surfaces into a single data.frame
				 model_data <- do.call("rbind", model_data)

				 # If the function has summarized the data into a single pixel value, then we assign row names
				 if(nrow(model_data) == length(surfaces) && all(surfaces %in% model_data$ID)) {
					 rownames(model_data) <- model_data$ID
				 }

				 # Computing the linear model based on this data
				 model <- stats::lm(temp ~ pixel, data = model_data)

				 # Returning the fitted model
				 model
			  }, pixels = pixel_values, surfaces = surfaces, functions = summary_functions, mc.cores = ncores)

	# Naming the elements of the model list
	names(models) <- names(pixel_values)

	# Updating the metadata with the model data and model parameters
	for(i in surfaces) metadata[[paste0(i, "pix")]] <- NA
	metadata$intercept <- NA
	metadata$slope <- NA

	# Filling in the values only for the metadata rows for which we have surfaces
	for(i in names(models)) {

		if(nrow(models[[i]]$model) == length(surfaces) && all(surfaces %in% rownames(models[[i]]$model))) {
			for(j in surfaces) {
				if(j %in% rownames(models[[i]]$model)) metadata[i, paste0(j, "pix")] <- models[[i]]$model[j, "pixel"]
			}
		}

		metadata[i, "intercept"] <- stats::coef(models[[i]])[1]
		metadata[i, "slope"] <- stats::coef(models[[i]])[2]
	}

	# Returning a list with the metadata and the models
	list(metadata = metadata, models = models, pixels = pixel_values)
}

#' Create a global radiometric calibration model from all images used in a flight
#'
#' Whereas \code{\link{thermal_lm}} creates one radiometric calibration model
#' for each image in a flight, this function combines data from all images in
#' the flight to create a global radiometric calibration model from all the
#' images.
#'
#' @param model_data A list of flight metadata, individual models for each
#' image, and pixel values, such as returned by \code{\link{thermal_lm}}.
#'
#' @return A list similar to that provided as input, but in which the metadata,
#' models, and pixels have been updated to reflect the fact that the model now
#' globally applies to all images in the flight. The names of the model and pixel
#' values elements will be called "global". This list is suitable for input to
#' \code{\link{plot_pixtemp}}.
#'
#' @export
#' @examples
#' NULL
global_lm <- function(model_data) {
	# Checking that the inputs are appropriate
	if(length(model_data) != 3 || ! identical(names(model_data), c("metadata", "models", "pixels"))) {
		stop("The input to global_lm must be an object returned by function thermal_lm")
	}

	# Joining the data used for all models in a single data.frame
	model_df <- lapply(model_data$models, function(x) x$model)
	model_df <- do.call("rbind", model_df)

	# Recomputing the model using that data
	new_model <- stats::lm(temp ~ pixel, data = model_df)

	# Setting the intercept and slope in the metadata to the new values
	model_data$metadata$intercept <- stats::coef(new_model)[1]
	model_data$metadata$slope <- stats::coef(new_model)[2]

	# Replacing the models element by the global model
	model_data$models <- list(global = new_model)

	# Also replacing the pixel values by a single data.frame with all value
	model_data$pixels <- list(global = do.call("rbind", model_data$pixels))

	return(model_data)
}

#' Identify the mode from the density of a distribution
#'
#' This function is mainly meant to be used in the argument summary_functions
#' of \code{\link{thermal_lm}} to obtain a single pixel value from a large
#' number of observed pixels.
#'
#' @param x A numeric vector.
#'
#' @return A single numeric value that represents the value of x at which the
#' density of the distribution is maximized.
#' 
#' @export
#' @examples
#' NULL
max_density <- function(x) {
	output <- stats::density(x)
	output$x[which.max(output$y)]
}

#' Compute temperature predictions from thermal model parameters
#'
#' This function can be used to compute temperature predictions for a given
#' digital number value based on linear model parameters that have been fit at
#' various timepoints during a flight.
#'
#' @param metadata A data.frame of metadata on a given flight augmented with
#' slope and intercept parameters of the relationship between digital number
#' (DN) and temperature, such as returned in the metadata object of
#' \code{\link{thermal_lm}}.
#' @param dn_value A numeric. The fixed digital number (DN) value for which to
#' generate predictions.
#'
#' @return A numeric vector of predicted temperature values for a fixed DN at
#' each timepoint of the flight based on the linear model parameters.
#'
#' @export
#' @examples
#' NULL
thermal_predict <- function(metadata, dn_value) {
	# Sanity check
	stopifnot(all(c("slope", "intercept") %in% colnames(metadata)))

	# Return a vector of the predictions
	metadata$intercept + metadata$slope * dn_value
}


#' Plot surface temperature as a function of thermal digital numbers
#'
#' This function can be used to assess the quality of a linear model relating
#' pixel values to temperature by representing the raw data, summarized data,
#' and model output in the same plot.
#'
#' @param model_data A list, such as returned by \code{\link{thermal_lm}},
#' containing flight metadata along with radiometric calibration models that
#' were fitted by combining surface temperature data with thermal pixel values.
#' @param id A character string indicating the ID of the picture for which the
#' model data should be displayed. Must correspond to one of the polygon IDs
#' used for extracting pixel values. If NULL (the default), then the data for
#' all pictures is displayed; this may not be what you want if you have a lot
#' of pictures.
#' @param lcol Either a single value that can be interpreted as a color or a
#' named character vector (with names corresponding to surface IDs) indicating
#' the colors to use for plotting the pixel values of each surface. Defaults to
#' "black".
#' @param col.model.points The color(s) to use for the points used to compute
#' the model.
#' @param cex.model.points The cex parameter to use for the points used to
#' compute the model.
#' @param cex.line A single numeric value. The cex parameter for the dotted
#' line showing the values predicted by the linear model.
#' @param ask A logical value indicating whether the function should wait
#' for user input before displaying the next plot. Defaults to FALSE but it can
#' be helpful to set it to TRUE for interactive analysis.
#' @param ... Other graphical parametres passed to the main plotting function
#'
#' @return NULL, invisibly. This function is invoked for its plotting
#' side-effect.
#'
#' @export
#' @examples
#' NULL
plot_pixtemp <- function(model_data, id = NULL, lcol = "black",
			 col.model.points = "red", cex.model.points = 2,
			 cex.line = 1, ask = FALSE,
			 ...) {

	# Setting ID to all images if ID is not specified
	if(is.null(id)) {
		id <- names(model_data$pixels)
	}

	stopifnot(all(id %in% names(model_data$pixels)))

	# Formatting the lcol vector if only one value was provided
	surface_ids <- unique(unlist(lapply(model_data$pixels, function(x) x$ID)))

	if(length(lcol) == 1) {
		lcol <- rep(lcol, length(surface_ids))
		names(lcol) <- surface_ids
	}

	stopifnot(length(lcol) == length(surface_ids) && all(names(surface_ids) %in% names(lcol)))

	# Setting the ask parameter and making sure it returns to original value when exiting function
	opar <- graphics::par(ask = ask)
	on.exit(graphics::par(opar))

	# Looping over all the IDs
	for(i in id) {

		# Extracting the model and pixel values for that image
		model <- model_data$models[[i]]
		pixel_values <- model_data$pixels[[i]]

		# Plotting a scatterplot of the pixel/temperature values
		plot(pixel_values$thermal, pixel_values$temp,
		     xlab = "Pixel value", ylab = "Surface temperature",
		     main = i,
		     col = lcol[pixel_values$ID],
		     ...)

		# Plotting a linear model of temperature as a function of thermal values
		graphics::abline(reg = model, lty = 3, cex = cex.line)

		# Also plotting the points used for the model
		graphics::points(x = model$model$pixel,
				 y = model$model$temp,
				 col = col.model.points,
				 cex = cex.model.points)
	}

	return(invisible(NULL))
}

