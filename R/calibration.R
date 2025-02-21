#' Load panel temperature data
#'
#' At the moment, this function is specific to output produced by a CR1000X
#' data logger in a very specific format. Support for other data formats
#' could be added upon request. Otherwise, temperature data simply needs
#' to be provided as a data.frame with "time" and "temp" columns.
#'
#' @param filename A character, the name of the file.
#' @param tz A character which can be interpreted as a valid timezone.
#' Must be supplied such that the proper time zone is used when matching
#' the temperature data and images.
#'
#' @return A data.frame of temperature values with the time column set
#'         to the proper data type.
#'
#' @export
#' @examples
#' NULL
read_temp <- function(filename, tz) {
	tempdata <- utils::read.table(filename, sep = ",", skip = 4, na.strings = c("NA", "NAN"))
	colnames(tempdata) <- c("time", "record", "battery", "devtemp", "temp")
	tempdata$time <- as.POSIXct(tempdata$time, format = "%Y-%m-%d %H:%M:%S", tz = tz)
	tempdata
}

#' Plot panel temperature data over a given time range
#'
#' A function used to simplify the plotting of temperature data for the three
#' reference panels over a given time range. Since this function uses base R
#' graphics functionality, further lines or points can be drawn on top of the
#' graph produced by this one after the function returns.
#'
#' @param tempdata A list with three elements named "black", "gray" and "white"
#' for each of the three panel colors. Each element is a data.frame of
#' temperature data as read with \code{\link{read_temp}}.
#' @param xrange A pair of dates identifying the minimum and maximum values used
#' for the x-axis. If NULL, the whole input range is plotted.
#' @param at A date object indicating the location of the x-axis labels.  If
#' NULL, the default is used.
#' @param main A character. The title of the plot.
#' @param xlab A character string. The title of the x-axis.
#' @param ylab A character string. The title of the y-axis.
#' @param lcol A named ("black", "gray", "white") character vector indicating
#' the colors to use for plotting the temperature values of each panel.
#' @param ... Further graphical arguments passed to the drawing functions.
#'
#' @return NULL, invisibly. This function is invoked for its plotting
#' side-effect.
#'
#' @export
#' @examples
#' NULL
plot_temp <- function(tempdata, xrange = NULL, at = NULL, main = NULL,
		      xlab = "Time", ylab = "Temperature (\u00b0C)",
		      lcol = c(black = "black", gray = "gray", white = "blue"),
		      ...) {
	
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
	plot(1, type = "n", xlim = as.numeric(xlim), ylim = ylim,
	     xlab = xlab, ylab = ylab, main = main,
	     xaxt = "n", ...)
	
	# Adding lines for each of the panels
	# By default the white panel is plotted in skyblue (this should be allowed to vary as a parameter)
	for(i in c("black", "gray", "white")) graphics::lines(tempdata[[i]]$time, tempdata[[i]]$temp, col = lcol[i], ...)

	# Adding an axis for the time
	if(is.null(at)) {
		graphics::axis.POSIXct(1, tempdata$black$time, ...)
	} else {
		graphics::axis.POSIXct(1, tempdata$black$time, at = at, ...)
	}

	invisible(NULL)
}

#' Extract the temperature at the nearest time point to a thermal image
#'
#' This function can be used to populate the columns in a flight metadata
#' data.frame with the temperatures associated to a given panel as measured
#' by a datalogger.
#'
#' @param metadata A data.frame containing metadata on a set of thermal images,
#' as returned by \code{\link{read_metadata}}. Must minimally contain a column
#' called "DateTimeOriginal" which indicates when the picture was taken.
#' @param temperature A data.frame of panel temperature data. Should contain a
#' column called "time" to allow matching the time stamps of both datasets and
#' a "temp" column for the temperature.
#' @param tolerance A difftime object of length 1 indicating the maximum time
#' difference acceptable between a picture and a temperature measurement to allow
#' both values to be matched.
#'
#' @return A numeric vector containing the temperature at the nearest time
#' point for every row in the metadata input.
#'
#' @export
#' @examples
#' NULL
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
	if(any(abs(metadata$DateTimeOriginal - temperature[indices, "time"]) > tolerance)) stop("No temperature available given tolerance value.")

	temperature[indices, "temp"]
}

#' Add reference panel temperature readings to a metadata dataset
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
#' difference acceptable between a picture and a temperature measurement to allow
#' both values to be matched.
#'
#' @return A metadata data.frame similar to the input one, but with added
#' columns containing the temperatures of each reference surface at the time
#' points when each picture was taken.
#'
#' @export
#' @examples
#' NULL
add_temp_metadata <- function(metadata, temperature_list, tolerance = as.difftime(10, units = "secs")) {

	# Adding the temperature of each of the panels to the metadata
	for(i in names(temperature_list)) {
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
#' temperatures of reference surfaces panels (see the argument \code{columns})
#' at the moment when the pictures were taken.
#' @param polygons A list of sf polygons representing the location of reference
#' surfaces. A "color" column should be used to identify each of the reference
#' surfaces Each element of the list must be named according to the row names
#' of the metadata.
#' @param columns A character vector of columns in the metadata that represent
#' names of reference surfaces which are found in the column "color" of the
#' \code{polygons} argument.
#' @param ncores An integer indicating the number of cores to use. Passed to
#' \code{\link[parallel]{mclapply}}. Defaults to 1 (no parallel processing).
#'
#' @return A list of data.frames suitable for plotting and modelling with one
#' row per pixel and the following columns:
#' ID: the name of the reference surface
#' thermal: the value of the pixel
#' temp: the temperature of the pixel (fixed for a given panel color)
#' Each element of the list returned is named after the ID of the picture,
#' which also matches the row names of the metadata.
#' 
#' @export
#' @examples
#' NULL
join_thermal <- function(metadata, polygons, columns, ncores = 1) {

	# Ensuring that the metadata contains temperature data
	if(!all(columns %in% colnames(metadata))) {
		stop("metadata input does not contain the columns ", paste(columns, collapse = " "))
	}

	# Also ensuring that all polygons correspond to rows in the metadata
	if(!all(names(polygons) %in% rownames(metadata))) {
		stop("All input polygon names must correspond to row names in metadata.")
	}

	# Extracting the values of the thermal pixels based on the polygons
	output <- parallel::mclapply(names(polygons), function(i, polygons, mdata, columns){
					     # Checking that the colors of the polygons match the columns
					     if(!all(polygons[[i]]$color %in% columns)) {
						     stop("The 'color' column of the polygons must match the 'columns' arguments.")
					     }

					     # Extracting the values from the thermal raster and naming the column
					     irast <- suppressWarnings(terra::rast(mdata[i, "SourceFile"]))
					     pix_values <- terra::extract(irast, polygons[[i]])
					     names(pix_values)[2] <- "thermal"

					     # Using the color as ID and using this ID to extract temperature values
					     pix_values$ID <- polygons[[i]]$color[pix_values$ID]
					     temperature_lookup <- as.numeric(mdata[i, columns])
					     names(temperature_lookup) <- columns
					     pix_values$temp <- temperature_lookup[pix_values$ID]

					     # A sanity check before returning
					     if(any(!stats::complete.cases(pix_values))) stop("NA values not allowed in temperature or thermal pixels")
					     pix_values
	     }, polygons = polygons, mdata = metadata, columns = columns, mc.cores = ncores)

	# Naming the list elements after the polygon IDs
	names(output) <- names(polygons)
	return(output)
}

#' Create linear models of temperature as a function of thermal digital numbers
#'
#' This function creates thermal models for each dataset of temperature and
#' thermal digital numbers provided. This may be a subset of the provided
#' metadata (if for example some pictures did not include reference panels),
#' but all elements of the pixel_values list must be found within the input
#' metadata. Computation of thermal models from multiple pixel values at once
#' could eventually be supported, but at the moment every panel is summarized
#' into a single value which is used for fitting the linear models.
#'
#' @param metadata A data.frame of metadata on a given flight, such
#' as returned by \code{\link{read_metadata}}. Row names of this data.frame
#' must correspond to the names of the elements in pixel_values.
#' @param pixel_values A list of data.frames linking thermal pixel values to
#' temperature, as returned by \code{\link{join_thermal}}. The names of the
#' elements in this list must correspond to the row names of the metadata input.
#' @param summary_functions A named list with functions used to summarize the
#' values for each given reference surface. The names of the list elements
#' should correspond to the IDs of the reference surfaces. The functions should
#' return a single value from an input vector with multiple values.
#' @param surfaces A character vector with the names of the reference
#' surfaces to use for computing linear models. Not all available surfaces
#' may be used if any should be excluded from computation.
#'
#' @return A list of two elements, the first one containing the metadata
#' updated with critical model information (the pixel values used to fit
#' the model and the model estimates) , and the second element containing a
#' list of data used to fit the linear models and the models themselves.
#' The metadata rows for which no model was fitted will have model-related
#' values set to NA.
#'
#' @export
#' @examples
#' NULL
thermal_lm <- function(metadata, pixel_values, summary_functions, surfaces) {

	# We check that all pixel_values names are in the metadata
	if(!all(names(pixel_values) %in% rownames(metadata))) {
		stop("The names of pixel_values must match row names in metadata.")
	}

	# Checking that all the names of the surfaces to use are in the summary_functions list
	if(!all(surfaces %in% names(summary_functions))) {
		stop("You need to provide a summary function for each reference surface.")
	}

	# Computing model data and the model itself for all elements in pixel_values
	models <- lapply(names(pixel_values), function(i, pixels) {

				 # Extracting the pixels for this picture
				 i_pixels <- pixels[[i]]

				 # Initializing the model data  and naming the rows according to panel color
				 model_data <- data.frame(ID = unique(i_pixels$ID),
							  pixel = NA,
							  temp = NA)

				 rownames(model_data) <- model_data$ID

				 # Looping over the reference surfaces
				 for(i in rownames(model_data)) {
					 # Summarizing the pixel values based on the provided functions
					 model_data[i, "pixel"] <- summary_functions[[i]](i_pixels[i_pixels$ID == i, "thermal"])
					 stopifnot(length(i_temp <- unique(i_pixels[i_pixels$ID == i, "temp"])) == 1)
					 model_data[i, "temp"] <- i_temp
				 }

				 # Keeping only the surfaces that we do want to use
				 model_data <- model_data[rownames(model_data) %in% surfaces, ]

				 # Computing the linear model based on this data
				 model <- stats::lm(temp ~ pixel, data = model_data)

				 # We return a list with the data used to fit the model and the model itself
				 list(data = model_data, model = model)
			  }, pixels = pixel_values)

	# Naming the elements of the model list
	names(models) <- names(pixel_values)

	# Updating the metadata with the model data and model parameters
	for(i in surfaces) metadata[[paste0(i, "pix")]] <- NA
	metadata$intercept <- NA
	metadata$slope <- NA

	# Filling in the values only for the metadata rows for which we have panels
	for(i in names(models)) {
		for(j in surfaces) metadata[i, paste0(j, "pix")] <- models[[i]]$data[j, "pixel"]
		metadata[i, "intercept"] <- stats::coef(models[[i]]$model)[1]
		metadata[i, "slope"] <- stats::coef(models[[i]]$model)[2]
	}

	# Returning a list with the metadata and the models
	list(metadata = metadata, models = models)
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
#' digital number value based on linear model parameters that have been
#' fit at various timepoints during a flight.
#'
#' @param metadata A data.frame of metadata on a given flight augmented
#' with slope and intercept parameters of the relationship between digital
#' number (DN) and temperature, such as returned in the metadata object
#' of \code{\link{thermal_lm}}.
#' @param dn_value A numeric. The fixed digital number (DN) value for which
#' to generate predictions.
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


#' Plot panel temperature as a function of thermal digital numbers
#'
#' This function can be used to assess the quality of a linear model relating
#' pixel values to temperature by representing the raw data, summarized data,
#' and model output in the same plot.
#'
#' @param pixel_values A data.frame linking thermal pixel values to
#' temperature, as returned by \code{\link{join_thermal}}.
#' @param lcol A named ("black", "gray", "white") character vector indicating
#' the colors to use for plotting the pixel values of each panel.
#' @param model A list containing data used for fitting a linear model that
#' links the temperature to the pixel values, and the linear model itself, used
#' for plotting the regression line and values used for fitting the model. Such
#' objects are returned by \code{\link{thermal_lm}}.  If NULL (the default),
#' then the output of the model is not plotted.
#' @param col.model.points The color(s) to use for the points used to compute
#' the model.
#' @param cex.model.points The cex parameter to use for the points used to
#' compute the model.
#' @param cex.line A single numeric value. The cex parameter for the dotted
#' line showing the values predicted by the linear model.
#' @param ... Other graphical parametres passed to the main plotting function
#'
#' @return NULL, invisibly. This function is invoked for its plotting
#' side-effect.
#'
#' @export
#' @examples
#' NULL
plot_pixtemp <- function(pixel_values,
			 lcol = c(black = "black", gray = "gray", white = "blue"),
			 model = NULL,
			 col.model.points = "red", cex.model.points = 2,
			 cex.line = 1,
			 ...) {


	# Plotting a scatterplot of the pixel/temperature values
	plot(pixel_values$thermal, pixel_values$temp, col = lcol[pixel_values$ID], ...)

	# Plotting a linear model of temperature as a function of thermal values if requested
	if(!is.null(model)) {
		graphics::abline(reg = model$model, lty = 3, cex = cex.line)
		# Also plotting the points used for the model
		graphics::points(x = model$data$pixel,
				 y = model$data$temp,
				 col = col.model.points,
				 cex = cex.model.points)
	}

	return(invisible(NULL))
}

