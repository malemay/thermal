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
#'
read_temp <- function(filename, tz) {
	tempdata <- read.table(filename, sep = ",", skip = 4, na.strings = c("NA", "NAN"))
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
#' @param xmin A pair of dates identifying the minimum and maximum values used
#' for the x-axis.  If NULL, the whole input range is plotted.
#' @param at A date object indicating the location of the x-axis labels.  If
#' NULL, the default is used.
#' @param main A character. The title of the plot.
#' @param lcol A named ("black", "gray", "white") character vector indicating
#' the colors to use for plotting the temperature values of each panel.
#'
#' @return NULL, invisibly. This function is invoked for its plotting
#' side-effect.
#'
#' @export
#' @examples
#' NULL
#'
plot_temp <- function(tempdata, xrange = NULL, at = NULL, main = NULL,
		      lcol = c(black = "black", gray = "gray", white = "blue")) {
	
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
	     xlab = "Time", ylab = "Temperature (°C)", main = main,
	     xaxt = "n")
	
	# Adding lines for each of the panels
	# By default the white panel is plotted in skyblue (this should be allowed to vary as a parameter)
	for(i in c("black", "gray", "white")) lines(tempdata[[i]]$time, tempdata[[i]]$temp, col = lcol[i])

	# Adding an axis for the time
	if(is.null(at)) {
		axis.POSIXct(1, tempdata$black$time)
	} else {
		axis.POSIXct(1, tempdata$black$time, at = at)
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
#' @param temperature A data.frame of panel temperature data, as returned by
#' the function \code{\link{read_temp}}. Should contain a column called "time"
#' to allow matching the time stamps of both datasets and a "temp" column for
#' the temperature.
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
#'
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
#' This is a convenience function that is a wrapper around
#' \code{\link{extract_temp}} when temperature values for three (black, gray, white)
#' reference panels are available. It allows to automatically add columns with the
#' temperature at given time points to the flight metadata.
#'
#' @inheritParams extract_temp
#' @param temperature_list A named list of temperature data.frames. The names must
#' be "black", "gray", and "white" and the corresponding data.frames contain temperature
#' data on their respective panels.
#'
#' @return A metadata data.frame similar to the input one, but with three
#' added columns ("black", "gray", "white") with the temperatures of each
#' of the panels at the time points when each picture was taken.
#'
#' @export
#' @examples
#' NULL
add_temp_metadata <- function(metadata, temperature_list, tolerance = as.difftime(10, units = "secs")) {
	# Some sanity checks that are conditions for this function
	stopifnot(length(temperature_list) == 3 && identical(sort(names(temperature_list)), c("black", "gray", "white")))

	# Adding the temperature of each of the panels to the metadata
	for(i in sort(names(temperature_list))) {
		metadata[[i]] <- extract_temp(metadata, temperature_list[[i]], tolerance = tolerance)
	}

	# We return the modified data.frame
	metadata
}

#' Join thermal values to panel temperatures
#'
#' To complete
#'
#' @param thermal A raster representing a thermal image
#' @param polygons Polygons representing the location of black, gray and white panels.
#'                 A "color" column should be used to identify each of the three polygons.
#' @param black_temp The temperature of the black panel at the moment when the picture was taken.
#' @param gray_temp The temperature of the gray panel at the moment when the picture was taken.
#' @param white_temp The temperature of the white panel at the moment when the picture was taken.
#'
#' @return A data.frame suitable for plotting and modelling with onw row per pixel
#'         and the following columns:
#'         ID: the color of the panel
#'         thermal: the value of the pixel
#'         temp: the temperature of the pixel (fixed for a given panel color)
#' 
#' @examples
#' NULL
#'
#' @export
join_thermal <- function(thermal, polygons, black_temp, gray_temp, white_temp) {
	# Extracting the values of the thermal pixels based on the polygons
	pixel_values <- extract(thermal, polygons)
	names(pixel_values)[2] <- "thermal"
	pixel_values$ID <- polygons$color[pixel_values$ID]

	# Adding the information on panel temperature
	pixel_values$temp <- NA
	pixel_values[pixel_values$ID == "black", "temp"] <- black_temp
	pixel_values[pixel_values$ID == "gray", "temp"] <- gray_temp
	pixel_values[pixel_values$ID == "white", "temp"] <- white_temp

	return(pixel_values)
}

#' Create a linear model of temperature as a function of thermal digital numbes
#'
#' @param pixel_values A data.frame linking thermal pixel values to temperature,
#'                     as returned by \code{\link{join_thermal}}
#' @param summary_functions A list with functions used to summarize the values
#'              for a given panel color. Defaults to max for black and gray
#'              panels and min for white panels.
#'
#' @return A list of two elements, the first one containing a data.frame of
#'         DN and temperature values for each panel, and the second element
#'         containing a linear model linking the two.
#'
#' @examples
#' NULL
#'
#' @export
thermal_model <- function(pixel_values,
			  summary_functions = list(black = max, gray = max, white = min)) {

	output <- data.frame(ID = unique(pixel_values$ID),
			     pixel = NA,
			     temp = NA)

	for(i in 1:nrow(output)) {
		i_value <- output[i, "ID"]
		output[i, "pixel"] <- summary_functions[[i_value]](pixel_values[pixel_values$ID == i_value, "thermal"])
		stopifnot(length(i_temp <- unique(pixel_values[pixel_values$ID == i_value, "temp"])) == 1)
		output[i, "temp"] <- i_temp
	}

	# Now computing the linear model
	output_lm <- lm(temp ~ pixel, data = output)

	# Returning a list with the data and the linear model
	list(data = output, model = output_lm)
}


#' Plot panel temperature as a function of thermal digital numbers
#'
#' To complete
#'
#' @param pixel_values A data.frame linking thermal pixel values to temperature,
#'                     as returned by \code{\link{join_thermal}}
#' @param summary_functions A list with functions used to summarize the values
#'              for a given panel color. Defaults to max for black and gray
#'              panels and min for white panels.
#' @param plot_lm A logical value indicating whether a linear model should be computed
#'                using the thermal_model function and plotted using abline.
#' @param col.model.points The color(s) to use for the points used to compute the model.
#' @param cex.model.points The cex parameter to use for the points used to compute the model.
#'
#' @return NULL, invisibly. This function is invoked for its plotting side-effect.
#'
#' @examples
#' NULL
#'
#' @export
plot_pixtemp <- function(pixel_values,
			 summary_functions = list(black = max, gray = max, white = min),
			 plot_lm = TRUE,
			 col.model.points = "red",
			 cex.model.points = 2) {


	plot(pixel_values$thermal, pixel_values$temp, col = ifelse(pixel_values$ID == "white", "blue", pixel_values$ID))

	# Plotting a linear model of temperature as a function of thermal values if requested
	if(plot_lm) {
		lmod <- thermal_model(pixel_values, summary_functions = summary_functions)
		abline(reg = lmod$model, lty = 3)
		# Also plotting the points used for the model
		points(x = lmod$data$pixel, y = lmod$data$temp,
		       col = col.model.points,
		       cex = cex.model.points)
	}

	return(invisible(NULL))
}

