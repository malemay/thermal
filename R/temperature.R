#' Load panel temperature data
#'
#' To complete
#'
#' @param filename A character, the name of the file
#'
#' @return A data.frame of temperature values with the time colum set
#'         to the proper data type.
#'
#' @examples
#' NULL
#'
#' @export
read_temp <- function(filename) {
	tempdata <- read.table(filename, sep = ",", skip = 4, na.strings = c("NA", "NAN"))
	colnames(tempdata) <- c("time", "record", "battery", "devtemp", "temp")
	tempdata$time <- strptime(tempdata$time, "%Y-%m-%d %H:%M:%S")
	tempdata
}

#' Plot panel temperature data over a given time range
#'
#' To complete
#'
#' @param tempdata A list with three elements named "black", "gray" and "white" for
#'                 each of the three panel colors. Each element is a data.frame of
#'                 temperature data as read with \code{\link{read_temp}}
#' @param xmin A date object identifying the minimum value used for the x-axis.
#'             If NULL, the whole data range is plotted.
#' @param xmax A date object identifying the maximum value used for the x-axis.
#'             If NULL, the whole data range is plotted.
#'
#' @return NULL, invisibly. This function is invoked for its plotting side-effect.
#'
#' @examples
#' NULL
#'
#' @export
plot_temp <- function(tempdata, xmin = NULL, xmax = NULL) {
	
	# Subsetting the input data to the range of interest
	# BUG: xmin and/or xmax should be set even if only one of them is NULL
	if(!is.null(xmin) && !is.null(xmax)) {
		tempdata <- lapply(tempdata, function(x, xmin, xmax) {
					   x <- x[x$time >= xmin & x$time <= xmax, ]
					   x
				   }, xmin = xmin, xmax = xmax)
	}

	# Getting the limits of the x- and y-axes and setting the plotting region accordingly
	xrange <- range(do.call("rbind", tempdata)$time)
	yrange <- range(do.call("rbind", tempdata)$temp, na.rm = TRUE)

	plot(1, type = "n", xlim = as.numeric(xrange), ylim = yrange,
	     xaxt = "n")
	
	# Adding lines for each of the panels
	# By default the white panel is plotted in skyblue (this should be allowed to vary as a parameter)
	for(i in c("black", "gray", "white")) {
		lines(tempdata[[i]]$time, tempdata[[i]]$temp, col = if(i == "white") "skyblue" else i)
	}

	# Adding an axis for the time
	axis.POSIXct(1, tempdata$black$time)

	invisible(NULL)
}

#' Extracts the temperature at the nearest time point to a thermal image
#'
#' To complete
#'
#' @param metadata A data.frame containing metadata on a set of thermal images,
#'                 as returned by \code{\link{read_metadata}}. Must minimally
#'                 contain a column called "CreateDate" which indicates when
#'                 the picture was taken.
#' @param temperature A data.frame of panel temperature data, as returned by the
#'                 function \code{\link{read_temp}}. Should contain a column called
#'                 "time" to allow matching the time stamps of both datasets.
#'
#' @return A numeric vector containing the temperature at the nearest time point
#'         for every row in the metadata input.
#'
#' @examples
#' NULL
#'
#' @export
extract_temp <- function(metadata, temperature) {
	# Creating the output vector
	output_temp <- numeric(nrow(metadata))

	# Looping over all the rows
	for(i in 1:nrow(metadata)) {
		i_time <- metadata[i, "CreateDate"]
		temperature$timediff <- abs(i_time - temperature$time)
		# Adding a sanity check to ensure that no time difference is
		# Greater than 10 seconds
		if(min(temperature$timediff) > 10) stop("No temperature taken less than 10 seconds from row ", i)
		output_temp[i] <- temperature[which.min(temperature$timediff), "temp"]
	}

	output_temp
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

#' Plot panel temperature as a function of thermal digital numbers
#'
#' To complete
#'
#' @param pixel_values A data.frame linking thermal pixel values to temperature,
#'                     as returned by \code{\link{join_thermal}}
#'
#' @return NULL, invisibly. This function is invoked for its plotting side-effect.
#'
#' @examples
#' NULL
#'
#' @export
plot_pixtemp <- function(pixel_values) {

	# Creating a linear model of temperature as a function of thermal values
	lmod <- lm(temp ~ thermal, data = pixel_values)

	plot(pixel_values$thermal, pixel_values$temp, col = ifelse(pixel_values$ID == "white", "blue", pixel_values$ID))
	abline(reg = lmod)

	return(invisible(NULL))
}

