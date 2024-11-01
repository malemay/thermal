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
#' @param metadata A data.frame of metadata on a flight, such as returned by
#' \code{\link{read_metadata}}.  Must minimally contain columns GPSLongitude and
#' GPSLatitude for longitude and latitude coordinates in the WGS84 datum.
#' @param coords A character vector of length two indicating which columns
#' in metadata are to be interpreted as X and Y coordinates (in that order).
#' @param base_crs An object that can be interpreted as a CRS description. Will
#' be used to set the coordinate system of the metadata. If NULL (the default),
#' then WGS84 (EPSG:4326) is assumed.
#' @param new_crs A coordinate reference system (such as returned by
#' \code{\link[sf]{st_crs}} to project the points to before display. If NULL
#' (the default), then data is displayed in the original coordinate system.
#' @param color A vector of character, numeric or factor data to specify how
#' to color the points. If a character value, this is interpreted as a color
#' value. If a factor, when a distinct color is used (from a predetermined set)
#' to represent each value of the factor. If a numeric, then a color scale
#' (see the color_scale argument) is used to translate numeric values to colors.
#' @param color_palette A character specifying the color palette to use when
#' mapping numeric values onto a color scale. Must correspond to a color
#' palette in the RColorBrewer package. See
#' \code{\link[RColorBrewer]{brewer.pal.info}} for more details.
#' @param n_colors The number of colors to use in the color scale when mapping
#' numeric values to colors. See \code{\link[RColorBrewer]{brewer.pal.info}}
#' for more details.
#' @param round_digits The decimal position to round the numeric values to
#' when mapping numeric values to a color scale.
#' @param features An optional sf object representing the position of a feature
#' of interest, whose coordinates will be transformed if new_crs is not NULL. Will
#' be coerced to an object that can be plotted using ­\code{\link[sf]{st_geometry}}
#' and \code{\link[sf]{st_as_grob}}.
#' @param feature_color A single character value: the color to use for the
#' features.
#' @param arrows A logical value specifying whether arrows should be plotted
#' to show the trajectory of the drone.
#' @param arrow_angle A numeric value controlling the angle of the arrows.
#' See \code{\link[grid]{arrow}}.
#' @param arrow_length A unit object controlling the length of the arrows.
#' See \code{\link[grid]{arrow}}.
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
plot_metadata <- function(metadata,
			  coords = c("GPSLongitude", "GPSLatitude"), base_crs = NULL, new_crs = NULL,
			  color = "black", color_palette = "RdBu", n_colors = 11, round_digits = 0,
			  features = NULL, feature_color = "black",
			  arrows = FALSE, arrow_angle = 30, arrow_length = grid::unit(0.25, "inches"),
			  title = NULL, cex.axis = 1, cex.points = 1) {
	# Checking the inputs
	stopifnot(all(coords %in% colnames(metadata)))

	# Checking if the base CRS was specified
	if(is.null(base_crs)) {
		warning("base_crs not set, assuming WGS84 (EPSG:4326)")
		base_crs <- sf::st_crs("EPSG:4326")
	}

	# Coercing to a spatial object
	metadata <- sf::st_as_sf(metadata, coords = coords, crs = base_crs, dim = "XY")

	# Transforming the CRS if needed
	if(!is.null(new_crs)) {
		metadata <- sf::st_transform(metadata, crs = new_crs)

		if(!is.null(features)) features <- sf::st_transform(features, crs = new_crs)
	}

	# Setting the colors to use
	factcol <- is.factor(color)
	numcol <- is.numeric(color)

	if(factcol) {
		factor_values <- color
		color <- gg_hue(length(unique(color)))[as.integer(color)]
	} else if(numcol) {
		color_map <- map_color(values = color, pal = color_palette, n_colors = n_colors)
		color <- color_map$mapped_colors
	}

	# Setting up the viewports
	# We reserve one-tenth of the viewport width for plotting a color scale
	grid::grid.newpage()
	
	# Setting the viewport for the plot
	grid::pushViewport(grid::viewport(layout = grid::grid.layout(ncol = 3, widths = grid::unit(c(9, 1, 1), "null"))))
	grid::pushViewport(grid::viewport(layout.pos.col = 1))

	grid::pushViewport(grid::plotViewport(margins = c(5.1, 4.1, 4.1, 0.6)))
	grid::pushViewport(grid::dataViewport(xData = sf::st_coordinates(metadata)[, "X"], yData = sf::st_coordinates(metadata)[, "Y"]))

	# Plotting the data and axes
	grid::grid.rect()
	grid::grid.points(sf::st_coordinates(metadata)[, "X"],
			  sf::st_coordinates(metadata)[, "Y"],
			  gp = grid::gpar(col = "black", fill = color, cex = cex.points),
			  pch = 21,
			  default.units = "native")

	# Plotting the arrows if required
	if(arrows) {
		for(i in 1:(nrow(metadata ) - 1)) {
			grid::grid.lines(x = sf::st_coordinates(metadata)[c(i, i + 1), "X"],
					 y = sf::st_coordinates(metadata)[c(i, i + 1), "Y"],
					 arrow = grid::arrow(angle = arrow_angle, length = arrow_length),
					 default.units  = "native")
		}
	}

	# Adding the position of the panels (if provided)
	if(!is.null(features)) {
		grid::grid.draw(sf::st_as_grob(sf::st_geometry(features),
					       gp = grid::gpar(col = feature_color),
					       pch = 16))
	}

	grid::grid.xaxis(gp = grid::gpar(cex = cex.axis))
	grid::grid.yaxis(gp = grid::gpar(cex = cex.axis))

	# Adding a title (if requested)
	if(!is.null(title)) {
		grid::grid.text(label = title, y = grid::unit(1, "npc") + grid::unit(1, "lines"))
	}

	# Going to the viewport in which the color scale will be plotted
	if(factcol || numcol) {
		grid::upViewport(3)
		grid::pushViewport(grid::viewport(layout.pos.col = 2))

		if(factcol) {
			grid.factorscale(factor_values)
		}

		if(numcol) {
			grid.colorscale(breaks = color_map$breaks,
					base_palette = color_map$base_palette,
					label_text = "",
					direction = "vertical",
					round_digits = round_digits)
		}
	}

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

#' Map a set of numeric values onto a color palette
#'
#' This function takes a set of numeric values and maps them on a color palette
#' provided by the RColorBrewer package. \code{NA} values are mapped to black
#' as a default color for missing values.
#'
#' @param values A numeric vector of values to be mapped onto the color scale.
#' \code{NA} values are accepted and mapped to the "black" color.
#' @param pal A character. The palette to use; must correspond to a palette in
#' the RColorBrewer package.
#' @param n_colors A numeric. The number of colors to use for mapping onto a
#' color scale. The maximum value allowed may differ depending on the
#' particular palette used.
#'
#' @return A list of length three comprising the following elements:
#'   \itemize{
#'     \item mapped_colors: a character vector of colors mapped from the values input vector
#'     \item breaks: a numeric vector of break positions used for mapping
#'     \item base_palette: the colors that were used for mapping the numeric values
#'   }
#'
#' @examples
#' NULL
map_color <- function(values, pal, n_colors) {
	base_palette <- RColorBrewer::brewer.pal(n = n_colors, name = pal)

	# Getting the breaks based on the minimum/maximum values and the scale_extend parameter
	breaks <- seq(min(values, na.rm = TRUE),
		      max(values, na.rm = TRUE),
		      length.out = n_colors + 1)

	# A vector of colors to use for plotting; NA values are changed to black
	mapped_colors <- base_palette[as.integer(cut(values, breaks = breaks, include.lowest = TRUE))]
	mapped_colors[is.na(mapped_colors)] <- "black"

	# We return a list with the mapped values, color levels and base palette colors
	list(mapped_colors = mapped_colors,
	     breaks = breaks,
	     base_palette = base_palette)
}

#' Plot the color scale used in a plot
#'
#' This function takes a set of numeric breaks and a set of colors and uses
#' this information to produce a graphical scale linking those breaks to the
#' colors used in a plot. This function is not exported from the package and is
#' therefore meant to be used only internally by the package.
#'
#' @param breaks A numeric vector of breaks used for the color scale.  There
#' must be one more break than the number of colors
#' @param base_palette A character vector of colors used in mapping the numeric
#' values onto the color scale.
#' @param label_text A character to use as a value for the legend label.
#' @param digits A numeric value indicating the number of digits that the
#' breaks in the scale should be rounded to.
#' @param direction A character indicating whether the direction in which the
#' color scale should be plotted. Supported values are "horizontal" and
#' "vertical".
#' @param fontsize A numeric. The font size to use for the axis label and tick
#' labels.
#'
#' @return \code{NULL}, invisibly. This function is invoked for its plotting
#' side-effects.
grid.colorscale <- function(breaks, base_palette, label_text, round_digits = 2,
			    direction = "horizontal", fontsize = 12) {

	# Checking that the number of breaks is one more than the number of colors
	if(length(breaks) != length(base_palette) + 1) {
		stop("Error in grid.colorscale(): length(breaks) should be equal to length(base_palette) + 1")
	}

	# Checking the value of direction
	if(!direction %in% c("horizontal", "vertical")) {
		stop("Error in grid.colorscale: direction must be one of 'horizontal' or 'vertical'")
	}

	horiz <- direction == "horizontal"

	# Getting the x-scale from the minimum and maximum values of the break labels
	xscale <- range(breaks)

	# Pushing a viewport with the appropriate scale
	grid::pushViewport(grid::viewport(width = 0.5,
					  height = 0.5,
					  xscale = xscale,
					  yscale = xscale))

	# Draw the scale colors as rectangles
	grid::grid.rect(x = if(horiz) breaks[-length(breaks)] else grid::unit(0, "npc"),
			width = if(horiz) diff(breaks) else grid::unit(1, "npc"),
			y = if(horiz) grid::unit(0, "npc") else breaks[-length(breaks)],
			height = if(horiz) grid::unit(1, "npc") else diff(breaks),
			default.units = "native",
			just = c(0, 0),
			gp = grid::gpar(fill = base_palette))

	# Adding the breaks for the scale
	if(horiz) {
		grid::grid.xaxis(at = breaks,
				 label = as.character(round(breaks, digits = round_digits)),
				 gp = grid::gpar(fontsize = fontsize))
	} else {
		grid::grid.yaxis(at = breaks,
				 label = as.character(round(breaks, digits = round_digits)),
				 gp = grid::gpar(fontsize = fontsize),
				 main = FALSE)
	}

	# Adding a name for the scale
	if(horiz) {
		grid::grid.text(label_text,
				y = grid::unit(-3, "lines"),
				gp = grid::gpar(fontsize = fontsize))
	} else {
		grid::grid.text(label_text,
				x = grid::unit(5, "lines"),
				gp = grid::gpar(fontsize = fontsize),
				rot = 90)
	}

	grid::upViewport()

	return(invisible(NULL))
}

#' Draw a color scale representing a factor
#'
#' Add a color scale to a plot where color was encoded as a factor
#'
#' @param values The factor vector that was used for plotting
#'
#' @return NULL, invisibly. This function is invoked for its plotting
#' side-effect.
grid.factorscale <- function(values) {
	n_colors <- length(unique(values))

	grid::grid.points(x = rep(0, n_colors),
			  y = seq(0.6, 0.4, length.out = n_colors),
			  gp = grid::gpar(col = "black", fill = gg_hue(n_colors)[1:n_colors]),
			  pch = 21)

	grid::grid.text(label = levels(values),
			x = rep(0.2, n_colors),
			y = seq(0.6, 0.4, length.out = n_colors),
			just = c(0, 0.5))

	invisible(NULL)
}

#' Generate a vector of colors comparable to those used by ggplot2
#'
#' See https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#'
#' @param n An integer of length one. The number of colors to generate
#'
#' @return A vector of colors to be used for plotting
gg_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}
