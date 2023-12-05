#' Get panel coordinates through interactive clicking
#'
#' Panel corners are to be clicked on clockwise.
#'
#' @param x A SpatRaster object to query for coordinates.
#'
#' @return An polygon object of class sfc that representes the coordinates of a
#' panel in a reference coordinate system.
#'
#' @examples
#' NULL
get_panel_coords <- function(x) {

	# Interactively getting the coordinates of the panel
	coords <- terra::click(x, n = 4, xy = TRUE)

	# Formatting the coordinates and closing the shape by repeating the first row
	coords <- rbind(coords[, c("x", "y")], coords[1, c("x", "y")])

	# Creating an object of the sfc class from the supplied coordinates
	output <- sf::st_as_sfc(sf::st_as_binary(sf::st_polygon(list(as.matrix(coords)))))

	# Setting the coordinate reference system (same as the raster)
	sf::st_crs(output) <- sf::st_crs(x)

	output
}

#' Interactively identify and write panel coordinates from a georeferenced raster
#'
#' This function will write the polygons representing the reference panels in
#' ESRI Shapefile format to the specified output directory. The names of the
#' files are hard-coded as "black.shp", "gray.shp" and "white.shp" and will not
#' be overwritten unless otherwise specified.
#'
#' @param image_file A character. The path to a georeferenced raster file.
#' @param output_dir A character representing a directory path to which the
#' files will be written. The directory will be created if it does not exist.
#' @param overwrite A logical value indicating whether existing files should
#' be overwritten by this function.
#'
#' @return A named list of polygon objects representing each of the three
#' panels.
#'
#' @export
#' @examples
#' NULL
write_panel_coords <- function(filepath, output_dir, overwrite = FALSE) {

	# Reading an georeferenced picture and interactively zooming into it
	ortho <- terra::rast(filepath)
	terra::RGB(ortho) <- 1:4
	terra::plot(ortho)

	message("Click on the image to zoom on an area of interest")
	panel_region <- terra::zoom(ortho)

	# Getting the coordinates for each of the panels
	message("Click clockwise on the four corners of the black panel")
	black <- get_panel_coords(terra::crop(ortho, panel_region))

	message("Click clockwise on the four corners of the gray panel")
	gray <- get_panel_coords(terra::crop(ortho, panel_region))

	message("Click clockwise on the four corners of the white panel")
	white <- get_panel_coords(terra::crop(ortho, panel_region))

	# Exporting the results as shapefiles
	dir.create(output_dir, recursive = TRUE)
	sf::st_write(black, paste0(output_dir, "/black.shp"), delete_layer = overwrite)
	sf::st_write(gray, paste0(output_dir, "/gray.shp"), delete_layer = overwrite)
	sf::st_write(white, paste0(output_dir, "/white.shp"), delete_layer = overwrite)

	return(list(black = black, gray = gray, white = white))
}

#' Compare panel area to expectations
#'
#' @param panels A list of polygons representing reference panel coordinates.
#' @param expected_area A numeric of lenght one indicating the expected area of
#' each panel. The function only supports cases where all panels are expected to
#' be of the same size.
#'
#' @return A numeric vector of the same length as the input list, showing the
#' difference between each panel's area and the expected area
#'
#' @export
#' @examples
#' NULL
check_panels <- function(panels, expected_area = 2.25) {
	sapply(panels, sf::st_area) - expected_area
}

#' Convert polygons from visible to thermal coordinates
#'
#' To complete
#'
#' @param x The polygons to which the conversion should be applied
#' @param optimout A model that describes the translation from visible to thermal
#'                 coordinates, as output by \code{\link{align_images}}
#' @param distortion_center A numeric vector of length two. The position of the center
#'                 of the image to be assumed for the distortion model.
#'
#' @return An polygon object of the same length as the input, with coordinates
#'         transformed to their equivalent on the thermal image.
#'
#' @examples
#' NULL
#'
#' @export
convert_polygons <- function(x, optimout, distortion_center = c(2000, 1500)) {
	# Converting the coordinates according to the optimized model
	xcoords <- sf::st_coordinates(x)
	new_coords <- convert_coordinates(x = xcoords[, 1], y = xcoords[, 2],
					  optimout = optimout,
					  distortion_center = distortion_center)
	xcoords[, 1] <- new_coords$x
	xcoords[, 2] <- new_coords$y

	# Creating new polygons with this geometry
	new_poly <- lapply(split(as.data.frame(xcoords[, 1:2]), xcoords[, 4]), function(x) sf::st_polygon(list(as.matrix(x))))
	new_poly <- do.call(sf::st_sfc, new_poly)

	# Giving that geometry to the input polygons
	sf::st_geometry(x) <- new_poly
	x
}

#' Remove polygons outside the bounds of an image
#'
#' To complete
#'
#' @param polygons A set of polygons with coordinates expressed as pixel positions in the image
#' @param xmin A numeric. The minimum allowed position on the x-axis.
#' @param xmax A numeric. The maximum allowed position on the x-axis.
#' @param ymin A numeric. The minimum allowed position on the y-axis.
#' @param ymax A numeric. The maximum allowed position on the y-axis.
#'
#' @return A list of polygons similar to that input, but with polygons located
#'         partially or fully outside the bounds of the image removed.
#'
#' @examples
#' NULL
#'
#' @export
filter_polygons <- function(polygons, xmin = 0, xmax = 640, ymin = 0, ymax = 512) {
	filtered <- sapply(polygons, function(x, xmin, xmax, ymin, ymax) {
				   xcoords <- sf::st_coordinates(x)
				   any(xcoords[, 1] < xmin | xcoords[, 1] > xmax | xcoords[, 2] < ymin | xcoords[, 2] > ymax) },
				   xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
	polygons[!filtered]
}

