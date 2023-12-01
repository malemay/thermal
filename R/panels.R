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

