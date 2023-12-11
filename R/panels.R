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

#' Parse panel coordinates returned by the panel_coordinates.py script
#'
#' This function parses the coordinates of a given panel (black, gray or white)
#' on different pictures as written out by the Python script that translates
#' panel coordinates from the oprthorectified images to the undistorted images.
#' The coordinates on the undistorted images will then be transferred to
#' thermal coordinates through the visible/thermal image registration
#' functions.
#'
#' @param filename A character. The name of the file containing the panel coordinates.
#' @param color A character indicating the color of the panel being read.  will
#' be stored as metadata in the output data.frames.
#' @param image_height The height of the image. Used for the purposes of
#' flipping y coordinates from image coordinates (origin is top-right corner)
#' to raster coordinates (origin is bottom-left corner) for proper use in
#' downstream applications.
#' @param flip_y A logical. Whether the y coordinates should be flipped to
#' convert from image to raster coordinates. Defaults to TRUE.
#'
#' @return A named list of data.frames, each with the x- and y-coordinates of
#' the panels and their color. The list is named according to the index of the
#' image that the coordinates refer to.
#'
#' @examples
#' NULL
parse_coords <- function(filename, color, image_height, flip_y = TRUE) {

	# Sanity checks
	stopifnot(is.character(color) && color %in% c("black", "gray", "white"))

	# Reading the file and properly splitting along images and newlines
	coords <- paste(readLines(filename), collapse = "\n")
	coords <- strsplit(coords, "ID:")[[1]]
	coords <- lapply(coords, function(x) strsplit(x, "\n")[[1]])
	coords <- Filter(length, coords)

	# Naming the list elements after image ID and formatting as matrices
	names(coords) <- sapply(coords, function(x) x[1])
	coords <- lapply(coords, function(x) x[-1])
	coords <- lapply(coords, function(x) do.call("rbind", strsplit(x, ",")))

	# Properly formatting each coordinates matrix as a data.frame with color as a column
	coords <- lapply(coords, function(x) {
				 x <- as.data.frame(x)
				 names(x) <- c("x", "y")
				 x$x <- as.numeric(x$x)
				 x$y <- as.numeric(x$y)
				 if(flip_y) x$y <- image_height - x$y
				 x$color <- color
				 x
	     })

	# Returning the list of data.frames
	coords
}

#' Parse all panel coordinates in a given directory
#'
#' This function is a wrapper for \code{\link{parse_coords}} which properly
#' reads all panel coordinates files (black, gray and white panels) in a
#' directory and formats them as sf polygon objects for use in downstream
#' processing. This function makes a lot of assumptions about the format of the
#' coordinates files (they must have been output by panel_coordinates.py) and
#' their names (black_image_coords.txt, gray_image_coords.txt,
#' white_image_coords.txt).
#'
#' @param coord_dir A character. The path to the directory where the
#' coordinates files are.
#' @inheritParams parse_coords
#'
#' @return A names list with as many elements as there are pictures in the
#' dataset, with the sf-formatted polygons describing the coordinates of the
#' black, gray and white panels.
#'
#' @export
#' @examples
#' NULL
parse_panels <- function(coord_dir, image_height, flip_y = TRUE) {

	# Read the coordinates of all three panels
	black_coords <- parse_coords(paste0(coord_dir, "/black_image_coords.txt"), "black", image_height, flip_y)
	gray_coords <-  parse_coords(paste0(coord_dir, "/gray_image_coords.txt"),  "gray", image_height, flip_y)
	white_coords <- parse_coords(paste0(coord_dir, "/white_image_coords.txt"), "white", image_height, flip_y)

	# Initializing the output list (will have one element per image)
	output <- list()

	# Sanity check
	stopifnot(identical(names(black_coords), names(gray_coords)) && identical(names(black_coords), names(white_coords)))

	# Looping over all the images
	for(i in names(black_coords)) {
		output[[i]] <- rbind(black_coords[[i]], gray_coords[[i]], white_coords[[i]])
		output[[i]] <- lapply(split(as.data.frame(output[[i]][, 1:2]), output[[i]][, 3]), function(x) sf::st_polygon(list(as.matrix(x))))
		output[[i]] <- do.call(sf::st_sfc, output[[i]])
		output[[i]] <- sf::st_as_sf(output[[i]])
		output[[i]]$color <- c("black", "gray", "white")
	}

	output
}

#' Translate polygons from visible to thermal coordinates
#'
#' This function computes the position of a given set of polygons on a thermal
#' image based on their positions on a visible image. This transformation
#' depends on a set of four transformation parameters that must be supplied
#' as input. More details on the transformation model used can be found in
#' the documentation for \code{\link{align_images}}.
#'
#' This function has been developed with the coordinates of reference thermal
#' panels in mind and may need to be improved for more general usage.
#'
#' @param x The sf polygons to which the translation should be applied.
#' @param params A numeric vector of length 4 that describes the translation
#' from visible to thermal coordinates, as computed by \code{\link{align_images}}.
#' @param distortion_center A numeric vector of length two. The position of the
#' center of the image to be used for the distortion model.
#'
#' @return A polygon object of the same length as the input, with coordinates
#' transformed to their equivalent on the thermal image.
#'
#' @export
#' @examples
#' NULL
#'
translate_polygons <- function(x, params, distortion_center = c(2000, 1500)) {
	# Converting the coordinates according to the optimized model
	xcoords <- sf::st_coordinates(x)
	new_coords <- convert_coordinates(as.matrix(xcoords[, 1:2]),
					  params = params,
					  distortion_center = distortion_center)
	xcoords[, 1] <- new_coords[, 1]
	xcoords[, 2] <- new_coords[, 2]

	# Creating new polygons with this geometry
	new_poly <- lapply(split(as.data.frame(xcoords[, 1:2]), xcoords[, 4]), function(x) sf::st_polygon(list(as.matrix(x))))
	new_poly <- do.call(sf::st_sfc, new_poly)

	# Giving that geometry to the input polygons
	sf::st_geometry(x) <- new_poly
	x
}

#' Remove polygons outside the bounds of an image
#'
#' @param polygons A set of polygons in image coordinates (pixel positions)
#' relative to the image that they were extracted from.
#' @param extent A SpatExtent object representing the extent of the raster on
#' which the polygons are located.
#' @param overlap A numeric of length 1 indicating the proportion (between 0 and 1)
#' of the input polygons that must overlap the input extent to be kept in the
#' output.
#'
#' @return A list of polygons similar to that input, but with polygons located
#' partially or fully outside the bounds of the image removed.
#'
#' @export
#' @examples
#' NULL
#'
filter_polygons <- function(polygons, extent, overlap = 1) {
	# Convert the extent region to an sf polygon object such that the intersection can be computed
	ext_poly <- sf::st_as_sfc(sf::st_bbox(extent))

	# Compute the proportion of overlap for each set of polygons
	overlaps <- sapply(polygons, polygon_overlap, ext_poly = ext_poly)

	# Return the polygons that meet the output criterion
	polygons[overlaps >= overlap]
}

#' Determine the proportion of a set of polygons that overlaps an image
#'
#' @param polygon An sf object representing polygons.
#' @param ext_poly An sf object representing an area whose intersection
#' with the input polygons will be computed.
#'
#' @return A single numeric value between 0 and 1 describing the proportion
#' of the input polygon object that overlaps the extent.
#'
#' @examples
#' NULL
polygon_overlap <- function(polygon, ext_poly) {
	# Find the total area of the input polygons
	total_area <- sum(sf::st_area(polygon))

	# Find the area that overlaps the image
	overlap_area <- sum(sf::st_area(suppressWarnings(sf::st_intersection(polygon, ext_poly))))

	# Return the ratio of the overlapping to the total area
	overlap_area / total_area
}

#' Transform polygon positions based on a translation-rotation model
#'
#' This function can be used to transform polygon coordinates given rotation
#' and translation parameters, as obtained from the function
#' \code{\link{optimize_transform}}. Not to be confused with the function
#' \code{\link{translate_polygons}}, which modifies polygon coordinates
#' according to a visible/thermal image registraition model.
#'
#' This function can be used, for example, to adjust panel coordinates
#' such that they align properly with panels on thermal images if
#' the coordinates following visible to thermal translation has not worked
#' perfectly.
#'
#' @param x An sf object representing polygons.
#' @param theta A numeric value representing the angle of rotation (in degrees).
#' @param center A numeric vector of length two. The center around which the
#' rotation should be performed.
#' @param htrans A numeric of length one. The horizontal translation parameter.
#' @param vtrans A numeric of length one. The vertical translation parameter.
#'
#' @return An sf polygon object similar to the input one, but with updated
#' coordinates.
#'
#' @export
#' @examples
#' NULL
transform_polygons <- function(x, theta, center, htrans, vtrans) {
	# Converting the coordinates according to the optimized model
	xcoords <- sf::st_coordinates(x)

	new_coords <- transform_coords(xcoords[, 1:2],
				       center = center,
				       theta = theta,
				       htrans = htrans,
				       vtrans = vtrans,
				       reverse = FALSE)

	# Putting the coordinates back in the original matrix
	xcoords[, 1] <- new_coords[, 1]
	xcoords[, 2] <- new_coords[, 2]

	# Creating new polygons with this geometry
	new_poly <- lapply(split(as.data.frame(xcoords[, 1:2]), xcoords[, 4]), function(x) sf::st_polygon(list(as.matrix(x))))
	new_poly <- do.call(sf::st_sfc, new_poly)

	# Giving that geometry to the input polygons
	sf::st_geometry(x) <- new_poly
	x
}

#' Create a raster representation from reference panel coordinates
#'
#' This function creates a raster from reference panel coordinates such
#' that this raster can be used to find the coordnates of the reference panels
#' that best align to the thermal image. This function is not exported and not
#' meant to be called onits own. See \code{\link{adjust_panels}} for more details.
#'
#' @param panels An sf object of polygons representing reference panel coordinates.
#' @param template A SpatRaster object to be used as a template (for dimensions and
#' resolution) for creating the new SpatRaster from the panel positions.
#' @param black_value A numeric of length one. The value that thermal pixels corresponding
#' to the black panel should be given.
#' @param gray_value Similar to black_value, but for the gray panel.
#' @param white_value Similar to black_value, but for the white panel.
#' @param avg_value A numeric of length one. The value given to background pixels that
#' do not correspond to any panel.
#'
#' @return A SpatRaster object that contains values consistent with the positions of
#' reference panels. This raster object can be used to find the best alignment between
#' panel coordinates and the target thermal image.
#'
#' @examples
#' NULL
rasterize_panels <- function(panels, template, black_value, gray_value, white_value, avg_value) {

	# Creating a raster from the template using the average background value
	panel_raster <- template
	panel_raster[] <- avg_value

	# Extracting the cell indices of each panel
	vect_panels <- terra::vect(panels)
	black_vect <- vect_panels[terra::values(vect_panels)$color == "black"]
	gray_vect <- vect_panels[terra::values(vect_panels)$color == "gray"]
	white_vect <- vect_panels[terra::values(vect_panels)$color == "white"]

	black_cells <- terra::cells(panel_raster, black_vect)[, 2]
	gray_cells <- terra::cells(panel_raster, gray_vect)[, 2]
	white_cells <- terra::cells(panel_raster, white_vect)[, 2]

	# Using these indices to
	panel_raster[black_cells] <- black_value
	panel_raster[gray_cells] <- gray_value
	panel_raster[white_cells] <- white_value

	panel_raster
}

#' Adjust the positions of reference panels on a thermal image
#'
#' This function aims to enable the use of \code{\link{optimize_transform}} for
#' finding the transformation of reference panel coordinates that best
#' matches the thermal image that they are associated with. The idea is to
#' create a raster that will best correlate to the thermal image when the panel
#' positions are aligned with their actual position with the image. To achieve this,
#' higher values should be attributed to the black and gray panels, and lower values
#' to the white panel. This rasterized form of the panel coordinates is then passed
#' to \code{\link{optimize_transform}} to find the rotation/translation parameters
#' that best align the panel coordinates to the target thermal image.
#'
#' @param target_raster A SpatRaster object representing a thermal image to which the
#' coordinates of the reference panels should be aligned.
#' @inheritParams rasterize_panels
#' @inheritParams optimize_transform
#'
#' @return An sf polygon object similar to the input one, but with updated
#' coordinates corresponding to the optimized panel positions.
#'
#' @export
#' @examples
#' NULL
adjust_panels <- function(panels, target_raster, black_value, gray_value, white_value, avg_value,
			  theta1 = 0, htrans1 = 0, vtrans1 = 0,
			  theta_range = 3, theta_length = 3,
			  htrans_range = 5, htrans_length = 5,
			  vtrans_range = 5, vtrans_length = 5,
			  min_overlap = 100,
			  fact = 1, cores = 1, min_cor = 0.5, maxiter = 3,
			  method = "Nelder-Mead", reltol = 10^-5, trace = 1) {

	# Creating a raster from the panel coordiantes, to used for optimization
	panel_raster <- rasterize_panels(panels = panels, template = target_raster,
					 black_value = black_value, gray_value = gray_value,
					 white_value = white_value, avg_value = avg_value)

	# Finding the optimized transform from the panel raster to the target raster
	panel_transform <- optimize_transform(x = target_raster, y = panel_raster,
					      theta1 = theta1, htrans1 = htrans1, vtrans1 = vtrans1,
					      theta_range = theta_range, theta_length = theta_length,
					      htrans_range = htrans_range, htrans_length = htrans_length,
					      vtrans_range = vtrans_range, vtrans_length = vtrans_length,
					      min_overlap = min_overlap,
					      fact = fact, cores = cores, min_cor = min_cor, maxiter = maxiter,
					      method = method, reltol = reltol, trace = trace)

	params <- panel_transform$optim$par

	# Transforming the coordinates of the input panels using the optimized parameters
	transform_polygons(panels,
			   center = c(terra::xmax(panel_raster) - terra::xmin(panel_raster),
				      terra::ymax(panel_raster) - terra::ymin(panel_raster)) / 2,
			   theta = params[1], htrans = params[2], vtrans = params[3])
}

