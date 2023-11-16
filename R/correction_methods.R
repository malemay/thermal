#' Compute corner positions of a raster on another given transformation parameters
#'
#' @param x A raster whose coordinates are to be transformed.
#' @param theta A numeric. A rotation parameter (in degrees) for the transform.
#' @param htrans A numeric. A horizontal translation parameter for the transform.
#' @param vtrans A numeric. A vertical translation parameter for the transform.
#' @param reverse A logical. Whether the transformation should be reversed relative
#' to the input parameters.
#'
#' @return An sf polygon that gives the coordinates of the raster in the transformed
#' coordinate system.
#'
#' @examples
#'NULL
get_corners <- function(x, theta, htrans, vtrans, reverse = FALSE) {
	corner_coords <- matrix(c(terra::xmin(x), terra::ymax(x),
				  terra::xmax(x), terra::ymax(x),
				  terra::xmax(x), terra::ymin(x),
				  terra::xmin(x), terra::ymin(x)),
				ncol = 2, byrow = TRUE)

	output <- transform_coords(coords = corner_coords,
				   theta = theta,
				   center = c(terra::xmax(x) - terra::xmin(x), terra::ymax(x) - terra::ymin(x)) / 2,
				   htrans = htrans,
				   vtrans = vtrans,
				   reverse = reverse)

	# Transforming the coordinates into an sf polygon
	sf::st_as_sf(sf::st_sfc(sf::st_polygon(list(rbind(output, output[1, ])))))
}

#' Compute the difference between the mean of the thermal regions in overlapping images
#'
#' The transformation parameters (theta, htrans, and vtrans) apply to the transformation
#' from y to x.
#'
#' @param x,y Two rasters that overlap each other with known transformation parameters
#' from raster y to raster x.
#' @inherit get_corners
#'
#' @return The difference of the mean of raster y to that of raster x in their shared extent.
#'
#' @examples
#' NULL
get_diff <- function(x, y, theta, htrans, vtrans) {

	x_corners <- get_corners(x = suppressWarnings(terra::rast(y)),
				 theta = theta,
				 htrans = htrans,
				 vtrans = vtrans)

	y_corners <- get_corners(x = suppressWarnings(terra::rast(x)),
				 theta = theta,
				 htrans = htrans,
				 vtrans = vtrans, 
				 reverse = TRUE)

	x_mean <- mean(terra::extract(x, x_corners)$FLIR_RAW_THERMAL_IMAGE)
	y_mean <- mean(terra::extract(y, y_corners)$FLIR_RAW_THERMAL_IMAGE)

	# Getting the difference of y relative to x
	y_mean - x_mean
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

#' Compute difference in thermal values for all thermal pictures in a dataset
#'
#' @param metadata A data.frame of metadata on a set of thermal pictures, including
#' transformation parameters (theta, htrans, vtrans).
#' @param ncores An integer. The number of processes to launch (passed to \code{\link[parallel]{mclapply}})
#' @param verbose A logical. Whether the function should output information on its progress.
#'
#' @return A numeric vector of the same length as the number of rows in the input
#' data.frame. The first value will be NA because the first picture has no previous
#' picture.
#'
#' @export
#' @examples
#' NULL
compute_diffs <- function(metadata, ncores = 1, verbose = TRUE) {
	# Sanity checks
	if(! all(c("theta", "htrans", "vtrans") %in% names(metadata))) {
		stop("Metadata columns must include transformation parameters theta, htrans and vtrans.")
	}

	# Using mclapply to compute the differences in parallel using multiple cores
	diffs <- parallel::mclapply(2:nrow(metadata), FUN = function(i, params, verbose) {
					     if(verbose) message("Processing row ", i, " out of ", nrow(params))

					     get_diff(x = terra::rast(params[i - 1, "SourceFile"]),
						      y = terra::rast(params[i, "SourceFile"]),
						      theta = params[i, "theta"],
						      htrans = params[i, "htrans"],
						      vtrans = params[i, "vtrans"])

				 }, params = metadata, verbose = verbose, mc.cores = ncores)

	# No diff is available for the first image so we use NA
	c(NA, unlist(diffs))
}

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

#' Correct thermal drift in a set of pictures
#'
#' By default this function will output files to the specified directory
#' and replace the file extension by ".tiff".
#'
#' @param metadata A data.frame of metadata on a set of thermal pictures, sorted
#' according to the time when the picture was taken ("CreateDate" column).
#' @param output_dir The directory to which the output images should be written.
#' The function will not allow source files to be overwritten, thus the output
#' directory should be different from the one containing the source data.
#' @param method A character. The method to use for drift correction. At the moment
#' only "overlap" is supported.
#' @param midpoint An integer. The index to start the correction from. Image values
#' will be adjusted to match the image at that index. If NULL (the default), then
#' the middle of the range is used.
#' @param overwrite_dst A logical. Whether the destination files should be overwritten if they already exist.
#' @param verbose A logical. Whether the function should output information on its progress (default = TRUE).
#'
#' @return A character vector of the files that were written to disk, invisibly.
#'
#' @export
#' @examples
#' NULL
correct_drift <- function(metadata, output_dir, method = "overlap", midpoint = NULL, overwrite_dst = FALSE, verbose = TRUE) {
	# Some sanity checks
	stopifnot("SourceFile" %in% names(metadata))

	# Should also check that we are not going to overwrite the original files
	src_files <- metadata$SourceFile
	dst_files <- paste0(output_dir, "/", sub("\\..*$", ".tiff", basename(src_files)))
	if(any(dst_files %in% src_files)) stop("correct_drift does not allow overwriting source files")

	# Also checking whether we are going to overwrite destination files
	if(!overwrite_dst && any(file.exists(dst_files))) stop("Overwriting destination files not allowed with overwrite_dst = FALSE")

	# Checking that all files are in increasing time order
	stopifnot(!is.unsorted(metadata$CreateDate))

	# Overlap method adjusts the mean value of images based on the shared overlap between successive images
	if(method == "overlap") {
		if(! "diff" %in% names(metadata)) stop("The difference between successive images must be precomputed to use method = 'overlap'")

		# Setting the midpoint if it is not already set
		if(is.null(midpoint)) midpoint <- floor(length(src_files) / 2)

		# Precomputing the adjustment factors based on the differences and midpoint
		adj_factors <- numeric(length(src_files))

		for(i in midpoint:1) {
			adj_factors[i] <- ifelse(i == midpoint, 0, adj_factors[i + 1] + metadata[i + 1, "diff"])
		}

		for(i in midpoint:length(src_files)) {
			adj_factors[i] <- ifelse(i == midpoint, 0, adj_factors[i - 1] - metadata[i, "diff"])
		}

		# Looping over all pictures in the dataset
		for(i in 1:length(src_files)) {
			src <- src_files[i]
			dst <- dst_files[i]

			raw_image <- suppressWarnings(terra::rast(src))

			if(verbose) message("Processing index ", i, " with adjustment = ", adj_factors[i])

			processed_image <- raw_image + adj_factors[i]

			writeRaster(processed_image,
				    filename = dst,
				    datatype = "INT2U",
				    overwrite = overwrite_dst)

			transfer_exif(src, dst)
		}

	} else {
		stop("Only method = 'overlap' is supported at the moment.")
	}

	invisible(dst_files)
}

#' Transfer exif metadata from the original to corrected files
#'
#' @param src_file The file from which the exif data should be taken.
#' @param dst_file The file to which the exif data will be transferred.
#'
#' @return NULL, invisibly
#'
#' @examples
#' NULL
transfer_exif <- function(src_file, dst_file) {
	stopifnot(all(file.exists(c(src_file, dst_file))))
	exifr::exiftool_call('-overwrite_original -tagsFromFile', c(src_file, dst_file))
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

#' Correct the effect of vignetting in a thermal picture dataset
#'
#' @param metadata A data.frame of metadata on a set of thermal pictures.
#' @param output_dir The directory to which the output images should be written.
#' The function will not allow source files to be overwritten, thus the output
#' directory should be different from the one containing the source data.
#' @param method A character. The method to use for drift correction. At the moment
#' only "overall" is supported.
#' @param overwrite_dst A logical. Whether the destination files should be overwritten if they already exist.
#' @param verbose A logical. Whether the function should output information on its progress (default = TRUE).
#'
#' @return A character vector of the files that were written to disk, invisibly.
#'
#' @export
#' @examples
#' NULL
correct_vignetting <- function(metadata, output_dir, method = "overall", overwrite_dst = FALSE, verbose = TRUE) {

	# Some sanity checks
	stopifnot("SourceFile" %in% names(metadata))

	# Should also check that we are not going to overwrite the original files
	src_files <- metadata$SourceFile
	dst_files <- paste0(output_dir, "/", sub("\\..*$", ".tiff", basename(src_files)))
	if(any(dst_files %in% src_files)) stop("correct_drift does not allow overwriting source files")

	# Also checking whether we are going to overwrite destination files
	if(!overwrite_dst && any(file.exists(dst_files))) stop("Overwriting destination files not allowed with overwrite_dst = FALSE")

	# Overall method adjusts the mean value of images based on the global vignetting pattern
	if(method == "overall") {

		# We compute the vignetting adjustment as the difference between the pattern of each pixel and the overall mean
		vignetting_pattern <- compute_vignetting(tmeta_params)
		vignetting_adjustment <- terra::global(vignetting_pattern, mean)$mean -  vignetting_pattern

		# Looping from the first to the last picture in the dataset
		for(i in 1:length(src_files)) {
			src <- src_files[i]
			dst <- dst_files[i]

			if(verbose) message("Processing index ", i)

			raw_image <- suppressWarnings(terra::rast(src))
			processed_image <- raw_image + vignetting_adjustment

			writeRaster(processed_image,
				    filename = dst,
				    datatype = "INT2U",
				    overwrite = overwrite_dst)

			transfer_exif(src, dst)
		}

	} else {
		stop("Only method = 'overall' is supported at the moment.")
	}

	invisible(dst_files)
}

