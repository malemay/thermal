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

	x_corners <- get_corners(x = y,
				 theta = theta,
				 htrans = htrans,
				 vtrans = vtrans)

	y_corners <- get_corners(x = x,
				 theta = theta,
				 htrans = htrans,
				 vtrans = vtrans, 
				 reverse = TRUE)

	x_mean <- terra::extract(x, x_corners, mean, ID = FALSE)$FLIR_RAW_THERMAL_IMAGE
	y_mean <- terra::extract(y, y_corners, mean, ID = FALSE)$FLIR_RAW_THERMAL_IMAGE

	# Getting the difference of y relative to x
	y_mean - x_mean
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

#' Correct thermal drift in a set of pictures
#'
#' By default this function will output files to the specified directory
#' and replace the file extension by ".tiff".
#'
#' @param metadata A data.frame of metadata on a set of thermal pictures, sorted
#' according to the time when the picture was taken ("DateTimeOriginal" column).
#' @param output_dir The directory to which the output images should be written.
#' The function will not allow source files to be overwritten, thus the output
#' directory should be different from the one containing the source data.
#' @param method A character. The method to use for drift correction. At the moment
#' only "overlap" is supported.
#' @param midpoint An integer. The index to start the correction from. Image values
#' will be adjusted to match the image at that index. If NULL (the default), then
#' the middle of the range is used.
#' @param overwrite_dst A logical. Whether the destination files should be overwritten if they already exist.
#' @param ncores An integer. The number of cores to use for parallel processing.
#' @param verbose A logical. Whether the function should output information on its progress (default = TRUE).
#'
#' @return A character vector of the files that were written to disk, invisibly.
#'
#' @export
#' @examples
#' NULL
correct_drift <- function(metadata, output_dir, method = "overlap", midpoint = NULL, overwrite_dst = FALSE, ncores = 1, verbose = TRUE) {
	# Some sanity checks
	stopifnot("SourceFile" %in% names(metadata))

	# Should also check that we are not going to overwrite the original files
	src_files <- metadata$SourceFile
	dst_files <- paste0(output_dir, "/", sub("\\..*$", ".tiff", basename(src_files)))
	if(any(dst_files %in% src_files)) stop("correct_drift does not allow overwriting source files")

	# Also extracting the directory part and the extension for later use in transfer_exif
	src_dir <- unique(dirname(src_files))
	src_ext <- paste0(".", unique(tools::file_ext(src_files)))

	stopifnot(length(src_dir) == 1 && length(src_ext) == 1)

	# Also checking whether we are going to overwrite destination files
	if(!overwrite_dst && any(file.exists(dst_files))) stop("Overwriting destination files not allowed with overwrite_dst = FALSE")

	# Checking that all files are in increasing time order
	stopifnot(!is.unsorted(metadata$DateTimeOriginal))

	# Overlap method adjusts the mean value of images based on the shared overlap between successive images
	if(method == "overlap") {
		if(! "diff" %in% names(metadata)) stop("The difference between successive images must be precomputed to use method = 'overlap'")

		# Setting the midpoint if it is not already set
		if(is.null(midpoint)) midpoint <- floor(length(src_files) / 2)

		# Computing the pre- and post-midpoint adjustment factors
		pre_midpoint  <- if(midpoint == 1) c() else rev(cumsum(metadata[midpoint:2, "diff"]))
		post_midpoint <- if(midpoint == nrow(metadata)) c() else cumsum(-metadata[(midpoint + 1):nrow(metadata), "diff"])

		# Precomputing the adjustment factors based on the differences and midpoint
		adj_factors <- c(pre_midpoint, 0, post_midpoint)

	} else {
		stop("Only method = 'overlap' is supported at the moment.")
	}

	# Looping over all pictures in the dataset
	parallel::mclapply(1:length(src_files), FUN = function(i, src_files, dst_files, adj_factors, verbose, overwrite_dst) {
				   src <- src_files[i]
				   dst <- dst_files[i]

				   raw_image <- suppressWarnings(terra::rast(src))

				   if(verbose) message("Processing index ", i, " with adjustment = ", adj_factors[i])

				   processed_image <- raw_image + adj_factors[i]

				   terra::writeRaster(processed_image,
						      filename = dst,
						      datatype = "INT2U",
						      overwrite = overwrite_dst)
				     },
				   src_files = src_files,
				   dst_files = dst_files,
				   adj_factors = adj_factors,
				   verbose = verbose,
				   overwrite_dst = overwrite_dst,
				   mc.cores = ncores)

	# Transfering the EXIF metadata from the source files to the destination files
	transfer_exif(src_dir, src_ext, output_dir)

	invisible(dst_files)
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
#' @param ncores An integer. The number of cores to use for parallel processing.
#' @param verbose A logical. Whether the function should output information on its progress (default = TRUE).
#'
#' @return A character vector of the files that were written to disk, invisibly.
#'
#' @export
#' @examples
#' NULL
correct_vignetting <- function(metadata, output_dir, method = "overall", overwrite_dst = FALSE, ncores = 1, verbose = TRUE) {

	# Some sanity checks
	stopifnot("SourceFile" %in% names(metadata))

	# Should also check that we are not going to overwrite the original files
	src_files <- metadata$SourceFile
	dst_files <- paste0(output_dir, "/", sub("\\..*$", ".tiff", basename(src_files)))
	if(any(dst_files %in% src_files)) stop("correct_drift does not allow overwriting source files")

	# Also extracting the directory part and the extension for later use in transfer_exif
	src_dir <- unique(dirname(src_files))
	src_ext <- paste0(".", unique(tools::file_ext(src_files)))

	stopifnot(length(src_dir) == 1 && length(src_ext) == 1)

	# Also checking whether we are going to overwrite destination files
	if(!overwrite_dst && any(file.exists(dst_files))) stop("Overwriting destination files not allowed with overwrite_dst = FALSE")

	# Overall method adjusts the mean value of images based on the global vignetting pattern
	if(method == "overall") {

		# We compute the vignetting adjustment as the difference between the pattern of each pixel and the overall mean
		vignetting_pattern <- compute_vignetting(metadata)
		vignetting_adjustment <- terra::global(vignetting_pattern, mean)$mean -  vignetting_pattern

		# Looping from the first to the last picture in the dataset
		parallel::mclapply(1:length(src_files), FUN = function(i, src_files, dst_files, verbose, vignetting_adjustment, overwrite_dst) {
					   src <- src_files[i]
					   dst <- dst_files[i]

					   if(verbose) message("Processing index ", i)

					   raw_image <- suppressWarnings(terra::rast(src))
					   processed_image <- raw_image + vignetting_adjustment

					   terra::writeRaster(processed_image,
							      filename = dst,
							      datatype = "INT2U",
							      overwrite = overwrite_dst)
			     },
			     src_files = src_files,
			     dst_files = dst_files,
			     verbose = verbose,
			     vignetting_adjustment = vignetting_adjustment,
			     overwrite_dst = overwrite_dst,
			     mc.cores = ncores)

	} else {
		stop("Only method = 'overall' is supported at the moment.")
	}

	transfer_exif(src_dir, src_ext, output_dir)

	invisible(dst_files)
}

#' Correct thermal data and fit linear models of temperature
#'
#' This function is a convenient wrapper for much of the functions
#' provided by the thermal package. From a metadata data.frame,
#' such as returned by \code{\link{read_metadata}} and other thermal
#' flight related data, it performs thermal drift or vignetting
#' correction in response to user-supplied parameters and returns
#' a data.frame of metadata on the corrected pictures. Optionally,
#' thermal models can also be fitted using the corrected data if
#' panel coordinates are supplied.
#'
#' More details on correction methods should be added here.
#'
#' @param base_data A metadata data.frame of the pictures to correct,
#' such as returned by \code{\link{read_metadata}}. Must have row names
#' corresponding to names of list elements in panels if thermal models
#' are to be fitted.
#' @param correction_type A chracter. The type of correction to apply to
#' the images. At the moment, the values "overall", "overlap" and
#' "vignetting_overlap" are supported.
#' @param output_dir A character. The directory to which the corrected
#' images should be output.
#' @param camera_tz A character that can be interpreted as a valid timezone.
#' This is the timezone of the original pictures taken by the camera.
#' @param display_tz A character that can be interpreted as a valid timezone.
#' This is the timezone that will be used for the output metadata, and therefore
#' the timezone that the user intends to use for downstream analyses.
#' @param tags A character vector of tags to read from the metadata or
#' a single character string that identifies a set of vectors. See
#' \code{\link{read_metadata}} and \code{\link{exif_tags}} for more details.
#' @param tparams A set of transformation parameters, as determined by
#' \code{\link{optimize_transform}}. See \code{\link{add_tparams}} for more
#' details. Only needs to be specify for methods based on overlap correction.
#' @param panels A list of sf polygons containing the coordinates of
#' reference (black, gray and white) thermal panels in the images taken
#' during the flight. The names of the list must correspond to row names
#' in the input data such that proper matching can be done. If NULL,
#' then thermal linear models will not be fitted.
#' @param temperature A list of three data.frames with the temperature
#' data for each of the panels, such as returned by \code{\link{read_temp}}.
#' @param use_panels A character vector with the names of the panels to
#' use for fitting the linear models. See \code{\link{thermal_lm}}.
#' @param midpoint The index of the image in metadata to use as a reference
#' for adjusting the mean of the thermal pictures if drift correction is
#' performed. See \code{\link{correct_drift}}. If NULL (the default), then
#' the median image is used.
#' @param overwrite_dst A logical. Whether overwriting destination files
#' (the corrected images) should be allowed.
#' @param rownames_tolerance A numeric interpreted as the difference
#' (in seconds) that is allowed between the timestamps of the original
#' and corrected pictures for the row names to be transferred between
#' them. Mostly used as a safety measure to prevent bad transfers of
#' row names.
#' @param verbose A logical. Whether the function should output information
#' on its progres.
#' @param ncores An single integer value specifying the number of cores
#' to use for specific parts of the workflow. Defaults to 1 (no multithreading).
#'
#' @return A data.frame of metadata describing the corrected pictures and
#' (optionally) thermal linear models fit using those pictures.
#'
#' @export
#' @examples
#' NULL
correct_thermal <- function(base_data, correction_type, output_dir,
			    camera_tz, display_tz, tags = "minimal",
			    tparams = NULL, panels = NULL,
			    temperature = NULL, use_panels = NULL,
			    midpoint = NULL, overwrite_dst = FALSE,
			    rownames_tolerance = 1.5,
			    verbose = TRUE, ncores = 1) {

	# base_data must be a data.frame of preprocessed metadata
	# it is the metadata on the source files to be used for correction
	stopifnot(is.data.frame(base_data))

	# We add the transformation parameters to the data.frame
	# but first we need to remove them if they already exist
	for(i in c("theta", "htrans", "vtrans")) base_data[[i]] <- NULL
	if(!is.null(tparams)) base_data <- add_tparams(base_data, tparams)

	# Creating the output directory
	dir.create(output_dir, recursive = TRUE)

	# Next we check for the correction that was requested
	if(correction_type == "overall") {
		# In this case we need to artificially set the transformation parameters to 0
		base_data$theta <- 0
		base_data$htrans <- 0
		base_data$vtrans <- 0

		# Then we compute the differences based on those parameters
		base_data$diff <- compute_diffs(base_data, ncores = ncores)

		# We then run the corretion routine, which returns the names of the modified files
		corrected_files <- correct_drift(metadata = base_data, output_dir = output_dir,
						 method = "overlap", midpoint = midpoint,
						 overwrite_dst = overwrite_dst, ncores = ncores,
						 verbose = verbose)

	} else if(correction_type == "overlap") {
		# Same as "overall" but we keep the parameters provided
		# Then we compute the differences based on those parameters
		base_data$diff <- compute_diffs(base_data, ncores = ncores)

		# We then run the corretion routine, which returns the names of the modified files
		corrected_files <- correct_drift(metadata = base_data, output_dir = output_dir,
						 method = "overlap", midpoint = midpoint,
						 overwrite_dst = overwrite_dst, ncores = ncores,
						 verbose = verbose)

	} else if(correction_type == "vignetting_overall") {

		# Vignetting is computed and corrected on the fly
		corrected_files <- correct_vignetting(base_data, output_dir, method = "overall",
						      overwrite_dst = overwrite_dst,
						      ncores = ncores, verbose = verbose)
	} else {
		stop("Unsupported correction type")
	}

	# Now we need to re-read the metadata from the files we just created
	corrected_meta <- read_metadata(corrected_files, camera_tz = camera_tz, display_tz = display_tz, tags = tags)

	# Compute the thermal models if panels are provided
	if(!is.null(panels)) {
		corrected_meta <- add_temp_metadata(metadata = corrected_meta,
						    temperature_list = temperature,
						    tolerance = as.difftime(10, units = "secs"))

		# Transferring the rownames of the original data to the corrected data
		stopifnot(all.equal(base_data$DateTimeOriginal, corrected_meta$DateTimeOriginal, tolerance = rownames_tolerance))
		rownames(corrected_meta) <- rownames(base_data)

		pixels <- join_thermal(corrected_meta, panels, ncores = ncores)
		corrected_meta <- thermal_lm(corrected_meta, pixels, use_panels = use_panels)$metadata
	}

	corrected_meta
}
