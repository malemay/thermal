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

