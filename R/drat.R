#' Compute representative digital numbers (RDN) according to the DRAT method
#'
#' @param metadata A data.frame of metadata on a set of thermal pictures.
#' @param peak A character indicating which peak to select from the local
#' maxima of the digital number probability density function. Supported values
#' are currently "minimum", "median", and "maximum", which respectively select
#' the peak with the lowest, median, and highest DN.
#' @param max_threshold A numeric value between 0 and 1, interpreted as a
#' proportion. Any local maximum must be at least that much relative to the
#' global maximum to be considered for selection as the RDN. Defaults to 0.1.
#' @param verbose A logical value. Whether or not progress and RDN values are
#' printed as part of the computation.
#' @param ncores An integer value, the number of cores to use for computing the
#' probability density functions in parallel. Defaults no 1 (no parallel
#' computation).
#'
#' @return A list of three elements:
#' - A vector of numeric values containing the representative DNs (RDN) for
#' each image.
#' - A vector of numeric values corresponding to the density for the selected
#' RDNs.
#' - The probability density function computed from the DNs, used for plotting
#' and diagnostics.
#' @noRd
compute_rdn <- function(metadata, peak, max_threshold = 0.1, verbose = TRUE, ncores = 1) {
	# Computing the densities for each image
	densities <- parallel::mclapply(1:nrow(metadata), function(i, meta) {
						if(verbose) message("compute_rdn: processing row ", i, " out of ", nrow(meta))
						values <- terra::values(terra::rast(meta[i, "SourceFile"], noflip = TRUE, subds = 1))
						stats::density(values)}, meta = metadata, mc.cores = ncores)

	# Compute the position and density value of all local maxima in the distribution
	maxima <- lapply(densities, function(x) {
				 # Getting the indexes of the local maxima 
				 indexes <- local_maxima(x)

				 # Extracting the global maximum and subsetting the indexes
				 # based on the relative importance of the local maxima
				 # A local maximum is kept if it is greater than some proportion
				 # max_threshold compared to the global maximum
				 yvals <- x$y[indexes]
				 global_maximum <- max(yvals)
				 indexes <- indexes[yvals >= max_threshold * global_maximum]

				 # Re-subsetting the x and y values prior to output
				 xvals <- x$x[indexes]
				 yvals <- x$y[indexes]
				 list(x = xvals, y = yvals)
				    })

	# Extract the RDN and density of the selected peak
	if(peak == "minimum") {
		rdn <- lapply(maxima, function(x) list(x = x$x[1], y = x$y[1]))
	} else if(peak == "maximum") {
		rdn <- lapply(maxima, function(x) list(x = x$x[length(x$x)], y = x$y[length(x$y)]))
	} else if(peak == "median") {
		rdn <- lapply(maxima, function(x) {
				      med <- floor(length(x$x) / 2)
				      list(x = x$x[med], y = x$y[med])
				    })
	}

	# Preparing the output list
	output <- list(rdn = sapply(rdn, function(x) x$x),
		       y   = sapply(rdn, function(x) x$y),
		       density = densities)

	# Printing the RDNs if verbose is TRUE
	if(verbose) {
		message("RDN values identified for this flight: ")
		verbose_df <- data.frame(file = basename(metadata$SourceFile), RDN = output$rdn)
		print(verbose_df)
	}

	return(output)
}

#' Find the index of local maxima in a probability density function
#'
#' This function accepts an object returned by stats::density and returns the
#' indexes of the local maxima in this distribution using a very simple
#' algorithm which compares each value in the previous and the next one.  Since
#' this algorithm is very sensitive to local variations it may be necessary to
#' filter the local maxima found for downstream applications, as compute_rdn
#' does.
#'
#' @param x An object returned by stats::density with element x as the values
#' considered and y being the density of that particular value.
#'
#' @return A vector of integer values representing the indexes of the local
#' maxima in the probability density function.
#'
#' @noRd
local_maxima <- function(x) {
	# Extacting the density values
	y <- x$y
	n <- length(y)

	# A vector which identifies whether a value is greater than the preceding one
	greater_than_previous <- y[2:(n - 1)] > y[1:(n - 2)]

	# Same thing when comparing to the next value
	greater_than_next <- y[2:(n - 1)] > y[3:n]

	# The index will be one more than the index in the vector because we did not consider the first value
	maxima <- which(greater_than_previous & greater_than_next) + 1

	return(maxima)
}

#' Plot the probability density functions and RDNs of a dataset
#' 
#' This function plots the probability density functions of a dataset
#' of thermal images and displays them all on the same plot with the
#' location of the representative digital numbers (RDNs) identified
#' according to the DRAT method of Wang et al. (2023b)
#'
#' @param metadata A data.frame of metadata on a flight, as returned by
#' \code{\link{read_metadata}}.
#' @param drat_peak A character indicating which peak to select from the local
#' maxima of the digital number probability density function. Supported values
#' are currently "minimum", "median", and "maximum", which respectively select
#' the peak with the lowest, median, and highest DN.
#' @param drat_threshold A numeric value between 0 and 1, interpreted as a
#' proportion. Any local maximum must be at least that much relative to the
#' global maximum to be considered for selection as the RDN. Defaults to 0.1.
#' @param alpha A numeric value between 0 and 1. The transparency value to use
#' for plotting. Values closer to 0 mean more transparency.
#' @param verbose A logical value. Whether or not progress and RDN values are
#' printed as part of the computation.
#' @param ncores An integer value, the number of cores to use for computing the
#' probability density functions in parallel. Defaults no 1 (no parallel
#' computation).
#' @param ... Other graphical parametres passed to the main plotting function
#'
#' @return Invisibly, a list of three elements:
#' \itemize{
#' \item rdn: a vector of representative digital numbers
#' \item y: a vector of probability density values for those representative
#' digital numbers
#' \item density: the probability density functions computed for each image in
#' the dataset
#' }
#' 
#' @references
#' Wang, Z., Zhou, J., Ma, J., Wang, Y., Liu, S., Ding, L., Tang, W., Pakezhamu,
#' N. and Meng, L., 2023. Removing temperature drift and temporal variation in
#' thermal infrared images of a UAV uncooled thermal infrared imager. ISPRS
#' Journal of Photogrammetry and Remote Sensing, 203, pp.392-411.
#' \url{https://doi.org/10.1016/j.isprsjprs.2023.08.011}
#'
#' @export
#' @examples
#' # Reading the metadata from the test dataset
#' tfiles <- dir(system.file("extdata/", package = "thermal"), 
#'               pattern = "thermal.tiff$", full.names = TRUE)
#' tmeta <- read_metadata(tfiles, camera_tz = "Etc/GMT+5", display_tz = "Etc/GMT+4", tags = "minimal")
#'
#' # Plotting the RDNs in the test dataset using the default threshold
#' plot_rdn(tmeta, drat_peak = "minimum", drat_threshold = 0.1)
#'
#' # Using a ridiculously low threshold to show that earlier peaks are detected
#' plot_rdn(tmeta, drat_peak = "minimum", drat_threshold = 0.01)
#'
#' # Using the same threshold with drat_peak = "maximum"
#' plot_rdn(tmeta, drat_peak = "maximum", drat_threshold = 0.01)
#'
#' # Using the same threshold with drat_peak = "median"
#' plot_rdn(tmeta, drat_peak = "median", drat_threshold = 0.01)
plot_rdn <- function(metadata, drat_peak, drat_threshold = 0.1, alpha = 0.5, verbose = FALSE, ncores = 1, ...) {
	# Computing the RDNs from the metadata
	rdn <- compute_rdn(metadata, peak = drat_peak, max_threshold = drat_threshold, verbose = verbose, ncores = ncores)

	# Extracting the range of x and y values to set up the plot limits
	xrange <- range(sapply(rdn$density, '[[', "x"))
	yrange <- range(sapply(rdn$density, '[[', "y"))

	# Creating a dummy plotting window with the right limits
	plot(0, ylab = "Density", xlab = "Digital number", type = "n", xlim = xrange, ylim = yrange, ...)

	# Computing the hex value for the color based on the alpha value
	alpha <- floor(alpha * 255)
	hexcol <- paste0("#000000", as.character(as.hexmode(alpha)))

	# Plotting all the lines from the density probability functions
	for(i in 1:length(rdn$density)) {
		graphics::lines(rdn$density[[i]]$x, rdn$density[[i]]$y, col = hexcol)
	}

	# Adding the points for the RDNs
	graphics::points(x = rdn$rdn, y = rdn$y, col = "red")

	return(invisible(rdn))
}

