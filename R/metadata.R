#' Read image metadata
#' 
#' Read EXIF metadta for a set of images
#'
#' @param files a character vector of files for which exif metadata
#' is to be read.
#' @param camera_tz A character that can be interpreted as a valid timezone
#' corresponding to that of the date/time values written by the camera.
#' It is mandatory to specify this value so that downstream functions
#' know how to interpret the time values of the camera.
#' @param display_tz A character that can be interpreted as a valid timezone
#' to use when displaying the data. If NULL (the default), then the values
#' are the same as camera_tz.
#' @param tags A character vector of tags to extract or a single character
#' string denoting a set of tags to extract. Supported values for single
#' string are: "default" (set of tags that may be interesting or useful),
#' "minimal" (minimal set of tags needed for downstream analyses). If NULL (the default),
#' then all tags are extracted.
#'
#' @return A data.frame with the metadata for each file.
#'
#' @export
#'
#' @examples
#' # Get the names of the thermal images in the test dataset
#' tfiles <- dir(system.file("extdata/", package = "thermal"), 
#'               pattern = "thermal.tiff$", full.names = TRUE)
#' # Reading the metadata ; we want to process the metadata in a different timezone
#' # as the camera was set to winter time but data was acquired under the summer time
#' tmeta <- read_metadata(tfiles, camera_tz = "Etc/GMT+5", display_tz = "Etc/GMT+4", tags = "minimal")
read_metadata <- function(files, camera_tz, display_tz = NULL, tags = NULL) {

	# Determining which EXIF tags to extract from the file
	if(!is.null(tags)) {
		if(length(tags) == 1 && tags %in% c("default", "minimal")) tags <- exif_tags(tags)
	} else {
		tags <- "all"
	}

	# Running exiftool to extract the tags of interest
	output <- exifr::read_exif(files, tags = tags)

	# Handling the date columns
	stopifnot(all(c("CreateDate", "DateTimeOriginal") %in% colnames(output)))

	# WARNING: this may fail if a column contains "Date" but is not a date
	for(i in grep("Date", colnames(output), value = TRUE)) {

		# The columns CreateDate and DateTimeOriginal need to be timezone-aware
		if(i %in% c("CreateDate", "DateTimeOriginal")) {

			output[[i]] <- as.POSIXct(output[[i]], format = "%Y:%m:%d %H:%M:%OS", tz = camera_tz)

			if(!is.null(display_tz) && display_tz != camera_tz) {
				attr(output[[i]], "tzone") <- display_tz
			}

		} else {
			# Time zone is not as critical for other date/time columns
			output[[i]] <- strptime(output[[i]], format = "%Y:%m:%d %H:%M:%S")
		}
	}

	# Making sure that CreateDate and DateTimeOriginal are not too far apart
	stopifnot(all(abs(difftime(output$CreateDate, output$DateTimeOriginal, units = "secs")) <= 2))

	# Sorting according to the moment when the picture was taken
	output <- output[order(output$DateTimeOriginal), ]

	as.data.frame(output)
}

#' Match thermal and visible images
#'
#' Takes metadata on thermal and visible images and returns a data.frame
#' with their correspondence.
#'
#' There may not be a universal way to match visible and thermal images
#' with 100% certainty in the quality of the output. Therefore, this
#' function provides a general framework for provding functions that
#' do the matching. It also provides sanity checks to make sure that
#' the output data makes sense. The function \code{\link{dji_filename_match}}
#' may provided an appropriate match_func for some camera models.
#' Otherwise, users will need to write their own function or file
#' a bug report to have the requested functionality added.
#'
#' @param visible A data.frame with metadata on visible files.
#' @param thermal A data.frame with metadata on thermal files.
#' @param match_func A function that takes both the visible and
#' thermal data.frames as input and returns an integer vector
#' of the indices in thermal that correspond to indices in
#' visible.
#' @param max_difftime A difftime object specifying the maximum
#' time difference between matching pictures. This function will
#' throw an error if any two matching images have timestamps whose
#' difference exceeds this parameter.
#'
#' @return A data.frame with the correspondence between the images and
#' the time difference between both pictures.
#'
#' @export
#' @examples
#' # Reading the metadata of test thermal images
#' tfiles <- dir(system.file("extdata/", package = "thermal"), 
#'               pattern = "thermal.tiff$", full.names = TRUE)
#' tmeta <- read_metadata(tfiles, camera_tz = "Etc/GMT+5",
#'                        display_tz = "Etc/GMT+4", tags = "minimal")
#'
#' # Reading the metadata of test visible images
#' vfiles <- dir(system.file("extdata/", package = "thermal"), 
#'               pattern = "visible.tiff$", full.names = TRUE)
#' vmeta <- read_metadata(vfiles, camera_tz = "Etc/GMT+5",
#'                        display_tz = "Etc/GMT+4", tags = "minimal")
#'
#' # Creating a function that matches the two files based on the fact
#' # that thermal files are numbered one less than the visible files.
#' # The function must take two metadata data.frames as input and return
#' # a vector of indices in the thermal data.frame that match the rows
#' # of the visible data.frame.
#' match_func <- function(visible, thermal) {
#'     indices <- numeric(nrow(visible))
#'     
#'     for(i in 1:nrow(visible)) {
#'         visible_num <- as.numeric(substr(basename(visible[i, "SourceFile"]), 5, 8))
#'         thermal_num <- formatC(visible_num - 1, width = 4, flag = 0)
#'         i_index <- which(grepl(thermal_num, basename(thermal$SourceFile)))
#'         stopifnot(length(i_index) == 1)
#'         indices[i] <- i_index
#'     }
#'
#'     return(indices)
#' }
#'
#' # Using the function that we just created to get the matches between the images
#' file_matches <- match_images(visible = vmeta, thermal = tmeta,
#'                              match_func = match_func,
#'                              max_difftime = as.difftime(1.5, units = "secs"))
match_images <- function(visible, thermal, match_func, max_difftime = as.difftime(1, units = "secs")) {

	# A first sanity check
	stopifnot(nrow(visible) == nrow(thermal))

	# Using the provided function to compute the matches
	indices <- match_func(visible, thermal)

	# Adding the names of the corresponding thermal images to the output
	output <- data.frame(visible_file = visible$SourceFile)
	output$thermal_file <- thermal[indices, "SourceFile"]

	# Doing some sanity checks
	stopifnot(all(file.exists(output$thermal_file)))
	stopifnot(length(output$thermal_file) == length(unique(output$thermal_file)))

	# Computing the time difference between the two files
	output$difftime <- difftime(visible[, "DateTimeOriginal"], thermal[indices, "DateTimeOriginal"], units = "secs")
	stopifnot(all(abs(output$difftime) < max_difftime))

	output
}

#' Find the thermal image paired with a given visible image based on DJI file names
#'
#' This function is meant to be provded as the match_func argument to
#' function match_images. It works for file names as written by the
#' DJI Zenmuse XT2 camera, for which the thermal images are written
#' with a number right before the one for the visible camera, but different
#' functions may need to be written for other cameras.
#'
#' @param visible A data.frame of metadata on a set of visible images,
#' such as returned by \code{\link{read_metadata}}.
#' @param thermal Similar to visible, but for a set of thermal images.
#'
#' @return An integer vector of indices in the thermal dataset that correspond
#' to rows in the visible dataset.
#'
#' @export
#' @examples
#' NULL
dji_filename_match <- function(visible, thermal) {
	# Checking that all files are in the same directory
	vdir <- unique(dirname(visible$SourceFile))
	tdir <- unique(dirname(thermal$SourceFile))
	stopifnot(length(vdir) == 1 && length(tdir) == 1 && vdir == tdir)

	# Extracting the numeric ID from the visible files
	vdji <- regmatches(basename(visible$SourceFile), regexpr("[0-9]{4}", basename(visible$SourceFile)))

	# The thermal ID is this number minus 1, which we must prefix with a 0
	tdji <- formatC(as.numeric(vdji) - 1, width = 4, flag = 0)

	# The value "0000" is a special case which must be wrapped back to "0999"
	tdji[tdji == "0000"] <- "0999"

	# Reassembling the thermal IDs into file names
	thermal_files <- paste0(tdir, "/", "DJI_", tdji, "_R.JPG")

	# At this point we must check that all files are indeed in the thermal data.frame (and the other way around)
	stopifnot(all(thermal_files %in% thermal$SourceFile) && all(thermal$SourceFile %in% thermal_files))

	# Returning the indices of thermal files in the visible data.frame
	match(thermal_files, thermal$SourceFile)
}


#' Set the EXIF tags to extract
#'
#' @param option A character indicating which set of options to choose
#'
#' @return A character vector of tags to extract.
#'
#' @examples
#' NULL
exif_tags <- function(option) {
	if(option == "default") {
		output <- c("Make", "Model", "Orientation", "XResolution", "YResolution", "ResolutionUnit",
			    "ModifyDate", "AbsoluteAltitude", "RelativeAltitude",
			    "GimbalRollDegree", "GimbalYawDegree", "GimbalPitchDegree",
			    "FlightRollDegree", "FlightYawDegree", "FlightPitchDegree",
			    "CreateDate", "FocalLength", "ExifImageWidth", "ExifImageHeight",
			    "FocalPlaneXResolution", "FocalPlaneYResolution", "FocalPlaneResolutionUnit",
			    "GPSLongitudeRef", "GPSLatitudeRef", "GPSAltitudeRef",
			    "DateTimeOriginal", "GPSMapDatum", "GPSAltitude", "GPSLatitude", "GPSLongitude")
	} else if(option == "minimal") {
		output <- c("GimbalYawDegree", "CreateDate", "DateTimeOriginal",
			    "GPSLongitudeRef", "GPSLatitudeRef", "GPSAltitudeRef",
			    "GPSLongitude", "GPSLatitude", "GPSAltitude",
			    "GPSMapDatum")
	}

	return(output)
}

#' Transfer exif metadata from the original to corrected files
#'
#' @param src_dir The directory where the source files are found.
#' @param src_ext The extension of the source files (including the dot prefix).
#' @param dst_dir The directory containing the files to write the metadata to.
#' @param tags A character vector of tags to transfer or a single character
#' string denoting a set of tags to transfer. Supported values for single string
#' are: "default" (set of tags that may be interesting or useful), "minimal"
#' (minimal set of tags needed for downstream analyses). If NULL (the default),
#' then all tags are transferred.
#' @param verbose A logical indicating whether informative messages should be
#' printed out to the console. This controls both the output of the command sent
#' through \code{\link[exifr]{exiftool_call}} and the messages output by exiftool
#' itself. It does not suppress warnings.
#'
#' @return NULL, invisibly
#'
#' @examples
#' NULL
transfer_exif <- function(src_dir, src_ext, dst_dir, tags = NULL, verbose = TRUE) {

	stopifnot(all(dir.exists(c(src_dir, dst_dir))))

	# Determining which EXIF tags to extract from the file
	if(!is.null(tags)) {
		if(length(tags) == 1 && tags %in% c("default", "minimal")) tags <- exif_tags(tags)
	} else {
		tags <- "all"
	}

	# Using the batch format syntax to apply the command to all files at once
	exifr::exiftool_call(args = paste0(ifelse(verbose, "", "-q "),
					   '-overwrite_original -tagsFromFile ',
					   src_dir, "/%f", src_ext, " ",
					   paste0("-", tags, collapse = " ")),
			     fnames = dst_dir, quiet = !verbose)
}

