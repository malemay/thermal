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

