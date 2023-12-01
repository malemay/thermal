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
#' NULL
#'
read_metadata <- function(files, camera_tz, display_tz = NULL, tags = NULL) {

	# Determining which EXIF tags to extract from the file
	if(!is.null(tags)) {
		if(length(tags) == 1) {
			if(tags == "default") {
				tags <- exif_tags("default")
			} else if(tags == "minimal") {
				tags <- exif_tags("minimal")
			}
		}
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
#' with their correspondence
#'
#' @param visible A data.frame with metadata on visible files
#' @param thermal A data.frame with metadata on thermal files
#'
#' @return A data.frame with the correspondence between the images
#'
#' @examples
#' NULL
#'
#' @export
match_images <- function(visible, thermal) {
	output <- data.frame(visible_file = visible$SourceFile)
	output$thermal_file <- NA
	output$timediff <- NA

	for(i in 1:nrow(visible)) {
		time_diff <- abs(visible[i, "DateTimeOriginal"] - thermal$DateTimeOriginal)
		output[i, "thermal_file"] <- thermal[which.min(time_diff), "SourceFile"]
		output[i, "timediff"] <- min(time_diff)
	}

	output
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
		output <- c("ExifToolVersion",
			    "FileModifyDate", "FileAccessDate", "FileType", "FileTypeExtension",
			    "Make", "Model", "Orientation", "XResolution", "YResolution", "ResolutionUnit",
			    "ModifyDate", "AbsoluteAltitude", "RelativeAltitude",
			    "GimbalRollDegree", "GimbalYawDegree", "GimbalPitchDegree",
			    "FlightRollDegree", "FlightYawDegree", "FlightPitchDegree",
			    "CreateDate", "FocalLength", "ExifImageWidth", "ExifImageHeight",
			    "FocalPlaneXResolution", "FocalPlaneYResolution", "FocalPlaneResolutionUnit",
			    "DateTimeOriginal", "GPSMapDatum", "GPSAltitude", "GPSLatitude", "GPSLongitude")
	} else if(option == "minimal") {
		output <- c("GimbalYawDegree", "CreateDate",
			    "DateTimeOriginal", "GPSAltitude", "GPSLatitude",
			    "GPSLongitude")
	}

	return(output)
}

#' Transfer exif metadata from the original to corrected files
#'
#' @param src_dir The directory where the source files are found.
#' @param src_ext The extension of the source files (including the dot prefix).
#' @param dst_dir The directory containing the files to write the metadata to.
#'
#' @return NULL, invisibly
#'
#' @examples
#' NULL
transfer_exif <- function(src_dir, src_ext, dst_dir) {

	stopifnot(all(dir.exists(c(src_dir, dst_dir))))

	# Using the batch format syntax to apply the command to all files at once
	exifr::exiftool_call(args = paste0('-overwrite_original -tagsFromFile ',
				    src_dir, "/%f", src_ext, " ",
				    paste0("-", thermal:::exif_tags("minimal"), collapse = " ")),
			     fnames = dst_dir)
}

