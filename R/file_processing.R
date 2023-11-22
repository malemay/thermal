#' Read image metadata
#' 
#' Read EXIF metadta for a set of images
#'
#' @param files a character vector of files for which exif metadata
#' is to be read
#' @param tags A character vector of tags to extract or a single character
#' value denoting a set of tags to extract. Supported values for single
#' character are: "default" (set of tags that may be interesting or useful),
#' "minimal" (minimal set of tags needed for downstream analyses). If NULL (the default),
#' then all tags are extracted.
#' @param return_df A logical, whether to return a data.frame instead of a
#' 	            tibble (defaults to TRUE)
#'
#' @return A data.frame with the metadata for each file
#'
#' @examples
#' NULL
#'
#' @export
read_metadata <- function(files, tags = NULL, return_df = TRUE) {

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

	output <- exifr::read_exif(files, tags = tags)

	# Format some of the columns
	if("ModifyDate" %in% colnames(output)) {
		output$ModifyDate <- strptime(output$ModifyDate, format = "%Y:%m:%d %H:%M:%S")
	}

	if("CreateDate" %in% colnames(output)) {
		output$CreateDate <- strptime(output$CreateDate, format = "%Y:%m:%d %H:%M:%S")
	}

	if("DateTimeOriginal" %in% colnames(output)) {
		output$DateTimeOriginal <- strptime(output$DateTimeOriginal, format = "%Y:%m:%d %H:%M:%S")
	}

	# Sanity checks
	stopifnot(max(abs(output$CreateDate - output$DateTimeOriginal)) < 2)

	# Sorting according to creation date
	output <- output[order(output$CreateDate), ]

	if(return_df) as.data.frame(output) else output
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

