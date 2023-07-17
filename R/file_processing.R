#' Read image metadata
#' 
#' Read EXIF metadta for a set of images
#'
#' @param files a character vector of files for which exif metadata
#'        is to be read
#' @param exclude a character of columns to remove from the output
#' @param return_df A logical, whether to return a data.frame instead of a
#' 	            tibble (defaults to TRUE)
#'
#' @return A data.frame with the metadata for each file
#'
#' @examples
#' NULL
#'
#' @export
read_metadata <- function(files, exclude = c("ThumbnailImage", "RawThermalImage"), return_df = TRUE) {
	output <- exifr::read_exif(files)
	output <- output[, ! colnames(output) %in% exclude]

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

