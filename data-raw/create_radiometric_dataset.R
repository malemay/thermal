# This script creates the datasets used for the vignette on radiometric
# calibration

# Loading the required packages
library(usethis)
library(sf)

# A function to read temperature datasets
read_temp <- function(filename, tz) {
	tempdata <- utils::read.table(filename, sep = ",", skip = 4, na.strings = c("NA", "NAN"))
	colnames(tempdata) <- c("time", "record", "battery", "devtemp", "temp")
	tempdata$time <- as.POSIXct(tempdata$time, format = "%Y-%m-%d %H:%M:%S", tz = tz)
	tempdata
}

# Reading and formatting the temperature data
black <- read_temp("~/Documents/gennaretti_lab/data/panel_temperature/07092023/BlackPanel_Table1.dat", tz = "Etc/GMT+4")[, c("time", "temp")]
gray <- read_temp("~/Documents/gennaretti_lab/data/panel_temperature/07092023/GrayPanel_Table1.dat", tz = "Etc/GMT+4")[, c("time", "temp")]
white <- read_temp("~/Documents/gennaretti_lab/data/panel_temperature/07092023/WhitePanel_Table1.dat", tz = "Etc/GMT+4")[, c("time", "temp")]

temperature <- list(black = black, gray = gray, white = white)
temperature <- lapply(temperature, function(x) x[x$time >= "2023-09-01 00:00:00" & x$time < "2023-09-02 00:00:00", ])

# Saving the dataset
usethis::use_data(temperature, overwrite = TRUE)



# Reading and formatting the panel data
panels <- readRDS("~/Documents/gennaretti_lab/analyses/thermal_correction/data/processed_undistorted/20230901_12/thermal_analysis/final_panels.rds")

# Reading the names of the images that are included in the package
thermal_files <- dir(system.file("extdata", package = "thermal"), pattern = "thermal\\.tiff$")
file_ids <- substr(thermal_files, 1, 8)

panels <- panels[names(panels) %in% file_ids]

# A function that divides the coordinates by a given values and returns an updated
# polygon object
update_polygon <- function(x, fact = 8) {
	xcoords <- st_coordinates(x)

	new_poly <- lapply(split(as.data.frame(xcoords[, 1:2] / fact), xcoords[, 4]),
	       function(x) sf::st_polygon(list(as.matrix(x))))

	new_poly <- do.call(sf::st_sfc, new_poly)

	sf::st_geometry(x) <- new_poly

	x
}

panels <- lapply(panels, update_polygon)

usethis::use_data(panels, overwrite = TRUE)
