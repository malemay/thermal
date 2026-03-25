#' Location of reference panels on test dataset pictures
#'
#' This dataset provides the spatial location of the panels used as surface
#' temperature references for the test images used in the vignette and
#' examples. It can be used by users to familiarize themselves with the
#' package's functions and to check that their input data matches the package's
#' requirements. Polygon datasets for use with \code{\link{thermal_lm}} should
#' be formatted similarly to this dataset. Function
#' \code{\link{create_surfaces}} may help with this.
#'
#' The reference surfaces "black", "gray" and "white" are aluminum panels that
#' were painted with matte paint for concrete of the respective color to allow
#' for temperature differences to arise in response to absorbed solar
#' radiation. The panels had their temperature monitored through type T
#' thermocouples (see the \code{\link{temperature}} dataset for more details).
#'
#' @format A list of sf objects named according to their identifier in the
#' dataset, which can be used to associate them with the corresponding picture.
#' Each element contains the coordinates of the polygons corresponding to the
#' black, gray, and white panels, according to column "ID". Coordinates are in
#' pixels as the images are not georeferenced.
#'
#' @seealso \code{\link{temperature}}, \code{\link{create_surfaces}}
#' 
"panels"

#' Surface temperature of the reference panels used in test dataset
#'
#' This dataset contains the surface temperature over the course of September
#' 1, 2023 of the three reference panels used in the vignettes and examples.
#' This time period covers that of the pictures in the test dataset.
#' Temperature was recorded every five seconds using type T thermocouples. This
#' dataset can be used by users to familiarize themselves with the package's
#' functions and to check that their input data matches the package's
#' requirements. Temperature datasets for use with
#' \code{\link{add_temp_metadata}} should be formatted similarly to this
#' dataset.
#'
#' The reference surfaces "black", "gray" and "white" are aluminum panels that
#' were painted with matte paint for concrete of the respective color to allow
#' for temperature differences to arise in response to absorbed solar
#' radiation. The coordinates of the panels on the test dataset pictures can be
#' found in the \code{\link{panels}} dataset.
#'
#' @format A list of three elements named according to the panel identifier,
#' each formatted as a data.frame with the following two columns:
#' \itemize{
#'     \item time: the time when the temperature was recorded (POSIXct, Etc/GMT+4 time zone)
#'     \item temp: surface temperature in degrees Celsius
#' }
#'
#' @seealso \code{\link{panels}}, \code{\link{add_temp_metadata}}
"temperature"

