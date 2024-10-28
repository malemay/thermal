# Tests for the code in R/metadata.R

# Tests of the read_metadata function
test_that("read_metadata() works properly", {
		  tmeta <- read_metadata(thermal_files, camera_tz = "Etc/GMT+5", tags = "minimal")

		  # Test that a data.frame is returned
		  expect_s3_class(tmeta, "data.frame")

		  # Test that the results are returned sorted by time
		  expect_true(!is.unsorted(tmeta$DateTimeOriginal))

		  # Test that tags = "minimal" returns the expected names
		  expect_named(tmeta, c("SourceFile", "GimbalYawDegree", "CreateDate", "DateTimeOriginal",
					"GPSLongitudeRef", "GPSLatitudeRef", "GPSAltitudeRef", "GPSLongitude",
					"GPSLatitude", "GPSAltitude", "GPSMapDatum"))

		  # Test that a camera time zone needs to be supplied
		  expect_error(read_metadata(thermal_files, tags = "minimal"))
})
