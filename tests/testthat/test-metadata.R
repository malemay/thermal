# Tests for the code in R/metadata.R

# The tmeta object is created in the helper script

# Tests of the read_metadata function
test_that("read_metadata() works properly", {

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

# Tests of the transfer_exif function
test_that("transfer_exif works properly", {
		  # We transfer the EXIF tags of the first thermal test files
		  # We need temporary directories to write the files to
		  tmpdir1 <- withr::local_tempdir("test")
		  tmpdir2 <- withr::local_tempdir("test")
		  tmpdir3 <- withr::local_tempdir("test")
		  tmpdir4 <- withr::local_tempdir("test")

		  # A data.frame of the original metadata
		  original_metadata <- read_metadata(thermal_files[1:5], camera_tz = "Etc/GMT+5")

		  # Copying only the image data of the files
		  for(i in thermal_files[1:5]) {
			  terra::writeRaster(terra::rast(i),
					     filename = paste0(tmpdir1, "/", basename(i)),
					     datatype = "INT2U")
		  }

		  # Also copying the files to other temporary directories on which other tests will be done
		  file.copy(dir(tmpdir1, full.names = TRUE), tmpdir2)
		  file.copy(dir(tmpdir1, full.names = TRUE), tmpdir3)
		  file.copy(dir(tmpdir1, full.names = TRUE), tmpdir4)

		  # Reading the metadata prior to tag update
		  intermediate_metadata <- as.data.frame(exifr::read_exif(dir(tmpdir1, full.names = TRUE)))

		  # --- FIRST CALL to transfer_exif: default parameters
		  transfer_exif(src_dir = unique(dirname(thermal_files)),
				src_ext = ".tiff",
				dst_dir = tmpdir1,
				verbose = FALSE)

		  # Reading the metadata after the tag update
		  updated_metadata <- read_metadata(dir(tmpdir1, full.names = TRUE), camera_tz = "Etc/GMT+5")

		  # We do not test the "GPSPosition" column because numerical inaccuracies
		  # can occur and since this is a string column, it cannot be numerically compared
		  # The GPSLatitude column should be able to assess whether these are similar enough
		  added_tags <- colnames(updated_metadata)[!colnames(updated_metadata) %in% colnames(intermediate_metadata)]
		  added_tags <- added_tags[added_tags != "GPSPosition"]

		  # --- SECOND CALL to transfer_exif: specifying a set of tags to transfer
		  time_tags <- c("CreateDate", "DateTimeOriginal")
		  transfer_exif(src_dir = unique(dirname(thermal_files)),
				src_ext = ".tiff",
				dst_dir = tmpdir2,
				tags = time_tags,
				verbose = FALSE)

		  # Reading the metadata after the second call to transfer_exif
		  time_metadata <- read_metadata(dir(tmpdir2, full.names = TRUE), camera_tz = "Etc/GMT+5")

		  # --- THIRD CALL to transfer_exif: transferring the default set of tags
		  transfer_exif(src_dir = unique(dirname(thermal_files)),
				src_ext = ".tiff",
				dst_dir = tmpdir3,
				tags = "default",
				verbose = FALSE)

		  default_metadata <- read_metadata(dir(tmpdir3, full.names = TRUE), camera_tz = "Etc/GMT+5")

		  # --- FOURTH CALL to transfer_exif: transferring the minimal set of tags
		  transfer_exif(src_dir = unique(dirname(thermal_files)),
				src_ext = ".tiff",
				dst_dir = tmpdir4,
				tags = "minimal",
				verbose = FALSE)

		  minimal_metadata <- read_metadata(dir(tmpdir4, full.names = TRUE), camera_tz = "Etc/GMT+5")

		  # EXPECTATIONS
		  # FIRST CALL
		  # We check that the number of tags in both original and updated files are the same
		  # when transferring all tags (the default)
		  expect_equal(ncol(original_metadata), ncol(updated_metadata))
		  # We check that the transferred metadata tags have the same value
		  expect_equal(original_metadata[, added_tags], updated_metadata[, added_tags])

		  # SECOND CALL
		  # Checking that the transferred time tags are the same
		  expect_identical(original_metadata[, time_tags], time_metadata[, time_tags])

		  # THIRD CALL
		  expect_equal(original_metadata[, exif_tags("default")], default_metadata[, exif_tags("default")])

		  # FOURTH CALL
		  minimal_tags <- exif_tags("minimal")[exif_tags("minimal") != "GPSPosition"]
		  expect_equal(original_metadata[, minimal_tags], minimal_metadata[, minimal_tags])

		  expect_lt(ncol(intermediate_metadata), ncol(time_metadata))
		  expect_lt(ncol(time_metadata), ncol(minimal_metadata))
		  expect_lt(ncol(minimal_metadata), ncol(default_metadata))
		  expect_lt(ncol(default_metadata), ncol(updated_metadata))

})

