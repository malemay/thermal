# Setup variables that will be used throughout testing
thermal_files <- dir(system.file("extdata", package = "thermal"),
                     pattern = "thermal\\.tiff$", full.names = TRUE)

tmeta <- read_metadata(thermal_files, camera_tz = "Etc/GMT+5", display_tz = "Etc/GMT+4", tags = "minimal")
