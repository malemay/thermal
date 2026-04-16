## ----options, include = FALSE-------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6,
  fig.align = "center")

## ----read-metadata------------------------------------------------------------
# Loading the package
library(thermal)

# The functions would work without explicitly loading the sf package
# However, this will enhance the display of spatial objects in this vignette
library(sf)

# Getting a vector of absolute paths to the test TIR images
thermal_files <- dir(system.file("extdata", package = "thermal"),
                     pattern = "thermal\\.tiff$", full.names = TRUE)

# Reading metadata from EXIF tags
# Setting the proper display time zone is critical here because we need to
# match the time zone of the separately recorded surface temperature data
tmeta_df <- read_metadata(thermal_files,
                          camera_tz = "Etc/GMT+5",
                          display_tz = "Etc/GMT+4",
                          tags = "minimal")

# Providing IDs for each picture as row names
# This is necessary to match the panel coordinates with each picture
rownames(tmeta_df) <- substr(basename(tmeta_df$SourceFile), 1, 8)

head(tmeta_df)

## ----panels-------------------------------------------------------------------
length(panels)
head(panels, n = 3)

## ----temperature--------------------------------------------------------------
str(temperature)

## ----plot-temp, fig.width = 8, fig.height = 6---------------------------------
plot_temp(temperature,
          lcol = c("black" = "black", "gray" = "gray", "white" = "blue"),
          legend.pos = "topright")

## ----plot-panels, fig.width = 8, fig.height = 6-------------------------------
# We can use identifiers to index the metadata because we have set row names
terra::plot(terra::rast(tmeta_df["DJI_0671", "SourceFile"], noflip = TRUE))

# The argument add = TRUE plots the shapes on top of the raster
plot(panels[["DJI_0671"]],
     col = "transparent",
     border = c("black", "gray", "white"),
     lwd = 2, add = TRUE)

## ----add-temp-metadata--------------------------------------------------------
# Will return an error if any picture was taken more than 5 seconds apart from
# the closest temperature measurement
tmeta_df <- add_temp_metadata(tmeta_df, temperature, tolerance = as.difftime(5, units = "secs"))

head(tmeta_df)

## ----thermal-lm---------------------------------------------------------------
# Parallel processing is supported through argument ncores
# However, this does not work on Windows as it uses function mclapply from the
# parallel package
thermal_models <- thermal_lm(metadata = tmeta_df,
                             polygons = panels,
                             summary_functions = median,
                             surfaces = c("black", "gray", "white"),
                             ncores = if(.Platform$OS.type == "unix") 2 else 1)

## ----thermal-lm-metadata------------------------------------------------------
head(thermal_models$metadata)

## ----thermal-model-output-----------------------------------------------------
# Looking at the linear model for the first image
thermal_models$models[["DJI_0671"]]

## ----thermal-model-pixels-----------------------------------------------------
head(thermal_models$pixels)

## ----plot-pixtemp, fig.width = 8, fig.height = 8------------------------------
# We will display data for the first 9 images
par(mfrow = c(3, 3))

plot_pixtemp(model_data = thermal_models,
             id = names(panels)[1:9],
             lcol = c(black = "black", gray = "gray", white = "blue"))

## ----plot-pixtemp-alternative, fig.width = 8, fig.height = 8------------------
alternative_models <- thermal_lm(metadata = tmeta_df,
                                 polygons = panels,
                                 summary_functions = c(white = mean,
                                                       gray = median,
                                                       black = identity),
                                 surfaces = c("black", "gray", "white"),
                                 ncores = if(.Platform$OS.type == "unix") 2 else 1)

par(mfrow = c(3, 3))

plot_pixtemp(model_data = alternative_models,
             id = names(panels)[1:9],
             lcol = c(black = "black", gray = "gray", white = "blue"))

## ----plot-pixtemp-white-black, fig.width = 8, fig.height = 8------------------
bw_models <- thermal_lm(metadata = tmeta_df,
                        polygons = panels,
                        summary_functions = median,
                        surfaces = c("black", "white"),
                        ncores = if(.Platform$OS.type == "unix") 2 else 1)

par(mfrow = c(3, 3))

plot_pixtemp(model_data = bw_models,
             id = names(panels)[1:9],
             lcol = c(black = "black", white = "blue"))

## ----global-pixtemp-----------------------------------------------------------
global_model <- global_lm(thermal_models)

plot_pixtemp(model_data = global_model,
             lcol = c(black = "black", gray = "gray", white = "blue"))

## ----lm-drift, fig.width = 6, fig.height = 8----------------------------------
# We rename tmeta_df to the metadata stored by the output of thermal_models
tmeta_df <- thermal_models$metadata

par(mfrow = c(3, 1))

# A plot of linear model intercept over time
plot(tmeta_df$DateTimeOriginal, tmeta_df$intercept,
     main = "Intercept", xlab = "Time", ylab = "Intercept")
lines(tmeta_df$DateTimeOriginal, tmeta_df$intercept, lty = 2)

# A plot of linear model slope over time
plot(tmeta_df$DateTimeOriginal, tmeta_df$slope,
     main = "Slope", xlab = "Time", ylab = "Slope")
lines(tmeta_df$DateTimeOriginal, tmeta_df$slope, lty = 2)

# A plot of temperature predicted from a DN of 3500 from models fitted at various time points
plot(tmeta_df$DateTimeOriginal,  3500 * tmeta_df$slope + tmeta_df$intercept,
     main = "Prediction", xlab = "Time", ylab = "Temp. predicted when DN = 3500 (°C)")
lines(tmeta_df$DateTimeOriginal, 3500 * tmeta_df$slope + tmeta_df$intercept, lty = 2)

## ----lm-correction-calibration, fig.width = 6, fig.height = 6-----------------
corrected_df <- read_metadata(thermal_files,
                              camera_tz = "Etc/GMT+5",
                              display_tz = "Etc/GMT+4",
                              tags = "minimal")

corrected_df <- correct_thermal(metadata = corrected_df,
                                correction_type = "lm",
                                output_dir = tempdir(),
                                camera_tz = "Etc/GMT+5",
                                display_tz = "Etc/GMT+4",
                                tags = "minimal",
                                nuc_threshold = 80,
                                overwrite_dst = TRUE,
                                verbose = FALSE)

rownames(corrected_df) <- substr(basename(corrected_df$SourceFile), 1, 8)

corrected_df <- add_temp_metadata(corrected_df, temperature,
                                  tolerance = as.difftime(5, units = "secs"))

corrected_models <- thermal_lm(metadata = corrected_df,
                               polygons = panels,
                               summary_functions = median,
                               surfaces = c("black", "gray", "white"),
                               ncores = if(.Platform$OS.type == "unix") 2 else 1)

# Extracting the metadata from the output of thermal_lm
corrected_df <- corrected_models$metadata

# A plot of temperature predicted from a DN of 3500 from models fitted at various time points
plot(corrected_df$DateTimeOriginal,  3500 * corrected_df$slope + corrected_df$intercept,
     main = "Prediction", xlab = "Time", ylab = "Temp. predicted when DN = 3500 (°C)")
lines(x = corrected_df$DateTimeOriginal,
      y = 3500 * corrected_df$slope + corrected_df$intercept,
      lty = 2)

## ----radiometric-calibration-pixtemp------------------------------------------
global_corrected <- global_lm(corrected_models)

plot_pixtemp(global_corrected,
             lcol = c(black = "black", gray = "gray", white = "blue"))


## ----radiometric-calibration-rast, fig.width = 8, fig.height = 8--------------
par(mfrow = c(3, 3))

for(i in rownames(corrected_df)[1:9]) {
        i_rast <- terra::rast(corrected_df[i, "SourceFile"], noflip = TRUE)
        i_rast <- thermal_predict(x = i_rast, model = global_corrected$models$global)
        terra::plot(i_rast, main = i)
}

## ----session-info-------------------------------------------------------------
sessionInfo()

