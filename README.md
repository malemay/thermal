# Overview

`thermal` is an R package for analyzing thermal infrared images acquired using
unmanned aerial vehicles (UAVs), or drones. It provides functionality for
reading such images and their metadata, correcting thermal drift using various
statistical or heuristic models, and conducting radiometric calibration of
images using the surface temperature of reference surfaces.

Thermal infrared images acquired using cameras mounted on UAVs are vulnerable
to drift, i.e. changes in sensor output over time. This drift is unrelated to
real changes in surface temperature of the objects being surveyed, which can
lead to erroneous interpretation of results if left uncorrected. The extent to
which camera firmware can compensate for these changes is limited, even when
using so-called radiometrically calibrated cameras which should enable direct
surface temperature readings

`thermal` provides functions for modelling this drift either statistically
(linear models, cubic splines) or heuristically by using the overlapping
fraction of successive pictures to adjust for drift. `thermal` implements these
models and corrections using a software framework and user interface that can
easily integrate new models for correcting drift.

Another issue with thermal infrared images is converting camera readings
(typically stored as a digital number, DN) into actual surface temperature,
which is needed for some applications or to compare surface temperature across
flights or images. `thermal` provides functionality for easily integrating
surface temperature measurements of reference surfaces used during *in situ*
flights to conduct radiometric calibration, i.e., to compute surface
temperature from the DNs stored in images.

# Installation

The package can be installed directly from GitHub. The R package `devtools`
provides convenient functionality to do so:

```r
# Simple package installation
devtools::install_github("malemay/thermal")

# Run this to also build the documentation vignettes, which we recommend
devtools::install_github("malemay/thermal", build_vignettes = TRUE)
```

# Quick start

`thermal` provides a test dataset of 30 images whose resolution was downgraded
to 80 x 64 pixels for the purpose of running tests or examples. The following
commands are based on this test dataset. First, we read the metadata on the
test images:

```r
# Loading the package
library(thermal)

# Reading the metadata of the test dataset provided with the package
tfiles <- dir(system.file("extdata/", package = "thermal"), 
              pattern = "thermal.tiff$", full.names = TRUE)

tmeta <- read_metadata(tfiles, camera_tz = "Etc/GMT+5",
                       display_tz = "Etc/GMT+4", tags = "minimal")
```

Next, we correct the images by assuming a linear model of drift. The corrected
images are stored under a temporary directory for the purposes of this example,
but they would typically be saved to permanent storage.

```r
# Creating a temporary directory for this example
tmpdir <- tempfile("corr_example")
dir.create(tmpdir)

# Computing the mean value of all images in the dataset
tmeta$mean <- thermal_mean(tmeta)

# Performing correction assuming a linear model of drift
# The object returned is a metadata data.frame with
# file paths pointing to the corrected data
lm_meta <- correct_thermal(metadata = tmeta,
                           correction_type = "lm", nuc_threshold = 80,
                           output_dir = tmpdir,
                           camera_tz = "Etc/GMT+5", display_tz = "Etc/GMT+4",
                           tags = "minimal",
                           overwrite_dst = TRUE,
                           verbose = TRUE)
```

We can fit radiometric calibration models by combining the images with
information on the surface temperature (dataset `temperature`) of reference
surfaces whose location on the images is stored in the dataset `panels`.

```r
# Adding the data on the temperature of the reference surfaces
lm_meta <- add_temp_metadata(lm_meta, temperature)

# Naming the rows of the metadata according to the image identifier
# This is necessary to match the images with the right reference surfaces
rownames(lm_meta) <- substr(basename(lm_meta$SourceFile), 1, 8)

# Computing the radiometric calibration model on all surfaces
# We use the median of the distribution of pixel values for each surface
thermal_models <- thermal_lm(metadata = lm_meta,
                             polygons = panels,
                             surfaces = c("black", "gray", "white"),
                             summary_functions = median)

# `thermal_models` stores one radiometric calibration model per image
# Here, we want to combine all this data in a single model
global_model <- global_lm(thermal_models)
```

Once we have a radiometric calibration model, we can use it to compute surface
temperature for any picture for which we believe that this model is accurate.

```r
# Reading the first picture of the dataset
image <- terra::rast(global_model$metadata[1, "SourceFile"], noflip = TRUE)

# Computing surface temperature
surface_temp <- thermal_predict(image, model = global_model$models$global)

# Visualizing the result
terra::plot(surface_temp)
```

# Documentation

The `thermal` package provides two vignettes. One of them shows how to correct
thermal drift, while the other shows how to conduct radiometric calibration. We
highly encourage users to build the vignettes when installing the package and
to follow along the analyses shown there to determine the best approach for
their analytical needs. Otherwise, every exported function is documented and
examples are provided for them directly in the function documentation. Let us
know if you feel like the documentation for any of the functions is unclear.

# Limitations

`thermal` has only been tested on images from a DJI Zenmuse XT2 camera and may
therefore fail to work using images acquired with other camera models. If you
need images from a specific camera to be supported, please post an issue on the
GitHub page along with an example image and we will do our best to add support
for it.

# Bugs

If you find any bug or would like to request support for a particular feature,
please let us know by filing an issue on the package's GitHub page.

# Citation

If you use the `thermal` package as part of your analyses, please cite it as
follows:
 
Lemay M, Olugbadieye G, Maylal B, Rosa É, Gennaretti F (2026).
  _thermal: Correct Thermal Infrared Images Acquired by Drones_. R
  package version 0.0.0.9000.
