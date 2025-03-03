## ---------------------------
##
## Script name: ENERSCAPE Test
##
## Author: Rhemi Toth
##
## Date Created: 2025-02-26
##
## Email: rhemitoth@g.harvard.edu
##
## ---------------------------
##
## Notes: Testing out the enerscape package
##  
##
## ---------------------------
 
rm(list = ls())

# Packages ----------------------------------------------------------------

library(terra)
library(enerscape)


# Prepare the data --------------------------------------------------------

# Load the DEM and study area boundary
dem_raw <- rast("/Users/rhemitoth/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/Rhemi/Cembra/GIS/Images/TINITALY/w51065_s10.tif")
cembra <- vect("/Users/rhemitoth/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/Rhemi/Cembra/GIS/SHP/Cembra_Boundary/Cembra_Boundary.shp")

# Crop the DEM to the study area boundary
dem <- crop(dem_raw, cembra)

# Plot the DEM
plot(dem)
plot(cembra, add = TRUE)

# Confirm that DEM CRS is in meters
crs(dem, proj = TRUE)


# Calculate the energy landscape ------------------------------------------

# specify the body mass
kg <- 30

# calculate the energy landscape
en <- enerscape(dem, kg, neigh = 8, unit = "kcal")

# visualize the landscape
en_log <- log(en)
plot(en_log)
terra::contour(dem, add = TRUE, nlevels = 5)

