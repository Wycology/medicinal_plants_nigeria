###############################################################
# Script name: medicine_nigeria.R
#
# Purpose:
#   This script implements a complete species distribution modeling (SDM)
#   workflow for selected wild medicinal plants of Nigeria.
#   The workflow includes data acquisition, cleaning, predictor preparation,
#   model fitting, ensemble forecasting, and visualization of habitat
#   suitability under current environmental conditions and future predictions.
#
# Study context:
#   The script supports analyses presented in:
#
#   "Towards conservation of threatened wild medicinal plants in Nigeria: Distribution and niche modelling"
#    Published in: ......
#    DOI: https://doi.org/10.....
#
# Authors:
#   Abdulwakeel Ajao¹, Dorcas Oluwayemisi Olaniyan², Wyclife Agumba Oluoch³, Yusuf Ola Mukaila¹,⁴
#
# ¹Department of Botany and Plant Biotechnology, University of Johannesburg,
#                         P.O. Box 524, Auckland Park, Johannesburg, 2006, South Africa
# ²Department of Biological and Environmental Sciences, Kampala International University,
#                          Kampala P.O. Box 20000, Uganda
# ³Land Economics Group, Institute for Food and Resource
#                             Economics, University of Bonn, Bonn, Germany
# ⁴Department of Botany, Obafemi Awolowo University, Ile-Ife, 220005, Osun State, Nigeria
#
# Last modified:
#   2026-02-19
###############################################################

# Load libraries and specify the versions, including R itself
# One can install using install.packages("package_name")
# Or directly from gitHub using devtools::install_github()

# R                         # Version 4.5.1
library(sf)                 # Version 1.0.24
library(kewr)               # Version 0.6.1
library(dismo)              # Version 1.3.16
library(dplyr)              # Version 1.1.4
library(tidyr)              # Version 1.3.2
library(readr)              # Version 2.1.6
library(purrr)              # Version 1.2.1
library(scales)             # Version 1.4.0
library(ggplot2)            # Version 4.0.2
library(blockCV)            # Version 3.2.0
library(ENMeval)            # Version 2.0.5.2
library(ecospat)            # Version 4.1.3
library(geodata)            # Version 0.6.6
library(viridis)            # Version 0.6.5
library(ggspatial)          # Version 1.1.10
library(rnaturalearth)      # Version 1.2.0
library(CoordinateCleaner)  # Version 3.0.1

# Read the species names
# This has genus and species columns separated

species_list <- read.csv("data/species_list.csv")

genus <- species_list$genus # Get the genus names
species <- species_list$species # Get the species names

# Obtain the region of interest (West African countries)

africa <- ne_countries(continent = "Africa") # From rnaturalearth package

countries <- c(
  "Nigeria",
  "Ghana",
  "Benin",
  "Senegal",
  "Gambia",
  "Guinea-Bissau",
  "Guinea",
  "Sierra Leone",
  "Liberia",
  "Côte d'Ivoire",
  "Burkina Faso",
  "Togo",
  "Cameroon",
  "Central African Rep.",
  "Chad",
  "Niger",
  "Mali",
  "Mauritania",
  "Congo",
  "Gabon",
  "Eq. Guinea"
) # Selected countries

roi <- vect(africa[africa$name %in% countries, ]) # Subset from africa
roi_agg <- aggregate(roi) # Merge all the selected countries into one polygon

sp_list <- list() # Create an empty list to carry the records

# Loop through the species to obtain their records from GBIF database
# Accessed February 8th 2026
for (i in seq_along(genus)) {
  sp_list[[i]] <- sp_occurrence( # This is from geodata package
    genus = genus[i], # The genus name
    species = species[i], # The species name
    ext = roi_agg, # Region of interest
    geo = TRUE, # Pick records with coordinates
    ntries = 5                   # Try the server 5 times in case of congestion
  )
}

# Count the records returned for each species

for (i in seq_along(sp_list)) {
  print(paste0(nrow(sp_list[[i]]), " records of ", 
               unique(sp_list[[i]]$species)))
}

saveRDS(object = sp_list, file = "data/species_occ.rds") # Save obtained records to disk

sp_list <- readRDS("data/species_occ.rds") # Load the saved records (list)

# Data cleaning
coord_certainty <- lapply(sp_list, function(df) {
  
  keep_uncertainty <- TRUE
  
  if ("coordinateUncertaintyInMeters" %in% names(df)) {
    keep_uncertainty <- is.na(df$coordinateUncertaintyInMeters) |
      df$coordinateUncertaintyInMeters < 2000
  }
  
  df <- df[
    !is.na(df$lon) &
      !is.na(df$lat) &
      df$lon != 0 &
      df$lat != 0 &
      df$year >= 1950 &
      df$basisOfRecord %in% c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN") &
      df$taxonomicStatus %in% c("ACCEPTED", "SYNONYM") &
      keep_uncertainty,
  ]
  
  df <- df[complete.cases(df[, c("lon", "lat")]), ]
  df[!duplicated(df[, c("lon", "lat")]), ]
})

# Count kept records per species

for (df in coord_certainty) {
  sp <- df$species[!is.na(df$species)][1]
  cat(nrow(df), "records of", sp, "\n")
}

# Further cleaning using the CoordinateCleaner package

sp_clean_list <- lapply(coord_certainty, function(df) {
  clean_coordinates(
    x = df,
    lon = "lon",
    lat = "lat",
    tests = c(
      "capitals",      # Drop records falling in capitals of countries
      "centroids",     # Drop records falling in centroids of countries
      "equal",         # Drop records where longitude and latitude have same value
      "gbif",          # Drop records falling within 1 degree of GBIF headquarters in Copenhagen, Denmark
      "institutions",  # Drop records within 100 m of institutions like national herbaria
      "seas",          # Drop records falling in seas and oceans
      "outliers",      # Drops records > 1000 km from the rest
      "urban",         # Drops records falling within urban centers
      "zeros"          # Drops points with both longitude and latitude values as zero or equal
    ),
    value = "clean"    # Keep only clean records, drop the flagged records
  )
})

# Count the returned records per species

for (df in sp_clean_list) {
  sp <- df$species[!is.na(df$species)][1]
  cat(nrow(df), "records of", sp, "\n")
}

# Now keep only species with at least 30 records for model building

final_dataset <- sp_clean_list[sapply(sp_clean_list, nrow) >= 30]

# Keep only the species names, longitude, and latitude values

export_data <- lapply(final_dataset, function(df) {
  df[, c("species", "lon", "lat")]
})

# Save the clean dataset, ready for model building

saveRDS(final_dataset, "data/clean_dataset.rds")

# Load it for modeling

final_dataset <- readRDS("data/clean_dataset.rds")

# Arrange the columns from lon, lat, to species

clean_dataset <- lapply(final_dataset, function(df) {
  df[, c("lon", "lat", "species")]
})

# Bind the list of dataframes into a single big dataframe

clean_dataset <- do.call(rbind, clean_dataset)

# Keep only where we have complete cases, further check.

dataset <- clean_dataset[complete.cases(clean_dataset$lon), ]
dataset <- clean_dataset[complete.cases(clean_dataset$lat), ]

# Rename lon and lat as x and y and ensure they are numeric

dataset$x <- as.numeric(dataset$lon)
dataset$y <- as.numeric(dataset$lat)
clean_dataset <- dataset

# Get native ranges of each species from POWO

# Native range retrieval and species-specific spatial filtering
#   Retrieve the native geographic ranges of each
#   modeled species from the Plants of the World Online (POWO)
#   database and use these ranges to spatially constrain
#   occurrence records.
#
#   Occurrence points are filtered to retain only those falling
#   within the documented native range of each species, based
#   on TDWG level-3 regions.
#
# Data sources:
#   - POWO (Plants of the World Online) via the kewr package
#   - World Geographical Scheme for Recording Plant Distributions
#     (WGSRPD / TDWG level 3 polygons)
#
# Key assumptions:
#   - The first matched POWO record represents the correct taxon
#   - Native ranges reported by POWO are authoritative
#   - TDWG level-3 regions adequately capture species’ native
#     distributions for SDM calibration
###############################################################

# Extract unique species names from the cleaned occurrence dataset
species_set <- levels(as.factor(dataset$species))

# Initialize an empty dataframe to store native TDWG codes per species
ranges <- data.frame()

# Loop through each species to retrieve native range information from POWO
for (i in 1:length(species_set)) {
  # Search POWO using the full species name
  s <- search_powo(query = species_set[i])
  
  # Proceed only if a match is found
  if (s[[1]] != 0) {
    # Extract the IPNI taxon identifier
    ipni <- s$results[[1]]$fqId
    ipni <- strsplit(ipni, split = "names:")[[1]][2]
    
    # Look up detailed species information, including distribution
    s <- lookup_powo(taxonid = ipni, distribution = TRUE)
    
    # Convert POWO output into a tidy structure
    tidied <- tidy(s)
  }
  
  # Store species name and its native TDWG level-3 codes
  ranges <- rbind(
    ranges,
    data.frame(
      species = species_set[i],
      tdwg_code = tidied$distribution[[1]]$natives[[1]]$tdwgCode,
      stringsAsFactors = FALSE
    )
  )
}

# Spatial filtering of occurrences using native range polygons
#   For each species, occurrence records are spatially filtered
#   to retain only those points falling within its native TDWG
#   level-3 distribution.

# Load TDWG level-3 polygons (WGSRPD)
# We got the polygons from https://github.com/tdwg/wgsrpd/tree/master, accessed on 01.02.2026
tdwg <- vect("data/wgsrpd-master/level3/level3.shp")

# Initialize an empty dataframe to store filtered occurrence records
result <- data.frame()

# Loop through each species in the cleaned dataset
for (i in 1:length(species_set)) {
  sp_name <- species_set[i]
  
  # Subset POWO-derived native range codes for the current species
  distributionspowo <- ranges[ranges$species == sp_name, ]
  
  # Skip species if no TDWG codes are found
  if (nrow(distributionspowo) == 0)
    next
  
  # Select TDWG polygons corresponding to native range codes
  powo_area <- tdwg[tdwg$LEVEL3_COD %in% distributionspowo$tdwg_code, ]
  
  # Skip if no polygons found
  if (nrow(powo_area) == 0)
    next
  
  # Extract occurrence records for the current species
  sp <- clean_dataset[clean_dataset$species == sp_name, ]
  
  # Convert occurrence data to spatial points
  points <- vect(sp, geom = c("lon", "lat"))
  
  # Spatially filter points to those within the native range
  filtered <- mask(points, powo_area)
  
  # Append valid points to the final dataset
  if (nrow(filtered) > 0) {
    result <- rbind.data.frame(result, as.data.frame(filtered))
  }
}


# Download bioclim variables

# We got the envidatS3paths.txt by selecting all the variables we need from CHELSA
# https://envicloud.wsl.ch/#/?bucket=https%3A%2F%2Fos.unil.cloud.switch.ch%2Fchelsa02%2F&prefix=chelsa%2Fglobal%2Fbioclim%2F
# Then brought the wget file to R to loop and download all 874 raster files:

# 19 + (19 * 3 * 3 * 5), 19 current bioclim variables + 19 future bioclim variables * 3 time periods *
# 3 SSPS and * 5 GCMs, hence the total of 874 files

19 + (19 * 3 * 3 * 5) # 874 files

urls <- readLines("data/envidatS3paths.txt")
urls <- trimws(urls) # The urls have some white space that need to be removed

# Create directory to hold the downloaded files
output_dir <- "data/chelsa_preds/"
dir.create(path = output_dir,
           showWarnings = FALSE,
           recursive = TRUE)

options(timeout = 300) # R to wait for up to 300 seconds (5 minutes) before giving up on a file, adjust based on your internet speed

# Download the bioclimatic variables and save to disk
# Note: The total data size will be 265 Gbs

for (url in urls) {
  file_name <- basename(url)
  dest_file <- file.path(output_dir, file_name)
  download.file(url, destfile = dest_file, mode = "wb")
}

# Remove records falling within the same grid cell

base <- rast("data/chelsa_preds/CHELSA_bio01_1981-2010_V.2.1.tif") # Read one of the rasters
names(base) <- "bio1" # Rename it to bio1

clean_dataset <- result[, c("species", "x", "y")]
records_points <- vect(clean_dataset, geom = c("x", "y"))
records_points$species <- as.factor(records_points$species)
filtered_no_cell_duplicates <- data.frame()

# Loop through and drop duplicate records per grid cell

for (i in 1:length(levels(records_points$species))) {
  print(i)
  v <- records_points[records_points$species == levels(records_points$species)[i], ]
  vals <- terra::extract(
    x = base$bio1,
    y = v,
    ID = FALSE,
    cells = T
  )
  ndup_recs <- cbind.data.frame(as.data.frame(x = v, geom = "XY"), vals$cell)
  ndup_recs <- ndup_recs[!duplicated(ndup_recs$`vals$cell`), ]
  filtered_no_cell_duplicates <- rbind.data.frame(filtered_no_cell_duplicates, ndup_recs)
}

write.csv(filtered_no_cell_duplicates,
          "data/filtered_no_cell_duplicates.csv",
          row.names = FALSE)

#############################
# Environmental layers setup #
#############################

# CHELSA current
# From the folder containing all bioclim variables, pick only those for present
present <- list.files(path = "data/chelsa_preds/",
                      full.names = TRUE,
                      pattern = "_1981-2010_")
present <- rast(present) # Read them using rast function
names(present) <- sub(".*(bio[0-9]{2}).*", "\\1", names(present)) # Change their names to simpler form
names(present) # Confirm the simple names

# Crop them to extent of Africa (Native for all the species)

africa <- aggregate(vect(ne_countries(continent = "Africa"))) #
present <- crop(present, africa)

# Soil (Obtained for Africa using geodata package)

bdod      <- soil_af(var = "BLKD", depth = 30, path = tempdir())
cec       <- soil_af(var = "CEC", depth = 30, path = tempdir())
coarse    <- soil_af(var = "coarse", depth = 30, path = tempdir())
clay      <- soil_af(var = "clay", depth = 30, path = tempdir())
nitrogen  <- soil_af(var = "Ntot", depth = 20, path = tempdir())
soc       <- soil_af(var = "SOC", depth = 30, path = tempdir())
ph        <- soil_af(var = "pH", depth = 30, path = tempdir())
sand      <- soil_af(var = "sand", depth = 30, path = tempdir())
silt      <- soil_af(var = "silt", depth = 30, path = tempdir())
bdr_path  <- file.path(tempdir(), "af_bdr_30s.tif") # For bedrock depth, a bit different

download.file(url = "https://geodata.ucdavis.edu/geodata/soil/afsis/af_bdr_30s.tif",
              destfile = bdr_path, mode = "wb")
bdr <- rast(bdr_path)

soils <- c(bdod, cec, coarse, clay, nitrogen, soc, ph, sand, silt, bdr) # Combine soil rasters together

writeRaster(soils, "data/soils.tif") # Write the soil rasters to disk

soils_af <- rast("data/soils.tif") # Load the soil rasters

soils_af <- crop(soils_af, africa) # Crop them to the extent of Africa like bioclim

names(soils_af) <- c("blkd",
                     "cec",
                     "coarse",
                     "clay",
                     "ntot",
                     "soc",
                     "pH",
                     "sand",
                     "silt",
                     "bdr") # Rename them, simpler names

elev <- elevation_global(res = 0.5, path = tempdir()) # Get global elevation data

writeRaster(elev, filename = "data/elevation.tif") # Write the elevation to disk

elev <- rast("data/elevation.tif") # Load the elevation layer

elev_af <- crop(elev, africa) # Crop elevation to the extent of Africa like bioclim and soils
names(elev_af) <- "elevation" # Rename elevation, simpler name.

slope   <- terrain(elev_af, "slope")      # Create slope layer
aspect  <- terrain(elev_af, "aspect")     # Create aspect layer
hnd     <- rast("data/hnd.tif")           # Load height above nearest drainage layer (Got from Google Earth Engine)

# Below is the JavaScript code we used in Google Earth Engine to obtain 
# height above the nearest drainage data

# var image = ee.Image("MERIT/Hydro/v1_0_1");
# var hnd = image.select("hnd");
# var poly = ee.Geometry.Polygon([
#   [
#     [-22.62838731962641, -36.12684573434243],
#     [55.06692518037359,  -36.12684573434243],
#     [55.06692518037359,   39.54847013492884],
#     [-22.62838731962641,  39.54847013492884],
#     [-22.62838731962641, -36.12684573434243]
#   ]
# ]);
# Export.image.toDrive({
#   image: hnd,
#   description: "hnd",
#   folder: "hnd",
#   fileNamePrefix: "hnd",
#   region: roi,
#   scale: 1000,
# });
# Then downloaded directly from GEE into the current working directory

hnd     <- crop(hnd, africa)              # crop to the extent of Africa
hnd     <- resample(hnd, slope)           # Resample to attain properties of slope, can use any other existing layer
topo    <- c(elev_af, slope, aspect, hnd) # Put all topographic variables together

writeRaster(x = topo,
            filename = "data/topography.tif",
            overwrite = T) # Save topographic variables to disk

topo <- rast("data/topography.tif") # Load the topographic variables

preds_all <- c(present, soils_af, topo) # Bring together all the predictor variables

# Drop bio08, bio09, bio18, and bio19 (seasonal variables) due to their poor qualities near equator
vars_selected <- subset(
  x = preds_all,
  subset = c("bio08", "bio09", "bio18", "bio19"),
  negate = TRUE
)

# Save the selected variables and specify datatype as FLT4S, otherwise bio03 will be saved as integer and fall to 0s.
writeRaster(
  vars_selected,
  "data/vars_selected.tif",
  datatype = "FLT4S",
  overwrite = TRUE
)

vars_selected <- rast("data/vars_selected.tif") # Load all the current selected variables for modeling

# Maxent tuning and calibration #
options(java.parameters = "-Xmx8g")     # Allocate JAVA more memory

records_points <- read.csv("data/filtered_no_cell_duplicates.csv") # Records with no cell duplicates

records_points <- vect(records_points, geom = c("x", "y"))

records_points$species <- as.factor(records_points$species)

# Start the model training loop

counts <- as.data.frame(table(records_points$species))
records_points <- records_points[records_points$species %in% counts[counts$Freq >=
                                                                      10, ]$Var1, ]
records_points$species <- droplevels(records_points$species)
df <- as.data.frame(records_points, geom = "xy")
nig <- geodata::gadm(country = "Nigeria", level = 0, path = tempdir)

for (s in 1:length(levels(records_points$species))) {
  species <- levels(records_points$species)[s]
  species_dir <- paste0("models/", species)
  if (!dir.exists(species_dir)) {
    dir.create(species_dir, recursive = TRUE)
  }
  # Subsets occurrences
  occurrences <- records_points[records_points$species == species, ]
  crs(occurrences) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  occurrences_df <- as.data.frame(occurrences, geom = "XY")
  
  occ_vect <- vect(occurrences_df, geom = c("x", "y"))
  tdwg_extract <- terra::extract(tdwg, occ_vect)
  all_codes <- unique(tdwg_extract$LEVEL3_COD[!is.na(tdwg_extract$LEVEL3_COD)])
  
  if (length(all_codes) > 0) {
    range <- tdwg[tdwg$LEVEL3_COD %in% all_codes, ]
  } else {
    stop(
      "ERROR: No valid TDWG codes found for species: ",
      species,
      ". Check if points overlap the tdwg shapefile."
    )
  }
  
  env <- crop(vars_selected, range)
  env <- mask(env, range)
  
  # Background points
  set.seed(248)
  background_points <- spatSample(
    env,
    10000,
    method = "random",
    xy = TRUE,
    exhaustive = TRUE,
    na.rm = TRUE
  )
  
  # Background points
  eA <- background_points[, -c(1, 2)]
  
  # Presence values
  eP <- terra::extract(x = env, y = occurrences)
  
  # Create data frames from values
  p <- rbind.data.frame(
    cbind.data.frame(occurrences_df[c("x", "y")], eP[, -1], "Pres" = 1),
    cbind.data.frame(background_points, "Pres" = 0)
  )
  xy <- p[, c(1, 2)]
  Pres.cov <- data.frame(eP[, -1], Pres = 1)
  Pres.cov <- Pres.cov[complete.cases(Pres.cov), ]
  Back.cov <- data.frame(eA, Pres = 0)
  Back.cov <- Back.cov[complete.cases(Back.cov), ]
  all.cov <- rbind(Pres.cov, Back.cov)
  all.cov <- all.cov[complete.cases(all.cov), ]
  
  # Setup the cross-validation
  eMAX <- list()
  varImp <- list()
  folds <- 5
  r <- env[[1]]
  p <- vect(p, geom = c("x", "y"))
  p <- p[!is.na(terra::extract(env[[1]], p)[, 1]), ]
  crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
  crs(p) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  #Creating spatial blocks
  scv1 <- cv_spatial(
    seed = 248,
    x = p,
    column = "Pres", # the response column (binary or multi-class)
    r = r,
    k = folds, # number of folds
    size = 250000, # arbitrary
    selection = "random", # random blocks-to-fold
    iteration = 50, # find evenly dispersed folds
    progress = TRUE, # turn off progress bar
    biomod2 = FALSE, # also create folds for biomod2
    plot = TRUE,
    raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
  )
  p$fold <- scv1$folds_ids
  p <- cbind(xy, p)
  xy <- p[, c(1, 2)]
  
  # spatial block partition groups for use in ENMevaluate
  j <- which(p$Pres == 1)
  occs.grp <- p$fold[j]
  bg.grp <- p$fold[-j]
  user_grp <- list(occs.grp, bg.grp)
  names(user_grp) <- c("occs.grp", "bg.grp")
  
  # prepare presence and bg dataframes for use in ENMevaluate
  j <- which(p$Pres == 1)
  k <- which(names(p) %in% c("Pres", "fold"))
  pres_df <- data.frame(p[j, -k])
  bg_df <- data.frame(p[-j, -k])
  
  # maxent model tuning with ENMeval
  tune.args <- list(fc = c("LQ", "H", "LQH", "LQHP"),
                    rm = c(1, 3, 5))
  maxent_tuning <- ENMevaluate(
    occs = pres_df,
    bg = bg_df,
    algorithm = "maxent.jar",
    tune.args = tune.args,
    user.grp = user_grp,
    partitions = "user",
    parallel = T,
    numCores = 12
  )
  
  # select the four models with highest testing AUC
  t_results <- maxent_tuning@results
  write.csv(t_results,
            paste0("models/", species, "/model_results.csv"),
            row.names = FALSE)
  t_results <- t_results[order(t_results$cbi.val.avg, decreasing = TRUE), ]
  t_results <- t_results[1:4, ]
  
  # among these, select the model with the lowest AUC_diff
  l <- which.min(t_results$auc.diff.avg)
  tune_args <- as.character(t_results$tune.args[l])
  bestmod <- which(maxent_tuning@results$tune.args == tune_args)
  
  maxent_mod <- maxent_tuning@models[[bestmod]]
  save(maxent_mod, file = paste0("models/", species, "/model.RData"))
  
  # get evaluation metrics of this model
  bestmod_eval <- maxent_tuning@results[bestmod, ]
  
  write.csv(
    bestmod_eval,
    paste0("models/", species, "/evaluation_best_model.csv"),
    row.names = FALSE
  )
  preds_nig <- crop(vars_selected, nig, mask = TRUE)
  pr_nig <- terra::predict(preds_nig,
                       maxent_tuning@models[[bestmod]],
                       type = 'cloglog',
                       na.rm = T)
  
  writeRaster(pr_nig,
              paste0("models/", species, "/suitability.tif"),
              overwrite = TRUE)
  
  # make predictions of the best model
  pr <- terra::predict(env,
                       maxent_tuning@models[[bestmod]],
                       type = 'cloglog',
                       na.rm = T)
  
  writeRaster(pr,
              paste0("models/", species, "/suit_all.tif"),
              overwrite = TRUE)
  
  # extract suitability for presence points
  pres <- p[which(p$Pres == 1), ]
  pres_coordinates <- data.frame(pres[, c("x", "y")])
  suit_pres <- terra::extract(pr, pres_coordinates)
  
  # extract suitability for background points
  bg <- p[which(p$Pres == 0), ]
  bg_coordinates <- data.frame(bg[, c("x", "y")])
  suit_bg <- raster::extract(pr, bg_coordinates)
  
  # suitability at presence locations
  spg <- occurrences
  #spg_coords <- data.frame(spg[, c("x", "y")])
  pres_vals <- terra::extract(pr, spg)[, 2]
  
  # suitability across the whole landscape (available environment)
  avail_vals <- terra::values(pr, na.rm = TRUE)
  
  # Boyce index (continuous)
  boyce <- ecospat.boyce(
    fit = avail_vals,
    obs = pres_vals[!is.na(pres_vals)],
    nclass = 0,
    # moving window
    window.w = "default",
    res = 100
  )
  
  # Save Boyce results
  write.csv(
    data.frame(HS = boyce$HS, F_ratio = boyce$F.ratio),
    paste0("models/", species, "/boyce_curve.csv"),
    row.names = FALSE
  )
  
  # overall Boyce correlation
  boyce_cor <- boyce$cor
  
  write.csv(
    data.frame(cBoyce = boyce_cor),
    paste0("models/", species, "/boyce_summary.csv"),
    row.names = FALSE
  )
  
  # ---- SINK → SOURCE THRESHOLD ----
  
  df <- data.frame(hs = boyce$HS, pe = boyce$F.ratio)
  
  df <- df[!is.na(df$pe), ]
  df <- df[df$pe >= 1, ]
  df <- df[order(df$hs), ]
  
  th <- df$hs[1]   # minimum HS where F.ratio >= 1
  
  # binarize suitability
  pr_thr <- pr >= th
  
  writeRaster(pr_thr,
              paste0("models/", species, "/presence.tif"),
              overwrite = TRUE)
}

# --------------------------------------------------------------
# Post-modeling cBoyce analysis for all species
# --------------------------------------------------------------

# Get list of the processed species

species_list <- list.dirs("models/", full.names = FALSE, recursive = FALSE)
species_list <- species_list[species_list != ""]

for (species in species_list) {
  cat("Processing:", species, "...")
  
  # Load suitability raster
  suit <- rast(paste0("models/", species, "/suit_all.tif"))
  
  # Get occurrence points
  occ_points <- records_points[records_points$species == species, ]
  occ_coords <- as.data.frame(occ_points, geom = "XY")[, c("x", "y")]
  
  # Extract values
  pres_vals <- terra::extract(suit, occ_coords)$lyr1
  pres_vals <- pres_vals[!is.na(pres_vals)]
  
  if (length(pres_vals) < 10)
    next
  
  avail_vals <- values(suit, na.rm = TRUE)
  
  # Calculate Boyce
  boyce <- ecospat.boyce(
    fit = avail_vals,
    obs = pres_vals,
    nclass = 0,
    window.w = "default",
    res = 100,
    PEplot = FALSE
  )
  
  # Calculate threshold (minimum suitability where F.ratio >= 1)
  if (length(boyce$HS) > 0 && length(boyce$F.ratio) > 0) {
    th <- min(boyce$HS[boyce$F.ratio >= 1], na.rm = TRUE)
    
    # Save binary map
    writeRaster(suit >= th,
                paste0("models/", species, "/presence_boyce.tif"),
                overwrite = TRUE)
    
    # Save simple results
    write.csv(
      data.frame(
        species = species,
        cBoyce = boyce$cor,
        threshold = th,
        n_presence = length(pres_vals)
      ),
      paste0("models/", species, "/boyce_summary.csv"),
      row.names = FALSE
    )
    
    cat(" cBoyce =",
        round(boyce$cor, 3),
        "threshold =",
        round(th, 4),
        "\n")
  } else {
    cat(" failed\n")
  }
}

# Get species list
species_list <- list.dirs("models/", full.names = FALSE, recursive = FALSE)
species_list <- species_list[species_list != ""]

for (species in species_list) {
  cat("Processing:", species, "...")
  
  # Load suitability raster
  suit <- rast(paste0("models/", species, "/suit_all.tif"))
  suit_nig <- rast(paste0("models/", species, "/suitability.tif"))
  # Get occurrence points
  occ_points <- records_points[records_points$species == species, ]
  occ_coords <- as.data.frame(occ_points, geom = "XY")[, c("x", "y")]
  
  # Extract values
  pres_vals <- terra::extract(suit, occ_coords)$lyr1
  pres_vals <- pres_vals[!is.na(pres_vals)]
  
  if (length(pres_vals) < 10) {
    cat(" skipping (too few points)\n")
    next
  }
  
  avail_vals <- values(suit, na.rm = TRUE)
  
  # Calculate Boyce
  boyce <- ecospat.boyce(
    fit = avail_vals,
    obs = pres_vals,
    nclass = 0,
    window.w = "default",
    res = 100,
    PEplot = FALSE
  )
  
  # Calculate threshold (minimum suitability where F.ratio >= 1)
  if (length(boyce$HS) > 0 && length(boyce$F.ratio) > 0) {
    th <- min(boyce$HS[boyce$F.ratio >= 1], na.rm = TRUE)
    
    # Save binary map
    writeRaster(suit >= th,
                paste0("models/", species, "/presence_boyce.tif"),
                overwrite = TRUE)
    
    writeRaster(suit_nig >= th,
                paste0("models/", species, "/presence_boyce_nig.tif"),
                overwrite = TRUE)
    # Save detailed Boyce data for plotting
    boyce_df <- data.frame(HS = boyce$HS, F.ratio = boyce$F.ratio)
    boyce_df <- boyce_df[complete.cases(boyce_df), ]
    
    write.csv(boyce_df,
              paste0("models/", species, "/boyce_detailed.csv"),
              row.names = FALSE)
    
    # Create and save the Boyce plot
    if (nrow(boyce_df) > 0) {
      p <- ggplot(boyce_df, aes(x = HS, y = F.ratio)) +
        geom_line(linewidth = 1.2, color = "darkblue") +
        geom_point(size = 2,
                   color = "#2E4053",
                   alpha = 0.7) +
        geom_hline(
          yintercept = 1,
          linetype = "dashed",
          color = "orange",
          linewidth = 0.8
        ) +
        geom_vline(
          xintercept = th,
          linetype = "dashed",
          color = "firebrick",
          linewidth = 1.2
        ) +
        geom_ribbon(aes(ymin = 1, ymax = F.ratio),
                    fill = "#5dade2",
                    alpha = 0.15) +
        annotate(
          "text",
          x = th / 2,
          y = max(boyce_df$F.ratio, na.rm = TRUE) * 0.6,
          label = "Sink",
          size = 3,
          fontface = "bold",
          color = "darkred"
        ) +
        annotate(
          "text",
          x = (th + 1) / 2,
          y = max(boyce_df$F.ratio, na.rm = TRUE) * 0.6,
          label = "Source",
          size = 3,
          fontface = "bold",
          color = "darkgreen"
        ) +
        labs(
          title = paste("Boyce Index:", species),
          x = "Habitat suitability",
          y = "Predicted/expected ratio",
          caption = paste(
            "cBoyce =",
            round(boyce$cor, 3),
            " | Threshold =",
            round(th, 4)
          )
        ) +
        scale_x_continuous(
          limits = c(0, 1),
          breaks = seq(0, 1, 0.2),
          expand = c(0, 0)
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_bw(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14),
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_blank(),
          legend.position = "none"
        )
      
      # Save plot
      ggsave(
        filename = paste0("models/", species, "/boyce_plot.png"),
        plot = p,
        width = 8,
        height = 6,
        dpi = 300
      )
    }
    
    # Save simple summary
    write.csv(
      data.frame(
        species = species,
        cBoyce = boyce$cor,
        threshold = th,
        n_presence = length(pres_vals)
      ),
      paste0("models/", species, "/boyce_summary.csv"),
      row.names = FALSE
    )
    
    cat(" cBoyce =",
        round(boyce$cor, 3),
        "threshold =",
        round(th, 4),
        "plot saved\n")
  } else {
    cat(" failed\n")
  }
}

# --------------------------------------------------------------
# Create summary plots
# --------------------------------------------------------------

# List species directories
species_dirs <- list.dirs("models", recursive = FALSE)

# Read and combine all csvs
boyce_df <- do.call(rbind, lapply(species_dirs, function(d) {
  read.csv(file.path(d, "boyce_summary.csv"))
}))

# Check and plot the result
boyce_df

ggplot(boyce_df, aes(
  x = reorder(species, cBoyce),
  y = cBoyce,
  fill = cBoyce
)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "orange", high = "darkblue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      face = "italic",
      hjust = 1,
      color = "black"
    ),
    text = element_text(size = 24, color = "black"),
    legend.position = "none"
  ) +
  labs(x = "Species", y = "cBoyce")

ggsave("figures/Figure_1_boyce_distribution2.png",
       width = 8,
       height = 6)


# Create and visualize current richness ---------------

species_dirs <- list.dirs("models", recursive = FALSE)

ras_list <- lapply(species_dirs, function(d) {
  
  f <- file.path(d, "presence_boyce_nig.tif")
  if (!file.exists(f)) return(NULL)
  
  r <- rast(f)
})

boyce_rich <- rast(ras_list)

richness <- sum(boyce_rich, na.rm = TRUE)

current_richness <- "current_richness"
dir.create(current_richness, showWarnings = FALSE)

writeRaster(richness,
            "current_richness/current_richness.tif",
            overwrite = TRUE)

# Convert richness raster to data.frame
rich_df <- as.data.frame(richness, xy = TRUE, na.rm = TRUE)
colnames(rich_df) <- c("lon", "lat", "richness")

# nigeria as sf (if not already)
nigeria_sf <- st_as_sf(nig)

rich_cols <- c(
  "gray",
  "#08306B", # dark blue
  "#2171B5", # blue
  "#6BAED6", # light blue
  "#A1D99B", # light green
  "#FEE08B", # yellow
  "#F46D43", # orange
  "#7F3B08"  # dark brown
)

# Plot
ggplot() +
  geom_raster(data = rich_df, aes(x = lon, y = lat, fill = richness)) +
  geom_sf(
    data = nigeria_sf,
    fill = NA,
    color = "white",
    linewidth = 0.4
  ) +
  scale_fill_gradientn(colours = rich_cols, name = "Species richness (per cell)") +
  coord_sf(expand = FALSE) +
  scale_x_continuous(
    labels = function(x)
      paste0(x, "°E")
  ) +
  scale_y_continuous(
    labels = function(y)
      paste0(y, "°N")
  ) +
  annotation_scale(
    location = "br",
    # bottom-right
    width_hint = 0.2,
    # adjust width
    bar_cols = c("black", "white"),
    text_cex = 1.2
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = 15,
    barheight = 1
  )) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.text = element_text(size = 12, color = "black")
  )

ggsave("figures/Figure_2_current_richness2.png",
       width = 8,
       height = 6)

##################
# Model analysis #
##################
# Load all evaluations
eval_files <- list.files(
  "models/",
  pattern = "evaluation_best_model.csv",
  full.names = T,
  recursive = T
)

result <- data.frame()
for (i in eval_files) {
  eval_table <- read.csv(i)
  species <- basename(dirname(i))
  result <- rbind(result, cbind(species = species, eval_table))
}

# Load the binaries, crop, and merge
bin_models <- list.files(
  "models/",
  pattern = "presence_boyce_nig.tif",
  full.names = T,
  recursive = T
)
results <- rast(bin_models)

# Sum everything
rich <- sum(results, na.rm = T)

######################################################################
# Future Layers and Predictions
######################################################################

# Load models and buffers
models <- list.files(
  "models/",
  pattern = "model.RData",
  full.names = T,
  recursive = T
)
evaluations <- list.files(
  "models/",
  pattern = "evaluation",
  full.names = T,
  recursive = T
)
thresholds <- list.files(
  "models/",
  pattern = "boyce_summary",
  full.names = T,
  recursive = T
)
presence <- list.files(
  "models/",
  pattern = "presence_boyce_nig.tif",
  full.names = T,
  recursive = T
)

clim_vars <- c(
  "bio01",
  "bio02",
  "bio03",
  "bio04",
  "bio05",
  "bio06",
  "bio07",
  "bio10",
  "bio11",
  "bio12",
  "bio13",
  "bio14",
  "bio15",
  "bio16",
  "bio17"
)

static_vars <- setdiff(names(vars_selected), clim_vars)

static_stack <- vars_selected[[static_vars]]
static_stack <- crop(static_stack, nig, mask = TRUE)

clean_bio_names <- function(x) {
  sub(".*_(bio[0-9]{2})_.*", "\\1", x)
}

gcms  <- c("gfdl-esm4",
           "ipsl-cm6a-lr",
           "mri-esm2-0",
           "ukesm1-0-ll",
           "mpi-esm1-2-hr")
ssps  <- c("ssp126", "ssp370", "ssp585")
times <- c("2011-2040", "2041-2070", "2071-2100")

fut <- list()

for (gcm in gcms) {
  for (ssp in ssps) {
    for (time in times) {
      pat <- paste0(gcm, "_", ssp, "_bio(0[1-7]|1[0-7])_", time)
      
      files <- list.files("data/chelsa_preds",
                          pattern = pat,
                          full.names = TRUE)
      
      if (length(files) == 0)
        next
      
      # read climate
      fut_clim <- rast(files)
      
      # crop to Nigeria
      fut_clim <- crop(fut_clim, nig, mask = TRUE)
      
      # clean layer names → bio01, bio02, ...
      names(fut_clim) <- clean_bio_names(names(fut_clim))
      
      # append static variables
      fut_full <- c(fut_clim, static_stack)
      
      # enforce identical layer order to vars_selected
      fut_full <- fut_full[[names(vars_selected)]]
      
      # store
      key <- paste(ssp, time, gcm, sep = "_")
      fut[[key]] <- fut_full
    }
  }
}

get_species <- function(path) {
  basename(dirname(path))
}

dir.create("projections", showWarnings = FALSE)

for (i in seq_along(models)) {
  # ---- load species-specific inputs ----
  species <- get_species(models[i])
  
  load(models[i])              # loads maxent_mod
  model <- maxent_mod
  
  thr <- read.csv(thresholds[i])   # Boyce summary
  th  <- thr$threshold             # <-- IMPORTANT: this must exist
  # (the sink→source threshold you computed)
  
  message("Projecting: ", species)
  
  # ---- project to all futures ----
  for (nm in names(fut)) {
    fut_stack <- fut[[nm]]
    
    # continuous prediction
    pr <- terra::predict(fut_stack, model, type = "cloglog", na.rm = TRUE)
    
    # Boyce-based binary (sink vs source)
    pr_bin <- pr >= th
    
    # ---- parse scenario info from name ----
    parts <- strsplit(nm, "_")[[1]]
    ssp   <- parts[1]
    time  <- parts[2]
    gcm   <- parts[3]
    
    outdir <- file.path("projections", ssp, time)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    outfile <- file.path(outdir, paste0(species, "_", gcm, "_boyce_binary.tif"))
    
    writeRaster(pr_bin, outfile, overwrite = TRUE)
  }
}

consensus_dir <- "projections_consensus"
dir.create(consensus_dir, showWarnings = FALSE)

ssps  <- c("ssp126", "ssp370", "ssp585")
times <- c("2011-2040", "2041-2070", "2071-2100")

for (ssp in ssps) {
  for (time in times) {
    indir <- file.path("projections", ssp, time)
    if (!dir.exists(indir))
      next
    
    files <- list.files(indir, pattern = "_boyce_binary.tif$", full.names = TRUE)
    
    if (length(files) == 0)
      next
    
    # extract species names
    spp <- unique(sub("_.*", "", basename(files)))
    
    for (sp in spp) {
      sp_files <- files[grepl(paste0("^", sp, "_"), basename(files))]
      
      # safety: need at least 3 GCMs
      if (length(sp_files) < 3)
        next
      
      r_sum <- sum(rast(sp_files), na.rm = TRUE)
      
      # majority rule (>=3 of 5)
      r_cons <- r_sum >= 3
      
      outdir <- file.path(consensus_dir, ssp, time)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      
      writeRaster(r_cons, file.path(outdir, paste0(sp, "_boyce_consensus.tif")), overwrite = TRUE)
    }
  }
}

richness_dir <- "richness"

dir.create(richness_dir, showWarnings = FALSE)

for (ssp in ssps) {
  for (time in times) {
    files <- list.files(file.path(consensus_dir, ssp, time),
                        pattern = "_boyce_consensus.tif$",
                        full.names = TRUE)
    
    if (length(files) == 0)
      next
    
    richness <- sum(rast(files), na.rm = TRUE)
    
    writeRaster(richness, file.path(richness_dir, paste0("richness_", ssp, "_", time, ".tif")), overwrite = TRUE)
  }
}


# Create richness change maps -------------------------

richness_present <- rast("current_richness/current_richness.tif")

change_dir <- "richness_change"
dir.create(change_dir, showWarnings = FALSE)

for (ssp in ssps) {
  for (time in times) {
    fut_file <- file.path(richness_dir, paste0("richness_", ssp, "_", time, ".tif"))
    
    if (!file.exists(fut_file))
      next
    
    richness_fut <- rast(fut_file)
    
    change <- richness_fut - richness_present
    
    writeRaster(change, file.path(
      change_dir,
      paste0("richness_change_", ssp, "_", time, ".tif")
    ), overwrite = TRUE)
  }
}

# Nigeria boundary (sf)
nga <- st_as_sf(nig)

# Load richness rasters
rich_files <- list.files("richness", pattern = "richness_.*tif$", full.names = TRUE)

# Load change rasters
change_files <- list.files("richness_change", pattern = "richness_change_.*tif$", full.names = TRUE)

rast_to_df <- function(r, value_name = "value") {
  as.data.frame(r, xy = TRUE, na.rm = TRUE) |>
    rename(!!value_name := 3)
}

rich_df <- lapply(rich_files, function(f) {
  r <- rast(f)
  
  meta <- strsplit(basename(f), "_")[[1]]
  ssp  <- meta[2]
  time <- gsub(".tif", "", meta[3])
  
  rast_to_df(r, "richness") |>
    mutate(ssp = ssp, time = time)
}) |>
  bind_rows()

rich_df$ssp  <- factor(rich_df$ssp, levels = c("ssp126", "ssp370", "ssp585"))
rich_df$time <- factor(rich_df$time, levels = c("2011-2040", "2041-2070", "2071-2100"))

p_richness <- ggplot() +
  geom_raster(data = rich_df, aes(x = x, y = y, fill = richness)) +
  geom_sf(
    data = nigeria_sf,
    fill = NA,
    color = "white",
    linewidth = 0.4
  ) +
  facet_grid(time ~ ssp) +
  scale_fill_gradientn(colours = rich_cols, name = "Species richness (per cell)") +
  coord_sf(crs = st_crs(4326), expand = FALSE) +
  scale_x_continuous(
    labels = function(x)
      paste0(x, "°E")
  ) +
  scale_y_continuous(
    labels = function(y)
      paste0(y, "°N")
  ) +
  annotation_scale(
    location = "br",
    width_hint = 0.2,
    bar_cols = c("black", "white"),
    text_cex = 1,
    unit_category = "metric",
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm")
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = 15,
    barheight = 1
  )) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

p_richness

ggsave("figures/Figure_3_future_richness.png",
       width = 10,
       height = 10)

change_df <- lapply(change_files, function(f) {
  r <- rast(f)
  
  meta <- strsplit(basename(f), "_")[[1]]
  ssp  <- meta[3]
  time <- gsub(".tif", "", meta[4])
  
  rast_to_df(r, "change") |>
    mutate(ssp = ssp, time = time)
}) |>
  bind_rows()

change_df$ssp  <- factor(change_df$ssp, levels = c("ssp126", "ssp370", "ssp585"))
change_df$time <- factor(change_df$time,
                         levels = c("2011-2040", "2041-2070", "2071-2100"))

change_cols <- c(
  "red", # strong loss (red)
  "#EF8A62",
  "#FDDBC7",
  "#F7F7F7", # zero
  "#D1E5F0",
  "#67A9CF",
  "blue" # strong gain (blue)
)

max_abs_change <- max(abs(change_df$change), na.rm = TRUE)

p_change <- ggplot() +
  geom_raster(data = change_df, aes(x = x, y = y, fill = change)) +
  geom_sf(
    data = nigeria_sf,
    fill = NA,
    color = "gray90",
    linewidth = 0.4
  ) +
  facet_grid(time ~ ssp) +
  scale_fill_gradientn(colours = change_cols,
                       limits = c(-6, 6),
                       name = "Change in species richness\n(future − present)") +
  coord_sf(crs = st_crs(4326), expand = FALSE) +
  scale_x_continuous(
    labels = function(x)
      paste0(x, "°E")
  ) +
  scale_y_continuous(
    labels = function(y)
      paste0(y, "°N")
  ) +
  annotation_scale(
    location = "br",
    width_hint = 0.2,
    bar_cols = c("black", "white"),
    text_cex = 1,
    unit_category = "metric",
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm")
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = 15,
    barheight = 1
  )) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

p_change

ggsave("figures/Figure_4_richness_change.png",
       width = 10,
       height = 10)

# Create richness per gcm -----------------------------

gcms  <- c("gfdl-esm4",
           "ipsl-cm6a-lr",
           "mri-esm2-0",
           "ukesm1-0-ll",
           "mpi-esm1-2-hr")
ssps  <- c("ssp126", "ssp370", "ssp585")
times <- c("2011-2040", "2041-2070", "2071-2100")

dir.create("richness_gcm",
           showWarnings = FALSE,
           recursive = TRUE)

for (ssp in ssps) {
  for (time in times) {
    for (gcm in gcms) {
      files <- list.files(
        file.path("projections", ssp, time),
        pattern = paste0("_", gcm, "_boyce_binary.tif$"),
        full.names = TRUE
      )
      
      if (length(files) == 0)
        next
      
      r_stack <- rast(files)
      r_rich  <- sum(r_stack, na.rm = TRUE)
      
      outdir <- file.path("richness_gcm", ssp, time)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      
      writeRaster(r_rich, file.path(outdir, paste0("richness_", gcm, ".tif")), overwrite = TRUE)
    }
  }
}

# Create uncertainty in richness across gcms ----------

dir.create("richness_uncertainty", showWarnings = FALSE)

uncertainty_rasters <- list()

for (ssp in ssps) {
  for (time in times) {
    files <- list.files(
      file.path("richness_gcm", ssp, time),
      pattern = "^richness_.*\\.tif$",
      full.names = TRUE
    )
    
    if (length(files) < 2)
      next  # need ≥2 GCMs
    
    r_stack <- rast(files)
    
    # inter-GCM uncertainty = SD of richness
    r_sd <- app(r_stack, sd, na.rm = TRUE)
    
    outname <- paste0("richness_uncertainty_", ssp, "_", time, ".tif")
    
    writeRaster(r_sd,
                file.path("richness_uncertainty", outname),
                overwrite = TRUE)
    
    uncertainty_rasters[[paste(ssp, time, sep = "_")]] <- r_sd
  }
}

unc_files <- list.files("richness_uncertainty",
                        pattern = "^richness_uncertainty_.*tif$",
                        full.names = TRUE)

unc_df <- lapply(unc_files, function(f) {
  r <- rast(f)
  meta <- strsplit(basename(f), "_")[[1]]
  
  rast_to_df(r, "uncertainty") |>
    mutate(ssp  = meta[3], time = gsub(".tif", "", meta[4]))
}) |>
  bind_rows()

unc_df$ssp  <- factor(unc_df$ssp, levels = c("ssp126", "ssp370", "ssp585"))
unc_df$time <- factor(unc_df$time, levels = c("2011-2040", "2041-2070", "2071-2100"))

unc_cols <- c("#F7FCF0", "#CCEBC5", "#7BCCC4", "#2B8CBE", "#084081")

p_uncertainty <- ggplot() +
  geom_raster(data = unc_df, aes(x = x, y = y, fill = uncertainty)) +
  geom_sf(
    data = nigeria_sf,
    fill = NA,
    color = "gray90",
    linewidth = 0.4
  ) +
  facet_grid(time ~ ssp) +
  scale_fill_gradientn(colours = unc_cols,
                       limits = c(0, 2),
                       name = "Richness uncertainty\n(SD across GCMs)") +
  coord_sf(crs = st_crs(4326), expand = FALSE) +
  scale_x_continuous(
    labels = function(x)
      paste0(x, "°E")
  ) +
  scale_y_continuous(
    labels = function(y)
      paste0(y, "°N")
  ) +
  annotation_scale(
    location = "br",
    width_hint = 0.2,
    bar_cols = c("black", "white"),
    text_cex = 1,
    unit_category = "metric",
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm")
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = 15,
    barheight = 1
  )) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

p_uncertainty

ggsave("figures/Figure_5_richness_uncertainty.png",
       width = 10,
       height = 10)

# Create cBoyce curves for the species ----------------

species_dirs <- list.dirs("models", recursive = FALSE, full.names = TRUE)

boyce_df <- map_dfr(species_dirs, function(sp_dir) {
  sp <- basename(sp_dir)
  
  det_file <- file.path(sp_dir, "boyce_detailed.csv")
  sum_file <- file.path(sp_dir, "boyce_summary.csv")
  
  if (!file.exists(det_file) || !file.exists(sum_file))
    return(NULL)
  
  det <- read_csv(det_file, show_col_types = FALSE)
  sum <- read_csv(sum_file, show_col_types = FALSE)
  
  det |>
    mutate(
      species   = sp,
      threshold = sum$threshold,
      cBoyce    = sum$cBoyce
    )
})

boyce_df$species <- factor(boyce_df$species, levels = sort(unique(boyce_df$species)))

p_boyce <- ggplot(boyce_df, aes(x = HS, y = F.ratio)) +
  
  # confidence / dominance area
  geom_ribbon(aes(ymin = 1, ymax = F.ratio),
              fill = "#5DADE2",
              alpha = 0.18) +
  
  # curve
  geom_line(linewidth = 0.9, color = "#1F4E79") +
  
  # reference lines
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "#E67E22",
    linewidth = 0.6
  ) +
  
  geom_vline(
    aes(xintercept = threshold),
    linetype = "dashed",
    color = "#B03A2E",
    linewidth = 0.7
  ) +
  
  facet_wrap( ~ species, ncol = 4, scales = "free_y") +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  
  labs(x = "Habitat suitability", y = "Predicted / expected ratio") +
  
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(size = 9)
  )

p_boyce

ggsave("figures/Figure_6_boyce_curve.png",
       width = 10,
       height = 10)

# Extract best model properties or metrics for --------

# Path to your models folder
models_dir <- "models"

# List all species folders
species_folders <- list.dirs(models_dir, full.names = TRUE, recursive = FALSE)

# Columns to keep from evaluation_best_model
keep_cols <- c("fc", "rm", "cbi.train", "cbi.val.avg", "cbi.val.sd")

# Initialize list to store results
all_results <- list()

for (sp_folder in species_folders) {
  species <- basename(sp_folder)
  
  # Paths to CSVs
  eval_path <- file.path(sp_folder, "evaluation_best_model.csv")
  boyce_path <- file.path(sp_folder, "boyce_summary.csv")
  
  if (!file.exists(eval_path) || !file.exists(boyce_path))
    next
  
  # Read evaluation_best_model and keep only relevant columns
  eval_df <- read_csv(eval_path, show_col_types = FALSE) %>%
    select(any_of(keep_cols))
  
  # Read Boyce summary to get n_presence
  boyce_df <- read_csv(boyce_path, show_col_types = FALSE) %>%
    select(n_presence)
  
  # Add species and n columns
  eval_df <- eval_df %>%
    mutate(
      species = species,
      n = boyce_df$n_presence[1]  # assuming one row per species
    )
  
  all_results[[species]] <- eval_df
}

# Combine all species results
combined_results <- bind_rows(all_results) %>%
  # Reorder columns
  select(species, n, fc, rm, cbi.train, cbi.val.avg, cbi.val.sd) %>%
  mutate(
    cbi.train = round(cbi.train, 2),
    cbi.val.avg = round(cbi.val.avg, 2),
    cbi.val.sd = round(cbi.val.sd, 2)
  ) %>%
  arrange(species, n, fc, rm)

# Write out combined CSV
write.csv(combined_results,
          "tables/all_species_best_model_results.csv",
          row.names = FALSE)

# Get the permutation importance for each vari --------

# Path to the models directory
models_dir <- "models"

# List all species folders
species_folders <- list.dirs(models_dir, full.names = TRUE, recursive = FALSE)

# Initialize an empty list to store results
importance_list <- list()

# Loop through all species folders
for (species_folder in species_folders) {
  # Path to the model.RData file
  model_path <- file.path(species_folder, "model.RData")
  
  if (file.exists(model_path)) {
    # Load the model
    load(model_path)
    
    # Ensure the model is a MaxEnt model (checking the class)
    if (inherits(maxent_mod, "MaxEnt_model")) {
      # Extract the results from the model
      results <- maxent_mod@results
      
      # Convert the results to a data frame
      results_df <- as.data.frame(results)
      
      # Get permutation importance columns by filtering for those containing "permutation.importance"
      importance_cols <- grep("permutation.importance", rownames(results_df), value = TRUE)
      
      # Check if any permutation importance columns exist
      if (length(importance_cols) > 0) {
        # Create a temporary data frame to store the variable and its permutation importance
        importance_df <- data.frame(
          variable = importance_cols,
          importance = as.numeric(results_df[importance_cols, "V1"]),
          species = gsub("\\..*", "", basename(species_folder))  # Clean up species name
        )
        
        # Remove "permutation.importance" from the variable names
        importance_df$variable <- gsub("\\.permutation.importance",
                                       "",
                                       importance_df$variable)
        
        # Append to the list
        importance_list[[basename(species_folder)]] <- importance_df
      }
    }
  }
}

# Combine all species data into one data frame
importance_data <- do.call(rbind, importance_list)

# Plot the permutation importance using ggplot2
ggplot(importance_data, aes(x = reorder(variable, importance), y = importance)) +
  geom_boxplot(fill = "#56B4E9",
               color = "black",
               alpha = 0.7) +
  coord_flip() +  # Flip axes for better readability
  labs(x = "Variable", y = "Permutation Importance (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(hjust = 1),
        strip.text = element_text(size = 12, face = "bold"))

# Save the plot
ggsave(
  "figures/Figure_7_current_richness2.png",
  width = 8,
  height = 9,
  dpi = 400
)

# Get the presence areas per species now and i --------

fs <- list.files("models")

for (i in seq_along(fs)) {
  r <- expanse(rast(paste0("models/", fs[i], 
                           "/presence_boyce_nig.tif")), 
               unit = "km", byValue = TRUE)
  print(paste(fs[i], round(r[r$value == 1, "area"])))
}

ssps <- c("ssp126", "ssp370", "ssp585")
times <- c("2011-2040", "2041-2070", "2071-2100")

areas_all <- list()
u <- 0
for (ssp in ssps) {
  for (time in times) {
    files <- list.files(
      file.path("projections_consensus", ssp, time),
      pattern = "consensus.tif$",
      full.names = TRUE
    )
    
    r_stack <- rast(files)
    names(r_stack) <- basename(files)
    # inter-GCM uncertainty = SD of richness
    r_areas <- expanse(r_stack, unit = "km", byValue = TRUE)
    r_areas$the_ssp <- ssp
    r_areas$the_time <- time
    r_areas <- r_areas[r_areas$value == 1, ]
    r_areas <- r_areas %>% 
      mutate(species = case_when(layer == 1 ~ "Afrofittonia silvestris", 
                                 layer == 2 ~ "Afzelia africana",
                                 layer == 3 ~ "Albizia ferruginea",
                                 layer == 4 ~ "Antrocaryon micraster",
                                 layer == 5 ~ "Baillonella toxisperma",
                                 layer == 6 ~ "Diospyros barteri",
                                 layer == 7 ~ "Entandrophragma angolense",
                                 layer == 8 ~ "Entandrophragma utile",
                                 layer == 9 ~ "Gambeya albida",
                                 layer == 10 ~ "Garcinia kola",
                                 layer == 11 ~ "Irvingia gabonensis",
                                 layer == 12 ~ "Khaya grandifoliola",
                                 layer == 13 ~ "Khaya ivorensis",
                                 layer == 14 ~ "Khaya senegalensis",
                                 layer == 15 ~ "Lophira alata",
                                 layer == 16 ~ "Nauclea diderrichii",
                                 layer == 17 ~ "Okoubaka aubrevillei",
                                 layer == 18 ~ "Terminalia ivorensis",
                                 layer == 19 ~ "Vitellaria paradoxa",
                                 TRUE ~ NA))
    u <- u + 1
    areas_all[[u]] <- r_areas
  }
}

bound_all <- bind_rows(areas_all)

bound_all <- bound_all %>% 
  mutate(present_area = case_when(
    species == "Afrofittonia silvestris" ~ 16121,
    species == "Afzelia africana" ~ 54873,
    species == "Albizia ferruginea" ~ 180045,
    species == "Antrocaryon micraster" ~ 89058,
    species == "Baillonella toxisperma" ~ 16566,
    species == "Diospyros barteri" ~ 65936,
    species == "Entandrophragma angolense" ~ 118125,
    species == "Entandrophragma utile" ~ 135521,
    species == "Gambeya albida" ~ 190285,
    species == "Garcinia kola" ~ 149510,
    species == "Irvingia gabonensis" ~ 87325,
    species == "Khaya grandifoliola" ~ 147043,
    species == "Khaya ivorensis" ~ 97575,
    species == "Khaya senegalensis" ~ 72097,
    species == "Lophira alata" ~ 102701,
    species == "Nauclea diderrichii" ~ 89309,
    species == "Okoubaka aubrevillei" ~ 144436,
    species == "Terminalia ivorensis" ~ 64296,
    species == "Vitellaria paradoxa" ~ 113382,
    TRUE ~ NA
  )) %>% 
  mutate(area = round(area)) %>% 
  mutate(area_change = area - present_area) %>% 
  mutate(change = (area_change / present_area)*100)

# Visualize the winners and losers --------------------

bound_all %>% 
  filter(the_time == "2071-2100") %>% 
  mutate(pos = ifelse(change >= 0, TRUE, FALSE)) %>% 
  ggplot(aes(reorder(species, change), change, fill = pos)) +
  geom_col() +
  scale_fill_manual(values = c("brown", "green"), guide = "none") +
  facet_grid(~the_ssp, scales = "free") +
  coord_flip() +
  theme_classic() +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.title = element_text(face = "bold")) +
  labs(y = "Presence area change relative to current (%)",
       x = "Species")

ggsave("figures/Figure_10_winers_losers_2071_2100.png",
       width = 10,
       height = 8)

#################  END  ####################################