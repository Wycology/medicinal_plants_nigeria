###############################################################
# Script: medicinal_plants_nigeria.R
#
# Purpose:  Build species distribution models for medicinal plants
#           of Nigeria, make predictions in current and future scenarios
#           and generate plots and figures
#
# Paper: "The fate of endemic medicinal plants of Nigeria in the face of climate change"
#         published in Conservation Biology, doi: https:/doi/10.....
#
# Workflow:
#   0. Create directories for the study
#   1. Obtain region of interest
#   2. Download occurrence records and clean
#   3. Obtain predictor variables
#   4. Generate background points
#   5. Build the models
#   6. Ensemble for each species
#   7. Make future predictions
#   8. Generate suitability maps
#   9. Generate richness maps
#  10. Write about the methodology
#
# Author: Abdulwakeel Ajao, Wyclife Agumba Oluoch, Yusuf Mukaila, Dorcas Oluwayemisi Olaniyan
# Last modified: 2026-01-21
###############################################################

# We recommend you run the script inside a new R Project

# ---------------------------------------------------------------------------
# 0. Create directories
# ---------------------------------------------------------------------------

create_project_dirs <- function(dirs = c("data", "figures", "models", "output", "script")) {
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      message("Created directory: ", d)
    } else {
      message("Directory already exists: ", d)
    }
  }
}

create_project_dirs() # This creates the five directories automatically

# ---------------------------------------------------------------------------
# Load libraries
# ---------------------------------------------------------------------------
#
# The following is the R version and list of libraries used in generating the code for this script.
# You can install any of them by running install.packages("package_name")
#
# R                         # Version 4.5.2:  The programming environment
library(sf)                 # Version 1.0.24: Spatial vector handling
library(sdm)                # Version 1.2.67: For building the models
library(usdm)               # Version 2.1.7:  Multicollinearity check
library(geodata)            # Version 0.6.6:  For downloading roi, occurrance, and environmental data.
library(ggplot2)            # Version 4.0.1:  Data visualization
library(rnaturalearth)      # Version 1.2.0:  Getting roi map
library(CoordinateCleaner)  # Version 3.0.1:  Further clean coordinate

# ---------------------------------------------------------------------------
# 1. Obtain the region of interest
# ---------------------------------------------------------------------------

africa <- ne_countries(continent = "Africa")

countries <- c("Nigeria", "Ghana", "Benin", "Senegal",
               "Gambia", "Guinea-Bissau", "Guinea",
               "Sierra Leone", "Liberia", "Côte d'Ivoire",
               "Burkina Faso", "Togo", "Cameroon",
               "Central African Rep.", "Chad", "Niger",
               "Mali", "Mauritania", "Congo", "Gabon",
               "Eq. Guinea")

roi <- vect(africa[africa$name %in% countries, ])
roi_agg <- aggregate(roi)

# ---------------------------------------------------------------------------
# 2. Download occurrence records and clean
# ---------------------------------------------------------------------------

species_list <- read.csv("data/species_list.csv")

genus <- species_list$genus
species <- species_list$species

sp_list <- list()

for (i in seq_along(genus)) {
  sp_list[[i]] <- geodata::sp_occurrence(
    genus = genus[i],
    species = species[i],
    ext = roi_agg,
    geo = TRUE,
    ntries = 5
  )
}

saveRDS(object = sp_list, file = "data/species_occ.rds")

for (i in seq_along(sp_list)) {
  print(paste0(nrow(sp_list[[i]]), " records of ", unique(sp_list[[i]]$species)))
}

sp_split <- readRDS("data/species_occ.rds")

# Clean coordinates

coord_certainty <- lapply(sp_list, function(df) {
  
  df <- df[
    !is.na(df$lon) &
    !is.na(df$lat) &
    df$lon != 0 &
    df$lat != 0 &
    df$year >= 1950 &
    df$basisOfRecord %in% c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN") &
    df$taxonomicStatus == "ACCEPTED" &
    (is.na(df$coordinateUncertaintyInMeters) |
        df$coordinateUncertaintyInMeters < 2000),
  ]
  
  df <- df[complete.cases(df[, c("lon", "lat")]), ]
  
  df[!duplicated(df[, c("lon", "lat")]), ]
})

for (df in coord_certainty) {
  sp <- df$species[!is.na(df$species)][1]
  cat(nrow(df), "records of", sp, "\n")
}

sp_clean_list <- lapply(coord_certainty, function(df) {
  clean_coordinates(
    x = df,
    lon = "lon",
    lat = "lat",
    tests = c(
      "capitals",
      "centroids",
      "equal",
      "gbif",
      "institutions",
      "seas",
      "outliers",
      "urban",
      "zeros"
    ),
    value = "clean"
  )
})

for (df in sp_clean_list) {
  sp <- df$species[!is.na(df$species)][1]
  cat(nrow(df), "records of", sp, "\n")
}

final_dataset <- sp_clean_list[sapply(sp_clean_list, nrow) >= 30]

export_data <- lapply(final_dataset, function(df){
  df[, c("species", "lon", "lat")]
})

saveRDS(final_dataset, "data/clean_dataset.rds")

final_dataset <- readRDS("data/clean_dataset.rds")

final_dataset_geo <- lapply(final_dataset, function(df) {
  v <- vect(df, geom = c("lon", "lat"))
  crop(v, roi_agg)
})

nigeria <- roi[roi$name == "Nigeria",]

pdf("figures/species_for_model.pdf", width = 8, height = 6)

for (i in seq_along(final_dataset_geo)) {
  species_name <- final_dataset_geo[[i]]$species
  plot(roi, lwd = 0.2, border = "gray80", main = parse(text = paste0('italic("', 
                                       species_name, '")')))  
  plot(final_dataset_geo[[i]], col = "red", add = TRUE)
  plot(nigeria, border = "blue", lwd = 2, add = TRUE)
}

dev.off()

# ---------------------------------------------------------------------------
# 3. Obtain predictor variables
# ---------------------------------------------------------------------------

preds <- worldclim_global(res = 0.5,
                           var = "bio",
                           path = "data/preds.tif")

names(preds) <- paste0("bio", 1:19)

preds <- crop(preds, roi_agg, mask = TRUE)

writeRaster(preds, "data/clim_preds.tif")

# ---------------------------------------------------------------------------
# 4. Generate background points
# ---------------------------------------------------------------------------

# Used background

v <- vifcor(preds)

preds_used <- exclude(preds, v)

bg <- background(x = preds_used, 
                 n = 10000,
                 method = 'gRandom')

bg_v <- vect(st_as_sf(bg, coords = c("x", "y"), crs = 4326))

plot(preds_used[[1]], main = "Roi, presence, background and Nigeria")
plot(roi, lwd = 1, border = "gray99", add = TRUE)
plot(bg_v, col = "red", cex = 0.1, add = TRUE)
plot(final_dataset_geo[[16]], col = "blue", add = TRUE)
plot(nigeria, border = "cyan", lwd = 4, add = TRUE)
# ---------------------------------------------------------------------------
# 5. Build the models
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 6. Ensemble for each species
# ---------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Spatial prediction
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Ensemble modelling
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Threshold selection (presence–absence conversion)
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Visualization of binary distribution maps
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Estimation of suitable habitat area
# -------------------------------------------------------------------------
# Calculate spatial extent (km²) of predicted suitable habitat


# -------------------------------------------------------------------------
# Niche similarity/overlap analysis
# -------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 9. Generate maps and figures
# ---------------------------------------------------------------------------




###### END OF SCRIPT ######################################################
