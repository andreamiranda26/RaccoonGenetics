#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(tigris)
  library(exactextractr)
  library(dplyr)
  library(tidyr)
})

options(tigris_use_cache = TRUE)

# ---- paths ----
ROOT_DIR <- "/scratch/Raccoon_andrea/gabe_trial"
NLCD_TIF <- file.path(ROOT_DIR, "data/nlcd/Annual_NLCD_LndCov_2024_CU_C1V1.tif")

out_dir  <- file.path(ROOT_DIR, "data/landcover")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_csv  <- file.path(out_dir, "al2_county_landcover_nlcd2024.csv")

# ---- terra temp dir (speed + avoids small /tmp) ----
tmp_dir <- file.path(ROOT_DIR, "tmp")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = tmp_dir, progress = 1)

# ---- checks ----
if (!file.exists(NLCD_TIF)) {
  stop("NLCD tif not found at: ", NLCD_TIF,
       "\nPut the NLCD .tif at that path (recommended), then rerun.")
}

message("Reading NLCD raster: ", NLCD_TIF)
nlcd <- terra::rast(NLCD_TIF)

message("Loading Alabama counties (tigris)...")
al_counties <- tigris::counties(state = "AL", year = 2020, class = "sf")

message("Transforming counties to NLCD projection...")
al_counties <- st_transform(al_counties, crs(nlcd))

message("Cropping and masking NLCD to Alabama (speeds extraction)...")
nlcd_al <- terra::crop(nlcd, terra::vect(al_counties))
nlcd_al <- terra::mask(nlcd_al, terra::vect(al_counties))

# ---- quick check: confirm land cover codes look like NLCD ----
message("Quick check: first few land cover class counts in cropped raster:")
fq <- terra::freq(nlcd_al, value = TRUE)
print(head(fq, 15))

message("If values look like 11, 21-24, 41-43, 81-82, 90, 95, recoding below will work.")

message("Extracting raster values with coverage fractions per county...")
# Returns list of data frames, one per county polygon, with 'value' and 'coverage_fraction'
raw_list <- exactextractr::exact_extract(nlcd_al, al_counties)

message("Summarizing to percent cover...")
lc_long <- bind_rows(lapply(seq_along(raw_list), function(i) {
  df <- raw_list[[i]]
  if (is.null(df) || nrow(df) == 0) return(NULL)

  if (!("coverage_fraction" %in% names(df))) {
    stop("Expected column 'coverage_fraction' not found in exact_extract output.")
  }

  df %>%
    filter(!is.na(value)) %>%
    group_by(value) %>%
    summarise(frac = sum(coverage_fraction, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      GEOID = al_counties$GEOID[i],
      percent = 100 * frac / sum(frac)
    ) %>%
    select(GEOID, value, percent)
}))

# Recode to categories meaningful for raccoons
# This second version recodes to finer NLCD classes (keep agriculture lumped)
lc_cat <- lc_long %>%
  mutate(category = case_when(
    # Developed intensity (separate)
    value == 21 ~ "dev_open",
    value == 22 ~ "dev_low",
    value == 23 ~ "dev_med",
    value == 24 ~ "dev_high",

    # Water
    value == 11 ~ "water",

    # Forest (separate; optional but included)
    value == 41 ~ "forest_decid",
    value == 42 ~ "forest_ever",
    value == 43 ~ "forest_mix",

    # Wetlands (separate; optional but included)
    value == 90 ~ "wet_woody",
    value == 95 ~ "wet_emerg",

    # Agriculture (lumped)
    value %in% c(81, 82) ~ "agriculture",

    TRUE ~ "other"
  )) %>%
  group_by(GEOID, category) %>%
  summarise(percent = sum(percent), .groups = "drop") %>%
  pivot_wider(names_from = category, values_from = percent, values_fill = 0) %>%
  # OPTIONAL: create convenience combined columns too (handy later)
  mutate(
    forest = forest_decid + forest_ever + forest_mix,
    developed = dev_open + dev_low + dev_med + dev_high,
    wetland = wet_woody + wet_emerg
  )
# Add county name for joining
county_key <- al_counties %>%
  st_drop_geometry() %>%
  transmute(GEOID = GEOID, county = NAME)

lc_final <- lc_cat %>%
  left_join(county_key, by = "GEOID") %>%
  relocate(county)

# Write output
write.csv(lc_final, out_csv, row.names = FALSE)
message("Wrote: ", out_csv)

message("Preview:")
print(head(lc_final, 10))