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

out_csv  <- file.path(out_dir, "al_county_landcover_nlcd2024.csv")

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
lc_cat <- lc_long %>%
  mutate(category = case_when(
    value %in% c(41, 42, 43) ~ "forest",
    value %in% c(21, 22, 23, 24) ~ "developed",
    value %in% c(81, 82) ~ "agriculture",
    value %in% c(90, 95) ~ "wetland",
    value == 11 ~ "water",
    TRUE ~ "other"
  )) %>%
  group_by(GEOID, category) %>%
  summarise(percent = sum(percent), .groups = "drop") %>%
  pivot_wider(names_from = category, values_from = percent, values_fill = 0)

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