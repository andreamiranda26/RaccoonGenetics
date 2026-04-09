#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(vegan)
})

# -------------------------
# USER SETTINGS / PATHS
# -------------------------
ROOT_DIR <- "/scratch/Raccoon_andrea/gabe_trial"

# sample metadata (must contain: sampleID, county, lat, lon)
LABELS_FILE <- file.path(ROOT_DIR, "scripts/PopulationGenetics/labels_with_landcover.csv")

# saved from your IBD pipeline
GEN_DIST_RDS <- file.path(ROOT_DIR, "scripts/PopulationGenetics/gen_dist.rds")
GEO_DIST_RDS <- file.path(ROOT_DIR, "scripts/PopulationGenetics/geo_dist.rds")

# NHD flowlines (Alabama) you downloaded/unzipped
NHD0 <- file.path(ROOT_DIR, "data/rivers/nhd_alabama/Shape/NHDFlowline_0.shp")
NHD1 <- file.path(ROOT_DIR, "data/rivers/nhd_alabama/Shape/NHDFlowline_1.shp")

# outputs
OUT_DIR <- file.path(ROOT_DIR, "scripts/PopulationGenetics/river_barrier")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

RIVER_LINE_GPKG <- file.path(OUT_DIR, "alabama_coosa_river_line.gpkg")
RIVER_CROSS_COUNTY_RDS <- file.path(OUT_DIR, "river_cross_county_matrix.rds")
RIVER_CROSS_SAMPLE_RDS <- file.path(OUT_DIR, "river_cross_sample_dist.rds")
PARTIAL_MANTEL_TXT <- file.path(OUT_DIR, "partial_mantel_river_crossing.txt")
RIVER_NAMES_TXT <- file.path(OUT_DIR, "river_name_values_used.txt")

# mantel settings
N_PERM <- 1000

# -------------------------
# HELPER FUNCTIONS
# -------------------------
as_dist <- function(x) {
  if (inherits(x, "dist")) return(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) return(as.dist(x))
  stop("Object must be a dist or a symmetric matrix/data.frame.")
}

get_labels <- function(x) {
  if (inherits(x, "dist")) return(attr(x, "Labels"))
  if (is.matrix(x) || is.data.frame(x)) return(rownames(x))
  return(NULL)
}

# -------------------------
# LOAD INPUTS
# -------------------------
message("Reading labels: ", LABELS_FILE)
labs <- read.csv(LABELS_FILE, stringsAsFactors = FALSE)

req_cols <- c("sampleID", "county", "lat", "lon")
missing_cols <- setdiff(req_cols, names(labs))
if (length(missing_cols) > 0) {
  stop("labels_with_landcover.csv is missing columns: ", paste(missing_cols, collapse = ", "))
}

if (!file.exists(GEN_DIST_RDS)) stop("Missing: ", GEN_DIST_RDS)
if (!file.exists(GEO_DIST_RDS)) stop("Missing: ", GEO_DIST_RDS)

message("Reading genetic distance: ", GEN_DIST_RDS)
gen_dist <- as_dist(readRDS(GEN_DIST_RDS))

message("Reading geographic distance: ", GEO_DIST_RDS)
geo_dist <- as_dist(readRDS(GEO_DIST_RDS))

gen_ids <- get_labels(gen_dist)
geo_ids <- get_labels(geo_dist)

if (is.null(gen_ids) || is.null(geo_ids)) {
  stop("gen_dist and/or geo_dist do not contain sample labels. Ensure dist objects have Labels or matrices have rownames.")
}
if (!setequal(gen_ids, geo_ids)) {
  stop("gen_dist and geo_dist do not contain the same sampleIDs.")
}

# reorder labels to match dist label order
labs2 <- labs %>% filter(sampleID %in% gen_ids)
labs2 <- labs2[match(gen_ids, labs2$sampleID), ]
if (any(is.na(labs2$sampleID))) stop("Some sampleIDs in distances were not found in labels file.")

message("Samples in analysis: ", nrow(labs2))
message("Unique counties in analysis: ", dplyr::n_distinct(labs2$county))

# -------------------------
# READ NHD FLOWLINES (AL + COOSA)
# -------------------------
if (!file.exists(NHD0)) stop("Missing NHD file: ", NHD0)
if (!file.exists(NHD1)) stop("Missing NHD file: ", NHD1)

message("Reading NHD flowlines...")
flow0 <- st_read(NHD0, quiet = TRUE)
flow1 <- st_read(NHD1, quiet = TRUE)

# Make sure fields match (some datasets vary slightly across tiles)
common_fields <- intersect(names(flow0), names(flow1))
flowlines <- rbind(flow0[, common_fields], flow1[, common_fields])

# ---- auto-detect a usable name field (handles lowercase NHD schemas) ----
name_candidates <- c(
  "GNIS_NAME", "GNIS_Name", "GNISNAME", "gnis_name", "gnisname",
  "FULLNAME", "FullName", "FULL_NAME", "fullname", "full_name",
  "NAME", "Name", "name", "RIVERNAME", "RiverName", "rivername"
)

name_field <- name_candidates[name_candidates %in% names(flowlines)][1]

if (is.na(name_field)) {
  stop(
    "No recognized name field found in NHD flowlines.\nAvailable fields:\n",
    paste(names(flowlines), collapse = ", ")
  )
}

message("Using river name field: ", name_field)

if (is.na(name_field)) {
  stop(
    "No recognized name field found in NHD flowlines.\nAvailable fields:\n",
    paste(names(flowlines), collapse = ", ")
  )
}

message("Using river name field: ", name_field)

# ---- filter Alabama River + Coosa River (mainstem only) ----
rivers <- flowlines %>%
  mutate(.river_name = as.character(.data[[name_field]])) %>%
  filter(!is.na(.river_name)) %>%
  filter(grepl("^(Alabama River|Coosa River)$", .river_name, ignore.case = TRUE))

# If that returns nothing (some datasets omit "River"), fall back to broader match
if (nrow(rivers) == 0) {
  message("No exact 'Alabama River'/'Coosa River' matches; falling back to broader name match.")
  rivers <- flowlines %>%
    mutate(.river_name = as.character(.data[[name_field]])) %>%
    filter(!is.na(.river_name) & grepl("Alabama|Coosa", .river_name, ignore.case = TRUE))
}

message("Unique river names matched:")
print(sort(unique(rivers$.river_name)))

# save the unique names used (helpful for debugging / writeup)
writeLines(sort(unique(rivers$.river_name)), con = RIVER_NAMES_TXT)
message("Wrote river name values used to: ", RIVER_NAMES_TXT)

river_line <- rivers %>%
  st_transform(4326) %>%
  st_union() %>%
  st_cast("MULTILINESTRING")

# Save river line for reproducibility
river_sf <- st_sf(name = "Alabama_Coosa", geometry = st_sfc(river_line, crs = 4326))
st_write(river_sf, RIVER_LINE_GPKG, append = FALSE, quiet = TRUE)
message("Saved river line to: ", RIVER_LINE_GPKG)

# -------------------------
# BUILD COUNTY-LEVEL CROSSING MATRIX
# -------------------------
county_pts <- labs2 %>%
  distinct(county, lon, lat) %>%
  mutate(county = trimws(county))

m <- nrow(county_pts)
message("Computing river crossings among counties (m = ", m, ") ...")

county_sf <- st_as_sf(county_pts, coords = c("lon", "lat"), crs = 4326)

river_cross_county <- matrix(0L, nrow = m, ncol = m,
                             dimnames = list(county_pts$county, county_pts$county))

for (i in 1:(m - 1)) {
  pi <- st_coordinates(county_sf[i, ])
  for (j in (i + 1):m) {
    pj <- st_coordinates(county_sf[j, ])
    line <- st_sfc(st_linestring(rbind(pi, pj)), crs = 4326)
    cross <- as.integer(st_intersects(line, river_line, sparse = FALSE)[1, 1])
    river_cross_county[i, j] <- cross
    river_cross_county[j, i] <- cross
  }
}

saveRDS(river_cross_county, RIVER_CROSS_COUNTY_RDS)
message("Saved county crossing matrix: ", RIVER_CROSS_COUNTY_RDS)

# -------------------------
# MAP COUNTY CROSSINGS TO SAMPLE ORDER
# -------------------------
county_to_idx <- setNames(seq_len(m), county_pts$county)
sample_county <- trimws(labs2$county)

if (any(!sample_county %in% names(county_to_idx))) {
  missing <- unique(sample_county[!sample_county %in% names(county_to_idx)])
  stop("Some sample counties were not in county centroid list: ", paste(missing, collapse = ", "))
}

idx <- county_to_idx[sample_county]
n <- length(idx)

river_cross_sample <- river_cross_county[idx, idx]
dimnames(river_cross_sample) <- list(labs2$sampleID, labs2$sampleID)

river_cross_dist <- as.dist(river_cross_sample)
attr(river_cross_dist, "Labels") <- labs2$sampleID

saveRDS(river_cross_dist, RIVER_CROSS_SAMPLE_RDS)
message("Saved sample-level river crossing dist: ", RIVER_CROSS_SAMPLE_RDS)

# -------------------------
# PARTIAL MANTEL TEST
# -------------------------
message("Running partial Mantel: genetic ~ river_cross | geographic_distance")
res <- mantel.partial(gen_dist, river_cross_dist, geo_dist,
                      permutations = N_PERM, method = "pearson")

sink(PARTIAL_MANTEL_TXT)
cat("Partial Mantel test (genetic ~ river_cross | geographic distance)\n")
cat("Permutations:", N_PERM, "\n\n")
cat("River name field used: ", name_field, "\n")
cat("Unique river names matched were written to: ", RIVER_NAMES_TXT, "\n\n")
print(res)
sink()

message("Wrote Mantel output: ", PARTIAL_MANTEL_TXT)
message("DONE.")