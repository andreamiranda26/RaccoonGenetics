#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(vegan)
  library(terra)   # only for point-to-point distance (no sf needed)
})

# -------------------------
# PATHS / SETTINGS
# -------------------------
ROOT_DIR <- "/scratch/Raccoon_andrea/gabe_trial"

LABELS_FILE   <- file.path(ROOT_DIR, "scripts/PopulationGenetics/labels_with_landcover.csv")
GEN_DIST_RDS  <- file.path(ROOT_DIR, "scripts/PopulationGenetics/gen_dist.rds") # this was from the Isolation by Distance pipeline.

OUT_DIR <- file.path(ROOT_DIR, "results/landscape_simple")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PM_CSV <- file.path(OUT_DIR, "partial_mantel_landcover_summary.csv")
OUT_INFO   <- file.path(OUT_DIR, "landcover_analysis_info.txt")

N_PERM <- 1000

# -------------------------
# HELPERS
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

# A safe way to build a dissimilarity matrix from a numeric vector
# (absolute differences between counties)
absdiff_dist <- function(x_named, county_order) {
  x <- x_named[county_order]
  if (any(is.na(x))) stop("NA values present in habitat variable for some counties.")
  mat <- abs(outer(x, x, "-"))
  rownames(mat) <- county_order
  colnames(mat) <- county_order
  as.dist(mat)
}

# -------------------------
# LOAD INPUTS
# -------------------------
if (!file.exists(LABELS_FILE)) stop("Missing: ", LABELS_FILE)
if (!file.exists(GEN_DIST_RDS)) stop("Missing: ", GEN_DIST_RDS)

message("Reading labels_with_landcover: ", LABELS_FILE)
labs <- read.csv(LABELS_FILE, stringsAsFactors = FALSE)

# required columns
req_cols <- c("sampleID", "county", "lat", "lon")
miss <- setdiff(req_cols, names(labs))
if (length(miss) > 0) stop("labels_with_landcover.csv missing columns: ", paste(miss, collapse = ", "))

# landcover columns expected (from your NLCD county extraction)
# If your column names differ, just change them here:
LC_COLS <- c("forest", "developed", "agriculture", "wetland")

miss_lc <- setdiff(LC_COLS, names(labs))
if (length(miss_lc) > 0) {
  stop("labels_with_landcover.csv is missing expected landcover columns: ",
       paste(miss_lc, collapse = ", "),
       "\nOpen the CSV and confirm column names (forest/developed/agriculture/wetland).")
}

message("Reading genetic distance: ", GEN_DIST_RDS)
gen_dist <- as_dist(readRDS(GEN_DIST_RDS))

gen_ids <- get_labels(gen_dist)
if (is.null(gen_ids)) stop("gen_dist has no Labels. Re-save gen_dist with sample labels from your IBD pipeline.")

# align sample->county
labs_s <- labs %>%
  transmute(sampleID = trimws(sampleID), county = trimws(county)) %>%
  filter(sampleID %in% gen_ids)

labs_s <- labs_s[match(gen_ids, labs_s$sampleID), ]
stopifnot(all(labs_s$sampleID == gen_ids))

# county centroid table (unique counties) + landcover composition
county_tbl <- labs %>%
  transmute(
    county = trimws(county),
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    forest = as.numeric(forest),
    developed = as.numeric(developed),
    agriculture = as.numeric(agriculture),
    wetland = as.numeric(wetland)
  ) %>%
  distinct()

# remove counties with NA landcover (if any)
county_tbl <- county_tbl %>% filter(!is.na(county))

counties <- sort(unique(county_tbl$county))
m <- length(counties)

message("Counties in analysis: ", m)

# -------------------------
# COUNTY-LEVEL GENETIC DISTANCE
# mean of sample-level genetic distances between counties
# -------------------------
message("Building COUNTY-level genetic distance matrix by averaging sample distances...")

gen_mat <- as.matrix(gen_dist)
rownames(gen_mat) <- gen_ids
colnames(gen_mat) <- gen_ids

gen_county <- matrix(NA_real_, nrow = m, ncol = m, dimnames = list(counties, counties))

for (i in seq_len(m)) {
  si <- labs_s$sampleID[labs_s$county == counties[i]]
  for (j in seq_len(m)) {
    sj <- labs_s$sampleID[labs_s$county == counties[j]]
    if (length(si) == 0 || length(sj) == 0) next
    gen_county[i, j] <- mean(gen_mat[si, sj], na.rm = TRUE)
  }
}

gen_county <- (gen_county + t(gen_county)) / 2
diag(gen_county) <- 0
gen_county_dist <- as.dist(gen_county)

# -------------------------
# COUNTY-LEVEL GEOGRAPHIC DISTANCE (km)
# -------------------------
message("Computing geographic distances among county centroids (km)...")

coords_df <- county_tbl[match(counties, county_tbl$county), c("lon", "lat")]
pts <- terra::vect(coords_df, geom = c("lon", "lat"), crs = "EPSG:4326")

geo_km <- as.matrix(terra::distance(pts)) / 1000
rownames(geo_km) <- counties
colnames(geo_km) <- counties
geo_county_dist <- as.dist(geo_km)

# -------------------------
# HABITAT DISSIMILARITY MATRICES (absolute differences)
# -------------------------
message("Building habitat dissimilarity matrices...")

# named vectors
forest_v    <- setNames(county_tbl$forest, county_tbl$county)
developed_v <- setNames(county_tbl$developed, county_tbl$county)
ag_v        <- setNames(county_tbl$agriculture, county_tbl$county)
wet_v       <- setNames(county_tbl$wetland, county_tbl$county)

forest_dist    <- absdiff_dist(forest_v, counties)
developed_dist <- absdiff_dist(developed_v, counties)
ag_dist        <- absdiff_dist(ag_v, counties)
wet_dist       <- absdiff_dist(wet_v, counties)

# -------------------------
# PARTIAL MANTEL TESTS
# genetic ~ habitat | geographic
# -------------------------
message("Running partial Mantel tests (", N_PERM, " permutations)...")

run_pm <- function(hab_dist, name) {
  res <- vegan::mantel.partial(gen_county_dist, hab_dist, geo_county_dist,
                              method = "pearson", permutations = N_PERM)
  data.frame(
    predictor = name,
    r = unname(res$statistic),
    p = res$signif,
    permutations = N_PERM,
    stringsAsFactors = FALSE
  )
}

tab <- bind_rows(
  run_pm(forest_dist, "forest_absdiff"),
  run_pm(developed_dist, "developed_absdiff"),
  run_pm(ag_dist, "agriculture_absdiff"),
  run_pm(wet_dist, "wetland_absdiff")
)

write.csv(tab, OUT_PM_CSV, row.names = FALSE)
message("Wrote: ", OUT_PM_CSV)
print(tab)

# -------------------------
# OPTIONAL: MRM (Multiple regression on distance matrices)
# -------------------------
mrm_note <- ""
if (requireNamespace("ecodist", quietly = TRUE)) {
  message("Running optional MRM (ecodist::MRM)...")
  # ecodist wants matrix form
  mrm_res <- ecodist::MRM(
    gen_county ~ geo_km + as.matrix(forest_dist) + as.matrix(developed_dist) + as.matrix(ag_dist) + as.matrix(wet_dist),
    nperm = N_PERM
  )
  sink(file.path(OUT_DIR, "MRM_output.txt"))
  print(mrm_res)
  sink()
  mrm_note <- "MRM_output.txt written (ecodist installed)."
} else {
  mrm_note <- "ecodist not installed; skipped MRM. Partial Mantels still complete."
  message(mrm_note)
}

# -------------------------
# WRITE INFO FILE
# -------------------------
sink(OUT_INFO)
cat("Landscape genetics (simple county-level) run info\n")
cat("Date:", format(Sys.time()), "\n\n")
cat("Inputs:\n")
cat("  labels_with_landcover:", LABELS_FILE, "\n")
cat("  gen_dist.rds:", GEN_DIST_RDS, "\n\n")
cat("Counties:", m, "\n")
cat("Permutations:", N_PERM, "\n")
cat("Landcover columns used:", paste(LC_COLS, collapse = ", "), "\n\n")
cat("Notes:\n")
cat("  Habitat dissimilarity matrices are absolute differences in % cover between counties.\n")
cat("  Partial Mantel tests evaluate genetic ~ habitat | geographic distance.\n")
cat("  ", mrm_note, "\n", sep = "")
sink()

message("Wrote: ", OUT_INFO)
message("DONE.")