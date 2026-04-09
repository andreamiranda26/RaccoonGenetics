#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
})

ROOT_DIR <- "/scratch/Raccoon_andrea/gabe_trial"

labels_file <- file.path(ROOT_DIR, "scripts/PopulationGenetics/labels_coords.csv")
lc_file     <- file.path(ROOT_DIR, "data/landcover/al2_county_landcover_nlcd2024.csv")

out_file    <- file.path(ROOT_DIR, "scripts/PopulationGenetics/labels_with_landcover2.csv")

labs <- read.csv(labels_file, stringsAsFactors = FALSE)
lc   <- read.csv(lc_file, stringsAsFactors = FALSE)

# normalize names for joining
labs$county <- trimws(labs$county)
lc$county   <- trimws(lc$county)

# fix common naming differences
labs$county <- gsub("StClair", "St. Clair", labs$county)

labs2 <- labs %>%
  left_join(lc %>% select(-GEOID), by = "county")
  
# quick check
missing <- sum(is.na(labs2$forest))
cat("Wrote:", out_file, "\n")
cat("Rows with missing forest %:", missing, "\n")
if (missing > 0) {
  cat("Counties missing (examples):\n")
  print(head(unique(labs2$county[is.na(labs2$forest)]), 20))
}

write.csv(labs2, out_file, row.names = FALSE)