#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(igraph)
  library(maps)
})

set.seed(123)

ROOT_DIR <- "/scratch/Raccoon_andrea/gabe_trial"
vcfFile  <- file.path(ROOT_DIR, "data/vcf/Raccoon.filtered.vcf")
labFile  <- file.path(ROOT_DIR, "scripts/PopulationGenetics/labels_coords.csv")

out_dir <- file.path(ROOT_DIR, "results/PCA_Moran_county")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading labels: ", labFile)
labs <- read.csv(labFile, stringsAsFactors = FALSE)
labs$sampleID <- trimws(as.character(labs$sampleID))
labs$county   <- trimws(as.character(labs$county))
labs$lon      <- as.numeric(labs$lon)
labs$lat      <- as.numeric(labs$lat)

message("Reading VCF: ", vcfFile)
vcf <- read.vcfR(vcfFile, verbose = FALSE)
gl  <- vcfR2genlight(vcf)  # biallelic-only; multiallelic loci are dropped

vcf_ids <- trimws(as.character(indNames(gl)))
keep_ids <- intersect(vcf_ids, labs$sampleID)
if (length(keep_ids) < 10) stop("Too few overlapping samples.")

gl <- gl[match(keep_ids, vcf_ids), ]
ids2 <- indNames(gl)

labs2 <- labs[match(ids2, labs$sampleID), ]
stopifnot(all(labs2$sampleID == ids2))

# --- County sample sizes ---
county_counts <- as.data.frame(table(labs2$county), stringsAsFactors = FALSE)
colnames(county_counts) <- c("county", "N")
county_counts$county <- trimws(county_counts$county)
county_counts <- county_counts[order(county_counts$county), ]

write.csv(
  county_counts,
  file.path(out_dir, "county_sample_sizes.csv"),
  row.names = FALSE
)

# --- County-level allele-frequency matrix ---
message("Creating county-level allele-frequency matrix...")
X_ind <- adegenet::tab(gl, NA.method = "mean")

counties <- sort(unique(labs2$county))
X_cty <- sapply(counties, function(cty) {
  idx <- which(labs2$county == cty)
  colMeans(X_ind[idx, , drop = FALSE], na.rm = TRUE)
})
X_cty <- t(X_cty)
rownames(X_cty) <- counties

# --- County centroid coords (one per county) ---
county_pts <- aggregate(cbind(lon, lat) ~ county, data = labs2, FUN = function(z) z[1])
county_pts$county <- trimws(county_pts$county)
county_pts <- county_pts[match(counties, county_pts$county), ]
stopifnot(all(county_pts$county == counties))

coords <- as.matrix(county_pts[, c("lon", "lat")])

# --- PCA on county allele-freq matrix (all counties) ---
message("Running PCA with all counties...")
pca <- prcomp(X_cty, center = TRUE, scale. = TRUE)
scores <- pca$x[, 1:2, drop = FALSE]

# Percent variance explained
pve <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
pc1_lab <- sprintf("PC1 (%.1f%%)", pve[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", pve[2])

write.csv(
  data.frame(
    county = counties,
    lon = coords[, 1],
    lat = coords[, 2],
    PC1 = scores[, 1],
    PC2 = scores[, 2]
  ),
  file.path(out_dir, "county_PCA_scores.csv"),
  row.names = FALSE
)

# --- Build kNN weights (all counties) ---
D <- as.matrix(dist(coords))

build_knn_W <- function(Dmat, k) {
  m_local <- nrow(Dmat)
  nn <- apply(Dmat, 1, function(d) order(d)[2:(k + 1)])
  W <- matrix(0, nrow = m_local, ncol = m_local)
  for (i in 1:m_local) W[i, nn[, i]] <- 1
  W <- (W + t(W)) > 0
  W <- 1 * W
  rs <- rowSums(W)
  rs[rs == 0] <- 1
  W / rs
}

moran_I <- function(x, W) {
  x <- as.numeric(x)
  x <- x - mean(x)
  m_local <- length(x)
  S0 <- sum(W)
  num <- sum(W * (x %o% x))
  den <- sum(x^2)
  (m_local / S0) * (num / den)
}

moran_perm <- function(x, W, nperm = 1000) {
  obs <- moran_I(x, W)
  perm <- replicate(nperm, moran_I(sample(x), W))
  p <- (sum(abs(perm) >= abs(obs)) + 1) / (nperm + 1)
  list(I = obs, p = p, perm = perm)
}

# choose k automatically to ensure connectivity
k <- 5
repeat {
  W <- build_knn_W(D, k)
  g <- igraph::graph_from_adjacency_matrix((W > 0) * 1, mode = "undirected", diag = FALSE)
  comp <- igraph::components(g)$no
  message("All counties: k = ", k, " -> components = ", comp)
  if (comp == 1 || k >= 12) break
  k <- k + 1
}

message("Testing spatial autocorrelation of PCA axes (all counties)...")
res1 <- moran_perm(scores[, 1], W, nperm = 1000)
res2 <- moran_perm(scores[, 2], W, nperm = 1000)

sink(file.path(out_dir, "moran_results.txt"))
cat("County-level PCA + Moran's I spatial autocorrelation\n")
cat("All counties included\n")
cat("kNN k =", k, "\n\n")
cat("PC1 Moran's I:", res1$I, "p =", res1$p, "\n")
cat("PC2 Moran's I:", res2$I, "p =", res2$p, "\n")
sink()

# --- Robustness check: remove singleton counties (N = 1) ---
message("Running PCA excluding counties with N = 1...")
keep_counties <- county_counts$county[county_counts$N > 1]

labs_filt <- labs2[labs2$county %in% keep_counties, ]
gl_filt <- gl[which(labs2$county %in% keep_counties), ]

stopifnot(nInd(gl_filt) == nrow(labs_filt))

X_ind_filt <- adegenet::tab(gl_filt, NA.method = "mean")
counties_filt <- sort(unique(labs_filt$county))

X_cty_filt <- sapply(counties_filt, function(cty) {
  idx <- which(labs_filt$county == cty)
  colMeans(X_ind_filt[idx, , drop = FALSE], na.rm = TRUE)
})
X_cty_filt <- t(X_cty_filt)
rownames(X_cty_filt) <- counties_filt

county_pts_filt <- aggregate(cbind(lon, lat) ~ county, data = labs_filt, FUN = function(z) z[1])
county_pts_filt$county <- trimws(county_pts_filt$county)
county_pts_filt <- county_pts_filt[match(counties_filt, county_pts_filt$county), ]
stopifnot(all(county_pts_filt$county == counties_filt))

coords_filt <- as.matrix(county_pts_filt[, c("lon", "lat")])

pca_filt <- prcomp(X_cty_filt, center = TRUE, scale. = TRUE)
scores_filt <- pca_filt$x[, 1:2, drop = FALSE]

pve_filt <- (pca_filt$sdev^2 / sum(pca_filt$sdev^2)) * 100
pc1_lab_filt <- sprintf("PC1 (%.1f%%)", pve_filt[1])
pc2_lab_filt <- sprintf("PC2 (%.1f%%)", pve_filt[2])

write.csv(
  data.frame(
    county = counties_filt,
    lon = coords_filt[, 1],
    lat = coords_filt[, 2],
    PC1 = scores_filt[, 1],
    PC2 = scores_filt[, 2]
  ),
  file.path(out_dir, "county_PCA_scores_noSingletons.csv"),
  row.names = FALSE
)

# --- Moran's I for filtered set ---
D_filt <- as.matrix(dist(coords_filt))
k_filt <- 5
repeat {
  W_filt <- build_knn_W(D_filt, k_filt)
  g_filt <- igraph::graph_from_adjacency_matrix((W_filt > 0) * 1, mode = "undirected", diag = FALSE)
  comp_filt <- igraph::components(g_filt)$no
  message("No singletons: k = ", k_filt, " -> components = ", comp_filt)
  if (comp_filt == 1 || k_filt >= 12) break
  k_filt <- k_filt + 1
}

res1_filt <- moran_perm(scores_filt[, 1], W_filt, nperm = 1000)
res2_filt <- moran_perm(scores_filt[, 2], W_filt, nperm = 1000)

sink(file.path(out_dir, "moran_results_noSingletons.txt"))
cat("County-level PCA + Moran's I spatial autocorrelation\n")
cat("Singleton counties removed (N > 1 only)\n")
cat("kNN k =", k_filt, "\n\n")
cat("PC1 Moran's I:", res1_filt$I, "p =", res1_filt$p, "\n")
cat("PC2 Moran's I:", res2_filt$I, "p =", res2_filt$p, "\n")
sink()

# =========================================================
# NEW FIGURE: PCA panels colored by latitude + Alabama map
# =========================================================

# Shared latitude scale based on ALL counties
lat_min <- min(county_pts$lat, na.rm = TRUE)
lat_max <- max(county_pts$lat, na.rm = TRUE)

# South -> North gradient
pal <- colorRampPalette(c("#FF8B00", "#C06832", "#4059CA"))
#pal <- colorRampPalette(c("#FF8B00", "#4059CA"))  # orange (south) → blue (north)
ncols <- 100

lat_to_col <- function(lat, min_lat, max_lat, palette_fun, n = 100) {
  scaled <- (lat - min_lat) / (max_lat - min_lat)
  scaled[scaled < 0] <- 0
  scaled[scaled > 1] <- 1
  cols <- palette_fun(n)
  cols[pmax(1, pmin(n, floor(scaled * (n - 1)) + 1))]
}

cols_all  <- lat_to_col(county_pts$lat,      lat_min, lat_max, pal, ncols)
cols_filt <- lat_to_col(county_pts_filt$lat, lat_min, lat_max, pal, ncols)

png(file.path(out_dir, "Figure2_PCA_latitude_gradient.png"),
    width = 2400, height = 2200, res = 250)

layout(
  matrix(c(1, 2,
           3, 3,
           4, 4),
         nrow = 3, byrow = TRUE),
  heights = c(1, 1, 0.25)
)

# --- Panel A: all counties PCA ---
par(mar = c(5, 5, 2, 1))
plot(
  scores[, 1], scores[, 2],
  pch = 21, bg = cols_all, col = "black", cex = 1.7,
  xlab = pc1_lab,
  ylab = pc2_lab,
  main = ""
)
abline(h = 0, lty = 2, col = "grey75")
abline(v = 0, lty = 2, col = "grey75")
mtext("A", side = 3, adj = -0.08, line = 0.5, font = 2, cex = 1.8)

# --- Panel B: no singletons PCA ---
par(mar = c(5, 5, 2, 1))
plot(
  scores_filt[, 1], scores_filt[, 2],
  pch = 21, bg = cols_filt, col = "black", cex = 1.7,
  xlab = pc1_lab_filt,
  ylab = pc2_lab_filt,
  main = ""
)
abline(h = 0, lty = 2, col = "grey75")
abline(v = 0, lty = 2, col = "grey75")
mtext("B", side = 3, adj = -0.08, line = 0.5, font = 2, cex = 1.8)

# --- Panel C: Alabama map with county centroid colors ---
par(mar = c(1.5, 1.5, 2.5, 1.5))  # tighter margins = bigger map

map(
  "state", "alabama",
  fill = TRUE,
  col = "grey95",
  border = "grey40",
  bg = "white",
  lwd = 1.5
)

points(
  county_pts$lon, county_pts$lat,
  pch = 21,
  bg = cols_all,
  col = "black",
  cex = 1.9   # slightly bigger points
)

title(
  "Counties colored by latitude (south–north gradient)",
  cex.main = 1.2,
  line = 0.5
)

# Move "C" label so it doesn’t overlap title
mtext("C", side = 3, adj = -0.03, line = 1.4, font = 2, cex = 1.8)

# --- Legend panel ---
par(mar = c(0, 2, 0, 2))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

text(0.08, 0.82, labels = "Latitude gradient", adj = 0, cex = 1.2, font = 2)

legend_cols <- pal(ncols)
legend_x <- seq(0.08, 0.78, length.out = ncols + 1)

for (i in 1:ncols) {
  rect(
    legend_x[i], 0.48, legend_x[i + 1], 0.66,
    col = legend_cols[i], border = legend_cols[i]
  )
}

text(
  0.08, 0.34,
  labels = paste0("South (", round(lat_min, 2), "°N)"),
  adj = 0, cex = 1
)
text(
  0.78, 0.34,
  labels = paste0("North (", round(lat_max, 2), "°N)"),
  adj = 1, cex = 1
)

text(
  0.08, 0.10,
  labels = paste(
    "Panels A and B: PCA points colored by county latitude.",
    "Panel C: Alabama county centroids shown using the same color scale.",
    sep = "\n"
  ),
  adj = 0, cex = 1
)

dev.off()

message("DONE. Outputs in: ", out_dir)