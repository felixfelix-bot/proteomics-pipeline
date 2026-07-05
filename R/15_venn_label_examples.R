###############################################################################
# 15_venn_label_examples.R
# Three label-placement variants for the TurboID vs Flag IP Venn diagram.
#
# PROBLEM: Category labels ("TurboID", "Flag IP") overlap the circles.
# This script generates 3 variants with different label positions so
# the researchers can choose their preferred style.
#
# VARIANT 1: Labels below circles (cat.pos = c(200, 160))
# VARIANT 2: Labels far outside circles (cat.dist = 0.12, cat.pos = c(-30, 30))
# VARIANT 3: No category labels on diagram — labels in title/subtitle only
#
# Uses synthetic data mimicking real TurboID vs Flag IP scenario.
#
# Usage:
#   make venn-label-examples
###############################################################################
cat("\n=========================================\n")
cat(" Venn Label Position Examples\n")
cat("=========================================\n\n")

library(VennDiagram)
library(grid)

commit_hash <- get_git_hash()

# ---- Synthetic data mimicking TurboID vs Flag IP ----
n_turbo   <- 500
n_flag    <- 350
n_overlap <- 120

title_str <- "TurboID vs Flag IP Overlap"

fill_a <- GLOBAL_COLORS[["venn_a_only"]]
fill_b <- GLOBAL_COLORS[["venn_b_only"]]

# Helper to save a Venn viewport to file
save_venn <- function(vp, file_prefix, title) {
  safe_name <- sanitize_filename(file_prefix)
  versioned_name <- paste0(safe_name, "_", commit_hash)
  png_path <- file.path(FIGURE_DIR, paste0(versioned_name, ".png"))
  pdf_path <- file.path(FIGURE_DIR, paste0(versioned_name, ".pdf"))

  grDevices::png(png_path, width = 7, height = 6, units = "in", res = 300)
  grid.draw(vp)
  pushViewport(viewport())
  grid.text(title, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  popViewport()
  grDevices::dev.off()

  grDevices::pdf(pdf_path, width = 7, height = 6)
  grid.draw(vp)
  pushViewport(viewport())
  grid.text(title, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  popViewport()
  grDevices::dev.off()

  cat(sprintf("  Saved: %s\n", basename(png_path)))
}

# =====================================================================
# VARIANT 1: Labels BELOW the circles
# =====================================================================
# cat.pos uses degrees: 0 = right, 90 = top, 180 = left, 270 = bottom
# c(220, 140) places labels at lower-left and lower-right
cat("[Variant 1] Labels below circles...\n")

vp1 <- draw.pairwise.venn(
  area1     = n_turbo,
  area2     = n_flag,
  cross.area = n_overlap,
  category  = c("TurboID", "Flag IP"),
  fill      = c(fill_a, fill_b),
  alpha     = rep(0.5, 2),
  cat.cex   = 1.6,
  cex       = 2.0,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  col       = "transparent",
  cat.pos   = c(220, 140),     # Below the circles
  cat.dist  = c(0.06, 0.06),
  margin    = 0.08,
  ind       = FALSE
)
save_venn(vp1, "venn_label_v1_below", title_str)

# =====================================================================
# VARIANT 2: Labels FAR OUTSIDE circles (pushed away)
# =====================================================================
# Keep labels at top but increase distance so they don't touch circles
cat("\n[Variant 2] Labels far outside circles...\n")

vp2 <- draw.pairwise.venn(
  area1     = n_turbo,
  area2     = n_flag,
  cross.area = n_overlap,
  category  = c("TurboID", "Flag IP"),
  fill      = c(fill_a, fill_b),
  alpha     = rep(0.5, 2),
  cat.cex   = 1.6,
  cex       = 2.0,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  col       = "transparent",
  cat.pos   = c(-30, 30),       # Upper edges but spread wider
  cat.dist  = c(0.12, 0.12),    # Pushed much further out
  margin    = 0.12,
  ind       = FALSE
)
save_venn(vp2, "venn_label_v2_far_out", title_str)

# =====================================================================
# VARIANT 3: No category labels on diagram — shown in title only
# =====================================================================
# Remove category labels entirely, put "TurboID vs Flag IP" in the title
cat("\n[Variant 3] No labels on circles, title only...\n")

vp3 <- draw.pairwise.venn(
  area1     = n_turbo,
  area2     = n_flag,
  cross.area = n_overlap,
  category  = c("", ""),        # Empty category labels
  fill      = c(fill_a, fill_b),
  alpha     = rep(0.5, 2),
  cat.cex   = 0,                # Zero-size category text
  cex       = 2.0,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  col       = "transparent",
  margin    = 0.08,
  ind       = FALSE
)

# Save with combined title that includes category names
safe_name <- sanitize_filename("venn_label_v3_title_only")
versioned_name <- paste0(safe_name, "_", commit_hash)
png_path <- file.path(FIGURE_DIR, paste0(versioned_name, ".png"))
pdf_path <- file.path(FIGURE_DIR, paste0(versioned_name, ".pdf"))

full_title <- "TurboID vs Flag IP Overlap"

grDevices::png(png_path, width = 7, height = 6, units = "in", res = 300)
grid.draw(vp3)
pushViewport(viewport())
grid.text(full_title, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("TurboID (left) vs Flag IP (right)", 0.5, 0.05,
          gp = gpar(fontsize = 10, fontface = "italic"))
popViewport()
grDevices::dev.off()

grDevices::pdf(pdf_path, width = 7, height = 6)
grid.draw(vp3)
pushViewport(viewport())
grid.text(full_title, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("TurboID (left) vs Flag IP (right)", 0.5, 0.05,
          gp = gpar(fontsize = 10, fontface = "italic"))
popViewport()
grDevices::dev.off()

cat(sprintf("  Saved: %s\n", basename(png_path)))

cat("\n=========================================\n")
cat(" Venn label examples complete!\n")
cat("=========================================\n")
cat("\nThree variants saved to output/figures/:\n")
cat("  v1_below     — Labels below circles\n")
cat("  v2_far_out   — Labels pushed far outside circles\n")
cat("  v3_title_only — No labels on circles, names in subtitle\n")
cat("\nRun: make open-venn-label-examples\n")
