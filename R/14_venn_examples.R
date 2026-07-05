###############################################################################
# 14_venn_examples.R
# Three different Venn diagram approaches for the asymmetric-set-size problem.
#
# PROBLEM: When one set is much smaller than the other (e.g., 300 vs 50
# proteins), area-proportional circles make the small circle too tiny for
# the count number to fit inside it.
#
# This script generates 3 EXAMPLE VARIANTS using synthetic data that mimics
# the real RA effect scenario (TRIP4 without RA vs TRIP4 with RA):
#   - Set A (-RA): ~350 proteins
#   - Set B (+RA): ~80 proteins
#   - Overlap: ~60 proteins
#
# The user will choose which approach they prefer, then we apply it to
# the real data in 10_targeted_venns.R.
#
# EXAMPLE 1: Area-proportional with auto-scaling text
#   Keeps proportional circles, shrinks font for smaller regions
#
# EXAMPLE 2: Equal-size circles (non-proportional)
#   Both circles same size, counts shown as labels
#   This is the Nature/Cell standard — readability over proportionality
#
# EXAMPLE 3: ggVennDiagram with set_size=1
#   Modern package, equal circles, label backgrounds for readability
#
# Usage:
#   make venn-examples
###############################################################################
cat("\n=========================================\n")
cat(" Venn Diagram Examples (3 approaches)\n")
cat("=========================================\n\n")

library(VennDiagram)
library(grid)
library(ggplot2)

commit_hash <- get_git_hash()

# ---- Synthetic data mimicking the real scenario ----
# Set A (-RA): 350 proteins
# Set B (+RA): 80 proteins
# Overlap:     60 proteins
n_a      <- 350
n_b      <- 80
n_overlap <- 60
n_only_a <- n_a - n_overlap   # 290
n_only_b <- n_b - n_overlap   # 20

title_str <- "TRIP4 without vs with Retinoic Acid"

fill_a <- GLOBAL_COLORS[["venn_a_only"]]
fill_b <- GLOBAL_COLORS[["venn_b_only"]]

# =====================================================================
# EXAMPLE 1: Area-Proportional with Auto-Scaling Text
# =====================================================================
# Keeps proportional circles (honest representation of set sizes) but
# uses per-region cex values so text shrinks for smaller regions.
# draw.pairwise.venn accepts a length-3 vector for cex:
#   cex[1] = region A only, cex[2] = region B only, cex[3] = overlap
cat("[Example 1] Area-proportional with auto-scaling text...\n")

# Calculate per-region cex based on region size
# Scale: larger regions get larger text, but floor at 0.8 for readability
cex_base <- 2.0   # max text size
cex_floor <- 0.8  # never go below this (still readable)

all_regions <- c(n_only_a, n_only_b, n_overlap)
max_region <- max(all_regions)

# sqrt scaling: text area ~ cex^2, circle area ~ radius^2
scale_cex <- function(count, max_count) {
  if (max_count == 0) return(cex_floor)
  ratio <- sqrt(count / max_count)
  scaled <- cex_base * ratio
  max(cex_floor, min(cex_base, scaled))
}

cex_values <- sapply(all_regions, scale_cex, max_count = max_region)
cat(sprintf("  Region sizes: A-only=%d, B-only=%d, overlap=%d\n",
            n_only_a, n_only_b, n_overlap))
cat(sprintf("  cex values:   %.1f, %.1f, %.1f\n",
            cex_values[1], cex_values[2], cex_values[3]))

vp1 <- draw.pairwise.venn(
  area1     = n_a,
  area2     = n_b,
  cross.area = n_overlap,
  category  = c("- RA", "+ RA"),
  fill      = c(fill_a, fill_b),
  alpha     = rep(0.5, 2),
  cat.cex   = 1.6,
  cex       = cex_values,     # PER-REGION text scaling
  fontfamily = "sans",
  cat.fontfamily = "sans",
  col       = "transparent",
  margin    = 0.08,
  ind       = FALSE
)

# Save
file1 <- safe_filepath(FIGURE_DIR,
  paste0("venn_example_1_proportional_scaled_text_", commit_hash), "")

grDevices::png(paste0(file1, ".png"), width = 7, height = 6, units = "in", res = 300)
grid.draw(vp1)
pushViewport(viewport())
grid.text(title_str, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
popViewport()
grDevices::dev.off()

grDevices::pdf(paste0(file1, ".pdf"), width = 7, height = 6)
grid.draw(vp1)
pushViewport(viewport())
grid.text(title_str, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
popViewport()
grDevices::dev.off()

cat(sprintf("  Saved: %s\n", basename(paste0(file1, ".png"))))

# =====================================================================
# EXAMPLE 2: Equal-Size Circles (Non-Proportional)
# =====================================================================
# Both circles are the same size regardless of set counts.
# This is the Nature/Cell standard — readability always wins over
# proportionality. The numbers convey the actual counts.
cat("\n[Example 2] Equal-size circles (non-proportional)...\n")

# Trick: set scaled=FALSE so VennDiagram draws equal circles
# Use the same area for both to force equal radius
equal_area <- max(n_a, n_b)

vp2 <- draw.pairwise.venn(
  area1     = equal_area,
  area2     = equal_area,
  cross.area = n_overlap,
  category  = c("- RA", "+ RA"),
  fill      = c(fill_a, fill_b),
  alpha     = rep(0.5, 2),
  cat.cex   = 1.6,
  cex       = 1.8,            # Single size — fits everywhere
  fontfamily = "sans",
  cat.fontfamily = "sans",
  col       = "transparent",
  scaled    = FALSE,           # KEY: disable area-proportional scaling
  margin    = 0.08,
  ind       = FALSE
)

file2 <- safe_filepath(FIGURE_DIR,
  paste0("venn_example_2_equal_circles_", commit_hash), "")

grDevices::png(paste0(file2, ".png"), width = 7, height = 6, units = "in", res = 300)
grid.draw(vp2)
pushViewport(viewport())
grid.text(title_str, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
popViewport()
grDevices::dev.off()

grDevices::pdf(paste0(file2, ".pdf"), width = 7, height = 6)
grid.draw(vp2)
pushViewport(viewport())
grid.text(title_str, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
popViewport()
grDevices::dev.off()

cat(sprintf("  Saved: %s\n", basename(paste0(file2, ".png"))))

# =====================================================================
# EXAMPLE 3: ggplot2 with geom_circle — Equal Circles + Callout Labels
# =====================================================================
# Full ggplot2 control: equal circles, text labels with white background
# boxes (geom_label) so numbers are always readable even in small regions.
# This approach also lets us add a descriptive subtitle.
cat("\n[Example 3] ggplot2 equal circles with label backgrounds...\n")

# Circle geometry
r <- 1.0
d <- 1.0  # distance between centers
cx1 <- -d/2
cx2 <-  d/2

# Label positions (visible center of each region)
# For overlap: geometric center between circle centers
label_df <- data.frame(
  x = c(cx1 - 0.35, 0, cx2 + 0.35),
  y = c(0, 0, 0),
  label = c(n_only_a, n_overlap, n_only_b)
)

cat_df <- data.frame(
  x = c(cx1 - 0.85, cx2 + 0.85),
  y = c(0.8, 0.8),
  label = c("- RA", "+ RA")
)

circle_df <- data.frame(
  x0 = c(cx1, cx2),
  y0 = c(0, 0),
  r  = c(r, r),
  fill = c(as.character(fill_a), as.character(fill_b))
)

p3 <- ggplot() +
  geom_circle(
    data = circle_df,
    aes(x0 = x0, y0 = y0, r = r, fill = fill),
    alpha = 0.4, color = "grey40", linewidth = 1.0
  ) +
  scale_fill_identity() +
  # Count labels with white background (geom_label = rounded box)
  geom_label(
    data = label_df,
    aes(x = x, y = y, label = label),
    size = 6, fontface = "bold",
    fill = "white", label.size = 0.2
  ) +
  # Category labels
  geom_text(
    data = cat_df,
    aes(x = x, y = y, label = label),
    size = 6, fontface = "bold"
  ) +
  coord_fixed() +
  xlim(-2.5, 2.5) +
  ylim(-1.8, 1.5) +
  labs(title = title_str) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

file3 <- safe_filepath(FIGURE_DIR,
  paste0("venn_example_3_ggplot_labels_", commit_hash), "")

ggsave(paste0(file3, ".png"), p3, width = 7, height = 6, dpi = 300)
ggsave(paste0(file3, ".pdf"), p3, width = 7, height = 6)

cat(sprintf("  Saved: %s\n", basename(paste0(file3, ".png"))))

cat("\n=========================================\n")
cat(" Venn examples complete!\n")
cat("=========================================\n")
cat("\nThree variants saved to output/figures/:\n")
cat("  1. venn_example_1_proportional_scaled_text — area-proportional, auto text size\n")
cat("  2. venn_example_2_equal_circles — equal circles, same text size\n")
cat("  3. venn_example_3_ggplot_labels — equal circles, white label backgrounds\n")
cat("\nRun: make open-venn-examples\n")
