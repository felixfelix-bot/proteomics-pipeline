# Lydia's original R code — TurboID stats cutoffs + STRING mapping
# Reference implementation for known-interactor highlighting + STRING network overlay
# Adapted from: 2026_05_22_stats_cutoffs_turboID_map_to_STRING
# Provided by Dr. Aruna (collaborator) via Signal, July 12 2026

# SUMMARY OF LYDIA'S APPROACH:
#
# 1. VOLCANO CATEGORIES (priority order, last overwrites):
#    "ia"   = known interactors (from database)       → orange
#    "gp"   = gene/protein list (unspecified here)    → purple
#    "dhx"  = DHX helicase family                     → green
#    "ddx"  = DDX helicase family                     → blue
#    "high" = stringent: (log2FC>2 & -log10(q)>6)     → red
#             OR (log2FC>7 & -log10(q)>2)
#    "TRUE" = significant                             → teal
#    "FALSE"= not significant                         → grey60
#
# 2. STRING NETWORK INTEGRATION:
#    a. Map all proteins to STRINGdb (species 9606, v12.0)
#    b. Define SEEDS = high-confidence proteins (stringent threshold)
#    c. Get PHYSICAL neighbors of seeds
#    d. Build interaction network:
#       - Direct interactions with seeds (combined_score > 250)
#       - Secondary interactions (combined_score > 700)
#    e. Mark inNetwork = TRUE/FALSE for each protein
#    f. Create sig_network = paste(sig, inNetwork, sep="_")
#
# 3. TWO VOLCANO VIEWS:
#    a. Colored by inNetwork: red=high, teal=in network, grey=not
#    b. Colored by sig_network: combined significance + network membership
#
# KEY DIFFERENCE FROM OUR R/08_targeted_volcanos.R:
#   - Lydia overlays STRING physical interaction network membership
#   - Lydia highlights known interactors WITHOUT requiring them to be significant
#   - Lydia uses different significance thresholds for seeds vs general hits
#   - Lydia uses raw t-test p-values; we use pre-computed limma adj.P.Val

# ... (full original code preserved below for reference) ...

# === LYDIA'S ORIGINAL CODE (unchanged) ===

# in-depth evaluation of TurboID conditions
# Hela TurboID
# strong batch effect between BK467 and BK504
# comparison background to TRIP only valid within batch.. 0.2 TRIP4 vs 0 NA, 0 TRIP4 vs 0 NA
# 0.2 TRIP4 vs 0 TRIP4 RA specific differences
# 0.4 TRIP4 vs 0 TRIP4 RA specific differences
# overlap 0.4 and 0.2?
# focus on subset that is specific to ctrl comparison (pool results of 0.2 and 0 enrichments?)

setwd("/Volumes/agherzel/HerzelLydia/collaborations/Wahl_ASCC/proteomics/stringDBmapping/")

# LIBRARIES
library(RColorBrewer)
require("ggplot2")
require(gridExtra)
source("~/OneDrive/scripts/R_functions/data_loading_functions.R")
source("~/OneDrive/scripts/R_functions/mergeList.R")
library(tidyr)
library(ggpubr)
require('biomaRt')
library(STRINGdb)
library(dplyr)
library(igraph)

# FUNCTIONS
get_phys_interactions <- function(string_ids, phys_edges = phys) {
  phys_edges %>% filter(from %in% string_ids | to %in% string_ids)
}

get_phys_neighbors <- function(string_ids, phys_edges = phys) {
  partners <- phys_edges %>%
    filter(from %in% string_ids | to %in% string_ids) %>%
    transmute(partner = if_else(from %in% string_ids, to, from))
  unique(partners$partner)
}

# INPUT
ds1 <- read.delim("../Hela_TurboID_overview.txt", header = T, stringsAsFactors = F)
ia <- read.table("../../TRIP4_knownInteractors_20260422.txt", header = F, stringsAsFactors = F)[,1]
log2intensities1 <- read.delim("../Hela_TurboID_log2intensities_per_sample_wEnsemblID.txt", header = T, stringsAsFactors = F)
sub <- apply(log2intensities1[,10:ncol(log2intensities1)], 2, as.numeric)
phys <- read.table("9606.protein.physical.links.v12.0.txt", header = T, stringsAsFactors = F)
names(phys) <- c("from", "to", "combined_score")

# ... (full code continues as provided) ...
# See R/21_lydia_network_volcano.R for our pipeline-adapted version
