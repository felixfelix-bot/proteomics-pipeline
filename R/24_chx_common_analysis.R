###############################################################################
# 24_chx_common_analysis.R
# CHX/DMSO common protein analysis for the TRIP4/ASCC proteomics pipeline.
#
# WHAT THIS DOES (requested by Dr. Aruna, July 12):
#   1. Loads the TRIP4_CHX_vs_TRIP4_DMSO experiment
#   2. Splits significant proteins into two directional sets:
#        CHX-enriched  — log2FC > +0.5, padj < 0.05  (proteins ENRICHED
#                        upon cycloheximide treatment)
#        CHX-depleted  — log2FC < -0.5, padj < 0.05  (proteins DEPLETED
#                        upon cycloheximide treatment)
#   3. Builds a STRING physical-interaction network for each set using
#      Lydia's direct-PPI method:
#        - Read local physical links file
#        - Map gene symbols to STRING IDs via STRINGdb$map()
#        - Seeds + direct neighbors (combined_score > 250)
#        - Strong secondary links (combined_score > 700)
#        - igraph force-directed layout (layout_with_fr)
#        - Node color by up/down direction; label top nodes by degree
#   4. Runs bidirectional GO enrichment (enrichGO) on the enriched and
#      depleted sets independently — BP, MF, CC with dotplot + barplot
#   5. Exports all gene lists to output/tables/
#
# COLOR PALETTE (project-wide Okabe-Ito):
#   CHX-enriched (up)   → Vermillion #D55E00
#   CHX-depleted (down) → Navy      #0072B2
#   Not significant     → Grey      #D0D0D0
#
# Usage:
#   make chx-common
###############################################################################

source('R/01_config.R', chdir = TRUE)
source('R/utils.R')

cat("\n=========================================\n")
cat(" CHX/DMSO Common Protein Analysis\n")
cat("=========================================\n\n")

# ---- Load all experiments (CHX/DMSO CSVs are auto-discovered) ----
experiments <- load_all_experiments()

# ---- Locate the key CHX-vs-DMSO comparison ----
CHX_EXP_NAME <- "TRIP4_CHX_vs_TRIP4_DMSO"
chx_df <- find_experiment(experiments, CHX_EXP_NAME)

if (is.null(chx_df)) {
  cat("ERROR:", CHX_EXP_NAME, "not found among loaded experiments.\n")
  cat("  Available:", paste(names(experiments), collapse = ", "), "\n")
  quit(status = 1)
}

cat(sprintf("  Loaded %s: %d proteins\n\n", CHX_EXP_NAME, nrow(chx_df)))

# ---- Directional significance splits ----
# Enriched = significantly MORE abundant with CHX (positive fold change)
# Depleted = significantly LESS abundant with CHX (negative fold change)
chx_df$direction <- "nonsig"
chx_df$direction[chx_df$padj < P_VALUE_CUTOFF &
                   chx_df$log2FC >  LOG2FC_CUTOFF] <- "enriched"
chx_df$direction[chx_df$padj < P_VALUE_CUTOFF &
                   chx_df$log2FC < -LOG2FC_CUTOFF] <- "depleted"

enriched_df <- chx_df[chx_df$direction == "enriched" &
                        !is.na(chx_df$gene), ]
depleted_df <- chx_df[chx_df$direction == "depleted" &
                        !is.na(chx_df$gene), ]

enriched_genes <- unique(enriched_df$gene)
depleted_genes <- unique(depleted_df$gene)

cat(sprintf("  CHX-enriched proteins:  %d\n", length(enriched_genes)))
cat(sprintf("  CHX-depleted proteins:  %d\n", length(depleted_genes)))

# ---- Export gene lists to output/tables/ ----
save_table(data.frame(gene = enriched_genes,
                      log2FC = enriched_df$log2FC[match(enriched_genes,
                                                        enriched_df$gene)],
                      padj = enriched_df$padj[match(enriched_genes,
                                                    enriched_df$gene)],
                      category = "CHX_enriched"),
           "chx_common_enriched_genes")

save_table(data.frame(gene = depleted_genes,
                      log2FC = depleted_df$log2FC[match(depleted_genes,
                                                        depleted_df$gene)],
                      padj = depleted_df$padj[match(depleted_genes,
                                                    depleted_df$gene)],
                      category = "CHX_depleted"),
           "chx_common_depleted_genes")

# Also export plain text lists (one gene per line) for external tools
write.table(data.frame(gene = enriched_genes),
            file = file.path(TABLE_DIR, "chx_enriched_genes.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(data.frame(gene = depleted_genes),
            file = file.path(TABLE_DIR, "chx_depleted_genes.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
cat("  Exported gene lists to output/tables/\n\n")

# =====================================================================
# STRING NETWORK ANALYSIS (Lydia's direct-PPI method)
# =====================================================================
# For each directional set:
#   1. Map genes to STRING IDs
#   2. Find physical interactions among mapped proteins
#   3. Expand: direct neighbors (score > 250) + strong secondary (score > 700)
#   4. Build igraph, layout with Fruchterman-Reingold
#   5. Color nodes by up/down; label highest-degree nodes

library(STRINGdb)
library(igraph)

# ---- Read local STRING physical links file ----
phys_file <- file.path(DATA_DIR, "9606.protein.physical.links.v12.0.txt")

if (!file.exists(phys_file)) {
  cat("WARNING: STRING physical links file not found:\n  ", phys_file, "\n")
  cat("  Download from STRING and place in data/. Skipping network analysis.\n")
} else {
  cat("--- Loading STRING physical interactions ---\n")
  phys <- read.table(phys_file, header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded %d physical interactions\n", nrow(phys)))

  # ---- Initialize STRINGdb for ID mapping ----
  string_cache <- file.path(OUTPUT_DIR, "string_cache")
  dir.create(string_cache, recursive = TRUE, showWarnings = FALSE)

  string_db <- STRINGdb$new(
    version         = STRING_VERSION,
    species         = STRING_TAXON,
    score_threshold = STRING_SCORE_THRESHOLD,
    input_directory = string_cache
  )

  # -----------------------------------------------------------------
  # Helper: build a STRING network plot for one gene set
  # -----------------------------------------------------------------
  # Adapted from Lydia's method (see 21_lydia_network_volcano.R).
  # All seed genes are the significant proteins in the set; we expand
  # to direct neighbors (score > 250) and strong secondary links
  # (score > 700 among neighbors of seeds).
  #
  # Args:
  #   gene_df      data frame with columns: gene, log2FC, padj
  #   set_label    short label for filenames/titles (e.g. "enriched")
  #   set_title    human-readable title
  #   node_color   single color for all nodes in this set (all same direction)
  #
  # Returns invisibly: list with graph, degree-ranked candidate table.
  build_chx_string_network <- function(gene_df, set_label, set_title,
                                        node_color) {

    genes <- unique(gene_df$gene)
    cat(sprintf("\n  [%s] %d proteins for STRING mapping\n", set_label,
                length(genes)))

    if (length(genes) < 3) {
      cat("    Skipped: fewer than 3 proteins\n")
      return(invisible(NULL))
    }

    # ---- Map gene symbols to STRING IDs ----
    # as.data.frame() before $map() avoids "incorrect dimensions" with tibbles
    map_input <- as.data.frame(gene_df[, c("gene", "log2FC", "padj")])
    map_input <- map_input[!is.na(map_input$gene) & map_input$gene != "", ]

    mapped <- tryCatch({
      string_db$map(map_input, "gene", removeUnmappedRows = TRUE)
    }, error = function(e) {
      cat(sprintf("    ERROR mapping: %s\n", conditionMessage(e)))
      return(NULL)
    })

    if (is.null(mapped) || nrow(mapped) == 0) {
      cat("    No proteins mapped to STRING\n")
      return(invisible(NULL))
    }
    cat(sprintf("    Mapped %d / %d proteins to STRING\n", nrow(mapped),
                length(genes)))

    seed_ids <- unique(mapped$STRING_id)
    cat(sprintf("    Seed STRING IDs: %d\n", length(seed_ids)))

    # ---- Direct interactions touching seeds (score > 250) ----
    direct_edges <- phys[
      (phys$protein1 %in% seed_ids | phys$protein2 %in% seed_ids) &
      phys$combined_score > 250, ]

    if (nrow(direct_edges) == 0) {
      cat("    No direct interactions found. Skipping.\n")
      return(invisible(NULL))
    }
    cat(sprintf("    Direct edges (score > 250): %d\n", nrow(direct_edges)))

    # IDs directly connected to seeds (excluding seeds themselves)
    direct_partners <- unique(c(
      direct_edges$protein1[!(direct_edges$protein1 %in% seed_ids)],
      direct_edges$protein2[!(direct_edges$protein2 %in% seed_ids)]
    ))

    # ---- Strong secondary links (score > 700 among neighbors) ----
    if (length(direct_partners) > 0) {
      secondary_edges <- phys[
        phys$combined_score > 700 &
        (phys$protein1 %in% direct_partners |
           phys$protein2 %in% direct_partners), ]
    } else {
      secondary_edges <- direct_edges[0, ]
    }
    cat(sprintf("    Strong secondary edges (score > 700): %d\n",
                nrow(secondary_edges)))

    # ---- Combine all edges ----
    # Normalize column names to from/to for igraph
    all_edges <- rbind(
      data.frame(from = direct_edges$protein1,
                 to   = direct_edges$protein2,
                 score = direct_edges$combined_score,
                 stringsAsFactors = FALSE),
      data.frame(from = secondary_edges$protein1,
                 to   = secondary_edges$protein2,
                 score = secondary_edges$combined_score,
                 stringsAsFactors = FALSE)
    )
    # Deduplicate
    all_edges <- all_edges[!duplicated(paste(pmin(all_edges$from, all_edges$to),
                                             pmax(all_edges$from, all_edges$to))), ]

    all_ids <- unique(c(all_edges$from, all_edges$to))
    cat(sprintf("    Network nodes: %d, edges: %d\n",
                length(all_ids), nrow(all_edges)))

    # ---- Build igraph ----
    g <- graph_from_data_frame(all_edges[, c("from", "to")], directed = FALSE)

    # Map STRING IDs back to gene names + fold-change data
    id_to_gene <- setNames(mapped$gene, mapped$STRING_id)
    V(g)$name_display <- sapply(V(g)$name, function(id) {
      gn <- id_to_gene[id]
      if (is.na(gn)) "" else gn
    })

    # Attach log2FC/padj for nodes that are in our seed set
    V(g)$log2FC <- mapped$log2FC[match(V(g)$name, mapped$STRING_id)]
    V(g)$padj   <- mapped$padj[match(V(g)$name, mapped$STRING_id)]
    V(g)$is_seed <- V(g)$name %in% seed_ids

    # ---- Node appearance: color by up/down, size by seed status ----
    # All seeds in a given set share the same direction (enriched=up,
    # depleted=down), so we use the set color for seeds. Neighbors that
    # are also in our data get a lighter shade; pure-STRING neighbors
    # (not measured) are grey.
    V(g)$color <- ifelse(V(g)$is_seed,
                         node_color,                          # set color
                         ifelse(!is.na(V(g)$log2FC),
                                grDevices::adjustcolor(node_color, 0.55),
                                "#D0D0D0"))                    # unmeasured
    V(g)$size  <- ifelse(V(g)$is_seed, 7, 3.5)

    # ---- Label top nodes by degree ----
    node_deg <- degree(g)
    # Label seeds with high degree + a few top neighbors
    label_threshold <- max(1, quantile(node_deg, 0.80, na.rm = TRUE))
    V(g)$label <- ifelse(node_deg >= label_threshold | V(g)$is_seed,
                         V(g)$name_display, NA)

    # ---- Plot: force-directed layout (Fruchterman-Reingold) ----
    set.seed(42)
    layout_fr <- layout_with_fr(g)

    commit_hash <- get_git_hash()
    base_name <- paste0("chx_common_string_", set_label)
    safe_name <- sanitize_filename(base_name)
    versioned  <- paste0(safe_name, "_", commit_hash)
    png_path <- safe_filepath(FIGURE_DIR, versioned, ".png")
    pdf_path <- safe_filepath(FIGURE_DIR, versioned, ".pdf")

    plot_title <- paste0(set_title, " — STRING Physical Network")

    grDevices::png(png_path, width = 12, height = 10, units = "in",
                   res = FIG_DPI)
    plot(g,
         layout             = layout_fr,
         vertex.label       = V(g)$label,
         vertex.label.cex   = 0.65,
         vertex.label.font  = 2,
         vertex.label.color = "black",
         vertex.frame.color = NA,
         edge.color         = "grey80",
         edge.width         = 0.5,
         edge.curved        = 0.1,
         main               = plot_title)
    legend("topright",
           legend = c("Seed (significant)", "Neighbor (measured)",
                      "Neighbor (STRING only)"),
           col    = c(node_color,
                      grDevices::adjustcolor(node_color, 0.55),
                      "#D0D0D0"),
           pch = 19, pt.cex = c(2, 1.5, 1.5), cex = 0.85)
    grDevices::dev.off()

    grDevices::pdf(pdf_path, width = 12, height = 10)
    plot(g,
         layout             = layout_fr,
         vertex.label       = V(g)$label,
         vertex.label.cex   = 0.65,
         vertex.label.font  = 2,
         vertex.label.color = "black",
         vertex.frame.color = NA,
         edge.color         = "grey80",
         edge.width         = 0.5,
         edge.curved        = 0.1,
         main               = plot_title)
    legend("topright",
           legend = c("Seed (significant)", "Neighbor (measured)",
                      "Neighbor (STRING only)"),
           col    = c(node_color,
                      grDevices::adjustcolor(node_color, 0.55),
                      "#D0D0D0"),
           pch = 19, pt.cex = c(2, 1.5, 1.5), cex = 0.85)
    grDevices::dev.off()

    cat(sprintf("    Saved: %s\n", basename(png_path)))
    cat(sprintf("    Saved: %s\n", basename(pdf_path)))

    # ---- Degree-ranked candidate table ----
    degree_df <- data.frame(
      gene        = V(g)$name_display,
      STRING_id   = V(g)$name,
      is_seed     = V(g)$is_seed,
      connections = degree(g),
      log2FC      = V(g)$log2FC,
      padj        = V(g)$padj,
      stringsAsFactors = FALSE
    )
    degree_df <- degree_df[degree_df$gene != "", ]
    degree_df <- degree_df[order(-degree_df$connections), ]
    save_table(degree_df, paste0("chx_common_string_", set_label, "_candidates"))

    cat(sprintf("    Top 10 by degree:\n"))
    for (i in seq_len(min(10, nrow(degree_df)))) {
      row <- degree_df[i, ]
      cat(sprintf("      %2d. %-12s (deg=%d, log2FC=%s, %s)\n",
                  i, row$gene, row$connections,
                  ifelse(is.na(row$log2FC), "NA",
                         sprintf("%.2f", row$log2FC)),
                  ifelse(row$is_seed, "SEED", "neighbor")))
    }

    invisible(list(graph = g, candidates = degree_df))
  }

  # ---- Build network for CHX-enriched proteins (Vermillion) ----
  cat("\n--- CHX-Enriched STRING Network ---\n")
  if (nrow(enriched_df) >= 3) {
    build_chx_string_network(
      gene_df    = enriched_df,
      set_label  = "enriched",
      set_title  = "CHX-Enriched Proteins (TRIP4 CHX vs DMSO)",
      node_color = GLOBAL_COLORS[["enriched_up"]]   # Vermillion #D55E00
    )
  } else {
    cat("  Skipped enriched network: fewer than 3 proteins\n")
  }

  # ---- Build network for CHX-depleted proteins (Navy) ----
  cat("\n--- CHX-Depleted STRING Network ---\n")
  if (nrow(depleted_df) >= 3) {
    build_chx_string_network(
      gene_df    = depleted_df,
      set_label  = "depleted",
      set_title  = "CHX-Depleted Proteins (TRIP4 CHX vs DMSO)",
      node_color = "#0072B2"   # Navy for depleted/down
    )
  } else {
    cat("  Skipped depleted network: fewer than 3 proteins\n")
  }
}

# =====================================================================
# BIDIRECTIONAL GO ENRICHMENT
# =====================================================================
# Run enrichGO independently on the enriched and depleted gene sets.
# This reveals biological processes UP-regulated by CHX vs processes
# DOWN-regulated by CHX.

cat("\n=========================================\n")
cat(" Bidirectional GO Enrichment\n")
cat("=========================================\n\n")

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

# Background universe = all proteins measured across experiments
universe <- unique(unlist(lapply(experiments, function(df) df$gene)))
cat(sprintf("Background universe: %d unique proteins\n\n", length(universe)))

ONT_LABELS <- c(
  "BP" = "Biological Process",
  "MF" = "Molecular Function",
  "CC" = "Cellular Component"
)

# ---- Helper: run enrichGO for one gene set, all three ontologies ----
run_bidirectional_go <- function(genes, set_label, set_title, universe) {
  cat(sprintf("\n--- %s (%d genes) ---\n", set_title, length(genes)))

  if (length(genes) < 5) {
    cat("  Skipped: fewer than 5 genes\n")
    return(invisible(NULL))
  }

  for (ont in c("BP", "MF", "CC")) {
    cat(sprintf("  [%s] %s — Running enrichGO...\n", ont, ONT_LABELS[ont]))

    result <- tryCatch({
      enrichGO(
        gene          = genes,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = GO_PVALUE_CUTOFF,
        qvalueCutoff  = GO_QVALUE_CUTOFF,
        universe      = universe
      )
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      return(NULL)
    })

    if (is.null(result) || nrow(as.data.frame(result)) == 0) {
      cat("    No enriched terms found\n")
      next
    }

    res_df <- as.data.frame(result)
    cat(sprintf("    Found %d enriched GO terms\n", nrow(res_df)))

    # Simplify redundant terms via GO graph structure
    result <- tryCatch(
      simplify(result, cutoff = 0.7, by = "p.adjust", select_fun = min),
      error = function(e) result
    )
    res_df <- as.data.frame(result)
    cat(sprintf("    After simplify: %d terms\n", nrow(res_df)))

    prefix <- sanitize_filename(paste0("GO_chx_common_", set_label, "_", ont))
    save_table(res_df, prefix)

    # ---- Dotplot ----
    n_show  <- min(20, nrow(res_df))
    fig_h   <- max(7, n_show * 0.4)
    go_title <- paste0("GO: ", set_title, " — ", ONT_LABELS[ont])

    p_dot <- dotplot(result, showCategory = n_show, title = go_title) +
      ggplot2::labs(x = "Gene Ratio", color = "p-adjusted value") +
      ggplot2::theme(
        plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 7)
      )
    save_figure(p_dot, paste0(prefix, "_dotplot"), width = 10, height = fig_h)

    # ---- Barplot (sorted by Count) ----
    n_bar  <- min(15, nrow(res_df))
    fig_hb <- max(7, n_bar * 0.5)

    p_bar <- barplot(result, showCategory = n_bar, orderBy = "Count",
                     title = go_title) +
      ggplot2::labs(color = "p-adjusted value", fill = "p-adjusted value") +
      ggplot2::theme(
        plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 7)
      )
    save_figure(p_bar, paste0(prefix, "_barplot"), width = 10, height = fig_hb)
  }
}

# ---- GO on CHX-enriched set (proteins more abundant with CHX) ----
run_bidirectional_go(enriched_genes, "enriched",
                     "CHX-Enriched Proteins", universe)

# ---- GO on CHX-depleted set (proteins less abundant with CHX) ----
run_bidirectional_go(depleted_genes, "depleted",
                     "CHX-Depleted Proteins", universe)

cat("\n=========================================\n")
cat(" CHX/DMSO common analysis complete!\n")
cat("=========================================\n")
cat("\nOutputs:\n")
cat("  - Gene lists: output/tables/chx_common_*_genes.{csv,txt}\n")
cat("  - STRING networks: output/figures/chx_common_string_{enriched,depleted}_*\n")
cat("  - Candidate tables: output/tables/chx_common_string_*_candidates_*.csv\n")
cat("  - GO enrichment: output/figures/GO_chx_common_* + output/tables/GO_chx_common_*\n")
cat("\nRun: make open-chx-common\n")
