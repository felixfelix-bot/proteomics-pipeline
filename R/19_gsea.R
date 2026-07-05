###############################################################################
# 19_gsea.R
# GSEA enrichment analysis with STRING network visualization.
#
# TWO PARTS:
#   PART A — GSEA on ALL experiments (ranked gene lists, all BP/MF/CC)
#   PART B — STRING protein-protein interaction network map for each experiment
#
# The user specifically requested network MAP visualizations (node-and-edge
# graphs like the STRING web interface shows) rather than just dotplots.
#
# Output per experiment:
#   - GSEA results table (CSV with NES, p-value, leading edge)
#   - STRING network map (PNG + PDF) — proteins as nodes, interactions as edges
#   - Ridgeplot (running enrichment score)
#
# Usage:
#   make gsea
###############################################################################
cat("\n=========================================\n")
cat(" GSEA Enrichment + STRING Network Maps\n")
cat(" (all experiments — Lydia style)\n")
cat("=========================================\n\n")

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

# Load all mass spec experiments from the data/ directory
experiments <- load_all_experiments()
cat(sprintf("Loaded %d experiments\n\n", length(experiments)))

ONT_LABELS <- c(
  "BP" = "Biological Process",
  "MF" = "Molecular Function",
  "CC" = "Cellular Component"
)

# -------------------------------------------------------------------
# Helper A: Run GSEA on one experiment
# -------------------------------------------------------------------
run_gsea <- function(df, exp_name, display_name) {
  cat(sprintf("\n======== %s (%d proteins) ========\n", display_name, nrow(df)))

  # Build ranked gene list: sign(logFC) * -log10(adj.P.Val)
  # This captures both direction (sign of logFC) AND significance (padj)
  df$rank_metric <- sign(df$logFC) * -log10(df$adj.P.Val)
  df <- df[is.finite(df$rank_metric), ]
  df <- df[order(df$rank_metric, decreasing = TRUE), ]

  gene_list <- df$rank_metric
  names(gene_list) <- df$gene
  cat(sprintf("  Ranked %d genes (%.2f to %.2f)\n",
              length(gene_list), min(gene_list), max(gene_list)))

  safe_exp <- sanitize_filename(gsub("[ _]", "_", display_name))

  for (ont in c("BP", "MF", "CC")) {
    cat(sprintf("\n  [%s] %s\n", ont, ONT_LABELS[ont]))

    result <- tryCatch({
      gseGO(
        geneList      = gene_list,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = ont,
        minGSSize     = 10,
        maxGSSize     = 500,
        pvalueCutoff  = 0.05,
        pAdjustMethod = "BH",
        verbose       = FALSE
      )
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      return(NULL)
    })

    if (is.null(result) || nrow(as.data.frame(result)) == 0) {
      cat("    No enriched gene sets\n")
      next
    }

    res_df <- as.data.frame(result)
    cat(sprintf("    Found %d enriched gene sets\n", nrow(res_df)))

    prefix <- sanitize_filename(paste0("GSEA_", safe_exp, "_", ont))
    save_table(res_df, prefix)

    # Ridgeplot — shows enrichment score distribution across ranked list
    n_show <- min(20, nrow(res_df))
    fig_height_ridge <- max(7, n_show * 0.45)
    ridge_title <- paste0("GSEA: ", display_name, " — ", ONT_LABELS[ont])

    p_ridge <- ridgeplot(result, showCategory = n_show,
                         label_format = 50) +
      ggplot2::labs(
        title = ridge_title,
        x = "Running Enrichment Score",
        fill = "p-adjusted value"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 7)
      )
    save_figure(p_ridge, paste0(prefix, "_ridgeplot"),
                width = 12, height = fig_height_ridge)
  }
  cat("\n")
}

# -------------------------------------------------------------------
# Helper B: Build STRING network map for significant proteins
# -------------------------------------------------------------------
build_string_network <- function(df, display_name) {
  # Get significant genes
  sig_genes <- get_significant_genes(df)
  n_sig <- length(sig_genes)
  cat(sprintf("  Significant: %d genes — ", n_sig))

  if (n_sig < 5) {
    cat("too few for network\n")
    return(NULL)
  }
  if (n_sig > 200) {
    # Cap at 200 to keep network readable and runtime reasonable
    cat(sprintf("capping at 200 (from %d)\n", n_sig))
    # Take top 200 by significance: sort by -log10(padj) descending
    df_sig <- df[df$gene %in% sig_genes, ]
    df_sig$score <- -log10(df_sig$padj) * abs(df_sig$log2FC)
    df_sig <- df_sig[order(df_sig$score, decreasing = TRUE), ]
    sig_genes <- unique(df_sig$gene)[1:200]
    n_sig <- length(sig_genes)
  } else {
    cat("OK\n")
  }

  # Subset the data frame to significant genes only
  sig_df <- df[df$gene %in% sig_genes, ]
  sig_df <- sig_df[!duplicated(sig_df$gene), ]

  # ---- Initialize STRING ----
  cat("  Querying STRING database...\n")
  string_cache_dir <- file.path(OUTPUT_DIR, "string_cache")
  dir.create(string_cache_dir, showWarnings = FALSE, recursive = TRUE)

  string_db <- tryCatch({
    STRINGdb::STRINGdb$new(
      version = STRING_VERSION,
      species = STRING_TAXON,
      score_threshold = STRING_SCORE_THRESHOLD,
      input_directory = string_cache_dir
    )
  }, error = function(e) {
    cat(sprintf("    STRING init error: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(string_db)) return(NULL)

  # Map to STRING
  mapped <- string_db$map(sig_df, "gene", removeUnmappedRows = TRUE)
  n_mapped <- sum(!is.na(mapped$STRING_id))
  cat(sprintf("    Mapped %d / %d to STRING\n", n_mapped, nrow(sig_df)))

  if (n_mapped < 3) {
    cat("    Too few mapped genes\n")
    return(NULL)
  }

  # Get interactions
  mapped_ids <- mapped$STRING_id[!is.na(mapped$STRING_id)]
  interactions <- string_db$get_interactions(mapped_ids)

  if (is.null(interactions) || nrow(interactions) == 0) {
    cat("    No interactions found\n")
    return(NULL)
  }
  cat(sprintf("    Found %d interactions\n", nrow(interactions)))

  # ---- Build igraph network ----
  library(igraph)
  edges <- data.frame(
    from = interactions$from,
    to   = interactions$to,
    weight = interactions$combined_score / 1000
  )
  g <- graph_from_data_frame(edges, directed = FALSE)

  # Label with gene symbols
  string_to_gene <- setNames(mapped$gene, mapped$STRING_id)
  V(g)$name_display <- sapply(V(g)$name, function(id) {
    gene <- string_to_gene[id]
    if (is.na(gene)) id else gene
  })

  # Classify: seeds (most significant) vs candidates
  # Use stricter Lydia criteria for seeds
  mapped$is_seed <- (mapped$padj < 0.000001 & abs(mapped$log2FC) > 2) |
                    (abs(mapped$log2FC) > 7)

  seed_string_ids <- mapped$STRING_id[mapped$is_seed &
                                        !is.na(mapped$STRING_id)] %||% character(0)

  # Apply to graph vertices
  V(g)$is_seed <- V(g)$name %in% seed_string_ids
  V(g)$color <- ifelse(V(g)$is_seed,
                       GLOBAL_COLORS[["enriched_up"]],  # Orange
                       GLOBAL_COLORS[["known_ia"]])      # Green
  V(g)$size <- ifelse(V(g)$is_seed, 10, 5)

  # Label seeds (not candidates — too many)
  label_vec <- ifelse(V(g)$is_seed, V(g)$name_display, NA)

  # ---- Render network ----
  set.seed(42)
  layout_fr <- layout_with_fr(g)
  commit_hash <- get_git_hash()
  safe_name <- sanitize_filename(gsub("[ _]", "_", display_name))

  png_path <- safe_filepath(FIGURE_DIR,
    paste0("network_", safe_name, "_", commit_hash), ".png")
  pdf_path <- safe_filepath(FIGURE_DIR,
    paste0("network_", safe_name, "_", commit_hash), ".pdf")

  title_text <- paste0(display_name, " — STRING Interaction Network",
                       " (", vcount(g), " nodes, ", ecount(g), " edges)")

  grDevices::png(png_path, width = 14, height = 12, units = "in", res = 300)
  plot(g,
       layout = layout_fr,
       vertex.label = label_vec,
       vertex.label.cex = 0.8,
       vertex.label.font = 2,
       vertex.label.color = "black",
       vertex.frame.color = "white",
       edge.color = "grey70",
       edge.width = 0.5,
       edge.curved = 0.1,
       main = title_text)
  legend("topright",
         legend = c("Seeds (top hits)", "Candidates"),
         col = c(GLOBAL_COLORS[["enriched_up"]],
                 GLOBAL_COLORS[["known_ia"]]),
         pch = 19, pt.cex = c(2.5, 1.5), cex = 1)
  grDevices::dev.off()

  grDevices::pdf(pdf_path, width = 14, height = 12)
  plot(g,
       layout = layout_fr,
       vertex.label = label_vec,
       vertex.label.cex = 0.8,
       vertex.label.font = 2,
       vertex.label.color = "black",
       vertex.frame.color = "white",
       edge.color = "grey70",
       edge.width = 0.5,
       edge.curved = 0.1,
       main = title_text)
  legend("topright",
         legend = c("Seeds (top hits)", "Candidates"),
         col = c(GLOBAL_COLORS[["enriched_up"]],
                 GLOBAL_COLORS[["known_ia"]]),
         pch = 19, pt.cex = c(2.5, 1.5), cex = 1)
  grDevices::dev.off()

  cat(sprintf("  Network map saved: %s\n", basename(png_path)))

  # Save candidate table
  degree_df <- data.frame(
    gene = V(g)$name_display,
    is_seed = V(g)$is_seed,
    connections = degree(g),
    stringsAsFactors = FALSE
  )
  degree_df <- degree_df[order(-degree_df$connections), ]
  degree_df <- merge(degree_df, sig_df[, c("gene", "log2FC", "padj")],
                     by = "gene", all.x = TRUE)
  degree_df <- degree_df[order(-degree_df$connections), ]
  save_table(degree_df, paste0("network_", safe_name, "_candidates"))
}

# ====================================================================
# RUN: Iterate over ALL experiments
# ====================================================================
# Define the experiments the user cares about (all TurboID, FlagIP, CHX, CRAC)
TARGET_NAMES <- c(
  # TurboID core
  "BK467_TRIP4_vs_BK467_WT",
  "BK467_TRIP4_RA02_vs_BK467_WT",
  "BK467_TRIP4_RA02_vs_BK467_TRIP4",
  "BK467_TRIP4_vs_BK504_TRIP4",
  "BK467_TRIP4_RA02_vs_BK504_TRIP4_RA04",
  # TurboID BK504
  "BK504_TRIP4_RA04_vs_BK504_TRIP4",
  "BK504_TRIP4_RA04_vs_BK467_WT",
  # Flag IP
  "BK516_Cflag_vs_BK516_Ctrl",
  "BK516_Nflag_vs_BK516_Ctrl",
  "BK516_Cflag_vs_BK516_Nflag",
  # Flag IP + RA
  "BK523_Cflag_RA04_vs_BK516_Cflag",
  "BK523_Cflag_RA04_vs_BK523_Ctrl_RA04",
  "BK523_Cflag_RA04_vs_BK523_Nflag_RA04",
  "BK523_Nflag_RA04_vs_BK516_Nflag",
  "BK523_Nflag_RA04_vs_BK523_Ctrl_RA04",
  # CHX/DMSO
  "TRIP4_CHX_vs_TRIP4_DMSO",
  "TRIP4_CHX_vs_WT",
  "TRIP4_DMSO_vs_WT"
)

cat(sprintf("Part A — GSEA Enrichment on all experiments\n"))
cat(sprintf("==========================================\n"))
n_processed <- 0
for (exp_name in TARGET_NAMES) {
  if (exp_name %in% names(experiments)) {
    # Build a display name without underscores
    display <- gsub("_vs_", " vs ", exp_name)
    display <- gsub("_", " ", display)
    run_gsea(experiments[[exp_name]], exp_name, display)
    n_processed <- n_processed + 1
  } else if (grepl("Cflag_vs_BK516_Ctrl", exp_name)) {
    # Handle duplicate names (__1, __2)
    for (suffix in c("__1", "__2")) {
      dup_name <- paste0(exp_name, suffix)
      if (dup_name %in% names(experiments)) {
        display <- gsub("_vs_", " vs ", dup_name)
        display <- gsub("_", " ", display)
        run_gsea(experiments[[dup_name]], dup_name, display)
        n_processed <- n_processed + 1
      }
    }
  }
}
cat(sprintf("\nGSEA complete: %d experiments processed\n\n", n_processed))

# ---- CRAC data (separate loading, different file format) ----
cat(sprintf("Part B — CRAC RNA Interactome\n"))
cat(sprintf("=============================\n"))
crac_file <- file.path(DATA_DIR, "FLAG-TRIP4_list_CRACdata.csv")
if (file.exists(crac_file)) {
  crac_df <- readr::read_csv(crac_file, show_col_types = FALSE)
  crac_df$gene <- crac_df[[CRAC_GENE_COL]]
  crac_df$logFC <- crac_df[[CRAC_LOG2FC_COL]]
  crac_df$padj <- crac_df[[CRAC_PADJ_COL]]
  crac_df <- crac_df[!is.na(crac_df$gene) & !is.na(crac_df$logFC) &
                      !is.na(crac_df$padj), ]
  run_gsea(crac_df, "CRAC", "CRAC RNA Interactome")
  cat(sprintf("\n"))
} else {
  cat("  CRAC file not found\n\n")
}

# ====================================================================
# PART C: STRING Network Maps for ALL experiments
# ====================================================================
cat(sprintf("Part C — STRING Network Maps\n"))
cat(sprintf("============================\n"))

library(STRINGdb)

for (exp_name in names(experiments)) {
  # Skip duplicates (only process main copy)
  if (grepl("__[0-9]+$", exp_name)) next

  # Build display name
  display <- gsub("_vs_", " vs ", exp_name)
  display <- gsub("_", " ", display)

  cat(sprintf("\n--- %s ---\n", display))
  build_string_network(experiments[[exp_name]], display)
}

# ---- CRAC network (separate loading) ----
if (file.exists(crac_file)) {
  cat(sprintf("\n--- CRAC RNA Interactome ---\n"))
  build_string_network(crac_df, "CRAC RNA Interactome")
}

cat(sprintf("\n=========================================\n"))
cat(sprintf(" GSEA + STRING network analysis complete!\n"))
cat(sprintf("=========================================\n"))
