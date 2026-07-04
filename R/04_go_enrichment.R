##############################################################################
# 04_go_enrichment.R
#
# WHAT THIS SCRIPT DOES (the big picture, for non-R users):
# --------------------
# This script answers a key biological question:
#   "Of the proteins that bind TRIP4, are certain biological functions or
#    pathways OVER-REPRESENTED compared to what we'd expect by chance?"
#
# This is called GO ENRICHMENT ANALYSIS (Gene Ontology Enrichment).
#
# What is Gene Ontology (GO)?
# ---------------------------
# GO is a standardized vocabulary for describing what genes/proteins do.
# Every gene/protein can be tagged with GO terms in 3 categories:
#   - BP (Biological Process): what the protein DOES (e.g. "DNA repair")
#   - MF (Molecular Function): what it IS (e.g. "kinase activity")
#   - CC (Cellular Component): WHERE it is (e.g. "nucleus", "ribosome")
#
# The enrichment question:
#   If 5000 proteins were detected in total (the "universe" / background),
#   and 50 of them are significant binders, and 10 of those 50 are involved
#   in "mRNA splicing" — is that more than expected by chance?
#   (Statistical test answers: yes/no, with a p-value.)
#
# Two methods are used here:
#
#   1. ORA (Over-Representation Analysis):
#      Takes a binary list (significant or not) and asks "which GO terms
#      appear more often in the significant group than expected?"
#      Uses a Fisher's exact test (hypergeometric distribution).
#      Like asking: "Is 'mRNA splicing' over-represented among my hits?"
#
#   2. GSEA (Gene Set Enrichment Analysis):
#      Uses a RANKED list of ALL proteins (sorted by significance × fold change)
#      and asks "do proteins from a given GO term cluster at the top or bottom
#      of the ranking?" It does NOT require a hard cutoff.
#      More sensitive than ORA because it uses the full gradient of evidence.
#
# This script also uses rrvgo to REDUCE redundant GO terms:
#   Enrichment often returns many similar terms (e.g. "mRNA splicing",
#   "RNA splicing", "mRNA processing" all overlap). rrvgo clusters
#   similar terms and picks one representative, making results readable.
#
# THIS IS THE MOST COMPUTATIONALLY INTENSIVE STEP in the pipeline.
# It may take several minutes because it tests thousands of GO terms
# against every experiment × every ontology.
#
# Usage:
#   source("R/01_config.R")
#   source("R/utils.R")
#   source("R/04_go_enrichment.R")
##############################################################################

# Print a header banner to the console.
cat("\n=========================================\n")
cat(" GO Enrichment Analysis\n")
cat("=========================================\n\n")

# ---- Load annotation database ----
# These library() calls load R packages (like "import" in Python).
#
# org.Hs.eg.db:
#   A database package containing gene/protein annotations for HUMANS
#   (Hs = Homo sapiens). It maps gene symbols to GO terms, etc.
#   This is our lookup table that tells us which GO terms each gene belongs to.
#
# clusterProfiler:
#   THE main R package for enrichment analysis. Provides enrichGO() for ORA
#   and gseGO() for GSEA. Highly cited in bioinformatics papers.
#
# enrichplot:
#   Provides visualization functions for enrichment results: dotplot(),
#   barplot(), cnetplot(), ridgeplot(), etc.
cat("Loading annotation database (org.Hs.eg.db)...\n")
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

# ---- Load all experiment data ----
# Same as in script 03: returns a named list of data frames, one per experiment.
experiments <- load_all_experiments()

# ---- Define universe (all detected proteins across experiments) ----
# THE "UNIVERSE" CONCEPT — critical for correct enrichment statistics:
# -------------------------------------------------------------------
# The universe is the BACKGROUND set of all proteins that COULD have been found.
# Enrichment is always RELATIVE to a background:
#   "Is GO term X more common in my significant list than in the universe?"
#
# Using the right universe matters enormously:
#   - If you use ALL ~20000 human genes as the universe, but your mass spec
#     only detects 3000, then proteins your instrument can't see will skew
#     the statistics (you'll get false enrichment signals).
#   - The correct universe is the set of proteins actually DETECTED/quantified
#     in your experiment — this is what we use here.
#
# lapply() iterates over experiments (like map in Python).
# Each df is a data frame; df$gene is the column of gene names.
# unlist() flattens all lists into one big vector.
# unique() removes duplicates (a protein detected in 3 experiments counts once).
universe <- unique(unlist(lapply(experiments, function(df) df$gene)))
cat(sprintf("Background universe: %d unique proteins\n", length(universe)))

# ---- Extract significant gene sets ----
# For each experiment, get the list of significant gene/protein names.
# lapply() applies get_significant_genes() to each data frame.
# Returns a named list of character vectors (gene name lists).
gene_sets <- lapply(experiments, get_significant_genes)

# =====================================================================
# Helper: Run ORA for one gene set and one ontology
# =====================================================================
# This is a FUNCTION DEFINITION in R — like "def run_ora(...):" in Python.
# It takes 4 arguments and returns an enrichment result object (or NULL).
#
# Arguments:
#   genes           = vector of significant gene SYMBOLS (e.g. "TP53", "MYC")
#   ontology        = one of "BP", "MF", "CC"
#   universe        = the background gene set (all detected proteins)
#   experiment_name = name of the experiment (for logging/labels)
run_ora <- function(genes, ontology, universe, experiment_name) {

  # WHY WE NEED >= 5 SIGNIFICANT GENES:
  # ----------------------------------
  # With fewer than 5 genes, enrichment statistics become meaningless:
  #   - The hypergeometric test has almost no power with tiny numbers.
  #   - You'd get wildly unstable p-values that don't reproduce.
  #   - A single protein belonging to a GO term could look "enriched" by chance.
  # 5 is a common minimum threshold in bioinformatics.
  #
  # If we have fewer than 5, skip this combination and return NULL (nothing).
  if (length(genes) < 5) {
    cat(sprintf("  [%s/%s] Skipped: only %d significant genes (need >=5)\n",
                experiment_name, ontology, length(genes)))
    return(NULL)
  }

  cat(sprintf("  [%s/%s] Running enrichGO with %d genes...\n",
              experiment_name, ontology, length(genes)))

  # tryCatch() is R's try/except (like Python's try/except).
  # It runs the code in the first block; if an error occurs, it runs the
  # error-handling block instead of crashing the whole script.
  #
  #   Python:  try:
  #                result = enrichGO(...)
  #            except Exception as e:
  #                print(e); result = None
  #
  result <- tryCatch({

    # ---- enrichGO(): the core ORA function from clusterProfiler ----
    # This performs Over-Representation Analysis using a hypergeometric test
    # (a generalization of Fisher's exact test).
    #
    # Parameters explained:
    #   gene:          the SIGNIFICANT gene list to test for enrichment.
    #   OrgDb:         the annotation database (org.Hs.eg.db for human).
    #   keyType:       our genes are given as SYMBOLS ("TP53"), not Entrez IDs.
    #   ont:           which GO ontology to test: "BP", "MF", or "CC".
    #   pAdjustMethod: how to correct for multiple testing. We test thousands
    #                  of GO terms, so raw p-values would give many false
    #                  positives. "BH" (Benjamini-Hochberg) controls the
    #                  false discovery rate (FDR). (Set in config.)
    #   pvalueCutoff:  only keep GO terms with adjusted p-value below this.
    #   qvalueCutoff:  also filter by q-value (another FDR measure).
    #   universe:      THE BACKGROUND — all proteins that could have been found.
    #                  Enrichment is calculated relative to this set.
    enrichGO(
      gene          = genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ontology,
      pAdjustMethod = GO_PADJUST_METHOD,
      pvalueCutoff  = GO_PVALUE_CUTOFF,
      qvalueCutoff  = GO_QVALUE_CUTOFF,
      universe      = universe
    )
  }, error = function(e) {
    # If enrichGO() throws an error (e.g. no genes map to this ontology),
    # print the error message and return NULL instead of crashing.
    cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })

  # Check whether any GO terms passed the significance cutoffs.
  # as.data.frame() converts the enrichment result object to a regular table.
  # If the table has 0 rows, nothing was significantly enriched.
  if (is.null(result) || nrow(as.data.frame(result)) == 0) {
    cat(sprintf("    No enriched terms found\n"))
    return(NULL)
  }

  # Success! Report how many enriched terms we found.
  cat(sprintf("    Found %d enriched GO terms\n", nrow(as.data.frame(result))))
  return(result)
}

# =====================================================================
# Helper: Visualize ORA results
# =====================================================================
# This function takes an enrichment result and generates standard plots:
#   - A dotplot (the most common enrichment figure in papers)
#   - A barplot
#   - A cnetplot (concept-gene network)
# It also saves the results table to CSV.
#
# VISUALIZATION TYPES EXPLAINED:
#
# dotplot:  Each row is a GO term, sorted by significance. Dot size = number
#           of genes; dot color = p-value (redder = more significant).
#           THE standard figure for enrichment results.
#
# barplot:  Similar to dotplot but uses bar LENGTH for gene count and color
#           for p-value. Simpler but less informative than dotplot.
#
# cnetplot: A NETWORK showing how genes (dots) connect to GO terms (larger
#           nodes). Useful for seeing which genes drive multiple enriched
#           terms. Helps spot hub genes that appear in many pathways.
#
visualize_ora <- function(ego, experiment_name, ontology) {
  # If there's no result, do nothing (early return — like "return" in Python).
  if (is.null(ego)) return()

  # paste0() concatenates strings with NO separator (like Python's "a" + "b").
  # paste() (without 0) would add a space between parts.
  # This builds a filename prefix like "GO_turboid_BP".
  prefix <- paste0("GO_", experiment_name, "_", ontology)

  # Convert the enrichment result object to a plain data frame for saving.
  res_df <- as.data.frame(ego)

  # Save results table (CSV). Our helper from utils.R.
  save_table(res_df, prefix)

  # ---- Dotplot (top 20) ----
  # dotplot() is from the enrichplot package.
  # showCategory = 20 means show only the top 20 most significant GO terms
  # (a full result might have hundreds — too many to display clearly).
  # The "+" adds a title and centers it, just like ggplot2 styling.
  p_dot <- dotplot(ego, showCategory = 20, title = paste(experiment_name, ontology)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  save_figure(p_dot, paste0(prefix, "_dotplot"), width = 9, height = 7)

  # ---- Barplot (top 20) ----
  # barplot() from enrichplot. Same data, bar-style visualization.
  # In this context, "barplot" is the enrichplot function, NOT base R's barplot().
  p_bar <- barplot(ego, showCategory = 20, title = paste(experiment_name, ontology)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  save_figure(p_bar, paste0(prefix, "_barplot"), width = 9, height = 7)

  # ---- Cnetplot (gene-concept network) — top 10 categories ----
  # Only make this if we have at least 3 enriched terms (otherwise the
  # network is trivial/uninformative).
  # categorySize = "pvalue" means node size reflects significance.
  # foldChange = NULL means we don't color genes by fold change (ORA doesn't
  #   inherently use fold change; that's more of a GSEA thing).
  if (nrow(res_df) >= 3) {
    p_cnet <- tryCatch({
      cnetplot(ego, showCategory = 10,
               categorySize = "pvalue",
               foldChange = NULL) +
        ggplot2::ggtitle(paste(experiment_name, ontology))
    }, error = function(e) NULL)

    # Only save if cnetplot() succeeded (it can fail with very small sets).
    if (!is.null(p_cnet)) {
      save_figure(p_cnet, paste0(prefix, "_cnetplot"), width = 12, height = 10)
    }
  }
}

# =====================================================================
# MAIN: Run ORA for each experiment and each ontology
# =====================================================================
# This is a NESTED LOOP: for each experiment, we test all 3 ontologies.
# If we have 3 experiments × 3 ontologies = 9 enrichment tests total.
# Each test can take 30+ seconds, so this is the slow part.
cat("\n[1/3] Running ORA (over-representation analysis)...\n")

# Create an EMPTY list to collect all ORA results.
# We'll fill it as we go. Lists in R grow dynamically.
# <- is R's assignment operator (like "=" in Python).
all_ora_results <- list()

# for loop over experiment names. Like Python's:
#   for exp_name in gene_sets.keys():
for (exp_name in names(gene_sets)) {
  cat(sprintf("\n  --- %s ---\n", exp_name))

  # [[exp_name]] extracts the element by name (like dict[key] in Python).
  # This gives us the vector of significant gene symbols for this experiment.
  sig_genes <- gene_sets[[exp_name]]
  cat(sprintf("  Significant genes: %d\n", length(sig_genes)))

  # Inner loop over the 3 GO ontologies: "BP", "MF", "CC".
  # GO_ONTOLOGIES is defined in 01_config.R.
  for (ont in GO_ONTOLOGIES) {
    # Run ORA for this experiment + this ontology.
    ego <- run_ora(sig_genes, ont, universe, exp_name)

    # Generate plots and save tables (if there's a result).
    visualize_ora(ego, exp_name, ont)

    # If we got a result, store it for later use (rrvgo reduction).
    if (!is.null(ego)) {
      # paste() with sep="_" joins strings with an underscore.
      # e.g. "turboid" + "_" + "BP" = "turboid_BP".
      result_key <- paste(exp_name, ont, sep = "_")

      # Store in our results list under this key.
      # [[result_key]] <- ego is like dict[key] = value in Python.
      all_ora_results[[result_key]] <- ego
    }
  }
}

# =====================================================================
# GSEA (Gene Set Enrichment Analysis) — uses ranked gene list
# =====================================================================
# GSEA is different from ORA in a key way:
#   - ORA needs a BINARY list: significant or not (a hard cutoff).
#   - GSEA uses a RANKED list of ALL genes, sorted by a metric that captures
#     both effect size and significance. No arbitrary cutoff needed.
#
# GSEA asks: "If I rank all proteins from most-up to most-down, do the members
# of a given GO term cluster unexpectedly high (or low) in that ranking?"
#
# This makes GSEA MORE SENSITIVE than ORA — it can detect pathways where
# many genes show modest but consistent changes (which ORA would miss because
# no single gene crosses the significance threshold).
cat("\n[2/3] Running GSEA (gene set enrichment analysis)...\n")

# Function definition: runs GSEA for one experiment.
run_gsea_for_experiment <- function(df, experiment_name) {
  cat(sprintf("\n  --- %s GSEA ---\n", experiment_name))

  # ---- Build the ranking metric ----
  # We need a single number that ranks each gene from "most worth attention"
  # to "least worth attention". A common choice is:
  #
  #   metric = sign(log2FC) × -log10(padj)
  #
  # Breaking this down:
  #   -log10(padj): converts p-values to a "significance score".
  #       padj = 0.05  ->  -log10 = 1.3
  #       padj = 0.001 ->  -log10 = 3.0
  #       Smaller p-value = bigger score = more significant.
  #   sign(log2FC): +1 for up-regulated, -1 for down-regulated, 0 for no change.
  #
  # Multiplying them gives a SIGNED significance score:
  #   - Large positive = significantly UP-regulated
  #   - Large negative = significantly DOWN-regulated
  #   - Near zero = not significant or no fold change
  #
  # The + 1e-10 prevents log10(0) which would be -Inf (infinity).
  # 1e-10 is scientific notation for 0.0000000001.
  df$metric <- -log10(df$padj + 1e-10) * sign(df$log2FC)

  # Remove rows where metric is NA (missing data).
  # is.na() checks for NA (R's version of "not available" / NaN).
  # The "!" means NOT, so !is.na() keeps rows that are NOT missing.
  df <- df[!is.na(df$metric), ]

  # Sort by metric, descending (largest first).
  # order() returns the sort order; decreasing = TRUE means largest first.
  # In Python this is like: df.sort_values("metric", ascending=False)
  df <- df[order(df$metric, decreasing = TRUE), ]

  # ---- Prepare the named ranked vector ----
  # gseGO() expects a NAMED NUMERIC VECTOR sorted in decreasing order.
  # The names are gene symbols; the values are the ranking metric.
  #
  # This is a specific R data structure: a vector where each element has a name.
  # Like a Python dict that's also ordered.
  gene_list <- df$metric
  names(gene_list) <- df$gene

  # Remove duplicates (GSEA requires unique gene names).
  # duplicated() flags the 2nd, 3rd, ... occurrence of each value as TRUE.
  gene_list <- gene_list[!duplicated(names(gene_list))]

  # Remove genes with metric exactly 0 (no fold change at all).
  # These sit in the middle of the ranking and add noise.
  gene_list <- gene_list[gene_list != 0]

  # ---- Run GSEA via gseGO() ----
  # gseGO() is clusterProfiler's function for GSEA against GO databases.
  #
  # Parameters:
  #   geneList:     the ranked, named numeric vector we just built.
  #   OrgDb:        human annotation database.
  #   keyType:      genes are SYMBOLS.
  #   ont = "BP":   GSEA is typically run on Biological Process only
  #                 (BP has the most biologically interpretable terms).
  #   minGSSize:    minimum gene-set size — ignore GO terms with fewer than
  #                 10 genes (too small to be reliable).
  #   maxGSSize:    maximum gene-set size — ignore GO terms with more than
  #                 500 genes (too broad/unspecific).
  #   pvalueCutoff: significance threshold for reporting enriched sets.
  #   verbose:      suppress extra console output.
  #
  # Note: GSEA does NOT use a "universe" argument — it uses the full ranked
  # list directly, which is why the universe concept applies to ORA but not GSEA.
  gsea_result <- tryCatch({
    gseGO(
      geneList     = gene_list,
      OrgDb        = org.Hs.eg.db,
      keyType      = "SYMBOL",
      ont          = "BP",       # GSEA typically on Biological Process
      minGSSize    = 10,
      maxGSSize    = 500,
      pvalueCutoff = GO_PVALUE_CUTOFF,
      verbose      = FALSE
    )
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })

  # Check for empty results.
  if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) == 0) {
    cat("    No enriched gene sets found\n")
    return(NULL)
  }

  cat(sprintf("    Found %d enriched gene sets\n", nrow(as.data.frame(gsea_result))))

  # Save results table (CSV).
  save_table(as.data.frame(gsea_result), paste0("GSEA_", experiment_name, "_BP"))

  # ---- Dotplot for GSEA ----
  # Same dotplot() function, but GSEA results also include a "NES" column
  # (Normalized Enrichment Score) indicating direction of enrichment.
  # Positive NES = enriched at top of ranking (up-regulated genes).
  # Negative NES = enriched at bottom (down-regulated genes).
  p_gsea <- dotplot(gsea_result, showCategory = 20,
                    title = paste("GSEA:", experiment_name, "BP")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  save_figure(p_gsea, paste0("GSEA_", experiment_name, "_BP_dotplot"), width = 10, height = 8)

  # ---- Ridgeplot (distribution of fold changes within enriched sets) ----
  # A ridge plot shows, for each enriched GO term, the DISTRIBUTION of fold
  # changes among its member genes. This tells you not just THAT a pathway
  # is enriched, but HOW: are all genes up? A mix? Bimodal?
  # Only generate if we have >= 3 enriched sets (otherwise not worth plotting).
  if (nrow(as.data.frame(gsea_result)) >= 3) {
    p_ridge <- tryCatch({
      ridgeplot(gsea_result) +
        ggplot2::ggtitle(paste("GSEA Ridge:", experiment_name))
    }, error = function(e) NULL)

    if (!is.null(p_ridge)) {
      save_figure(p_ridge, paste0("GSEA_", experiment_name, "_BP_ridge"), width = 10, height = 8)
    }
  }

  return(gsea_result)
}

# Run GSEA for each experiment.
# This loop calls our GSEA function once per experiment.
for (exp_name in names(experiments)) {
  run_gsea_for_experiment(experiments[[exp_name]], exp_name)
}

# =====================================================================
# GO Term Reduction (rrvgo) — collapse redundant enriched terms
# =====================================================================
# WHY REDUCE GO TERMS?
# -------------------
# GO enrichment often returns MANY highly similar terms because GO is
# structured hierarchically. For example, you might get:
#   - "mRNA splicing, via spliceosome"   (p = 1e-6)
#   - "RNA splicing, via transesterification" (p = 2e-6)
#   - "mRNA processing"                  (p = 5e-5)
#   - "nuclear mRNA splicing, via spliceosome" (p = 8e-5)
#
# All of these describe essentially the same biology. Reporting 20 nearly
# identical terms is confusing and unhelpful for interpretation.
#
# rrvgo solves this by:
#   1. Calculating semantic SIMILARITY between all pairs of GO terms
#      (how much they overlap in meaning, based on the GO graph structure).
#   2. CLUSTERING similar terms together.
#   3. Picking one REPRESENTATIVE term per cluster (the most significant one).
#
# This typically reduces 50+ terms to ~5-15 readable, distinct terms.
cat("\n[3/3] Reducing GO term redundancy (rrvgo)...\n")

# Load the rrvgo package (specialized for GO similarity + reduction).
library(rrvgo)

# Loop over every ORA result we collected earlier.
# for (key in names(...)) iterates over the keys of our results list.
for (key in names(all_ora_results)) {
  # Get the enrichment result for this experiment+ontology combination.
  ego <- all_ora_results[[key]]
  res_df <- as.data.frame(ego)

  # Only reduce if we have at least 5 enriched terms.
  # With fewer terms, reduction isn't meaningful (nothing to collapse).
  # "next" skips to the next iteration of the loop (like "continue" in Python).
  if (nrow(res_df) < 5) {
    cat(sprintf("  [%s] Skipped: only %d terms (need >=5)\n", key, nrow(res_df)))
    next
  }

  # Extract the ontology type from the key.
  # The key looks like "turboid_BP", so we split on "_" and take the last part.
  # strsplit() splits a string — returns a LIST (hence [[1]] to get the first element).
  # ont[length(ont)] gets the LAST element (the ontology: BP, MF, or CC).
  ont <- strsplit(key, "_")[[1]]
  ontology <- ont[length(ont)]

  cat(sprintf("  [%s] Reducing %d GO terms...\n", key, nrow(res_df)))

  # ---- Step 1: Calculate semantic similarity matrix ----
  # calculateSimMatrix() computes pairwise similarity between all GO term IDs.
  # method = "Rel" is the Rel similarity measure (one of several options;
  #   "Rel" accounts for information content and works well for GO).
  #
  # The result is a square MATRIX: rows and columns are GO term IDs,
  # and each cell holds the similarity score (0 = unrelated, 1 = identical).
  sim_matrix <- tryCatch({
    calculateSimMatrix(
      res_df$ID,
      orgdb = org.Hs.eg.db,
      ont = ontology,
      method = "Rel"
    )
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })

  # If similarity calculation failed, skip to the next result.
  if (is.null(sim_matrix)) next

  # ---- Step 2: Assign scores and reduce ----
  # Each GO term gets a SCORE used to pick the representative.
  # We use -log10(p.adjust) as the score: more significant = higher score.
  # setNames() attaches the GO term IDs as names to the score vector.
  scores <- setNames(-log10(res_df$p.adjust), res_df$ID)

  # reduceSimMatrix() clusters similar terms and picks representatives.
  # threshold = 0.7: terms with similarity > 0.7 get merged into one cluster.
  #   (Higher threshold = more terms kept; lower = more aggressive merging.)
  # The term with the highest score in each cluster becomes the representative.
  reduced <- reduceSimMatrix(sim_matrix, scores, threshold = 0.7,
                              orgdb = org.Hs.eg.db)

  # Save the reduced term list to CSV.
  save_table(as.data.frame(reduced), paste0("reduced_", key))

  # ---- Treemap visualization ----
  # A treemap shows the reduced terms as rectangles sized by their score
  # (significance). Bigger rectangles = more important representative terms.
  # This gives a quick visual summary of the main biological themes.
  p_tree <- tryCatch({
    treemapPlot(reduced) +
      ggplot2::ggtitle(paste("Reduced GO:", key))
  }, error = function(e) NULL)

  if (!is.null(p_tree)) {
    save_figure(p_tree, paste0("reduced_", key, "_treemap"), width = 10, height = 8)
  }

  # Report how many terms we collapsed down to.
  cat(sprintf("    Reduced to %d representative terms\n", nrow(reduced)))
}

# Final summary banner.
cat("\n=========================================\n")
cat(" GO enrichment analysis complete!\n")
cat(sprintf(" Figures: %s/\n", FIGURE_DIR))
cat(sprintf(" Tables:  %s/\n", TABLE_DIR))
cat("=========================================\n")
