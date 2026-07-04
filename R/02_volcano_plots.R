###############################################################################
# 02_volcano_plots.R
# Lydia-style volcano plots with multi-category highlighting.
#
# ─────────────────────────────────────────────────────────────────────────────
# WHAT IS A VOLCANO PLOT?  (Biological background for beginners)
# ─────────────────────────────────────────────────────────────────────────────
# A volcano plot is the most common visualization in proteomics/genomics.
# Imagine ~3000 proteins measured in a mass spectrometer. For each protein
# we know two things:
#
#   1. FOLD CHANGE (log2FC): how MUCH the protein's abundance changed
#      between conditions (e.g., TRIP4-expressing cells vs. control).
#      - log2FC = 0  → no change
#      - log2FC = +1 → protein DOUBLED in abundance (2-fold increase)
#      - log2FC = +2 → protein went up 4-fold
#      - log2FC = -1 → protein HALVED (2-fold decrease)
#      We use log2 so that a 2-fold up and 2-fold down are symmetric.
#
#   2. ADJUSTED P-VALUE (padj): how CONFIDENT we are the change is real
#      (not just noise). We transform it with -log10() so that very
#      significant p-values (0.00001) become large numbers (5) at the
#      top of the plot. The p-value is "adjusted" to correct for testing
#      thousands of proteins at once (False Discovery Rate correction).
#
# The plot looks like a volcano (hence the name):
#   - X-axis = log2 fold change  (left = down, right = up)
#   - Y-axis = -log10(padj)      (higher = more confident)
#   - Each DOT = one protein
#   - The most interesting proteins are in the TOP-RIGHT and TOP-LEFT corners:
#     they changed a lot AND we're very confident about it.
#
# This script creates an ENHANCED volcano plot (Lydia's style) that doesn't
# just show "significant / not significant" — it overlays MULTIPLE biological
# categories so you can see at a glance which significant proteins are:
#   - Known interactors (already published as TRIP4 partners)
#   - Flag IP hits (confirmed by a separate co-IP experiment)
#   - Members of interesting gene families (RNA helicases, GPATCH, LARPs)
#   - High-confidence hits (extremely strong signal)
#
# ─────────────────────────────────────────────────────────────────────────────
# THE CATEGORY PRIORITY SYSTEM  (Read this carefully!)
# ─────────────────────────────────────────────────────────────────────────────
# Every protein gets ONE color on the plot. But a single protein could belong
# to several categories at once (e.g., a gene could be a known interactor
# AND a Flag IP hit AND a DHX family member). Since a dot can only have one
# color, we need a priority system.
#
# The categories are applied in LAYERS, like painting on a canvas:
#
#   Layer 0 (base):    Every protein starts as "TRUE" (significant) or
#                      "FALSE" (not significant).
#
#   Layer 1:           Significant proteins in gene families (DHX, DDX, LARP,
#                      GPATCH) get re-colored with their family color.
#                      This OVERWRITES the "TRUE" base layer.
#
#   Layer 2:           Flag IP hits overwrite gene-family coloring.
#                      A protein found in 2+ Flag IP conditions → "flagMulti".
#                      A protein found in exactly 1 Flag IP condition → "flagOnce".
#                      This OVERWRITES layer 1.
#
#   Layer 3:           Known interactors overwrite Flag IP coloring.
#                      These are the most important biology, so they get
#                      the highest visual priority (orange). This OVERWRITES
#                      layer 2.
#
#   Layer 4:           High-confidence hits overwrite EVERYTHING.
#                      If a protein has extreme fold-change AND extreme
#                      significance, it gets painted RED regardless of any
#                      other category. This has the HIGHEST priority.
#
#   RULE: Later layers always overwrite earlier ones.
#   So the priority order (lowest to highest) is:
#     FALSE/TRUE → gene families → Flag IP → known interactors → high-conf
#
# ─────────────────────────────────────────────────────────────────────────────
# Categories overlaid on each volcano:
#   - ia: Known interactors (from known_interactors.txt)
#   - flagMulti / flagOnce: Flag IP hits (from Flag IP data)
#   - CRAC: RNA interactome hits (if CRAC data available)
#   - gp / dhx / ddx / LARPs: Gene families
#   - high: High-confidence hits (top combined significance)
#
# Usage:
#   source("R/01_config.R")
#   source("R/utils.R")
#   source("R/02_volcano_plots.R")
###############################################################################

# ── cat() prints text to the console (like print() in Python). ──────────────
# \n means "new line" — it moves the cursor to the next line.
# We use cat() with multiple \n to create a visual separator banner.
cat("\n=========================================\n")
cat(" Volcano Plot Generation (Lydia-style)\n")
cat("=========================================\n\n")

# ── LOAD KNOWN INTERACTORS ──────────────────────────────────────────────────
# Known interactors are proteins already published or experimentally confirmed
# as binding partners of TRIP4. We load them from a text file (one gene name
# per line) so we can highlight them in orange on the volcano plot.
# This confirms that our mass spec data "rediscovers" known biology (validation).

# file.path() joins path pieces together correctly for the OS.
# It joins DATA_DIR (defined in 01_config.R, e.g., "project/data") with the
# filename. On Linux this produces "project/data/known_interactors.txt".
# The result is stored (assigned) to a variable.
#
# <- is R's ASSIGNMENT operator (like = in Python or most other languages).
# In R, <- and = do almost the same thing for simple assignments.
# The R community convention is to use <- for variables. We use it throughout.
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")

# load_known_interactors() is a function defined in utils.R.
# It reads the text file and returns a CHARACTER VECTOR (a list of gene-name
# strings). A "character vector" in R is like a Python list of strings:
#   c("TRIP4", "ASCC1", "ASCC2", ...)
# We store the result in known_interactors so we can use it later.
known_interactors <- load_known_interactors(interactors_file)

# ── LOAD ALL EXPERIMENT DATA ────────────────────────────────────────────────
# This loads every mass spec CSV file found in the data/ directory.
# Each file becomes a DATA FRAME (a table — think of an Excel spreadsheet):
#   - One row per protein
#   - Columns: gene (name), log2FC (fold change), padj (adjusted p-value)
# The result is a NAMED LIST: a collection where you access each experiment
# by name. For example:
#   experiments$turbo_trip4_vs_wt   → a data frame of ~3000 proteins
#   experiments$flag_cflag_vs_ctrl  → another data frame
#   (and so on for every experiment)
cat("Loading experiment data...\n")
experiments <- load_all_experiments()

# ── BUILD FLAG IP HIT LISTS FOR OVERLAY ─────────────────────────────────────
# In the Flag IP (co-immunoprecipitation) experiments, proteins that
# significantly associate with the Flag-tagged bait are "Flag IP hits."
# If a protein shows up as significant in MULTIPLE Flag IP conditions,
# that's stronger evidence it's a real interactor.
#
# We split Flag IP hits into two categories:
#   flagMulti = hit in 2 or more conditions (stronger evidence → purple)
#   flagOnce  = hit in exactly 1 condition   (weaker evidence → blue)

# grep() searches for a PATTERN in a vector of strings.
# Here we search the NAMES of the experiments list for those starting with
# "flag_" (the ^ means "start of string").
#   - grep(pattern, x) normally returns the INDEX (position number) of matches.
#   - value = TRUE changes the behavior: instead of returning positions,
#     it returns the actual matched STRING VALUES.
# So this gives us a character vector like:
#   c("flag_cflag_vs_ctrl", "flag_nflag_vs_ctrl", ...)
flag_exp_names <- grep("^flag_", names(experiments), value = TRUE)

# if () { ... } else { ... } is a CONDITIONAL, just like in Python/other languages.
# The condition inside the parentheses must be TRUE or FALSE.
# length() returns how many elements are in the vector.
# So "length(flag_exp_names) >= 2" means "do we have at least 2 Flag IP experiments?"
if (length(flag_exp_names) >= 2) {
  cat("\nBuilding Flag IP hit categories...\n")

  # lapply() applies a function to each element of a list and returns a new list.
  # "l" in lapply stands for "list" — it always returns a list.
  # Here:
  #   experiments[flag_exp_names]   → a sub-list of only the Flag IP experiments
  #   get_significant_genes         → a function (from 01_config.R) that returns
  #                                   the names of genes passing the significance
  #                                   thresholds (padj < 0.05 AND |log2FC| > 0.5)
  # So flag_sig_lists becomes a LIST of character vectors, one per Flag IP
  # experiment, each containing the significant gene names from that experiment:
  #   flag_sig_lists[[1]] = c("TRIP4", "DHX9", "DDX3", ...)
  #   flag_sig_lists[[2]] = c("TRIP4", "GPATCH1", ...)
  flag_sig_lists <- lapply(experiments[flag_exp_names], get_significant_genes)

  # unlist() flattens a list into a single vector.
  # If two experiments each found "TRIP4", the result has "TRIP4" twice.
  # This is intentional — we want to COUNT how many times each gene appears.
  flag_all_genes <- unlist(flag_sig_lists)

  # table() creates a frequency table: it counts how many times each unique
  # value appears. Example:
  #   table(c("A", "A", "B"))  →  A:2, B:1
  # So flag_counts tells us, for each gene, in how many Flag IP conditions
  # it was found significant.
  flag_counts <- table(flag_all_genes)

  # names() gets the names (the gene symbols) from the frequency table.
  # We use [condition] to keep only genes where the count meets a threshold:
  #   flag_counts >= 2  → keep genes found in 2+ conditions (flagMulti)
  #   flag_counts == 1  → keep genes found in exactly 1 condition (flagOnce)
  # These produce CHARACTER VECTORS of gene names.
  flag_multi <- names(flag_counts)[flag_counts >= 2]
  flag_once <- names(flag_counts)[flag_counts == 1]

  # sprintf() formats a string with placeholders (%d = integer).
  # length() returns the count of elements.
  # This prints something like:
  #   "  Flag IP hits: 45 multi-condition, 120 single-condition"
  cat(sprintf("  Flag IP hits: %d multi-condition, %d single-condition\n",
              length(flag_multi), length(flag_once)))
} else {
  # The else branch runs when there are fewer than 2 Flag IP experiments.
  # We can't compute "multi-condition" hits with just one experiment,
  # so we set both lists to EMPTY.
  # character(0) is an empty character vector (a list of zero strings).
  # This is R's way of saying "no elements" for text data.
  flag_multi <- character(0)
  flag_once <- character(0)
}

# ── ASSIGN GENE FAMILY CATEGORIES ───────────────────────────────────────────
# This function takes a vector of gene names and returns a vector identifying
# which gene FAMILY each gene belongs to (or NA if it's not in any family).
#
# Why gene families? TRIP4 and its partners are heavily involved in RNA
# processing. The families we highlight are all RNA-binding proteins:
#   DHX  = DEAD/H-box helicases (unwind RNA/DNA, e.g., DHX9, DHX15)
#   DDX  = DEAD-box helicases  (RNA unwinding, e.g., DDX3, DDX5)
#   LARP = La-related proteins (RNA binding/stabilization, e.g., LARP1)
#   GPATCH = G-patch domain proteins (RNA processing, e.g., GPATCH1)
# Highlighting these families helps spot patterns: if many helicases appear
# as hits, it supports TRIP4's role in RNA metabolism.

# function() defines a NEW function. The syntax is:
#   function(parameters) { ...body... }
# Here, the function takes one parameter called "genes" (a character vector
# of gene names). We assign the whole function to the name assign_gene_family.
assign_gene_family <- function(genes) {
  # rep() creates a vector by REPEATING a value a certain number of times.
  # NA_character_ is R's way of saying "missing text value" (NA for character
  # type). We create one NA per gene, as the starting default — meaning
  # "this gene is not in any family yet."
  # length(genes) = how many genes we received.
  fam <- rep(NA_character_, length(genes))

  # grepl() checks whether each string MATCHES a REGULAR EXPRESSION pattern.
  # It returns TRUE/FALSE for each gene (a logical vector).
  # "^DHX" means "starts with DHX" (^ = start of string).
  # We then use [mask] to ASSIGN a family label to just those positions.
  #
  # The assignment fam[grepl("^DHX", genes)] <- "dhx" means:
  #   "For every position where the gene starts with DHX, set fam to 'dhx'."
  # This OVERWRITES the NA default for matching genes.
  fam[grepl("^DHX", genes)] <- "dhx"
  fam[grepl("^DDX", genes)] <- "ddx"
  fam[grepl("^LARP", genes)] <- "LARPs"

  # GPATCH genes don't have a simple naming pattern, so we can't use a prefix
  # match like ^GPATCH. Instead, GENE_FAMILIES$GPATCH is a hard-coded list
  # of all known GPATCH-family gene names (defined in 01_config.R).
  # %in% checks membership: "is this gene IN the GPATCH list?"
  # It returns TRUE/FALSE for each gene. We use that to set "gp".
  fam[genes %in% GENE_FAMILIES$GPATCH] <- "gp"

  # return() sends the result back to whoever called this function.
  # The result is a character vector the same length as the input, with either
  # a family name ("dhx", "ddx", "LARPs", "gp") or NA for each gene.
  return(fam)
}

# =====================================================================
# Helper: Lydia-style volcano with category highlighting
# =====================================================================
# This is the MAIN plotting function. It takes a data frame of protein results
# and produces a beautiful, multi-layered volcano plot with category colors.
#
# Parameters (inputs to the function):
#   df                   → the data frame (gene, log2FC, padj columns)
#   title                → text to show at the top of the plot
#   known_ia = NULL      → known interactor gene names (default: none)
#   flag_m = NULL        → flagMulti gene names (default: none)
#   flag_o = NULL        → flagOnce gene names (default: none)
#   show_gene_families = TRUE  → highlight gene-family proteins? (default: yes)
#   highlight_high = TRUE      → highlight ultra-significant proteins? (default: yes)
#
# The "= NULL" and "= TRUE" syntax provides DEFAULT VALUES. If the caller
# doesn't supply these arguments, the defaults are used. This makes the
# function flexible — you can call it with just df and title, or override
# any of the optional parameters.
plot_lydia_volcano <- function(df, title, known_ia = NULL,
                               flag_m = NULL, flag_o = NULL,
                               show_gene_families = TRUE,
                               highlight_high = TRUE) {

  # Make a working copy of the input data frame (so we don't modify the
  # original). In R, <- creates a new variable; here toPlot references the
  # same data frame, and we'll add a new column to it.
  toPlot <- df

  # ── LAYER 0 (BASE): Significance classification ───────────────────────────
  # Start with significance as category.
  #
  # ifelse() is R's VECTORIZED conditional. It's like a Python list
  # comprehension that applies to every element at once:
  #
  #   ifelse(condition, value_if_TRUE, value_if_FALSE)
  #
  # For each ROW of the data frame, it evaluates the condition:
  #   - padj < P_VALUE_CUTOFF: is the adjusted p-value below 0.05?
  #   - & means logical AND (both conditions must be true)
  #   - abs(log2FC) > LOG2FC_CUTOFF: is the fold change (absolute value)
  #     bigger than 0.5? (abs() ignores the sign, so both up and down count)
  #
  # If BOTH are true → the protein is significant → category = "TRUE" (a string)
  # Otherwise → not significant → category = "FALSE"
  #
  # The $ operator accesses a COLUMN by name from the data frame:
  #   toPlot$padj   = "the padj column"
  #   toPlot$log2FC = "the log2FC column"
  #
  # We store the result in a NEW column called "category".
  toPlot$category <- ifelse(
    toPlot$padj < P_VALUE_CUTOFF & abs(toPlot$log2FC) > LOG2FC_CUTOFF,
    "TRUE", "FALSE"
  )

  # ── LAYER 1: Gene families (only for significant proteins) ────────────────
  # Now we start OVERWRITING. Only significant proteins (category == "TRUE")
  # can get a gene-family label. Non-significant proteins stay "FALSE."
  #
  # Why only significant ones? A gene might be in the DHX family, but if its
  # abundance didn't change significantly, it's not interesting for this plot.
  if (show_gene_families) {
    # Call our helper function to get a family label for every gene:
    fam <- assign_gene_family(toPlot$gene)

    # Create a LOGICAL MASK (a TRUE/FALSE vector, one value per protein):
    #   !is.na(fam)         → the gene IS in a family (fam is not missing)
    #   &                   → AND
    #   category == "TRUE"  → the protein is also significant
    # Only proteins meeting BOTH conditions get re-categorized.
    fam_mask <- !is.na(fam) & toPlot$category == "TRUE"

    # Apply the mask: for positions where fam_mask is TRUE, set the category
    # to the family name (e.g., "dhx", "ddx", "LARPs", "gp").
    # This OVERWRITES the "TRUE" base layer for these proteins.
    toPlot$category[fam_mask] <- fam[fam_mask]
  }

  # ── LAYER 2: Flag IP hits ─────────────────────────────────────────────────
  # Flag IP hits overwrite gene-family labels. These are proteins confirmed
  # by a separate co-immunoprecipitation experiment, so they get higher
  # visual priority than family membership.
  #
  # First check that the inputs exist and aren't empty:
  #   !is.null(flag_m) → flag_m was provided (not NULL)
  #   length(flag_m) > 0 → it has at least one gene name
  # Both conditions must be true (& is logical AND).
  if (!is.null(flag_m) && length(flag_m) > 0) {
    # %in% checks membership: which of toPlot's genes are in the flag_m list?
    # toPlot$category[...] <- "flagMulti" overwrites the category for those
    # proteins. This OVERWRITES layer 1 (gene families).
    toPlot$category[toPlot$gene %in% flag_m] <- "flagMulti"
  }
  if (!is.null(flag_o) && length(flag_o) > 0) {
    # Same logic, for single-condition Flag IP hits.
    toPlot$category[toPlot$gene %in% flag_o] <- "flagOnce"
  }

  # ── LAYER 3: Known interactors (highest priority so far) ──────────────────
  # Known interactors are the most biologically important category — these are
  # proteins already confirmed by previous research to bind TRIP4. They get
  # ORANGE and overwrite everything from layers 0–2.
  if (!is.null(known_ia) && length(known_ia) > 0) {
    toPlot$category[toPlot$gene %in% known_ia] <- "ia"
  }

  # ── LAYER 4: High-confidence hits (HIGHEST priority overall) ──────────────
  # Proteins with EXTREME fold-change AND EXTREME significance are the most
  # striking findings. They get RED and overwrite ALL other categories.
  #
  # The condition has two parts joined by | (logical OR):
  #   (A) log2FC > 2 AND -log10(padj) > 6
  #       → very strong change (4-fold+) at extreme confidence
  #   (B) log2FC > 7 AND -log10(padj) > 2
  #       → enormous change (128-fold+) at moderate confidence
  # Either condition qualifies a protein as "high."
  #
  # -log10(padj): converts tiny p-values to positive numbers for the Y-axis.
  #   padj = 0.000001 → -log10(0.000001) = 6 (very significant)
  #
  # !is.na(toPlot$padj): exclude proteins with a missing p-value
  # (some proteins have NA padj if they couldn't be tested).
  if (highlight_high) {
    toPlot$category[!is.na(toPlot$padj) &
                    ((toPlot$log2FC > 2 & -log10(toPlot$padj) > 6) |
                     (toPlot$log2FC > 7 & -log10(toPlot$padj) > 2))] <- "high"
  }

  # ── CONVERT CATEGORY TO A FACTOR ──────────────────────────────────────────
  # factor() converts a character vector into a CATEGORICAL variable with
  # a FIXED SET of possible values called "levels," and the levels have an
  # ORDER. This matters for the plot because:
  #   1. The legend will show ALL categories in a fixed order (even if some
  #      don't appear in this particular experiment).
  #   2. Colors are assigned based on level order.
  #
  # levels = names(CATEGORY_COLORS) sets the levels to the names of the
  # color scheme defined in 01_config.R (e.g., "ia", "flagMulti", "FALSE"...).
  # The ORDER of names determines the order in the legend.
  #
  # If a protein's category string isn't in the levels, it becomes NA
  # (which would be invisible on the plot) — but our categories always match.
  toPlot$category <- factor(toPlot$category, levels = names(CATEGORY_COLORS))

  # ── DECIDE WHICH POINTS TO LABEL WITH GENE NAMES ──────────────────────────
  # Labeling ALL ~3000 proteins would make the plot unreadable. We only label
  # proteins in "interesting" categories (not plain significant or not-significant).
  #
  # c() creates a CHARACTER VECTOR by combining values.
  # These are the category names that should get text labels.
  label_cats <- c("ia", "gp", "dhx", "ddx", "LARPs", "flagMulti", "flagOnce", "high")

  # Filter toPlot to keep only rows whose category is in label_cats.
  # toPlot[condition, ] means "select ROWS where condition is TRUE, keep ALL
  # columns" (the comma followed by nothing means "all columns").
  # %in% checks if each category value is in the label_cats vector.
  label_data <- toPlot[toPlot$category %in% label_cats, ]

  # ── BUILD THE GGPLOT (the actual visualization) ───────────────────────────
  # ggplot2 is R's premier plotting package. The key idea: you build a plot
  # by LAYERING components with the + operator. Think of it like stacking
  # transparent sheets, each adding one visual element.
  #
  # The :: syntax means "use this function from this package."
  #   ggplot2::ggplot = "use ggplot from the ggplot2 package"
  # This is explicit and avoids name conflicts with other packages.

  # ggplot(data, aes(...)) initializes the plot:
  #   - data = toPlot: the data frame to plot
  #   - aes() = "aesthetic MAPPINGS" — it maps DATA COLUMNS to VISUAL
  #     PROPERTIES of the plot:
  #       x = log2FC    → the X-axis position comes from the log2FC column
  #       y = -log10(padj) → the Y-axis comes from the negative log10 of padj
  #       color = category → dot COLOR is determined by the category column
  #     aes() doesn't draw anything — it just sets up the mapping. The actual
  #     drawing is done by geom_* functions added with +.
  #
  # The + operator ADDS layers to the plot. It's NOT arithmetic addition —
  # it's a special ggplot2 syntax for composing plots layer by layer.
  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = -log10(padj), color = category)) +

    # geom_point() draws one DOT (point) for each row in the data.
    # alpha = 0.3: dots are 30% opaque (70% transparent), so overlapping
    #   points create darker regions — you can see density.
    # size = 1.2: each dot is 1.2 mm in radius (small, since there are thousands).
    ggplot2::geom_point(alpha = 0.3, size = 1.2) +

    # scale_color_manual() lets you specify EXACT colors for each category
    # using HEX CODES (like HTML colors: #d95f02 = orange).
    #   values = CATEGORY_COLORS: the named vector from 01_config.R that maps
    #     each category name to a color hex code.
    #   drop = FALSE: keep ALL categories in the legend, even if some don't
    #     appear in this particular data set. This keeps legends consistent
    #     across all plots.
    ggplot2::scale_color_manual(
      values = CATEGORY_COLORS,
      drop = FALSE
    ) +

    # labs() sets the AXIS LABELS, TITLE, and LEGEND TITLE.
    # expression() creates MATHEMATICAL NOTATION for the axis labels:
    #   Log[2]~Fold~Change → renders as "Log₂ Fold Change" with a subscript 2
    #   -Log[10]~(...)     → renders as "-Log₁₀ (adj. p value)" with italic p
    # The ~ symbol means "space" in expression notation.
    # italic() makes the "p" italic (standard scientific notation).
    # color = "Category" sets the LEGEND TITLE.
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adj.~italic(p)~value)),
      title = title,
      color = "Category"
    )

  # ── ADD TEXT LABELS FOR HIGHLIGHTED POINTS ────────────────────────────────
  # We add gene-name labels ONLY for the highlighted points (label_data),
  # not for all 3000 proteins. Otherwise the plot would be unreadable.
  #
  # nrow() returns the number of rows in a data frame.
  # We check that there's at least one row to label.
  if (nrow(label_data) > 0) {
    # ggrepel::geom_text_repel() is from the "ggrepel" package.
    # "repel" means it uses a physics-like algorithm to PUSH labels apart
    # so they DON'T OVERLAP. Standard text labels would pile on top of each
    # other; ggrepel solves this by drawing thin connecting lines from each
    # dot to its displaced label.
    #
    # Parameters:
    #   data = label_data: use ONLY the highlighted points (not all data)
    #   aes(label = gene): map the "gene" column to the text shown
    #   size = 2.5: font size (in points)
    #   fontface = "bold": bold text
    #   max.overlaps = 30: if a label would overlap more than 30 others,
    #     skip it (otherwise the algorithm struggles with dense regions)
    #   show.legend = FALSE: don't add a legend entry for the text layer
    #
    # We add this layer to the plot with + and re-assign the result to p.
    p <- p + ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.5, fontface = "bold",
      max.overlaps = 30,
      show.legend = FALSE
    )
  }

  # ── THRESHOLD LINES ───────────────────────────────────────────────────────
  # Draw dashed reference lines so the viewer can see where the significance
  # thresholds are. Points beyond these lines are "significant."
  #
  # geom_hline(): draws a HORIZONTAL line across the plot.
  #   yintercept = -log10(P_VALUE_CUTOFF): the Y position of the line.
  #     Since P_VALUE_CUTOFF = 0.05, this is -log10(0.05) ≈ 1.3.
  #     Points ABOVE this line have padj < 0.05 (statistically significant).
  #   linetype = "dashed": dotted/dashed line style
  #   color = "grey50": medium grey
  #
  # geom_vline(): draws a VERTICAL line.
  #   xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF): TWO lines (the c() makes
  #     a vector with two values) at x = -0.5 and x = +0.5.
  #     Points to the RIGHT of +0.5 are up-regulated (more abundant).
  #     Points to the LEFT of -0.5 are down-regulated (less abundant).
  p <- p +
    ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50")

  # ── THEME (visual styling — Lydia's compact style) ────────────────────────
  # theme_bw() applies a clean "black and white" theme: white background,
  # black axis lines, grey gridlines. This is a good starting point.
  #
  # theme() allows FINE-TUNING of every visual element. Each argument
  # controls one piece of the plot:
  #
  #   element_text() configures text appearance (size, color, angle, etc.).
  #
  #   axis.text.x/y: the tick labels on the X and Y axes
  #     - colour = "black": black text (default is grey; black is more readable)
  #     - size = 8: font size in points
  #     - angle = 90, vjust = 0.5, hjust = 1: rotate X labels 90° and align
  #       them so they don't overlap (only matters if labels are long)
  #
  #   axis.title.x/y: the axis titles ("Log2 Fold Change", etc.)
  #   legend.title/text: the legend heading and entries
  #   plot.title: the main title at top
  #     - hjust = 0.5: center the title (0 = left, 0.5 = center, 1 = right)
  #     - face = "bold": bold font
  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(colour = "black", size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(colour = "black", size = 8),
      axis.title.x = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10)
    )

  # return() sends the finished plot object back to the caller.
  # The plot isn't drawn yet — it's a "p" object that can be saved or displayed.
  return(p)
}

# =====================================================================
# Helper: Overlay two experiments on one volcano
# =====================================================================
# This function takes TWO data frames (from two different experiments) and
# plots them on the SAME axes, colored by which experiment they came from.
# This lets you visually compare: "do the same proteins come up in both
# experiments?" Overlapping clusters suggest consistent findings.
#
# Parameters:
#   df1, df2     → the two data frames to overlay
#   name1, name2 → labels for each (shown in legend)
#   label_genes  → which genes to label with text (e.g., known interactors)
#   title        → plot title
plot_volcano_overlay <- function(df1, name1, df2, name2, label_genes, title) {
  # Add a new column "experiment" to each data frame so we can tell them
  # apart after combining. The $ operator creates/sets a column.
  # This is like adding a column in a spreadsheet.
  df1$experiment <- name1
  df2$experiment <- name2

  # Calculate significance for each data frame. We create a "sig" column
  # that is TRUE/FALSE for each protein. This uses the same thresholds as
  # elsewhere: padj < 0.05 AND |log2FC| > 0.5.
  # NOTE: Line 2 has a known quirk — it uses df1$log2FC for df2's condition
  # (likely a copy-paste artifact). We preserve the original code exactly.
  # & is logical AND. abs() is absolute value.
  df1$sig <- df1$padj < P_VALUE_CUTOFF & abs(df1$log2FC) > LOG2FC_CUTOFF
  df2$sig <- df2$padj < P_VALUE_CUTOFF & abs(df1$log2FC) > LOG2FC_CUTOFF

  # rbind() ("row bind") stacks two data frames VERTICALLY — it appends the
  # rows of the second below the first, creating one combined data frame.
  # But first, we select only the columns we need from each:
  #   df1[, c("gene", "log2FC", "padj", "experiment", "sig")]
  # means "all rows, only these 5 columns."
  # The two data frames must have the same columns for rbind() to work,
  # which is why we added the "experiment" and "sig" columns above.
  combined <- rbind(
    df1[, c("gene", "log2FC", "padj", "experiment", "sig")],
    df2[, c("gene", "log2FC", "padj", "experiment", "sig")]
  )

  # Calculate the Y-axis value: -log10 of the adjusted p-value.
  # This transforms small p-values (0.00001) into larger numbers (5) so
  # significant proteins appear near the top of the plot.
  # We add this as a new column to the combined data frame.
  combined$neglog10p <- -log10(combined$padj)

  # We only want to label each gene ONCE (not once per experiment), to avoid
  # duplicate labels. So we filter label_data to include only genes from the
  # FIRST experiment. This deduplicates the labels.
  #   combined$gene %in% label_genes: gene is in our "interesting" list
  #   & combined$experiment == name1: AND it's from the first experiment
  label_data <- combined[combined$gene %in% label_genes & combined$experiment == name1, ]

  # Build the overlay plot. Same ggplot2 layering approach as before:
  #
  # ggplot() initializes with the combined data. This time:
  #   color = experiment → dots are colored by WHICH EXPERIMENT they came from
  #     (TurboID vs. Flag IP), not by significance category.
  #
  # The + layers follow the same pattern as plot_lydia_volcano.
  p <- ggplot2::ggplot(combined, ggplot2::aes(x = log2FC, y = neglog10p, color = experiment)) +
    ggplot2::geom_point(alpha = 0.3, size = 1) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.5, fontface = "bold",
      max.overlaps = 25, show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = EXPERIMENT_COLORS) +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adj.~italic(p)~value)),
      title = title, color = "Experiment"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  return(p)
}

# =====================================================================
# MAIN: Generate all volcano plots
# =====================================================================
# This section actually CREATES and SAVES the plots. Everything above was
# function definitions and data preparation; now we put it all together.
# We generate 4 sets of plots:
#   [1/4] One Lydia-style volcano per TurboID experiment
#   [2/4] One Lydia-style volcano per Flag IP experiment
#   [3/4] An overlay comparing TurboID TRIP4 vs. Flag IP C-Flag
#   [4/4] An overlay comparing TRIP4 with vs. without Retinoic Acid

# ── 1. Lydia-style volcano for each TurboID experiment ─────────────────────
# TurboID is a proximity-labeling technique: TRIP4 is fused to an enzyme
# (TurboID) that biotinylates (tags) all nearby proteins. We then pull down
# the tagged proteins and measure them by mass spectrometry.
# This reveals TRIP4's NEIGHBORHOOD — not just direct binders, but proteins
# in the same complex or vicinity.

cat("\n[1/4] TurboID volcano plots (Lydia-style)...\n")

# Find all experiment names starting with "turbo_".
# grep() with value = TRUE returns the matching names (not positions).
# ^turbo_ means "starts with turbo_".
turbo_names <- grep("^turbo_", names(experiments), value = TRUE)

# for (variable in collection) { ... } is R's FOR LOOP.
# It iterates over each element in turbo_names, running the body once per name.
# Each iteration, "name" takes on the next experiment name string.
for (name in turbo_names) {
  # Call our main plotting function for this experiment.
  # experiments[[name]] uses DOUBLE BRACKETS to extract a single element
  # from the list by name. (Single [ would return a sub-list of length 1;
  # double [[ extracts the data frame itself.)
  #
  # We pass the known interactors and Flag IP hit lists so they can be
  # highlighted on the plot. show_gene_families = TRUE turns on family coloring.
  p <- plot_lydia_volcano(
    experiments[[name]], name,
    known_ia = known_interactors,
    flag_m = flag_multi, flag_o = flag_once,
    show_gene_families = TRUE
  )
  # save_figure() (defined in utils.R) saves the plot as both PNG and PDF.
  # paste0() concatenates strings WITHOUT spaces: paste0("volcano_", name)
  # produces "volcano_turbo_trip4_vs_wt". width and height are in inches.
  save_figure(p, paste0("volcano_", name), width = 7, height = 5)
}

# ── 2. Lydia-style volcano for each Flag IP experiment ─────────────────────
# Flag IP = Flag-tag immunoprecipitation. TRIP4 (or its fragments) is tagged
# with a Flag peptide. We use an anti-Flag antibody to pull down TRIP4 and
# anything bound to it, then identify those proteins by mass spectrometry.
# This is a more DIRECT method than TurboID — it captures actual physical
# complexes, but may miss transient or weak interactions.
#
# Note: We DON'T pass flag_m/flag_o here because these ARE the Flag IP
# experiments — highlighting them on their own plot would be circular.
# We DO pass known_ia (known interactors) for validation.

cat("\n[2/4] Flag IP volcano plots (Lydia-style)...\n")
flag_names <- grep("^flag_", names(experiments), value = TRUE)
for (name in flag_names) {
  p <- plot_lydia_volcano(
    experiments[[name]], name,
    known_ia = known_interactors,
    show_gene_families = TRUE
  )
  save_figure(p, paste0("volcano_", name), width = 7, height = 5)
}

# ── 3. Overlay: TurboID TRIP4 vs Flag IP C-Flag ────────────────────────────
# This overlay puts TWO experiments on the same axes:
#   - TurboID TRIP4 (proximity labeling — captures the neighborhood)
#   - Flag IP C-Flag (co-IP — captures direct/strong binders)
# Comparing them side-by-side reveals which proteins are found by BOTH
# methods (high confidence) vs. only one method (possibly method-specific).
# Proteins appearing in both the top-right of both colors are the strongest
# candidates for real TRIP4 interactors.

cat("\n[3/4] TurboID vs Flag IP overlay...\n")

# Define the canonical experiment names we want to overlay.
turbo_main <- "turbo_trip4_vs_wt"
flag_main <- "flag_cflag_vs_ctrl"

# Check that BOTH experiments exist before trying to plot.
# %in% checks membership. && is logical AND (returns a single TRUE/FALSE,
# unlike & which is vectorized — here we have single values so either works).
if (turbo_main %in% names(experiments) && flag_main %in% names(experiments)) {
  # Call the overlay function with both data frames and descriptive labels.
  # The label_genes are the known interactors — we want to see if these
  # known partners show up in both methods.
  p <- plot_volcano_overlay(
    experiments[[turbo_main]], "TurboID_TRIP4",
    experiments[[flag_main]], "FlagIP_C",
    label_genes = known_interactors,
    title = "TurboID (TRIP4 vs WT) vs Flag IP (C-Flag vs Ctrl)"
  )
  save_figure(p, "volcano_overlay_turboid_flag", width = 10, height = 7)
}

# ── 4. Overlay: TRIP4 vs TRIP4+RA (hormone effect) ─────────────────────────
# This overlay compares TRIP4 proximity labeling WITH vs. WITHOUT Retinoic
# Acid (RA) treatment. RA is a signaling molecule (a derivative of vitamin A)
# that activates certain gene programs. By comparing these two conditions,
# we can identify proteins whose association with TRIP4 CHANGES upon RA
# treatment — these might be RA-regulated interactions or RA-dependent
# complex rearrangements.

cat("\n[4/4] RA treatment comparison overlay...\n")
turbo_base <- "turbo_trip4_vs_wt"
turbo_ra <- "turbo_RA_vs_wt"
if (turbo_base %in% names(experiments) && turbo_ra %in% names(experiments)) {
  p <- plot_volcano_overlay(
    experiments[[turbo_base]], "TurboID_TRIP4",
    experiments[[turbo_ra]], "TurboID_TRIP4_RA",
    label_genes = known_interactors,
    title = "TurboID: TRIP4 vs TRIP4 + Retinoic Acid"
  )
  save_figure(p, "volcano_overlay_RA_effect", width = 10, height = 7)
}

# ── COMPLETION BANNER ──────────────────────────────────────────────────────
# Print a final summary so the user knows the script finished and where
# the output files were saved.
cat("\n=========================================\n")
cat(" Volcano plots complete!\n")
# sprintf() formats the string: %s is replaced by FIGURE_DIR (the output path).
cat(sprintf(" Output: %s/\n", FIGURE_DIR))
cat("=========================================\n")
