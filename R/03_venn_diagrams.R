##############################################################################
# 03_venn_diagrams.R
#
# WHAT THIS SCRIPT DOES (the big picture, for non-R users):
# --------------------
# A Venn diagram is a picture showing how much two (or more) lists of proteins
# OVERLAP with each other. In our proteomics study, we ran several experiments
# to find proteins that interact with TRIP4:
#   - TurboID:   a proximity-labeling experiment (catches nearby proteins)
#   - Flag IP:   an antibody pull-down experiment (catches direct binders)
#   - CRAC RNA:  proteins found associated with RNA
#
# If a protein shows up in BOTH TurboID AND Flag IP, that is strong evidence
# it is a real TRIP4 interaction partner. Proteins unique to only one
# experiment might be experiment-specific (artifacts or context-dependent).
#
# The Venn diagram shows:
#   - The overlap (shared proteins) in the middle
#   - Proteins unique to each experiment on the sides
#
# This script:
#   1. Builds 2-set and 3-set Venn diagrams (using the ggVennDiagram package)
#   2. Extracts the shared/unique protein lists and saves them to CSV files
#      (so we can feed them into GO enrichment analysis in script 04)
#   3. Builds an "UpSet plot" — a more powerful alternative for >4 sets,
#      because traditional Venn diagrams get unreadable with many circles
#
# Usage:
#   source("R/01_config.R")
#   source("R/utils.R")
#   source("R/03_venn_diagrams.R")
##############################################################################

# cat() prints text to the screen — like print() in Python.
# \n means "new line" so this draws a header banner.
cat("\n=========================================\n")
cat(" Venn Diagram Analysis\n")
cat("=========================================\n\n")

# ---- Load all experiment data ----
# load_all_experiments() is a custom helper defined in utils.R.
# It reads all our CSV data files and returns a "list" in R.
# A "list" is like a dictionary in Python: a collection of named items.
# Here each item is a data frame (a table — like a pandas DataFrame in Python)
# for one experiment.
experiments <- load_all_experiments()

# ---- Extract significant gene sets ----
# We need to pull out just the SIGNIFICANT proteins from each experiment
# (significant = passed statistical thresholds, defined in utils.R).
cat("Extracting significant gene sets...\n")

# lapply() is one of the most important R functions to know.
# Think of it EXACTLY like map() in Python or a for-each loop:
#
#   Python:  gene_sets = [get_significant_genes(df) for df in experiments]
#   R:       gene_sets  = lapply(experiments, function(df) { ... })
#
# lapply() takes two arguments:
#   1. A list to iterate over (here: experiments)
#   2. A function to apply to EACH element of the list
# It returns a NEW list where each element is the result of the function.
#
# Here, the function(df) takes one data frame (df) and:
#   - Calls get_significant_genes(df) to get the protein names
#   - Prints how many it found (cat() = print to screen)
#   - Returns those gene names
gene_sets <- lapply(experiments, function(df) {
  genes <- get_significant_genes(df)
  cat(sprintf("  %d significant genes\n", length(genes)))
  return(genes)
})

# Give friendly names for display
# names(gene_sets) returns the names of the list items, e.g. "turboid", "flag_ip".
set_names <- names(gene_sets)

# gsub() is "global substitute" — it replaces text using pattern matching.
# Here we replace every underscore "_" with a space " ".
# So "flag_ip" becomes "flag ip" for nicer-looking labels.
# gsub() is like Python's .replace() or re.sub().
display_names <- gsub("_", " ", set_names)

# =====================================================================
# 1. Venn diagram: TurboID vs Flag IP
# =====================================================================
# This is the most important comparison: do the two protein-finding methods
# agree on which proteins bind TRIP4?
cat("\n[1/2] TurboID vs Flag IP Venn diagram...\n")

# Check that BOTH "turboid" AND "flag_ip" exist in our list of gene sets.
# %in% checks membership (like "in" in Python).
# && means "logical AND" — both conditions must be true.
if ("turboid" %in% names(gene_sets) && "flag_ip" %in% names(gene_sets)) {

  # Build a named list for the Venn diagram.
  # The names ("TurboID", "Flag IP") will appear as labels on the diagram.
  # The $ operator accesses a named element of a list:
  #   gene_sets$turboid  is the same as  gene_sets[["turboid"]]
  venn_list_1 <- list(
    "TurboID"  = gene_sets$turboid,
    "Flag IP"  = gene_sets$flag_ip
  )

  # ggVennDiagram is an R package for drawing Venn diagrams.
  # The "::" syntax means "use the ggVennDiagram function FROM the
  # ggVennDiagram package" — it makes the package name explicit.
  #
  # Parameters:
  #   label = "count"    -> show protein counts in each region of the diagram
  #   set_size = 5       -> font size for the set names (TurboID, Flag IP)
  #   label_size = 4     -> font size for the count numbers
  #
  # The "+" operator adds layers/styling, like in Python's matplotlib.
  # This is the grammar-of-graphics style used by ggplot2:
  #   - scale_fill_gradient: color the overlap regions with a blue gradient
  #       (low = light blue, high = dark blue = more proteins overlap)
  #   - ggtitle: add a title
  #   - theme(...): center the title and make it bold
  p_venn1 <- ggVennDiagram::ggVennDiagram(
    venn_list_1,
    label = "count",
    set_size = 5,
    label_size = 4
  ) +
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggplot2::ggtitle("TurboID vs Flag Co-IP Overlap") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))

  # save_figure() is our custom helper (from utils.R) that saves the plot
  # as both PNG and PDF into the FIGURE_DIR folder.
  save_figure(p_venn1, "venn_turboid_vs_flagip", width = 6, height = 6)

  # ---- Extract the actual overlap sets ----
  # These are R's set operations — they work on vectors of strings (gene names).
  # They are identical to Python's set operations:
  #
  #   intersect(A, B)  ->  A ∩ B   (proteins in BOTH lists)  -> Python: set(A) & set(B)
  #   setdiff(A, B)    ->  A − B   (proteins in A but NOT B) -> Python: set(A) - set(B)
  #   union(A, B)      ->  A ∪ B   (proteins in either list) -> Python: set(A) | set(B)
  #
  # Biologically:
  #   - turbo_flag_shared: proteins found by BOTH methods (high confidence!)
  #   - turbo_only:        proteins only TurboID caught (maybe weak/transient)
  #   - flag_only:         proteins only Flag IP caught (maybe IP-specific)
  turbo_flag_shared <- intersect(gene_sets$turboid, gene_sets$flag_ip)
  turbo_only <- setdiff(gene_sets$turboid, gene_sets$flag_ip)
  flag_only <- setdiff(gene_sets$flag_ip, gene_sets$turboid)

  # sprintf() formats text (like Python's f-strings or % formatting).
  # %d is a placeholder for an integer.
  cat(sprintf("  Shared (TurboID ∩ Flag IP): %d proteins\n", length(turbo_flag_shared)))
  cat(sprintf("  Unique to TurboID: %d proteins\n", length(turbo_only)))
  cat(sprintf("  Unique to Flag IP: %d proteins\n", length(flag_only)))

  # ---- Save extracted sets to CSV ----
  # data.frame() creates a table (like a spreadsheet with columns).
  # We make a 2-column table: gene name + which category it belongs to.
  #
  # save_table() is our custom helper (from utils.R) that writes the table
  # to a CSV file in the TABLE_DIR folder.
  # CSV = "Comma-Separated Values" — a universal spreadsheet format you can
  # open in Excel, Google Sheets, or any text editor.
  #
  # These saved lists are useful for:
  #   - Sharing which proteins were unique vs shared
  #   - Feeding into GO enrichment analysis (script 04)
  save_table(data.frame(gene = turbo_flag_shared, category = "TurboID_AND_FlagIP"),
             "overlap_turboid_flagip")
  save_table(data.frame(gene = turbo_only, category = "TurboID_only"),
             "unique_turboid")
  save_table(data.frame(gene = flag_only, category = "FlagIP_only"),
             "unique_flagip")
}

# =====================================================================
# 2. Venn diagram: Protein Interactome vs CRAC RNA Interactome
# =====================================================================
# Now we add a third experiment (CRAC RNA) to the comparison.
# This shows whether the same proteins that bind TRIP4 also associate with RNA.
cat("\n[2/2] Protein Interactome vs CRAC RNA Venn diagram...\n")

# For 3+ sets, we can do a combined Venn
# Check whether the "crac_rna" experiment exists in our data.
if ("crac_rna" %in% names(gene_sets)) {
  # Figure out which of the 3 possible experiments we actually have.
  # intersect(c(...), names(...)) returns only the names that exist in BOTH lists.
  # So this gives us a subset of {turboid, flag_ip, crac_rna} that we have data for.
  sets_available <- intersect(c("turboid", "flag_ip", "crac_rna"), names(gene_sets))

  if (length(sets_available) >= 3) {
    # We have all 3 experiments — build a 3-circle Venn diagram.
    # A 3-set Venn has 7 regions: 3 unique-only, 3 pairwise-overlaps, and 1 all-three.
    venn_list_2 <- list(
      "TurboID"   = gene_sets$turboid,
      "Flag IP"   = gene_sets$flag_ip,
      "CRAC RNA"  = gene_sets$crac_rna
    )

    # Same ggVennDiagram approach as before, just with 3 sets.
    # Set sizes and label sizes are slightly smaller to fit 3 circles.
    p_venn2 <- ggVennDiagram::ggVennDiagram(
      venn_list_2,
      label = "count",
      set_size = 4,
      label_size = 3.5
    ) +
      ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
      ggplot2::ggtitle("Protein Interactome vs CRAC RNA Interactome") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))

    save_figure(p_venn2, "venn_turboid_flagip_crac", width = 7, height = 7)

    # ---- Extract overlaps for the 3-set case ----
    # Reduce() is a powerful R function that "reduces" a list to a single value
    # by applying a binary function (a function that takes 2 arguments) repeatedly.
    #
    # Think of it like this in Python:
    #   from functools import reduce
    #   all_three = reduce(lambda a, b: set(a) & set(b), [list1, list2, list3])
    #
    # Reduce(intersect, list(A, B, C)) does:
    #   Step 1:  intersect(A, B)        -> temp
    #   Step 2:  intersect(temp, C)     -> final result
    #
    # So it finds proteins present in ALL THREE experiments at once.
    # Biologically: these are the highest-confidence TRIP4 partners.
    all_three <- Reduce(intersect, gene_sets[sets_available])

    # Pairwise overlaps for completeness.
    turbo_crac <- intersect(gene_sets$turboid, gene_sets$crac_rna)
    flag_crac <- intersect(gene_sets$flag_ip, gene_sets$crac_rna)

    cat(sprintf("  All three overlap: %d proteins\n", length(all_three)))
    cat(sprintf("  TurboID ∩ CRAC: %d proteins\n", length(turbo_crac)))
    cat(sprintf("  Flag IP ∩ CRAC: %d proteins\n", length(flag_crac)))

    # Save each overlap set to its own CSV file.
    save_table(data.frame(gene = all_three, category = "ALL_THREE_overlap"),
               "overlap_all_three")
    save_table(data.frame(gene = turbo_crac, category = "TurboID_AND_CRAC"),
               "overlap_turboid_crac")
    save_table(data.frame(gene = flag_crac, category = "FlagIP_AND_CRAC"),
               "overlap_flagip_crac")

  } else if (length(sets_available) == 2) {
    # Only 2 of the 3 experiments are available — fall back to a 2-set Venn.
    # This is the "graceful fallback" branch.
    names_friendly <- gsub("_", " ", sets_available)

    # setNames() renames the elements of a list/vector.
    # Here we take the 2 gene sets we have and give them human-friendly names.
    venn_list_2 <- setNames(gene_sets[sets_available], names_friendly)

    p_venn2 <- ggVennDiagram::ggVennDiagram(
      venn_list_2,
      label = "count",
      set_size = 5,
      label_size = 4
    ) +
      ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
      ggplot2::ggtitle("Protein Interactome vs CRAC RNA Interactome") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))

    save_figure(p_venn2, "venn_protein_vs_rna", width = 6, height = 6)
  }
}

# =====================================================================
# 3. UpSet plot for all experiments (handles >4 sets gracefully)
# =====================================================================
# Why an UpSet plot?
# -----------------
# Traditional Venn diagrams work fine for 2 or 3 sets, but with 4+ sets
# the circles overlap in confusing, hard-to-read ways. An "UpSet plot"
# (from the UpSetR package) is a bar-chart-style alternative that scales
# to any number of sets.
#
# An UpSet plot has two parts:
#   - Left/bottom: a bar for each set showing its total size
#   - Main area:    bars showing the size of each INTERSECTION combination
#                   (dots connected by lines indicate which sets are included)
# For example, a bar with dots under "TurboID" and "Flag IP" connected
# shows how many proteins are in BOTH of those sets.
if (length(gene_sets) >= 2) {
  cat("\n[BONUS] Generating UpSet plot for all experiments...\n")

  # ---- Build a binary membership matrix ----
  # unlist() flattens a list of vectors into one big vector.
  # unique() removes duplicates.
  # So all_genes = every protein seen in ANY experiment, listed once.
  all_genes <- unique(unlist(gene_sets))

  # Start with a data frame that has one column: the gene name.
  # We'll add a 0/1 column for each experiment.
  binary_matrix <- data.frame(gene = all_genes)

  # A for loop iterates over the names of all experiments.
  # for (name in names(gene_sets)) is like Python's:
  #   for name in gene_sets.keys():
  for (name in names(gene_sets)) {
    # %in% checks membership: is each gene in this experiment's gene set?
    # as.integer() converts TRUE/FALSE to 1/0.
    # So this adds a column named after the experiment, with 1 if the gene
    # was significant in that experiment, 0 otherwise.
    # [[name]] accesses the element by name (like dict["key"] in Python).
    binary_matrix[[name]] <- as.integer(all_genes %in% gene_sets[[name]])
  }

  # Convert to the format UpSetR expects
  # Remove the gene-name column (column 1) since UpSetR only wants the 0/1 matrix.
  # [, -1] means "all rows, all columns except the first".
  upset_data <- binary_matrix[, -1]  # remove gene column

  # rownames() sets the row names (like an index in pandas).
  # We use the gene names so rows are labeled.
  rownames(upset_data) <- binary_matrix$gene

  # ---- Draw and save the UpSet plot ----
  # png() opens a PNG graphics device — everything we plot gets written to
  # this file until we call dev.off() (device off) to close it.
  # file.path() joins folder + filename in a cross-platform way.
  #   width/height/res control the output dimensions and quality.
  png(safe_filepath(FIGURE_DIR, "upset_all_experiments", ".png"),
      width = 10, height = 6, units = "in", res = FIG_DPI)

  # print() is needed here because the UpSetR plot is drawn to the active
  # graphics device (the PNG we just opened).
  # print() forces the plot to render.
  print(UpSetR::upset(
    upset_data,
    sets = names(gene_sets),               # which sets (experiments) to show
    order.by = "freq",                     # sort intersection bars by size
    main.bar.color = "grey30",             # color of the main bars
    sets.bar.color = unname(EXPERIMENT_COLORS[names(gene_sets)]),  # per-experiment colors
    point.size = 3.5,                      # dot size in the matrix
    line.size = 1.2,                       # line thickness connecting dots
    text.scale = 1.1                       # font scale factor
  ))

  # dev.off() closes the PNG file and saves it.
  # ALWAYS pair png() with dev.off() — forgetting it leaves the file empty.
  dev.off()
  cat("  Saved: upset_all_experiments.png\n")
}

# Final summary banner showing where outputs were saved.
cat("\n=========================================\n")
cat(" Venn diagrams complete!\n")
cat(sprintf(" Figures: %s/\n", FIGURE_DIR))
cat(sprintf(" Tables:  %s/\n", TABLE_DIR))
cat("=========================================\n")
