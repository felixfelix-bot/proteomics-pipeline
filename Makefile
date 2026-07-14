# Makefile for TRIP4/ASCC Proteomics Analysis Pipeline
#
# Usage:
#   make install     - install all R packages
#   make all         - run full pipeline on real data
#   make test        - run pipeline on synthetic test data
#   make volcano     - generate volcano plots only
#   make venn        - generate Venn diagrams only
#   make go          - GO enrichment analysis only
#   make string      - STRING network analysis only
#   make families    - gene family highlighting only
#   make overlap     - cross-experiment overlap only
#   make clean       - remove all generated output
#   make pull        - git pull latest changes
#   make push        - git add, commit, and push (usage: make push m="message")
#   make status      - show git status

# ---- Cross-platform file open command ----
# Uses R's browseURL() which works on Windows (start), macOS (open), and Linux (xdg-open)
# since R is guaranteed to be installed for this pipeline.
OPEN_CMD := $(RSCRIPT) -e "for(f in commandArgs(T)) browseURL(f)"

# ---- Configuration ----
# Rscript path - auto-detected, or override via: make RSCRIPT=/path/to/Rscript all
# On Windows: R is typically at "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
# On Linux/macOS: usually on PATH, or in conda/mamba env

# Windows: hardcode the known path (no Unix shell probing needed)
ifeq ($(OS),Windows_NT)
	RSCRIPT := "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
else
# Linux/macOS: try to find Rscript automatically
	RSCRIPT ?= $(shell \
		if command -v Rscript >/dev/null 2>&1; then \
			command -v Rscript; \
		elif [ -f "/usr/lib/R/bin/Rscript" ]; then \
			echo "/usr/lib/R/bin/Rscript"; \
		elif ls /home/*/.local/bin/Rscript >/dev/null 2>&1; then \
			ls /home/*/.local/bin/Rscript | head -1; \
		else \
			echo "Rscript"; \
		fi)
endif

# Pipeline scripts
SCRIPTS    := R/00_install_packages.R R/01_config.R R/utils.R R/02_volcano_plots.R \
              R/03_venn_diagrams.R R/04_go_enrichment.R R/05_string_network.R \
              R/06_gene_families.R R/07_overlap_analysis.R

FIGURE_DIR := output/figures
TABLE_DIR  := output/tables

# ---- Default target ----
.DEFAULT_GOAL := help

.PHONY: help install all test volcano venn go string families overlap clean pull push status check \
        targeted-volcano flagip-volcano targeted-venn targeted-go go-network-volcano context list-data chx-crac-analysis venn-examples venn-label-examples \
        string-network crac-network venn-overflow-examples gsea pathway-network lydia-volcano chx-common \
        all-volcano all-venn all-targeted targeted-plots \
        volcano-plots \
        open-targeted-volcano open-flagip-volcano open-targeted-venn open-targeted-go open-chx-crac open-venn-examples open-venn-label-examples \
        open-string-network open-crac-network open-venn-overflow-examples open-gsea open-pathway-network \
        open-targeted-plots open-all \
        chx-kegg-crac

# ---- Help ----
help:
	@echo "TRIP4/ASCC Proteomics Analysis Pipeline"
	@echo ""
	@echo "First time setup:"
	@echo "  make install     Install all required R packages"
	@echo "  make check       Verify all packages are installed"
	@echo ""
	@echo "Running the pipeline:"
	@echo "  make all         Run full pipeline on real data in data/"
	@echo "  make test        Run pipeline on synthetic test data"
	@echo ""
	@echo "Individual steps:"
	@echo "  make volcano     Volcano plots with interactor highlighting"
	@echo "  make venn        Venn diagrams and set overlap extraction"
	@echo "  make go          GO enrichment analysis (clusterProfiler)"
	@echo "  make string      STRING protein interaction network"
	@echo "  make families    Gene family highlighting (GPATCH/DHX/DDX/LARP)"
	@echo "  make overlap     Cross-experiment overlap analysis"
	@echo ""
	@echo "Targeted plots (researcher-specified):"
	@echo "  make targeted-volcano  3 custom volcanos: TRIP4 vs WT, RA effects"
	@echo "  make flagip-volcano    TurboID volcano with common C/N-Flag hits"
	@echo "  make flagip-validated-go  GO+KEGG on proteins validated by BOTH C+N Flag"
	@echo "  make targeted-venn     RA effect Venn + TurboID vs C/N-Flag Venn"
	@echo "  make targeted-go       GO enrichment on 3 targeted gene sets"
	@echo "  make context          Print data context (file info, overlaps, gene names)"
	@echo "  make chx-crac-analysis  CHX/DMSO volcanos + GO, CRAC volcano + GO"
	@echo "  make list-data        List ALL CSV files (paths, sizes, headers — no data values)"
	@echo "  make venn-examples    Generate 3 Venn diagram style variants (pick your favorite)"
	@echo "  make venn-label-examples  3 Venn label position variants"
	@echo "  make string-network   STRING network + GO by membership (Lydia style)"
	@echo "  make crac-network     STRING network for CRAC RNA interactome"
	@echo "  make gsea             GSEA (ranked gene list, Lydia style)"
	@echo "  make pathway-network  STRING network maps for enriched pathways"
	@echo "  make lydia-volcano   Lydia-style volcano w/ STRING network overlay"
	@echo "  make chx-common      CHX/DMSO common analysis (enriched/depleted + STRING + GO)"
	@echo "  make ra-common        RA common analysis across both concentrations"
	@echo "  make bidirectional-go Bidirectional GO dot plot (up=right, down=left)"
	@echo "  make shinygo-compare  Compare ShinyGO export with our STRING results"
	@echo "  make chx-kegg-crac   KEGG on CHX + CRAC/TurboID overlap + labeled volcano"
	@echo ""
	@echo "Group targets:"
	@echo "  make volcano-plots      All volcano plots (targeted + flagIP + Lydia)"
	@echo "  make all-volcano         All volcano plots (incl. Lydia network)"
	@echo "  make aruna-all         Run EVERYTHING (fast then slow)"
	@echo "  make all-plots         Run everything + poster figures (ONE COMMAND)"
	@echo "  make poster            Generate ALL figures + compile LaTeX poster"
	@echo "  make poster-figures    Generate poster-styled figures only (no LaTeX)"
	@echo "  make aruna-fast        Quick targets only (volcanos, Venns, GO)"
	@echo "  make aruna-slow        Slow targets only (STRING networks, GSEA)"
	@echo "  make all-venn          All Venn diagrams"
	@echo "  make all-targeted      All targeted analysis"
	@echo "  make targeted-plots    Alias for make all-targeted"
	@echo ""
	@echo "Open output PDFs (after running the corresponding make target):"
	@echo "  make open-targeted-volcano    Open targeted volcano PDFs"
	@echo "  make open-flagip-volcano      Open Flag IP volcano PDF"
	@echo "  make open-targeted-venn       Open Venn diagram PDFs"
	@echo "  make open-targeted-go         Open all GO enrichment PDFs"
	@echo "  make open-venn-examples       Open 3 Venn diagram example PDFs"
	@echo "  make open-targeted-plots      Open ALL targeted plot PDFs"
	@echo "  make open-string-network      Open STRING network PDFs"
	@echo "  make open-crac-network        Open CRAC STRING network PDFs"
	@echo "  make open-gsea                Open GSEA enrichment plots"
	@echo "  make open-pathway-network      Open pathway STRING network maps"
	@echo "  make open-all                 Open every output PDF"
	@echo ""
	@echo "Git operations:"
	@echo "  make pull        Pull latest changes from GitHub"
	@echo "  make push m=msg  Stage, commit with message, and push"
	@echo "  make status      Show git working tree status"
	@echo ""
	@echo "Cleanup:"
	@echo "  make clean       Remove all generated figures and tables"

# ---- Install ----
install: ## Install all required R packages
	$(RSCRIPT) R/00_install_packages.R

check: ## Verify all packages can be loaded
	$(RSCRIPT) R/check_packages.R

# ---- Clean old output files ----
# DISABLED: output files now persist across commits (no hash-based deletion).
# To re-enable, uncomment the Rscript line below.
clean-old:
	@echo "  clean-old: skipped (output persistence enabled)"

# Most analysis targets depend on both check + clean-old
# (clean-old removes stale files from previous runs before generating fresh ones)
all: check clean-old ## Run full pipeline on real data
	$(RSCRIPT) run_all.R

test: ## Run pipeline on synthetic test data
	$(RSCRIPT) run_all.R --test

# ---- Individual analysis steps ----
# Each runs that single module with full logging
volcano: check clean-old ## Generate volcano plots
	$(RSCRIPT) R/run_step.R volcano

venn: check clean-old ## Generate Venn diagrams
	$(RSCRIPT) R/run_step.R venn

go: check clean-old ## GO enrichment analysis
	$(RSCRIPT) R/run_step.R go

string: check clean-old ## STRING network analysis
	$(RSCRIPT) R/run_step.R string

families: check clean-old ## Gene family highlighting
	$(RSCRIPT) R/run_step.R families

overlap: check clean-old ## Cross-experiment overlap
	$(RSCRIPT) R/run_step.R overlap

# ---- Targeted plots (researcher-specified) ----
targeted-volcano: check clean-old ## 3 custom volcano plots: BK467 TRIP4 vs WT, RA effects
	$(RSCRIPT) R/run_step.R targeted_volcanos

flagip-volcano: check clean-old ## TurboID volcano with common C-Flag/N-Flag hits labeled
	$(RSCRIPT) R/run_step.R flagip_volcano

targeted-venn: check clean-old ## 2 Venn diagrams: RA effect + TurboID vs C/N-Flag
	$(RSCRIPT) R/run_step.R targeted_venns

targeted-go: check clean-old ## GO enrichment on 3 targeted gene sets
	$(RSCRIPT) R/run_step.R targeted_go

go-network-volcano: check clean-old ## Volcano highlighting GO network proteins (no labels)
	$(RSCRIPT) R/run_step.R go_network_volcano

context: check clean-old ## Print data context (file structure, overlaps, gene names)
	$(RSCRIPT) R/print_context.R

chx-crac-analysis: check clean-old ## CHX/DMSO volcanos + GO, CRAC volcano + GO
	$(RSCRIPT) R/run_step.R chx_crac_analysis

list-data: ## List ALL CSV files (paths, sizes, headers — no data values)
	$(RSCRIPT) R/list_data.R

headers: ## Print column names of every CSV file
	$(RSCRIPT) R/print_headers.R

venn-examples: check clean-old ## Generate 3 Venn diagram style examples (choose your favorite)
	$(RSCRIPT) R/run_step.R venn_examples

venn-label-examples: check clean-old ## Generate 3 Venn label position variants
	$(RSCRIPT) R/run_step.R venn_label_examples

string-network: check clean-old ## STRING network analysis + GO by network membership
	$(RSCRIPT) R/run_step.R string_network

crac-network: check clean-old ## STRING network for CRAC RNA interactome data
	$(RSCRIPT) R/run_step.R crac_string_network

venn-overflow-examples: check clean-old ## 3 Venn overflow solutions (ext.text, equal circles, boxed labels)
	$(RSCRIPT) R/run_step.R venn_overflow_examples

gsea: check clean-old ## GSEA enrichment analysis (ranked gene list via clusterProfiler)
	$(RSCRIPT) R/run_step.R gsea

pathway-network: check clean-old ## STRING network maps for enriched GO pathways
	$(RSCRIPT) R/run_step.R pathway_network

diagnostics: ## Print structural data summary (counts/ranges only, safe to share)
	@echo "Running diagnostics..."
	@Rscript R/run_step.R diagnostics

chx-volcano-venn: check clean-old ## CHX vs DMSO volcano plot + Venn diagram
	$(RSCRIPT) R/run_step.R chx_volcano_venn

flagip-validated-go: check clean-old ## GO + KEGG on proteins validated by BOTH C-Flag and N-Flag IP
	$(RSCRIPT) R/run_step.R flagip_validated_go

poster-figures: ## Generate poster-styled figures (standardized fonts) in poster/figures/
	$(RSCRIPT) R/run_step.R poster_figures

poster: aruna-fast poster-figures ## Generate ALL poster figures + compile LaTeX poster
	@echo ""
	@echo "========================================="
	@echo " Figures generated. Building poster..."
	@echo "========================================="
	@# Copy any non-ggplot figures (STRING networks) from output/ to poster/
	@for f in output/figures/*string_style*.pdf output/figures/*lydia_volcano*.pdf; do \
		if [ -f "$$f" ]; then \
			cp "$$f" poster/figures/ 2>/dev/null || true; \
			echo "  Copied: $$f"; \
		fi; \
	done
	@# Try to compile LaTeX poster (simple version first — no beamerposter needed)
	@if command -v pdflatex >/dev/null 2>&1; then \
		echo "  Compiling poster_simple.tex..."; \
		cd poster && pdflatex -interaction=nonstopmode poster_simple.tex && \
			pdflatex -interaction=nonstopmode poster_simple.tex; \
		echo ""; \
		echo "========================================="; \
		echo " POSTER COMPLETE: poster/poster_simple.pdf"; \
		echo "========================================="; \
	else \
		echo ""; \
		echo "  NOTE: pdflatex not found. Install MiKTeX (Windows) or TeX Live (Linux/Mac)."; \
		echo "  Figures are ready in poster/figures/."; \
		echo "  To compile manually: cd poster && pdflatex poster_simple.tex"; \
	fi

chx-kegg-crac: check clean-old ## KEGG on CHX data + CRAC/TurboID overlap + volcano
	$(RSCRIPT) R/run_step.R chx_kegg_crac_overlap

string-style-network: check clean-old ## STRING website-style bubble network for RA common proteins
	$(RSCRIPT) R/run_step.R string_style_network

bidirectional-go-ra: check clean-old ## Bidirectional GO: -RA (left) vs +RA (right)
	$(RSCRIPT) R/run_step.R bidirectional_go_ra

network-go: check clean-old ## GO enrichment split by in-network vs not-in-network
	$(RSCRIPT) R/run_step.R network_go

lydia-volcano: check clean-old ## Lydia-style volcano with STRING network overlay
	$(RSCRIPT) R/run_step.R lydia_network_volcano

chx-common: check clean-old ## CHX/DMSO common analysis — enriched/depleted + STRING + GO
	$(RSCRIPT) R/run_step.R chx_common_analysis

ra-common: check clean-old ## RA common protein analysis — enriched/depleted across RA02+RA04
	$(RSCRIPT) R/run_step.R ra_common

ra-common-network: ra-common ## Alias — same as ra-common
	@echo "RA common analysis complete (includes STRING networks)."

bidirectional-go: check clean-old ## Bidirectional GO dot plot — up right, down left
	$(RSCRIPT) R/run_step.R bidirectional_go

shinygo-compare: check clean-old ## Compare ShinyGO export with STRING pipeline
	$(RSCRIPT) R/run_step.R shinygo_comparison

# ---- Group targets (run multiple steps at once) ----
all-volcano: targeted-volcano flagip-volcano lydia-volcano ## All volcano plot targets
	@echo ""
	@echo "All volcano plots complete."

volcano-plots: all-volcano ## Alias: all volcano plots (targeted + flagip + lydia)
	@echo ""

all-venn: venn targeted-venn ## All Venn diagram targets
	@echo ""
	@echo "All Venn diagrams complete."

all-targeted: targeted-volcano flagip-volcano targeted-venn targeted-go ## All targeted plots
	@echo ""
	@echo "All targeted analysis complete."

# ---- Timing wrapper ----
# Prints elapsed time for each target so Aruna can see how long things take.
# Usage: make timed-target STEP=volcano
timed-target:
	@START=$$(date +%s) ; \
	echo "  >>> Running: $(STEP)" ; \
	$(RSCRIPT) R/run_step.R $(STEP) ; \
	RC=$$? ; \
	END=$$(date +%s) ; \
	ELAPSED=$$(($$END - $$START)) ; \
	echo "  <<< $(STEP) finished in $${ELAPSED}s (exit $$RC)" ; \
	exit $$RC

# ---- FAST targets (seconds, no STRINGdb/network downloads) ----
# These run in under a minute each. Bundle them so Aruna can review
# outputs quickly while slower network analyses are still running.
aruna-fast: targeted-volcano flagip-volcano targeted-venn targeted-go \
            venn bidirectional-go flagip-validated-go poster-figures lydia-volcano
	@echo ""
	@echo "================================================"
	@echo " FAST ANALYSIS + POSTER FIGURES COMPLETE"
	@echo " output/figures/  — pipeline figures (Arial 14)"
	@echo " poster/figures/  — poster-styled figures (PDF+PNG)"
	@echo " Now run: make aruna-slow"
	@echo "================================================"

# ---- SLOW targets (minutes — STRINGdb + network building) ----
aruna-slow: lydia-volcano string-network crac-network \
            ra-common chx-common gsea
	@echo ""
	@echo "================================================"
	@echo " SLOW ANALYSIS COMPLETE — check output/figures/"
	@echo "================================================"

# ---- Run EVERYTHING (fast first, then slow) ----
aruna-all: aruna-fast aruna-slow chx-kegg-crac chx-volcano-venn
	@echo ""
	@echo "================================================"
	@echo " ALL ANALYSIS COMPLETE — check output/figures/"
	@echo " Now run: make poster-figures"
	@echo "================================================"

# ---- Generate ALL plots including poster figures (one command) ----
all-plots: aruna-all poster-figures
	@echo ""
	@echo "================================================"
	@echo " ALL PLOTS GENERATED"
	@echo " output/figures/  — pipeline figures (PNG+PDF)"
	@echo " poster/figures/  — poster figures (Arial 24, PDF+PNG)"
	@echo " output/tables/   — all result tables (CSV)"
	@echo "================================================"

# ---- Convenience alias ----
targeted-plots: all-targeted ## Run all targeted plots (volcano + flagIP + venn + GO)
	@echo ""
	@echo "All targeted plots complete."

# ======================================================================
# OPEN TARGETS — Open output PDFs with default system viewer
# Uses R's browseURL() for cross-platform support (Windows/macOS/Linux)

open-targeted-volcano: ## Open targeted volcano PDFs
	$(OPEN_CMD) $(FIGURE_DIR)/targeted_volcano_BK467_TRIP4_vs_WT.pdf \
	            $(FIGURE_DIR)/targeted_volcano_BK467_RA_effect.pdf \
	            $(FIGURE_DIR)/targeted_volcano_BK504_RA_effect.pdf

open-flagip-volcano: ## Open Flag IP overlap volcano PDF
	$(OPEN_CMD) $(FIGURE_DIR)/flagip_overlap_volcano_BK467_TRIP4_vs_WT.pdf

open-targeted-venn: ## Open targeted Venn diagram PDFs
	$(OPEN_CMD) $(FIGURE_DIR)/targeted_venn_RA_effect_BK467.pdf \
	            $(FIGURE_DIR)/targeted_venn_turboid_flagip.pdf

open-venn-examples: ## Open the 3 Venn diagram example PDFs
	$(OPEN_CMD) $(wildcard $(FIGURE_DIR)/venn_example_*.pdf)

open-venn-label-examples: ## Open Venn label position variant PDFs
	$(OPEN_CMD) $(wildcard $(FIGURE_DIR)/venn_label_v*.pdf)

open-string-network: ## Open STRING network plot + GO by network category
	$(OPEN_CMD) $(wildcard $(FIGURE_DIR)/string_network_*.pdf) \
	            $(wildcard $(FIGURE_DIR)/GO_STRING_*.pdf)

open-crac-network: ## Open CRAC STRING network plot
	$(OPEN_CMD) $(wildcard $(FIGURE_DIR)/crac_string_network_*.pdf)

open-gsea: ## Open GSEA enrichment plots (dotplot, ridgeplot, enrichment)
	$(OPEN_CMD) $(wildcard $(FIGURE_DIR)/GSEA_*.pdf)

open-pathway-network: ## Open pathway STRING network maps
	$(OPEN_CMD) $(wildcard $(FIGURE_DIR)/pathway_net_*.pdf)

open-venn-overflow-examples: ## Open Venn overflow example PDFs
	$(OPEN_CMD) $(wildcard $(FIGURE_DIR)/example*_venn_*.pdf)

open-targeted-go: ## Open all GO enrichment plot PDFs
	$(OPEN_CMD) $(wildcard $(FIGURE_DIR)/targeted_GO_*.pdf)

open-chx-crac: ## Open CHX/DMSO + CRAC volcano and GO PDFs
	$(OPEN_CMD) $(wildcard $(FIGURE_DIR)/volcano_TRIP4_CHX*.pdf) \
	            $(wildcard $(FIGURE_DIR)/volcano_TRIP4_DMSO*.pdf) \
	            $(wildcard $(FIGURE_DIR)/volcano_CRAC*.pdf) \
	            $(wildcard $(FIGURE_DIR)/GO_TRIP4_CHX*.pdf) \
	            $(wildcard $(FIGURE_DIR)/GO_TRIP4_DMSO*.pdf) \
	            $(wildcard $(FIGURE_DIR)/GO_CRAC*.pdf)

open-targeted-plots: ## Open ALL targeted plot PDFs
	$(OPEN_CMD) $(FIGURE_DIR)/targeted_volcano_BK467_TRIP4_vs_WT.pdf \
	            $(FIGURE_DIR)/targeted_volcano_BK467_RA_effect.pdf \
	            $(FIGURE_DIR)/targeted_volcano_BK504_RA_effect.pdf \
	            $(FIGURE_DIR)/flagip_overlap_volcano_BK467_TRIP4_vs_WT.pdf \
	            $(FIGURE_DIR)/targeted_venn_RA_effect_BK467.pdf \
	            $(FIGURE_DIR)/targeted_venn_turboid_flagip.pdf \
	            $(wildcard $(FIGURE_DIR)/targeted_GO_*.pdf)

open-all: ## Open ALL output PDFs (all steps)
	$(OPEN_CMD) $(wildcard $(FIGURE_DIR)/*.pdf)

# ---- Cleanup ----
clean: ## Remove all generated output
	rm -rf $(FIGURE_DIR) $(TABLE_DIR) output/logs
	@echo "Removed: output/ (figures, tables, logs)"

# ---- Git operations ----
pull: ## Pull latest from GitHub
	git pull
	@echo ""
	@echo "  Up to date. Run: make all"

push: ## Stage, commit, and push (usage: make push m="commit message")
	@if [ -z "$(m)" ]; then \
		echo "Usage: make push m=\"commit message\""; \
		exit 1; \
	fi
	git add -A
	git commit -m "$(m)"
	git push
	@echo ""
	@echo "  Pushed to GitHub."

status: ## Show git status
	@git status -s
