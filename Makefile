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
        open-targeted-volcano open-flagip-volcano open-targeted-venn open-targeted-go open-chx-crac open-venn-examples open-venn-label-examples \
        open-string-network open-crac-network open-venn-overflow-examples open-gsea open-pathway-network \
        open-targeted-plots open-all

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
	@echo ""
	@echo "Group targets:"
	@echo "  make all-volcano       All volcano plots (incl. Lydia network)"
	@echo "  make aruna-all         Run EVERYTHING (all volcanos, GO, networks, RA, CHX, GSEA)"
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

# ---- Full pipeline ----
all: check ## Run full pipeline on real data
	$(RSCRIPT) run_all.R

test: ## Run pipeline on synthetic test data
	$(RSCRIPT) run_all.R --test

# ---- Individual analysis steps ----
# Each runs that single module with full logging
volcano: check ## Generate volcano plots
	$(RSCRIPT) R/run_step.R volcano

venn: check ## Generate Venn diagrams
	$(RSCRIPT) R/run_step.R venn

go: check ## GO enrichment analysis
	$(RSCRIPT) R/run_step.R go

string: check ## STRING network analysis
	$(RSCRIPT) R/run_step.R string

families: check ## Gene family highlighting
	$(RSCRIPT) R/run_step.R families

overlap: check ## Cross-experiment overlap
	$(RSCRIPT) R/run_step.R overlap

# ---- Targeted plots (researcher-specified) ----
targeted-volcano: check ## 3 custom volcano plots: BK467 TRIP4 vs WT, RA effects
	$(RSCRIPT) R/run_step.R targeted_volcanos

flagip-volcano: check ## TurboID volcano with common C-Flag/N-Flag hits labeled
	$(RSCRIPT) R/run_step.R flagip_volcano

targeted-venn: check ## 2 Venn diagrams: RA effect + TurboID vs C/N-Flag
	$(RSCRIPT) R/run_step.R targeted_venns

targeted-go: check ## GO enrichment on 3 targeted gene sets
	$(RSCRIPT) R/run_step.R targeted_go

go-network-volcano: check ## Volcano highlighting GO network proteins (no labels)
	$(RSCRIPT) R/run_step.R go_network_volcano

context: check ## Print data context (file structure, overlaps, gene names)
	$(RSCRIPT) R/print_context.R

chx-crac-analysis: check ## CHX/DMSO volcanos + GO, CRAC volcano + GO
	$(RSCRIPT) R/run_step.R chx_crac_analysis

list-data: ## List ALL CSV files (paths, sizes, headers — no data values)
	$(RSCRIPT) R/list_data.R

headers: ## Print column names of every CSV file
	$(RSCRIPT) R/print_headers.R

venn-examples: check ## Generate 3 Venn diagram style examples (choose your favorite)
	$(RSCRIPT) R/run_step.R venn_examples

venn-label-examples: check ## Generate 3 Venn label position variants
	$(RSCRIPT) R/run_step.R venn_label_examples

string-network: check ## STRING network analysis + GO by network membership
	$(RSCRIPT) R/run_step.R string_network

crac-network: check ## STRING network for CRAC RNA interactome data
	$(RSCRIPT) R/run_step.R crac_string_network

venn-overflow-examples: check ## 3 Venn overflow solutions (ext.text, equal circles, boxed labels)
	$(RSCRIPT) R/run_step.R venn_overflow_examples

gsea: check ## GSEA enrichment analysis (ranked gene list via clusterProfiler)
	$(RSCRIPT) R/run_step.R gsea

pathway-network: check ## STRING network maps for enriched GO pathways
	$(RSCRIPT) R/run_step.R pathway_network

lydia-volcano: check ## Lydia-style volcano with STRING network overlay
	$(RSCRIPT) R/run_step.R lydia_network_volcano

chx-common: check ## CHX/DMSO common analysis — enriched/depleted + STRING + GO
	$(RSCRIPT) R/run_step.R chx_common_analysis

ra-common: check ## RA common protein analysis — enriched/depleted across RA02+RA04
	$(RSCRIPT) R/run_step.R ra_common

ra-common-network: ra-common ## Alias — same as ra-common
	@echo "RA common analysis complete (includes STRING networks)."

bidirectional-go: check ## Bidirectional GO dot plot — up right, down left
	$(RSCRIPT) R/run_step.R bidirectional_go

shinygo-compare: check ## Compare ShinyGO export with STRING pipeline
	$(RSCRIPT) R/run_step.R shinygo_comparison

# ---- Group targets (run multiple steps at once) ----
all-volcano: targeted-volcano flagip-volcano lydia-volcano ## All volcano plot targets
	@echo ""
	@echo "All volcano plots complete."

all-venn: venn targeted-venn ## All Venn diagram targets
	@echo ""
	@echo "All Venn diagrams complete."

all-targeted: targeted-volcano flagip-volcano targeted-venn targeted-go ## All targeted plots
	@echo ""
	@echo "All targeted analysis complete."

# ---- Run EVERYTHING from Aruna's July 12 feedback ----
aruna-all: all-volcano all-venn targeted-go string-network crac-network \
           ra-common bidirectional-go chx-common gsea
	@echo ""
	@echo "================================================"
	@echo " ALL ANALYSIS COMPLETE — check output/figures/"
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
