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

# ---- Configuration ----
# Rscript path - auto-detected, or override via: make RSCRIPT=/path/to/Rscript all
# On Windows: R is typically at "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
# On Linux/macOS: usually on PATH, or in conda/mamba env

# Try to find Rscript automatically
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

# Windows fallback (checked at runtime via make)
ifeq ($(OS),Windows_NT)
	ifeq ($(shell command -v Rscript 2>/dev/null),)
		RSCRIPT := "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
	endif
endif

# Pipeline scripts
SCRIPTS    := R/00_install_packages.R R/01_config.R R/utils.R R/02_volcano_plots.R \
              R/03_venn_diagrams.R R/04_go_enrichment.R R/05_string_network.R \
              R/06_gene_families.R R/07_overlap_analysis.R

FIGURE_DIR := output/figures
TABLE_DIR  := output/tables

# ---- Default target ----
.DEFAULT_GOAL := help

.PHONY: help install all test volcano venn go string families overlap clean pull push status check

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
