# LaTeX Poster Guide — Uniform Font Sizes

## The Problem

When you resize images for a poster (making some bigger, some smaller), the
font sizes inside each image become inconsistent. A gene label readable in one
panel becomes tiny in another.

## The Solution

**Two-part approach: R standardizes fonts → LaTeX controls layout.**

### Part 1: Standardize R Figures

All R plots now use `theme_poster()` from `R/00_theme.R`. This enforces
identical font sizes across every figure:

```r
# In any R script, after creating a plot:
my_plot + theme_poster()
save_poster_figure(my_plot, "figure_name", width=8, height=6)
```

This saves BOTH `.pdf` (vector, for poster) and `.png` (backup).

**Why PDF?** PDF is vector format — text stays sharp at ANY size. When LaTeX
includes a PDF figure and resizes it, the text scales cleanly with no
pixelation.

### Part 2: LaTeX Template

Two templates provided:

| File | Use When | Dependency |
|------|----------|------------|
| `poster/poster.tex` | Professional conference poster (120x90cm) | beamerposter package |
| `poster/poster_simple.tex` | Quick A4 poster (any printer) | Base LaTeX only |

### Key Design Principle

In both templates, ALL figures use the **same width**:

```latex
\includegraphics[width=\linewidth]{figures/my_figure.pdf}
```

This means every figure is scaled to the column width. Since all R figures
use the same base font size (from `theme_poster()`), the fonts appear at
the same physical size in every panel.

**If you need a figure smaller:** reduce the width proportionally:
```latex
\includegraphics[width=0.8\linewidth]{figures/smaller_panel.pdf}
```
The fonts in this panel will be 80% of the others. This is unavoidable when
making one panel physically smaller — but it's proportional and predictable.

### Workflow

1. **Run the pipeline:**
   ```
   make aruna-all
   ```

2. **Copy figures to poster folder:**
   ```
   # On Aruna's machine (Windows):
   mkdir poster\figures
   copy output\figures\*.pdf poster\figures\
   ```

3. **Edit poster content:**
   - Open `poster/poster_simple.tex` in any text editor
   - Change figure paths, text, section titles
   - Add/remove sections as needed

4. **Compile:**
   ```
   # If you have LaTeX installed (MiKTeX/TeX Live):
   cd poster
   pdflatex poster_simple.tex
   pdflatex poster_simple.tex   # run twice for layout
   ```

5. **Result:** `poster/poster_simple.pdf` — ready to print.

### Which Figures to Include

Recommended figure set for a TRIP4/ASCC poster:

| Section | Figure | Makefile target |
|---------|--------|----------------|
| TurboID | `lydia_volcano_network.pdf` | `make lydia-volcano` |
| GO | `targeted_GO_turbo_BP_dotplot.pdf` | `make targeted-go` |
| KEGG | `targeted_KEGG_turbo_dotplot.pdf` | `make targeted-go` |
| Flag IP | `flagip_volcano_overlap.pdf` | `make flagip-volcano` |
| Venn | `targeted_venn_turbo_vs_flag.pdf` | `make targeted-venn` |
| STRING | `ra_string_style_network.pdf` | `make string-style-network` |
| Validated | `flagip_GO_validated_both_BP_dotplot.pdf` | `make flagip-validated-go` |

### Font Size Cheat Sheet

All sizes defined in `R/00_theme.R` → `POSTER_FONTS`:

| Element | Size (pt) |
|---------|-----------|
| Title | 14 |
| Axis title | 12 |
| Axis text | 10 |
| Legend title | 11 |
| Legend text | 9 |
| Annotations | 8 |

To change ALL fonts at once, edit the numbers in `POSTER_FONTS` in `R/00_theme.R`.
Every figure updates automatically.
