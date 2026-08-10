# PNP Diagram Fix Plan v2 — Consultant-Approved Approach

## Problem Analysis
The core issue: chemfig's `*5()` ring syntax inside a branch from a ring atom
breaks the parent ring's closure. This is a known chemfig limitation when nesting
ring commands inside ring branches.

## Technical Solution: Two-Phase Approach

### Phase 1: Draw MESG as separate molecules with manual bond (SIMPLEST)
Draw purine base and ribose as two separate `\chemfig{}` calls inside a single
tabular cell, with a horizontal rule/bond drawn between them using `\hspace`
and a rule, or using TikZ overlay.

The purine base renders correctly standalone:
```latex
\chemfig{*6(-\Natom=C(-NH_2)-\Natom=C(-[6]*5(-\Natom(-CH_3)-C=\Natom-))-C(=\Satom^{-})=)}
```

The ribose renders correctly standalone:
```latex
\chemfig{[:-18]*5(-\Oatom-C(-[::36]CH_2OH)-C(-[::36]OH)-C(-[::36]OH)-)}
```

Connection: place them side by side in a tabular with a drawn bond:
```latex
\begin{tabular}{cc}
  \chemfig{purine} &
  \chemfig{ribose}
\end{tabular}
```

### Phase 2: If Phase 1 insufficient, use full TikZ drawing
Draw the entire diagram in TikZ with manual node placement and edge drawing.
This gives full control over bond angles and ring closure. Use chemfig only for
individual molecule labels if needed.

## Worker Task Specification

### Task: Fix PNP diagram and verify with vision models

**File:** `/home/c03rad0r/repos/proteomics-pipeline/molecular_diagrams/pnp_reaction.tex`
**Original image:** `/home/c03rad0r/.hermes/profiles/manager/image_cache/img_048b8b2dc649.png`

**Approach:**
1. Rewrite pnp_reaction.tex using the tabular approach (separate chemfig calls for base and ribose)
2. Compile: `pdflatex -output-directory=molecular_diagrams molecular_diagrams/pnp_reaction.tex`
3. Render to PNG: `pdftoppm -png -r 200 molecular_diagrams/pnp_reaction.pdf /tmp/pnp_check`
4. Vision verify with Ollama Cloud API (NOT localhost — returns 401):
   ```python
   import requests, base64, os
   key = os.environ['OLLAMA_CLOUD_API_KEY']
   # Send BOTH original + render to each model
   resp = requests.post("https://api.ollama.com/api/generate", json={
       "model": "qwen3.5:cloud",  # or gemma4:cloud, minimax-m3:cloud
       "prompt": "Compare these two diagrams. Rate similarity 0-100%.",
       "images": [orig_b64, render_b64],
       "stream": False
   }, timeout=300, headers={"Authorization": f"Bearer {key}"})
   ```
5. If similarity < 85%, try TikZ approach
6. If TikZ approach needed, draw molecules as TikZ nodes with manual edges
7. Generate SVG: `dvisvgm --pdf --no-fonts molecular_diagrams/pnp_reaction.pdf --output=molecular_diagrams/pnp_reaction.svg`
8. Git commit + push

**Quality Gates:**
- PDF compiles without errors
- SVG generated
- At least 2 of 3 vision models report ≥85% similarity
- Git committed and pushed to origin

**Important Notes:**
- Ollama Cloud API can be slow/rate-limited. Use 300s timeout. Space calls by 10s.
- If a model times out, try the next one. Don't wait forever.
- The purine base syntax `*6(-N=C(-NH2)-N=C(-[6]*5(-N(-CH3)-C=N-))-C(=S)=)` is CONFIRMED working by vision models.
- The ribose-1-phosphate syntax is CONFIRMED working.
- ONLY the MESG (base+ribose connection) needs fixing.
- The product base (7-methyl-6-thioguanine) renders correctly.
- Document class: standalone (no \begin{center} — use \centering or tabular)