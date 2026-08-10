# PNP Reaction Diagram — Fix Plan

## Goal
Recreate the original PNP enzyme reaction diagram in LaTeX chemfig.
All vision models must agree it matches the original before showing to Felix.

## Original Diagram (img_048b8b2dc649.png)
- 3D ball-and-stick style
- Top: MESG (purine base LEFT + ribose RIGHT at N9) + Pi
- Arrow: vertical down, labeled "PNP"
- Bottom: 7-methyl-6-thioguanine (thione C=S) + ribose-1-phosphate
- Text: "Absorbance at 360nm"
- Colors: blue N, red O, yellow S, orange P

## Technical Approach

### Strategy: Branch bond from N9 (test_v5 syntax)
The working approach: attach ribose as a branch from N9 using `(-[::0]*5(...))` 
nested inside the N9 atom's substituent list.

**Working purine base (confirmed by vision):**
```latex
\chemfig{*6(-\Natom=C(-NH_2)-\Natom=C(-[6]*5(-\Natom(-CH_3)-C=\Natom-))-C(=\Satom)=)}
```

**V5 approach (compiles, needs vision verification):**
```latex
\chemfig{*6(-\Natom=C(-NH_2)-\Natom=C(-[6]*5(-\Natom(-CH_3)-C=\Natom(-[::0]*5(-\Oatom-C(-OH)-C(-OH)-C(-CH_2OH)-))-))-C(=\Satom^{-})=)}
```

The ribose branches from N9 at angle ::0 (rightward) as a nested *5() ring.

### Alternative approaches if V5 doesn't work:
1. **Two separate chemfig calls in tabular** — base and ribose side by side, with a manual bond drawn between them using Tikz overlay
2. **\schemestart with \chemmove** — draw base and ribose as separate schemes, connect with chemmove
3. **TikZ + chemfig hybrid** — use TikZ nodes to position molecules, chemfig to draw each one

### Full diagram structure:
```latex
\documentclass[border=15pt]{standalone}
\usepackage{chemfig,xcolor,amsmath,array}
% Color definitions
% Atom macros
\begin{document}
\begin{tabular}{c}
  % Reactants: MESG + Pi
  \begin{tabular}{ccc}
    \chemfig{MESG with ribose} & $+$ & $P_i$ \\
    \textbf{MESG} & & \textbf{Phosphate}
  \end{tabular} \\
  $\downarrow$ $\mathrm{PNP}$ \\
  % Products: Free Base + Ribose-1-P
  \begin{tabular}{ccc}
    \chemfig{free base} & $+$ & \chemfig{ribose-1-phosphate} \\
    \textbf{7-Methyl-6-thioguanine} & & \textbf{Ribose-1-phosphate}
  \end{tabular} \\
  \textit{Absorbance at 360\,nm}
\end{tabular}
\end{document}
```

## Worker Tasks

### Task 1: Fix chemfig syntax (worker-balloon or delegate_task leaf)
1. Start with test_v5 syntax (compiles, may need angle adjustments)
2. Compile: `cd /home/c03rad0r/repos/proteomics-pipeline && pdflatex -output-directory=molecular_diagrams molecular_diagrams/pnp_reaction.tex`
3. Render to PNG: `pdftoppm -png -r 200 molecular_diagrams/pnp_reaction.pdf /tmp/pnp_check`
4. Vision check with ALL models:
   ```python
   import requests, base64, os
   key = os.environ['OLLAMA_CLOUD_API_KEY']
   with open("/tmp/pnp_check-1.png", "rb") as f:
       img_b64 = base64.b64encode(f.read()).decode()
   # Test image: /home/c03rad0r/.hermes/profiles/manager/image_cache/img_048b8b2dc649.png
   # Send BOTH images to model, ask: "Do these show the same molecular diagram?"
   resp = requests.post("https://api.ollama.com/api/generate", json={
       "model": "qwen3.5:cloud",
       "prompt": "Compare these two images. Do they show the same PNP reaction diagram with the same molecular structures? List differences.",
       "images": [original_b64, render_b64],
       "stream": False
   }, timeout=300, headers={"Authorization": f"Bearer {key}"})
   ```
5. Run comparison with: qwen3.5:cloud, gemma4:cloud, minimax-m3:cloud
6. If ANY model finds issues, fix and recompile
7. ALL THREE must agree the diagrams match before proceeding

### Task 2: SVG generation
1. Check if pdf2svg installed: `which pdf2svg`
2. If not: `sudo apt install pdf2svg` or `pip install pdf2svg`
3. Alternative: `dvisvgm --pdf --no-fonts molecular_diagrams/pnp_reaction.pdf`
4. Generate: `pdf2svg molecular_diagrams/pnp_reaction.pdf molecular_diagrams/pnp_reaction.svg`
5. Verify SVG opens and shows the diagram

### Task 3: Git commit + push
1. `cd /home/c03rad0r/repos/proteomics-pipeline`
2. `git add molecular_diagrams/pnp_reaction.tex molecular_diagrams/pnp_reaction.pdf molecular_diagrams/pnp_reaction.svg`
3. `git commit -m "fix: PNP diagram matches original — all vision models verified"`
4. `git push origin main`

## Quality Gates
1. PDF compiles without errors
2. SVG generated successfully
3. ALL 3 vision models (qwen3.5:cloud, gemma4:cloud, minimax-m3:cloud) confirm:
   - Purine bicyclic rings closed (both reactant and product)
   - Ribose connected to base at N9
   - No overlapping or broken bonds
   - All labels present
   - Matches original diagram layout
4. Git pushed to remote

## Vision Verification Protocol
For each model, send BOTH the original and our render:
```
Prompt: "Compare these two biochemistry diagrams. Image 1 is the original. Image 2 is the recreation. 
Do they show the same PNP reaction? List specific differences in:
1. Molecular structures (ring closure, atom positions)
2. Layout (reactants, arrow, products)
3. Labels and text
4. Overall similarity (0-100%)"
```
Pass criteria: ALL models report ≥85% similarity AND no structural errors.

## Files
- Original: /home/c03rad0r/.hermes/profiles/manager/image_cache/img_048b8b2dc649.png
- LaTeX: /home/c03rad0r/repos/proteomics-pipeline/molecular_diagrams/pnp_reaction.tex
- Repo: https://github.com/felixfelix-bot/proteomics-pipeline
- Ollama Cloud API: https://api.ollama.com/api/generate
- API key: $OLLAMA_CLOUD_API_KEY
- Available models: qwen3.5:cloud, gemma4:cloud, minimax-m3:cloud, kimi-k2.6:cloud (if not rate-limited)