# PNP Assay Diagram Fix Plan

## Goal
Reproduce the original PNP phosphorolysis assay diagram in LaTeX/chemfig so that all Ollama Cloud vision models agree it matches the original.

## Original Diagram Description (from 3 vision models)

### Layout
- **Vertical flow**, top to bottom
- **Top**: MESG molecule (left) + text "MESG + Pi" (right). Pi is TEXT, not drawn structure.
- **Middle**: Vertical downward arrow with "PNP" label to the right
- **Bottom**: 7-methyl-6-thioguanine (left) + "+" + ribose-1-phosphate (right)
- **Very bottom**: "Absorbance at 360nm" (no space in 360nm)

### Molecular Structures

#### MESG (top reactant)
- **Purine base** (bicyclic: 6-membered pyrimidine fused to 5-membered imidazole)
  - 6-ring: NH2 at C2, S- (thiolate, yellow) at C6
  - 5-ring: CH3 at N7, N9 connects to sugar
  - Nitrogen atoms colored blue
  - Sulfur colored yellow
- **Ribose sugar** (5-membered furanose ring)
  - Oxygen (red) in ring
  - CH2OH at C5 (wedge bond, up)
  - OH at C2 and C3 (wedge bonds)
  - N-glycosidic bond from C1 to N9 of purine (hashed/dashed wedge)

#### 7-Methyl-6-thioguanine (bottom left product)
- Same purine base but:
  - C6 has C=S (thione, double bond to yellow S) instead of S-
  - N9 has H instead of ribose
  - CH3 still at N7

#### Ribose-1-phosphate (bottom right product)
- Ribose furanose ring (same as MESG sugar)
  - C1 connected via O to phosphate group
  - CH2OH at C5, OH at C2 and C3
  - Stereochemistry: wedges and dashes
- Phosphate group: P (orange) with =O, -O-, -OH, -O-

### Colors
- Nitrogen = Blue (#0066CC)
- Oxygen = Red (#CC0000)
- Sulfur = Yellow/Gold (#D4A017)
- Phosphorus = Orange (#E07000)
- Carbon = Black (implicit)
- Bonds = Black

### Text Labels (MINIMAL)
- "MESG + Pi" (top right, next to molecule)
- "PNP" (right of arrow)
- "Absorbance at 360nm" (bottom)
- NO long chemical names, NO IUPAC names

## Known chemfig Issues to Fix

### Problem 1: Fused rings not actually fused
- `*6(-*5(...))` syntax creates spiro connection (single bond), NOT fused rings
- Need to use proper fusion: branch from C4 of 6-ring, draw 5-ring that shares C4-C5 edge
- Test syntax: `\chemfig{*6(-\Natom=-(-[:30]*5(-\Natom=-\Natom=-))-\Natom=-)}`

### Problem 2: Ribose furanose ring not closing
- `*5(...)` ring syntax not producing visible closed pentagon
- Need proper 5-membered ring with O in ring: `*5(-C-C-C-C-O-)` style

### Problem 3: Isotope labels appearing as text
- `[::0]` and `[30]` showing as literal text in PDF
- These are chemfig angle syntax — need to ensure they're inside bond declarations, not atom labels

### Problem 4: Missing stereochemistry
- Need wedge bonds: `>:[]` (solid wedge, toward viewer)
- Need dash bonds: `:<[]` (hashed wedge, away from viewer)

## Implementation Plan

### Phase 1: Fix chemfig syntax (worker)
1. Test purine bicyclic ring fusion in isolation
2. Test ribose furanose ring closure in isolation
3. Test phosphate group drawing
4. Test wedge/dash stereochemistry bonds
5. Save test file as `test_structures.tex` in molecular_diagrams/
6. Compile and verify each structure renders correctly

### Phase 2: Assemble full diagram (worker)
1. Combine all tested structures into `pnp_reaction.tex`
2. Layout: vertical flow using tabular or tikz
3. Top: MESG molecule + "MESG + Pi" text
4. Middle: tikz vertical arrow with "PNP" label
5. Bottom: two products with "+" between
6. Bottom label: "Absorbance at 360nm"
7. Apply colors: blue N, red O, yellow S, orange P
8. Compile to PDF

### Phase 3: Vision verification (worker)  
1. Render PDF to PNG at 200 DPI
2. Send to Ollama Cloud API (https://ollama.com/v1/chat/completions)
3. Compare original vs render using ALL vision models:
   - qwen3.5:397b
   - gemma4:31b
   - minimax-m3:cloud
4. Each model must confirm: structures match, layout matches, colors match
5. If any model reports differences, go back to Phase 2 and fix

### Phase 4: Export and commit (worker)
1. Generate SVG: `pdftocairo -svg pnp_reaction.pdf pnp_reaction.svg`
2. Generate PNG: `pdftoppm -png -r 300 pnp_reaction.pdf pnp_reaction_preview`
3. Git add .tex, .pdf, .svg
4. Commit: "fix: PNP diagram matches original per vision model verification"
5. Push to origin/main

## Quality Gates

1. **Compiles clean**: `pdflatex -interaction=nonstopmode` with zero errors
2. **Structures correct**: Purine is fused bicyclic, ribose is closed furanose ring
3. **Colors correct**: Blue N, Red O, Yellow S, Orange P
4. **Labels minimal**: Only "MESG + Pi", "PNP", "Absorbance at 360nm"
5. **Vision consensus**: ALL 3 vision models agree render matches original
6. **No artifacts**: No `[::0]`, `[30]`, or other chemfig syntax leaking as text
7. **Pushed**: Git commit + push to felixfelix-bot/proteomics-pipeline main

## API Details for Workers

### Ollama Cloud Vision API
```
Endpoint: https://ollama.com/v1/chat/completions
Auth: Bearer $OLLAMA_CLOUD_API_KEY
Header: User-Agent: Ollama/0.30.10
Image format: {type: "image_url", image_url: {url: "data:image/png;base64,..."}}
Models: qwen3.5:397b, gemma4:31b, minimax-m3:cloud
```

### Files
- Original image: `/home/c03rad0r/.hermes/profiles/manager/image_cache/img_9aba242f8849.png`
- Working dir: `/home/c03rad0r/repos/proteomics-pipeline/molecular_diagrams/`
- Current .tex: `pnp_reaction.tex` (needs rewrite)
- Test file: `test_structures.tex` (new, for Phase 1)

## Reference: Original Diagram Text (from OCR)
- "PNP" (center, confirmed)
- "Absorbance at 360nm" (bottom, confirmed)
- Top has molecular structure with colored atoms