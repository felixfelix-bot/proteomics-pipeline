# Plan: Integrate STRING Network Analysis into GO Enrichment Pipeline

## Problem
Our current GO analysis (scripts 11, 13) runs enrichment on ALL significant
genes as one group. Lydia's original pipeline (scripts 04, 05) first mapped
proteins to the STRING interaction network, then split GO analysis by network
membership. This produces biologically more meaningful results.

## What Lydia Did (That We Don't)

1. **STRING network mapping**: Map significant proteins to STRING PPI database
2. **Seed identification**: Define high-confidence "seed" proteins using strict
   criteria (log2FC > 2 AND p < 1e-6, OR log2FC > 7)
3. **Network expansion**: Find all neighbors of seeds in STRING
4. **Categorization**: Split significant proteins into:
   - Seeds (highly enriched)
   - In-network (connected to seeds via STRING)
   - Not in network (significant but isolated)
5. **Split GO enrichment**: Run GO separately on each category
6. **Network visualization**: igraph plot showing protein interactions
7. **GSEA**: Ranked gene list analysis (in addition to ORA)

## Implementation Steps

### Step 1: Create 16_string_network_targeted.R
- Port Lydia's STRING analysis from 05_string_network.R
- Apply to TRIP4 vs WT experiment
- Output: network gene lists, candidate interactors, igraph visualization
- Uses Lydia's exact seed criteria

### Step 2: Update 11_targeted_go.R
- After STRING analysis, split GO enrichment:
  a) All significant genes (current behavior)
  b) Genes in STRING network of highly enriched
  c) Enriched genes NOT assigned to interaction network
- Add GSEA alongside ORA

### Step 3: Add make targets
- make string-network
- make open-string-network

### Step 4: Test and push
