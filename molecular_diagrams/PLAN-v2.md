# PNP Assay Diagram — Execution Plan v2
# Created: 2026-08-11
# Goal: TikZ diagram matching original at 7+/10 from 2+ vision models

## PROBLEM ANALYSIS

Previous failures:
- 3× worker timeout at 300s (compile+render+vision loop too long)
- gemma4:31b rates 3/10 but has poor spatial reasoning — misidentifies molecules
- z.ai vision 503 (exhausted)
- kimi-k3:cloud signin blocked on Ollama

Root causes:
1. Trying to do everything in one 300s delegate_task
2. Using gemma4 for both construction guidance AND verification — it's fast but dumb
3. No structured comparison — just "rate 1-10" with no reference description

## STRATEGY

Split into 4 independent stages, each fits in 300s:
- Stage 1: CONSTRUCT (no vision needed — pure code)
- Stage 2: VERIFY (vision comparison, one model at a time)
- Stage 3: FIX (code changes based on verification feedback)
- Stage 4: COMMIT (git push)

Use kimi-k2.7-code for construction (great spatial reasoning, no vision needed — 
give it the gemma4 description of the original + current .tex, ask it to fix).
Use qwen3.5:397b for verification (slow but smart, good at detailed comparison).
Use gemma4:31b for quick sanity checks between full qwen3.5 runs.

## STAGES

### Stage 1: CONSTRUCT (delegate_task, leaf, kimi-k2.7-code, 300s)
- Input: gemma4's detailed description of original + current pnp_reaction.tex
- Task: Rewrite pnp_reaction.tex to match the original description
- kimi-k2.7-code has excellent code/spatial reasoning, no vision needed
- It reads the description of what the original looks like and fixes the TikZ
- Success: .tex compiles, PDF renders
- Cron: NO — one-shot delegate_task

### Stage 2: VERIFY (direct execution, qwen3.5:397b, 120s)
- Compile + render the .tex
- Send both images (original + render) to qwen3.5:397b via Ollama Cloud
- Ask: "Rate 1-10. List top 5 specific differences with coordinates."
- 90s per call but worth it for quality
- Success: rating >= 7/10
- Cron: NO — direct terminal call

### Stage 3: FIX LOOP (repeat stages 1-2, max 3 iterations)
- If qwen3.5 rates < 7/10:
  - Take its specific feedback
  - Feed back to kimi-k2.7-code delegate_task with the feedback
  - Re-verify with qwen3.5
- If >= 7/10: proceed to Stage 4
- Cron: NO — manual loop

### Stage 4: CROSS-VERIFY + COMMIT (direct, 120s)
- Verify with gemma4:31b (quick sanity check)
- If both qwen3.5 AND gemma4 rate 7+/10:
  - pdftocairo -svg pnp_reaction.pdf pnp_reaction.svg
  - git add pnp_reaction.tex pnp_reaction.pdf pnp_reaction.svg
  - git commit -m "feat: PNP assay diagram verified by vision models (7+/10)"
  - git push origin main
- Cron: NO — direct execution

## SCHEDULING

This is NOT a cron job — it's a sequential pipeline of delegate_tasks 
and direct terminal calls. The manager (me) orchestrates:

1. Dispatch kimi-k2.7-code worker to fix .tex (300s)
2. Compile + render (30s direct)
3. Verify with qwen3.5:397b (90s direct)  
4. If < 7/10: goto 1 with feedback
5. If >= 7/10: verify with gemma4 (30s direct)
6. If both 7+: commit + push (10s direct)

Total estimated time: 3-5 iterations × (300s + 30s + 90s) = ~20-35 min

## KEY INPUTS FOR STAGE 1

Gemma4's description of the ORIGINAL (confirmed accurate):
- Top: MESG (purine base + ribose) with S⁻ at C6, NH2 at C2, CH3 at N7
- Pi text to the RIGHT of MESG
- "MESG + Pi" text label
- Middle: vertical arrow with "PNP" label
- Bottom LEFT: 7-methyl-6-thioguanine (C=S thione, N9-H, N1-H, NH2, CH3)
- Bottom CENTER: "+" sign
- Bottom RIGHT: ribose-1-phosphate (ribose + phosphate via C1'-O-P)
- Very bottom: "Absorbance at 360nm" (italic)

Colors: N=#0066CC, O=#CC0000, S=#D4A017, P=#E07000