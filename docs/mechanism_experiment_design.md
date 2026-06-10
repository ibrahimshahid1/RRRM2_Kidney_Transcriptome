# Resolving the upstream driver of spaceflight WNK–SPAK–NCC dephosphorylation

**Purpose.** A pre-registered, data-grounded experimental design to move from the v11 *constraint* result ("the transporter lesion sits at regulatory phosphorylation, but public data cannot identify the upstream cause") to a *mechanism*. Two experiments, sequenced cheapest-to-definitive. Intended as a grant-aim / protocol skeleton, not a claim of results.

---

## 1. The question, sharpened

Observed endpoint (OSD-462/RR-10, whole kidney): the NCC activating phosphocluster and SPAK's activation-loop site are suppressed in flight while protein abundance is flat. Two forks must be resolved, in order:

- **Fork A — regulation vs remodeling.** Is phosphorylation lower *per DCT cell* (active regulation), or does whole-kidney signal fall because DCT cells are lost/dedifferentiated (composition)? Whole-kidney data cannot separate these.
- **Fork B — which driver (if regulation).** (i) ionic/osmotic suppression of WNK kinase activity (WNKs are Cl⁻-inhibited); (ii) phosphatase activation (PP1/PP4 on NCC/SPAK); (iii) enhanced WNK turnover (KLHL3/CUL3); (iv) hormonal/volume signaling (aldosterone/SGK1).

## 2. What OSD-462 already fixes — and what it is blind to

This is the design's anchor. Direct readout of the cascade in the existing data:

| Cascade node | Site (OSD-462) | Flight effect (log₂) | q | Interpretation |
|---|---|---|---|---|
| NCC activating cluster | Slc12a3 pS65 | −0.79 | 8.8×10⁻⁵ | output suppressed |
| NCC activating cluster | Slc12a3 p65;68 | −1.56 | 4.7×10⁻⁵ | output suppressed |
| NCC activating cluster | Slc12a3 p89 | −0.93 | 8.5×10⁻⁴ | output suppressed |
| SPAK activation loop | Stk39 pS383 | −0.79 | 2.0×10⁻³ | mid-cascade input suppressed |
| SPAK S-motif | Stk39 pS366 | −0.52 | 3.0×10⁻⁴ | mid-cascade suppressed |

**The blind spots (these define the experiment's required new readouts):**

- **WNK1/WNK4 kinase-activation autophosphorylation** — the Cl⁻-sensor output that says whether WNK *activity* is actually down — is **not quantified** (the measured WNK1/WNK4 sites are C-terminal/regulatory, not the activation loop).
- **OSR1 activation loop** — 1 site, flat; effectively uncovered.
- **Phosphatases / turnover** — Ppp4c 0 sites, Klhl3 0 sites, Cul3 1 site; the phosphatase-activation and WNK-degradation hypotheses are untestable in this data.

**So the crux is precise:** the data shows the cascade is down at SPAK and NCC, but cannot say whether that is because **WNK activity fell** (ionic/turnover, *upstream*) or because **NCC/SPAK are being actively dephosphorylated** at normal WNK activity (phosphatase, *downstream*). Every experiment below is built to read exactly those missing nodes.

---

## 3. Experiment A — mpkDCT ionic-perturbation cascade readout (cell-intrinsic, cheap, first)

**Rationale.** Isolates Fork B(i) from remodeling entirely: a pure DCT-lineage cell line (mpkDCT; the same line behind PXD001729) lets you impose the candidate spaceflight cell state and read the whole cascade with no composition confound.

**Conditions (factorial):** (1) simulated microgravity / fluid-shift (clinostat or RWV bioreactor) vs static; crossed with (2) ionic perturbation along the Cl⁻/K⁺ axis (low-Cl⁻ → WNK-active control; high-Cl⁻/raised [K⁺] → predicted WNK-inhibiting); plus (3) a phosphatase-inhibitor arm (calyculin A / tautomycetin) to test whether dephosphorylation is phosphatase-driven.

**Readouts (targeted phospho-MS or phos-specific WB), in cascade order — including the OSD-462 blind spots:**
1. **WNK1/WNK4 activation-loop autophosphorylation** (the missing chloride-sensor readout) — *primary discriminator.*
2. **SPAK pS383 / OSR1 pS325** (T-loop) — activation of the middle.
3. **NCC pT53/pS65/pS68** — the endpoint; must recapitulate the flight pattern as the assay's positive anchor.
4. PP1/PP4 activity (phosphatase assay) in the inhibitor arm.

**Predicted patterns (the discriminating logic):**
- *Ionic/osmotic WNK suppression (hypothesis i):* high-Cl⁻/K⁺ → **WNK autophospho ↓ → SPAK pS383 ↓ → NCC ↓**, top-down. Phosphatase inhibitor does *not* rescue.
- *Phosphatase activation (hypothesis ii):* **WNK autophospho normal, SPAK/NCC ↓**, and phosphatase inhibitor **rescues** NCC/SPAK phosphorylation.
- *Neither reproduces flight in vitro:* points to a tissue-level/composition driver → go to Experiment B as primary.

**Controls/anchor:** dDAVP arm (ties to PXD001729); confirm the low-Cl⁻ condition gives the canonical WNK-active NCC pattern; require the flight-mimicking arm to reproduce the OSD-462 NCC pS65/p65;68 suppression before any mechanism is claimed.

---

## 4. Experiment B — per-cell localization (settles Fork A; needs tissue)

**Rationale.** Resolves regulation vs remodeling, which Experiment A cannot (cell line has no remodeling). Requires flight or hindlimb-unloaded analog kidney.

**Approach:** sorted/enriched DCT cells (e.g., NCC⁺ or parvalbumin-DCT) → targeted phospho-MS, **and/or** spatial phosphoproteomics / spatial transcriptomics + multiplex phospho-IF on sections, paired with snRNA-seq to quantify DCT fraction and state.

**Pivotal readout:** **NCC phospho *per DCT cell*.** If still ↓ after normalizing to DCT cell fraction → active regulation (proceed on Experiment A's mechanism). If whole-tissue ↓ but per-cell normal → composition/remodeling, and the "transporter suppression" reframes as cell loss. This is also the spatial prediction already stated in the v11 manuscript (DCT-adjacent vs DCT-marker-high spots).

---

## 5. Decision matrix

| Per-cell NCC phospho (Exp B) | WNK autophospho (Exp A) | Phosphatase-inhibitor rescue (Exp A) | Implicated driver |
|---|---|---|---|
| ↓ (regulation) | ↓ | no | **Ionic/osmotic WNK suppression** |
| ↓ (regulation) | normal | yes | **Phosphatase activation** |
| ↓ (regulation) | normal | no | Trafficking / WNK turnover → ubiquitinomics (TUBE) follow-up |
| normal per cell | — | — | **Remodeling / DCT cell loss**, not active suppression |

Add: matched plasma/urine electrolytes, [K⁺], aldosterone, and volume status in the Exp-B animals — the physiology the public data entirely lacks — to constrain ionic vs hormonal inputs.

---

## 6. Sequence & feasibility

1. **Exp A (mpkDCT)** — weeks–months, standard cell + phospho-MS budget; nominates the mechanism and is publishable on its own as the cell-intrinsic test.
2. **Hindlimb-unloading analog** — adds in-vivo regulation + rescue genetics/pharmacology at moderate cost without a flight slot.
3. **Exp B targeted add-on to the next rodent flight (RR-X)** — the definitive in-vivo confirmation; expensive, slot-limited, so it should ride an existing mission with the readout list above pre-specified.

**Bottom line:** Experiment A is the highest-value first move — it directly reads the one node (WNK activation-loop autophosphorylation) that the public data is blind to and that separates "WNK is off" from "something is dephosphorylating NCC." That single readout is the hinge of the whole mechanism question.
