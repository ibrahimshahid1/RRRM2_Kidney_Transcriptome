# Biology-First Framing

## Core Thesis

RRRM-2 should be framed as a rigor-remediated cross-cohort analysis of kidney remodeling after spaceflight. The biological question is not "which individual genes are significant everywhere?" but "which kidney programs recur, reverse, or disappear across mission, strain, sex, and recovery contexts?"

The fixed external panel is:

- PPAR/fatty-acid metabolism.
- Cholesterol biosynthesis.
- ECM remodeling.
- EMT/fibrosis.
- Tubular ion transport.
- TGF-beta/Wnt signaling.
- Oxidative stress.
- Translation machinery.

## Cohort Roles

OSD-771 is the discovery and network-rewiring cohort. It supplies the factorial Age x Arm x Environment structure and the ISS-T/LAR contrast logic.

OSD-102 is the closest directional replication cohort for LAR-Young-like biology: female, young, C57BL/6J kidney, FLT vs GC. It cannot validate ISS-T, old-arm, classifier, or global network claims.

OSD-513 is the secondary sex-robustness cohort: male, young, C57BL/6J kidney, FLT vs HGC. It asks whether OSD-771/OSD-102 pathway signals persist in males.

OSD-163 and OSD-253 are context-mapping cohorts. They are used to ask whether lipid, ECM, tubular, stress, and translation axes recur across broader strain/mission settings. They are not strict replication cohorts unless a nonzero directional discovery effect is explicitly registered before analysis.

## Claim Ladder

Use "replicated" only for pre-registered directional hypotheses with nonzero discovery effects, same external direction, and q passing the registered threshold.

Use "context_detected" for OSD-163/253 rows with discovery_effect = 0 and q passing the registered threshold. This means the pathway is active in that external cohort; it does not mean the RRRM-2 direction replicated. The preferred external pathway statistic is now preranked GSEA; the older mean pathway t-statistic is retained only as a sensitivity check.

Use "candidate" for edge-sum genes, top-decile rewiring genes, pathway enrichments from exploratory ORA, and silent-shifter-like genes that do not pass their calibration benchmarks.

Use "unsupported" for broad regression-level gene discoveries, strict silent-shifter enrichment, classifier superiority over expression baselines, ISS-T external validation, old-arm external validation, and global cohort expansion.

## Manuscript Shape

A defensible biology-first manuscript would center one figure on the OSD-771 remodeling map, one figure on the cross-cohort pathway panel across OSD-102/513/163/253, and one supplement on the negative statistical stress tests. The abstract should emphasize context-dependent lipid/ECM/tubular remodeling rather than per-gene discovery counts.
