# Isoform-resolved proteogenomics of CAD-associated splicing in VSMCs — analysis code

Analysis and figure-generation code for:

> Korchak JA\*, Perry RN\*, et al. *Isoform-resolved proteogenomics reveals the functional
> consequences of genetically regulated splicing in coronary artery disease.*

The scripts take the deposited and publicly available inputs listed in `data/README.md`
through colocalization integration, junction–isoform mapping, risk/nonrisk pair definition,
protein feature annotation, splice-event classification, peptide-evidence mapping, AP-MS scoring
and network annotation, and the generation of every R-made figure panel and Supplementary Table.

Upstream processing that lives in other repositories:

- Long-read transcriptome and predicted proteome (Iso-Seq3, SQANTI3, CPAT, SQANTI Protein):
  https://github.com/sheynkman-lab/Long-Read-Proteogenomics
- IS-PRM (Tomahto) target selection and post-acquisition filtering:
  https://github.com/sheynkman-lab/Alternative-Proteome-Detection-Project
- sQTL mapping and colocalization statistics: Aherrahrou et al. (2023), GEO GSE193817.

## Layout

```
config.R                                  all input/output paths, palette, shared plot helpers
01_colocalization_integration.Rmd
02_junction_isoform_mapping_and_pairs.Rmd
03_protein_feature_annotation.Rmd
04_splice_event_classification.Rmd
05_peptide_evidence_and_browser_panels.Rmd
06_pair_property_profiling.Rmd
07_apms_scoring_and_network_annotation.Rmd
08_apms_go_enrichment.Rmd
09_bait_prior_evidence_biogrid_intact.Rmd
10_figure_and_table_generation.Rmd
python/compute_af3_metrics.py             AlphaFold3 structural metrics (runs between 02 and 03)
data/                                     inputs (see data/README.md for every file and its source)
results/<script>/                         outputs of each script (created on first run)
figures/                                  figure panels, SVG + PNG (written by script 10)
```

Every script starts with `source("config.R")`; no script contains an absolute path. Files that
one script writes and another reads are named once in `config.R` (`OUT_*` variables), so the
chain can be followed from that file alone.

## Scripts

| script | what it does | manuscript outputs |
|---|---|---|
| 01 colocalization integration | Assigns each colocalized LeafCutter junction an effect direction relative to the CAD risk allele (Aragam et al. 2022 GWAS, with LD proxies), matches junctions by coordinates to the long-read and GENCODE v38 intron sets, and derives the colocalization / LD-expansion / SpliceAI counts | Fig. 1g source; numbers in Results (355 junctions, 168 sGenes, 304 lead variants, 1,534 shared variants, 115 with splicing evidence) |
| 02 junction–isoform mapping and pairs | Junction coding classification, LeafCutter cluster matching by coordinates, sibling PSI, risk/nonrisk junction assignment, candidate isoform pools, best-pair selection (highest-CPM non-NMD isoform), salvage of one-sided pairs, AP-MS pairs, inclusion flags (94 computed pairs), junction → amino-acid windows | Methods "Mapping CAD sQTL junctions…" and "Definition of risk- and nonrisk-associated junctions and isoform pairs"; inputs of 03–06 and 10 |
| *python/compute_af3_metrics.py* | Whole-isoform and pair-level AlphaFold3 metrics (pLDDT, radius of gyration, surface charge, DSSP secondary structure, TM-align, global alignment, junction-window metrics) | `data/protein_features/af3_metrics_per_{isoform,pair}.csv` read by 03 |
| 03 protein feature annotation | Per-pair property tables (length, AAindex, IUPred3, NetSurfP-3.0, PTM-POSE, DeepLoc 2.0, Pfam, UniProt feature overlap, AlphaFold3, ±10 aa junction windows), paired risk-vs-nonrisk tests, direction-matched and aggregate CAD-vs-non-CAD tests, NMD majority vote, named domain changes | Supplementary Table 15; Fig. 4f source; Extended Data Fig. 5b source; "43 of 94 pairs" domain statement |
| 04 splice-event classification | Three classifiers (junction pair, cluster topology, transcript pair) reconciled to a consensus event class per cluster; CAD vs non-CAD class distribution; predicted protein consequence per junction | Supplementary Table 14; Extended Data Fig. 5a; Fig. 4e source |
| 05 peptide evidence and browser panels | Sequence-first mapping of discovery and IS-PRM peptides to risk/nonrisk isoform pools, junction-level validation and exclusivity ceiling, Supplementary Table 12, and every genome-browser panel | Supplementary Tables 12 and 13 sources; Fig. 4g source; Extended Data Fig. 4a/5c sources; browser panels of Figs. 1e, 2c, 2f, 3e, 3f, 4i, 6a, 6f, 6k and Extended Data Figs. 1b, 2d/f/h, 3g/h, 4b/d, 6a/c/e/g |
| 06 pair property profiling | SD-scaled risk − nonrisk property differences, property families, dominant change mode per pair, oriented signature axes, sequence/structure similarity | Fig. 4h; Extended Data Figs. 5d and 7a sources; Supplementary Table 16 |
| 07 AP-MS scoring and network annotation | CompPASS ∩ SAINTexpress high-confidence interactions, risk/nonrisk/shared specificity, log2 fold change of spectral counts, VSMC expression support, BioGRID annotation, Aragam CAD-gene flag, control checks (RAB11B, GFP, beads, known-partner recovery) | Fig. 5b–e sources; Extended Data Fig. 7 (Cytoscape input); Supplementary Table 17; Results numbers (601 interactions, 83 % isoform-specific, 36 in BioGRID) |
| 08 AP-MS GO enrichment | clusterProfiler GO enrichment of risk-specific, nonrisk-specific and shared interactor sets per bait; top-term labels for the network figures | Fig. 5f; Fig. 5g/h labels; Supplementary Table 18 |
| 09 bait prior evidence | Gene-level BioGRID coverage and isoform-level IntAct records for the 13 bait genes / 26 bait isoforms | Methods and Results statements on prior interaction evidence |
| 10 figure and table generation | Reads only the exports above plus the experimental source data and draws every R-made panel at final size (SVG + 400 dpi PNG); assembles Supplementary Tables 8 and 10–20 | all R-made panels of Figs. 1–6 and Extended Data Figs. 1–8; Supplementary Tables |

Panels not made by these scripts: schematics (Illustrator/BioRender), LocusCompare plots,
AlphaFold3 models (ChimeraX), gel images, Cytoscape networks (from the 07/08 exports),
microscopy and spheroid images.

## Running

1. Assemble `data/` as described in `data/README.md`.
2. Knit the scripts in numerical order (`rmarkdown::render("01_colocalization_integration.Rmd")`,
   …). Run `python/compute_af3_metrics.py` after 02 and before 03 (or place the two deposited
   AF3 metric tables in `data/protein_features/`). Scripts 07–09 depend only on `data/` and on
   each other (08 and 09 read the 07 export); 10 depends on everything else.
3. If the working directory is not the repository root, set `CAD_REPO_ROOT` (used by
   `config.R`) to the repository path.

Each script checks that its inputs exist before doing any work (`require_files()`), asserts the
cardinalities the manuscript relies on (355 junctions, 146 long-read-captured junctions, 133
probeable junctions, 94 computed pairs) and ends with `sessionInfo()`.

### Software

R ≥ 4.3 with tidyverse, data.table, readxl, writexl, scales, patchwork, ggtext, ggrepel,
svglite, eulerr; Bioconductor GenomicRanges, rtracklayer, Biostrings, pwalign (≥ 3.19),
clusterProfiler, org.Hs.eg.db. Python ≥ 3.9 with numpy, pandas, gemmi, biopython, plus TM-align
and DSSP (`mkdssp`) for the AlphaFold3 metrics. Package versions used for the manuscript are
recorded in the `sessionInfo()` output at the end of each knitted script.

## Notes on reproducibility

- The IntAct check in script 09 is a live web-service query. The repository runs with
  `QUERY_INTACT <- FALSE` and reads the stored query result
  (`data/apms/intact_isoform_records.csv`); set it to `TRUE` to re-query and refresh that file.
- Fig. 2k P values are computed in script 10 (one-tailed paired t-test on the 12 measurements
  per condition); Fig. 2j adjusted P values are read from the deposited DESeq2 results of the
  PARP12 overexpression RNA-seq.
- Spheroid permutation tests in script 10 are seeded (`set.seed(1)`); the Wilcoxon statistics
  reported in the manuscript are deterministic.
- Supplementary Table 18 contains every GO term passing the enrichment cutoffs; the ten most
  significant terms per set are used only to pick the labels shown in Fig. 5g/h.
- Gene symbols are harmonised to HGNC for display (`RP11-77K12.1` → `CFDP1`, `NHERF2` →
  `SLC9A3R2`); input tables keep their original identifiers.

## Acknowledgement

The per-bait GO enrichment and network-label code in script 08 was adapted from
https://github.com/SIDR73/CADNetworks.

## Citation

Please cite the manuscript above. Processed data are available as described in its Data
availability statement (GEO GSE193817; NCBI SRA; ProteomeXchange/PRIDE; GWAS Catalog
GCST90132314).
