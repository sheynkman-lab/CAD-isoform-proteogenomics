# =============================================================================
# config.R — shared configuration for the CAD isoform manuscript analysis code
# =============================================================================
# Every script in this repository starts with `source("config.R")`. This file is
# the single place where input locations, output locations, the colour palette
# and the plotting helpers are defined, so that no script carries an absolute
# path and the whole chain can be re-run from a fresh checkout.
#
# Layout (relative to the repository root):
#   data/       inputs — deposited data, public resources and experimental source
#               data; see data/README.md for every file and where it comes from
#   results/    outputs — one sub-directory per script (created on demand)
#   figures/    figure panels (SVG + PNG) written by 10_figure_and_table_generation
#
# Set REPO_ROOT if the scripts are run from a working directory other than the
# repository root (for example when knitting from RStudio with a different
# project directory).
# -----------------------------------------------------------------------------

REPO_ROOT <- Sys.getenv("CAD_REPO_ROOT", unset = ".")
DATA_DIR  <- file.path(REPO_ROOT, "data")
RES_DIR   <- file.path(REPO_ROOT, "results")
FIG_DIR   <- file.path(REPO_ROOT, "figures")

# One results directory per script ---------------------------------------------
RES_01 <- file.path(RES_DIR, "01_colocalization_integration")
RES_02 <- file.path(RES_DIR, "02_junction_isoform_mapping_and_pairs")
RES_03 <- file.path(RES_DIR, "03_protein_feature_annotation")
RES_04 <- file.path(RES_DIR, "04_splice_event_classification")
RES_05 <- file.path(RES_DIR, "05_peptide_evidence_and_browser_panels")
RES_06 <- file.path(RES_DIR, "06_pair_property_profiling")
RES_07 <- file.path(RES_DIR, "07_apms_scoring_and_network_annotation")
RES_08 <- file.path(RES_DIR, "08_apms_go_enrichment")
RES_09 <- file.path(RES_DIR, "09_bait_prior_evidence_biogrid_intact")
RES_10 <- file.path(RES_DIR, "10_figure_and_table_generation")
TABLE_DIR <- file.path(RES_10, "tables")               # per-panel source-data tables
SUPP_DIR  <- file.path(RES_10, "supplementary_tables") # Supplementary Tables
BROWSER_DIR <- file.path(RES_05, "browser_panels")      # genome-browser panels (SVG + PNG)

ensure_dir <- function(...) {
  for (d in c(...)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  invisible(NULL)
}

# =============================================================================
# INPUT FILES (F_*) — all under data/; see data/README.md
# =============================================================================

# ---- sQTL study (Aherrahrou et al. 2023; GEO GSE193817) ----------------------
D_SQTL <- file.path(DATA_DIR, "sqtl")
F_LEAFCUTTER_BED_PROLIF   <- file.path(D_SQTL, "WithFBS.bed")              # LeafCutter clusters, proliferative
F_LEAFCUTTER_BED_QUIESC   <- file.path(D_SQTL, "WithoutFBS.bed")           # LeafCutter clusters, quiescent
F_LEAFCUTTER_COUNTS_PROLIF <- file.path(D_SQTL, "WithFBS_cleaned.counts")  # numerator/denominator counts
F_LEAFCUTTER_COUNTS_QUIESC <- file.path(D_SQTL, "WithoutFBS_cleaned.counts")
F_GENERAL_SQTL_LR          <- file.path(D_SQTL, "general_sqtls_mapped_to_long_read.tsv") # all VSMC sQTL junctions (non-CAD background)

# ---- CAD GWAS colocalization and variant annotation -------------------------
D_COLOC <- file.path(DATA_DIR, "colocalization")
F_COLOC_WITH_SLOPES  <- file.path(D_COLOC, "2026_coloc_with_slopes.csv")        # colocalized sQTLs with slopes, ref/alt, lead rsID
F_COLOC_NOVEL_ONLY   <- file.path(D_COLOC, "2026_novel_only_updated_coloc.csv") # colocalizations by novelty status (Fig. 1b)
F_CAD_SQTL_SNPS_LD   <- file.path(D_COLOC, "2026_CAD_sQTL_SNPs.csv")            # LD-expanded (r2 >= 0.8) sQTL/CAD SNPs (HaploReg v4.2)
F_SQTL_PSI_LONG      <- file.path(D_COLOC, "sqtl_plots_long.csv")               # per-donor PSI by genotype for plotted sQTLs
F_SQTL_SNP_JUNCTION  <- file.path(D_COLOC, "sqtl_snp_junction_lookup.csv")      # cluster/junction per plotted SNP
F_RAB5C_MST1_PSI     <- file.path(D_COLOC, "RAB5C_MST1_sQTL_plot_data.xlsx")    # PSI for LD-proxy SNPs (Extended Data Fig. 1d, 1f)
F_SUPP_TABLES_XLSX   <- file.path(D_COLOC, "Supplementary_Tables.xlsx")         # the manuscript workbook; sheets S1–S5 are inputs here
F_RT_PCR_PRIMERS     <- file.path(D_COLOC, "RT-PCR_primers.xlsx")               # amplicon coordinates (Supplementary Table 6)

# ---- CAD GWAS ----------------------------------------------------------------
D_GWAS <- file.path(DATA_DIR, "gwas")
F_ARAGAM_SNP_SUBSET  <- file.path(D_GWAS, "aragam_cad_snps.tsv")             # GCST90132314 rows for lead + proxy rsIDs
F_ARAGAM_S31_GENES   <- file.path(D_GWAS, "Aragam_2022_S31_CAD_GWAS.csv")    # most likely causal CAD genes (Aragam 2022, Supp. Table 31)

# ---- VSMC long-read transcriptome and predicted proteome ---------------------
D_LR <- file.path(DATA_DIR, "longread")
F_SMC_GTF             <- file.path(D_LR, "SMC_CORRECTED_with_cds.gtf")               # long-read transcript models with CDS
F_SMC_INTRON_COORDS   <- file.path(D_LR, "SMC_CORRECTED_with_cds_intron_coords.txt")  # intron coordinates per transcript
F_SMC_EXON_COORDS     <- file.path(D_LR, "SMC_CORRECTED_with_cds_exon_coords.txt")    # exon coordinates per transcript
F_SMC_PSL             <- file.path(D_LR, "SMC_corrected.psl")                          # exon block structure per transcript
F_SMC_PROTEIN_FASTA   <- file.path(D_LR, "SMC_CORRECTED_orf_refined.fasta")           # predicted protein sequences
F_SQANTI_PROTEIN_CPM  <- file.path(D_LR, "260115_sqanti_protein_w_cpm.csv")           # SQANTI Protein classes + CPM + NMD + ORF coords
F_SQANTI_PROTEIN_TSV  <- file.path(D_LR, "SMC.sqanti_protein_classification.tsv")     # SQANTI Protein classification (raw)
F_SQANTI3_CLASS       <- file.path(D_LR, "SMC_classification.txt")                    # SQANTI3 transcript classification
F_MAPPED_FL_COUNTS    <- file.path(D_LR, "SMC.mapped_fl_count.txt")                   # full-length read counts per transcript per sample
F_DEMUX_SAMPLE_COUNTS <- file.path(D_LR, "SMC.demuxed_with_sample_data.tsv")          # per-sample demultiplexed counts with donor
F_BEST_ORF_CPM        <- file.path(D_LR, "SMC_CORRECTED_best_orf.tsv")                # best ORF per transcript with CPM and gene

# ---- Public reference annotation --------------------------------------------
D_REF <- file.path(DATA_DIR, "reference")
F_GENCODE_V38_GTF        <- file.path(D_REF, "gencode.v38.annotation.gtf")
F_GENCODE_V38_INTRONS    <- file.path(D_REF, "gencode_v38_intron_coords_to_isoforms.txt") # intron -> ENST list, derived from the GTF
F_GENCODE_V38_EXONS      <- file.path(D_REF, "gencode_v38_exon_coords.txt")               # exon table, derived from the GTF
F_UNIPROT_ANNOT          <- file.path(D_REF, "uniprot_combined_annot.tsv")                # UniProt feature intervals (hg38)
D_UNIPROT_BED            <- file.path(D_REF, "uniprot_bed")                               # UCSC UniProt track BEDs (unipDomain.bed, ...)
F_BIOGRID_TAB3           <- file.path(D_REF, "BIOGRID-ORGANISM-Homo_sapiens-5.0.252.tab3.txt")
F_HPA_SMC_GENES          <- file.path(D_REF, "all_SMC_genes_expressed_from_HPA.tsv")      # Human Protein Atlas smooth-muscle-cell expressed genes

# ---- Per-isoform protein feature predictions ---------------------------------
D_FEAT <- file.path(DATA_DIR, "protein_features")
F_IUPRED_SUMMARY     <- file.path(D_FEAT, "iupred_smc_lr_summary.csv")
F_IUPRED_RESIDUES    <- file.path(D_FEAT, "iupred_smc_lr_all_residues.tsv")
F_HMMSCAN_DOMTBLOUT  <- file.path(D_FEAT, "hmmscan_pfam_SMC_lrRNA_hmmscan_domtblout.txt")  # Pfam-A v35.0
F_DEEPLOC_SMC        <- file.path(D_FEAT, "smc_deeploc_master.csv")                       # DeepLoc 2.0
F_NETSURFP_FEATURES  <- file.path(D_FEAT, "netsurfp3_smc_protein_features_slim.csv")       # NetSurfP-3.0 per isoform
F_NETSURFP_RESIDUES  <- file.path(D_FEAT, "netsurfp3_smc_master_residue_predictions.csv")  # NetSurfP-3.0 per residue
F_AAINDEX            <- file.path(D_FEAT, "AAIndex_flexibility_granthamvolume_kytedoolittlehydropathy_zimmermanpolarity.csv")
F_PTM_EXONS          <- file.path(D_FEAT, "smc_lr_exons_with_ptms.tsv")                   # PTM-POSE
F_PTM_SUMMARY        <- file.path(D_FEAT, "smc_lr_ptms_summary.tsv")                      # PTM-POSE
F_AF3_METRICS_PAIR   <- file.path(D_FEAT, "af3_metrics_per_pair.csv")                     # AlphaFold3 pair alignment metrics
F_AF3_METRICS_ISO    <- file.path(D_FEAT, "af3_metrics_per_isoform.csv")                  # AlphaFold3 per-isoform metrics
# Both AF3 tables are produced by python/compute_af3_metrics.py from the AlphaFold3 models.

# ---- VSMC proteomics (discovery and IS-PRM) ----------------------------------
D_PROT <- file.path(DATA_DIR, "proteomics")
F_CLEAVER_DIGEST    <- file.path(D_PROT, "231023_SMC_hybrid_gencode_v43_PEPTIDES.csv")          # in-silico tryptic digest (Cleaver)
F_CAD173_TARGETS    <- file.path(D_PROT, "CAD173_info.csv")                                      # the 173 IS-PRM target peptides
F_ISPRM_DETECTED    <- file.path(D_PROT, "240710_R_CAD173_rep_TOM_OUTPUT_include_NA_filtered_out.csv") # IS-PRM detections (Tomahto, filtered)
F_POGO_ISPRM_BED    <- file.path(D_PROT, "240710_R_CAD173_rep_TOM_POGO.bed")                     # PoGo coordinates, IS-PRM peptides
F_POGO_RTS_BED      <- file.path(D_PROT, "231030_rts_EDJ_hits_combined_POGO.bed")                 # PoGo coordinates, RTS (B94)
F_POGO_DDA_CAD_BED  <- file.path(D_PROT, "231031_dda_cad_peptides_POGO.bed")                      # PoGo coordinates, DDA CAD-gene peptides
F_POGO_DDA_ALL_BED  <- file.path(D_PROT, "231031_L101_B94_Z99_all_combined_Pogo.bed")             # PoGo coordinates, all DDA peptides

# ---- AP-MS --------------------------------------------------------------------
D_APMS <- file.path(DATA_DIR, "apms")
F_SAINT_SET          <- file.path(D_APMS, sprintf("saint_express_output_10_CAD_set_%d_isoforms.csv", 1:3))
F_COMPPASS_SET       <- file.path(D_APMS, sprintf("comppass_output_10_CAD_set_%d.csv", 1:3))
F_APMS_BAIT_PAIRS    <- file.path(D_APMS, "CAD_APMS_Baits_risk_v_nonrisk.csv")   # the 13 bait pairs (gene, risk PB, nonrisk PB)
F_APMS_BAIT_FASTA    <- file.path(D_APMS, "Bait_Protein_Sequences_CAD.fa")       # bait ORF protein sequences (Supplementary Table 8)
F_SMC_SPECIFICITY_MANUAL <- file.path(D_APMS, "SMC_specificity_table.csv")       # manual VSMC-expression check for preys absent from HPA
F_SRP54_SPECTRAL     <- file.path(D_APMS, "srp54_investigation.csv")             # SRP54 pulldown spectral counts for SRP subunits
F_INTACT_RECORDS     <- file.path(D_APMS, "intact_isoform_records.csv")          # stored IntAct PSICQUIC query result for the bait accessions

# ---- Experimental source data ------------------------------------------------
D_EXP <- file.path(DATA_DIR, "experimental")
F_GEL_DENSITOMETRY   <- file.path(D_EXP, "Isoform_gel_quant_plot_data.xlsx")     # RT-PCR gel densitometry (PARP12, WIPI1, DCLRE1B, ITGA1_raw)
F_PARP12_OE_DATA     <- file.path(D_EXP, "PARP12_OE_plots_raw_data.xlsx")        # qPCR QC, Hallmark GSEA, NF-kB targets, p65 ELISA
F_ITGA1_SPHEROID_IMAGING <- file.path(D_EXP, "ITGA1_spheroid_raw_data.xlsx")     # actin/membrane, SHG per spheroid
F_ITGA1_SPHEROID_MECH    <- file.path(D_EXP, "ITGA1_spheroid_compression_raw.xlsx") # per-spheroid moduli, circularity, batch
F_PARP12_DESEQ2          <- file.path(D_EXP, "PARP12_OE_DESeq2_results.csv")       # DESeq2 results, risk vs nonrisk PARP12 isoform OE (gene, log2FoldChange, pvalue, padj)

# =============================================================================
# FILES EXCHANGED BETWEEN SCRIPTS (OUT_*) — written by one script, read by others
# =============================================================================
# 01 — colocalization integration
OUT_LFC_INTEGRATED        <- file.path(RES_01, "lfc_sQTLs_integrated_SMC_corrected_intron_coords.txt")
OUT_VARIANT_VEP_SPLICEAI  <- file.path(RES_01, "variant_level_vep_spliceai.csv")
OUT_COLOC_NUMBER_SUMMARY  <- file.path(RES_01, "colocalization_number_summary.csv")

# 02 — junction–isoform mapping and risk/nonrisk pairs
OUT_CAD_JUNCTION_CODING     <- file.path(RES_02, "cad_junction_coding_classification.csv")
OUT_NONCAD_JUNCTION_CODING  <- file.path(RES_02, "noncad_junction_coding_classification.csv")
OUT_JUNCTION_UNIVERSE       <- file.path(RES_02, "full_junction_universe_for_AS.csv")
OUT_JUNCTION_ISOFORM_MAP    <- file.path(RES_02, "full_junction_isoform_map_for_AS.csv")
OUT_CAD_JUNCTION_PAIRS      <- file.path(RES_02, "cad_risk_nonrisk_junction_pairs.csv")
OUT_CAD_CLUSTER_SIBLINGS    <- file.path(RES_02, "cad_cluster_siblings_with_psi.csv")
OUT_CAD_CLUSTER_MAP         <- file.path(RES_02, "cad_junction_cluster_map.csv")
OUT_PAIRS_ALL               <- file.path(RES_02, "leafcutter_isoform_pairs_all.csv")
OUT_PAIRS_FULL_POOLS        <- file.path(RES_02, "leafcutter_isoform_pairs_full_pools.csv")
OUT_PAIRS_UNIFIED           <- file.path(RES_02, "leafcutter_isoform_pairs_unified.csv")
OUT_PAIRS_UNIFIED_FLAGGED   <- file.path(RES_02, "leafcutter_isoform_pairs_unified_flagged.csv")
OUT_JUNCTION_AA_MAPPING     <- file.path(RES_02, "junction_aa_mapping.csv")

# 03 — protein feature annotation and property statistics
OUT_PROPERTY_ISOFORM_WIDE   <- file.path(RES_03, "property_table_isoform_wide.csv")
OUT_PROPERTY_JUNCTION_WIDE  <- file.path(RES_03, "property_table_junction_wide.csv")
OUT_JUNCTION_AA_MAPPING_NONCAD <- file.path(RES_03, "junction_aa_mapping_NONCAD.csv")
OUT_AGGREGATE_JUNCTION      <- file.path(RES_03, "aggregate_cad_vs_noncad_junction_results.csv")
OUT_NMD_CAD_MAJORITY        <- file.path(RES_03, "nmd_cad_junction_majority_vote.csv")
OUT_NMD_NONCAD_MAJORITY     <- file.path(RES_03, "nmd_noncad_junction_majority_vote.csv")
OUT_PROPERTY_STATS_MASTER   <- file.path(RES_03, "master_results_summary.csv")   # Supplementary Table 15
OUT_DOMAIN_CHANGES_NAMED    <- file.path(RES_03, "domain_changes_named.csv")

# 04 — splice-event classification
OUT_AS_CONSENSUS            <- file.path(RES_04, "as_event_classification_consensus.csv")
OUT_AS_DISTRIBUTION         <- file.path(RES_04, "as_type_distribution_cad_vs_noncad.csv")
OUT_PROTEIN_EFFECTS         <- file.path(RES_04, "protein_effects_per_junction.csv")

# 05 — peptide evidence
OUT_PEPTIDE_EVIDENCE_DETAIL   <- file.path(RES_05, "peptide_isoform_evidence_detail.csv")
OUT_JUNCTION_VALIDATION       <- file.path(RES_05, "junction_validation_seqfirst.csv")
OUT_JUNCTION_EXCLUSIVE_PEPS   <- file.path(RES_05, "exclusive_peptide_possibility_per_junction.csv")
OUT_JUNCTION_NOT_USABLE       <- file.path(RES_05, "junction_not_usable_reasons.csv")
OUT_EVIDENCE_CHAIN            <- file.path(RES_05, "evidence_chain.csv")
OUT_PEPTIDE_SUPP_TABLE        <- file.path(RES_05, "SuppTable_CAD_sGene_peptide_evidence.csv")
OUT_INTERESTING_ISOFORMS      <- file.path(RES_05, "validated_interesting_isoform_candidates_bonafide.csv")

# 06 — pair property profiling
OUT_PAIR_DOMINANT_MODE      <- file.path(RES_06, "per_pair_dominant_mode.csv")
OUT_PAIR_GROUPINGS          <- file.path(RES_06, "pair_groupings_table.csv")
OUT_PAIR_SIGNATURES_LONG    <- file.path(RES_06, "pair_property_signatures_long.csv")
OUT_PROPERTY_AXIS_ORDER     <- file.path(RES_06, "property_axis_order.csv")
OUT_PAIR_SIMILARITY         <- file.path(RES_06, "pair_structure_vs_sequence_similarity.csv")
OUT_APMS_PAIR_SIGNATURES    <- file.path(RES_06, "apms_pair_property_signatures_long.csv")
OUT_APMS_BAIT_PAIRS         <- file.path(RES_06, "apms_bait_pairs.csv")

# 07 — AP-MS scoring and network annotation
OUT_SAINT_COMBINED          <- file.path(RES_07, "saint_combined.tsv")
OUT_MASTER_NETWORK_LFC      <- file.path(RES_07, "MASTER_network_LFC.csv")
OUT_MASTER_NETWORK_ANNOT    <- file.path(RES_07, "master_network_full_annotated.csv")
OUT_PREY_VSMC_SUPPORT       <- file.path(RES_07, "combined_shared_full_annot_summary.csv")
OUT_APMS_LONG_SUMMARY       <- file.path(RES_07, "CAD_APMS_long_summary.csv")
OUT_APMS_WIDE_SUMMARY       <- file.path(RES_07, "CAD_APMS_wide_summary.csv")

# 08 — AP-MS GO enrichment
D_GO_PAIR_TABLES            <- file.path(RES_08, "pair_by_pair")                 # <BAIT>_pathway_details.csv
OUT_GO_MASTER_SUMMARY       <- file.path(RES_08, "MASTER_pathway_summary.csv")
OUT_GO_NETWORK_LABELS       <- file.path(RES_08, "figure5_GO_labels.csv")

# 09 — bait prior evidence
OUT_BAIT_NOVELTY            <- file.path(RES_09, "bait_isoform_novelty_table.csv")

# =============================================================================
# SHARED PALETTE AND PLOT HELPERS
# =============================================================================
# Colours are fixed across the manuscript: risk red, nonrisk blue, CAD purple,
# non-CAD green, novel orange, annotated grey.
COL <- list(
  risk      = "#C42B2B", nonrisk   = "#2166AC",
  cad       = "#7B2D8E", noncad    = "#1B9E77",
  novel     = "#E8700A", annotated = "#9E9E9E",
  neutral   = "#9E9E9E", light     = "#BDBDBD",
  carried   = "#6E8B7B"
)
# Lighter risk/nonrisk fills used by the figure script for bars and heatmaps.
COL_FIG <- list(risk = "#D66B6B", nonrisk = "#6494C5")

# Publication constraints: text >= 7 pt; export width <= 7.5 in.
MAX_W     <- 7.5
BASE_SIZE <- 7
sz <- function(pt) pt / 2.845276          # ggplot text sizes are in mm; sz(7) draws 7 pt

theme_cell <- function(base_size = BASE_SIZE) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text        = ggplot2::element_text(colour = "black"),
      axis.line   = ggplot2::element_line(linewidth = 0.3),
      axis.ticks  = ggplot2::element_line(linewidth = 0.3),
      axis.text   = ggplot2::element_text(colour = "black", size = base_size),
      axis.title  = ggplot2::element_text(size = base_size),
      legend.key.size = grid::unit(0.3, "cm"),
      legend.title = ggplot2::element_text(size = base_size),
      legend.text  = ggplot2::element_text(size = base_size),
      plot.title   = ggplot2::element_text(size = base_size, face = "plain"),
      panel.grid   = ggplot2::element_blank()
    )
}

# Save a ggplot as SVG and 400-dpi PNG at its final size (never PDF).
save_fig <- function(plot, name, w, h, dir = FIG_DIR) {
  if (w > MAX_W) stop(sprintf("%s: width %.1f in exceeds the %.1f in page width", name, w, MAX_W))
  ensure_dir(dir)
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), plot, width = w, height = h,
                  units = "in", dpi = 400, bg = "white")
  ggplot2::ggsave(file.path(dir, paste0(name, ".svg")), plot, width = w, height = h,
                  units = "in", bg = "white")
  invisible(plot)
}

# Grid-based graphics (e.g. eulerr) do not go through ggsave.
save_grid <- function(obj, name, w, h, dir = FIG_DIR) {
  ensure_dir(dir)
  grDevices::png(file.path(dir, paste0(name, ".png")), width = w, height = h, units = "in",
                 res = 400, bg = "white"); print(obj); grDevices::dev.off()
  grDevices::svg(file.path(dir, paste0(name, ".svg")), width = w, height = h, bg = "white")
  print(obj); grDevices::dev.off()
  invisible(obj)
}

# Display names: the sQTL tables and the AP-MS bait list use RP11-77K12.1 for
# CFDP1 and NHERF2 for SLC9A3R2; all displayed output uses the HGNC symbol.
display_gene <- function(x) {
  x <- sub(",.*$", "", trimws(as.character(x)))
  dplyr::recode(x,
                "RP11-77K12.1"       = "CFDP1",
                "RP11-77K12-1"       = "CFDP1",
                "CFDP1,RP11-77K12.1" = "CFDP1",
                "NHERF2"             = "SLC9A3R2",
                "C1orf112"           = "FIRRM",
                .default = x)
}

# Significance labels used in tables and panels.
sig_label <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

# Assert that every required input exists before any work starts.
require_files <- function(...) {
  f <- unlist(list(...))
  missing <- f[!file.exists(f)]
  if (length(missing)) stop("Missing input file(s):\n  ", paste(missing, collapse = "\n  "))
  invisible(TRUE)
}
