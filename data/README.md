# Input data manifest

Every input read by the scripts, keyed by the `config.R` variable that names it. Paths are
relative to the repository root; `used by` lists the scripts that read the file. Provenance
states where the file comes from: *deposited* files are processed data generated for this study
and provided with the paper (SRA, PRIDE, GEO and the Supplementary Tables), *public* files are
downloaded from the named resource at the stated version, *experimental source data* are the
measurement tables behind the plotted panels, and *generated* files are produced by a script in
this repository.

Large public references (GENCODE v38 GTF, BioGRID 5.0.252 TAB3) are not redistributed; download
them to the listed path. The LeafCutter BED and count files are those of the VSMC sQTL study
(GEO GSE193817).

## sqtl/ — VSMC sQTL study

| file | variable | used by | description | provenance |
|---|---|---|---|---|
| `data/sqtl/WithFBS.bed` | `F_LEAFCUTTER_BED_PROLIF` | 02, 04 | LeafCutter junction BED (cluster identifier and coordinates) for proliferative VSMCs | deposited (this paper) — LeafCutter output of the Aherrahrou et al. 2023 sQTL study |
| `data/sqtl/WithFBS_cleaned.counts` | `F_LEAFCUTTER_COUNTS_PROLIF` | 02 | LeafCutter numerator/denominator junction counts per donor (proliferative) | deposited (this paper) — LeafCutter output of the Aherrahrou et al. 2023 sQTL study |
| `data/sqtl/WithoutFBS.bed` | `F_LEAFCUTTER_BED_QUIESC` | 02, 04 | LeafCutter cluster BED (all junctions per cluster) for quiescent VSMCs | deposited (this paper) — LeafCutter output of the Aherrahrou et al. 2023 sQTL study |
| `data/sqtl/WithoutFBS_cleaned.counts` | `F_LEAFCUTTER_COUNTS_QUIESC` | 02 | LeafCutter numerator/denominator junction counts per donor (quiescent) | deposited (this paper) — LeafCutter output of the Aherrahrou et al. 2023 sQTL study |
| `data/sqtl/general_sqtls_mapped_to_long_read.tsv` | `F_GENERAL_SQTL_LR` | 02, 03, 04 | All VSMC sQTL junctions mapped to long-read isoforms (chr, intron_start, intron_end, cluster_id); non-CAD background | deposited (this paper) — derived from the Aherrahrou et al. 2023 sQTL summary statistics (GEO GSE193817) and the long-read models |

## colocalization/ — CAD GWAS colocalization and variant annotation

| file | variable | used by | description | provenance |
|---|---|---|---|---|
| `data/colocalization/2026_CAD_sQTL_SNPs.csv` | `F_CAD_SQTL_SNPS_LD` | 01, 05 | LD expansion of the lead sQTL variants: proxy SNP (column 1) with its query lead variant (query_lead_snp_rsid), r2 >= 0.8 in 1000 Genomes EUR | deposited (this paper) — HaploReg v4.2 LD expansion of lead variants |
| `data/colocalization/2026_coloc_with_slopes.csv` | `F_COLOC_WITH_SLOPES` | 01, 05 | CAD GWAS-VSMC sQTL colocalization results, one row per colocalized LeafCutter junction x condition: sGene, cluster_id (chr:start:end:clu_N_strand:ENSG), lead sQTL rsid, sQTL slope and slope_se, ref/alt alleles, MAF, nominal/threshold/q-value statistics, Condition (Proliferative/Quiescent), eQTL_at_sGene, tensorqtl statistics | deposited (this paper) — Supplementary Tables S1/S2 source data; sQTL summary statistics from GEO GSE193817 (Aherrahrou et al. 2023) |
| `data/colocalization/2026_novel_only_updated_coloc.csv` | `F_COLOC_NOVEL_ONLY` | 10 | CAD GWAS-VSMC sQTL colocalizations with PPH4, cluster_id, Condition, sGene and a novelty flag (newly reported vs replicated); one row per colocalized junction x condition (Fig. 1b) | deposited (this paper) - Supplementary Tables S1/S2 source data |
| `data/colocalization/RAB5C_MST1_sQTL_plot_data.xlsx` | `F_RAB5C_MST1_PSI` | 10 | Donor x genotype x junction PSI for the LD-proxy SNPs rs2240010 (RAB5C) and rs13085791 (MST1); sheets RAB5C and MST1 (Extended Data Fig. 1d, 1f) | deposited (this paper) - sQTL genotypes and LeafCutter PSI from GEO GSE193817 |
| `data/colocalization/RT-PCR_primers.xlsx` | `F_RT_PCR_PRIMERS` | 05 | RT-PCR primer design with hg38 amplicon coordinates (Supplementary Table 6; sheet name containing Primer) | experimental source data — primer design workbook |
| `data/colocalization/Supplementary_Tables.xlsx` | `F_SUPP_TABLES_XLSX` | 01 | The manuscript Supplementary Tables workbook; this script reads sheets S1-new coloc results, S2-Aggregate coloc results, S3-LD snp list, S4-CAD_sQTL SNP overlap and S5-SpliceAI & VEP all data | deposited (this paper) — Supplementary Tables |
| `data/colocalization/sqtl_plots_long.csv` | `F_SQTL_PSI_LONG` | 10 | Per-donor scaled PSI by genotype for the plotted sQTLs (gene, rsid, psi_measure, genotype, psi) used for the PSI boxplots of Fig. 1d, 2b, 2d, 2e and Extended Data Fig. 1a/c/e, 2a-c | deposited (this paper) - sQTL genotypes and LeafCutter PSI from GEO GSE193817 |
| `data/colocalization/sqtl_snp_junction_lookup.csv` | `F_SQTL_SNP_JUNCTION` | 10 | Lookup of LeafCutter cluster and junction_key per plotted SNP and PSI measure (sqtl_snp_plotted, psi_measure, cluster, junction_key) for the panel subtitles | deposited (this paper) - derived from the colocalization results |

## gwas/ — CAD GWAS

| file | variable | used by | description | provenance |
|---|---|---|---|---|
| `data/gwas/Aragam_2022_S31_CAD_GWAS.csv` | `F_ARAGAM_S31_GENES` | 07 | Most likely causal genes at CAD GWAS loci (column most_likely_causal_gene); flags AP-MS preys encoded by independent CAD loci (Extended Data Fig. 7, Supplementary Table 17) | public — Aragam et al. 2022 Nat. Genet., Supplementary Table 31 |
| `data/gwas/aragam_cad_snps.tsv` | `F_ARAGAM_SNP_SUBSET` | 01 | Rows of the Aragam et al. (2022) multi-ancestry CAD GWAS summary statistics for the lead sQTL variants and their LD proxies; columns rsid, effect_allele, other_allele, beta, odds_ratio, p_value. Used to define the CAD risk allele per variant | public — GWAS Catalog GCST90132314 (GCST90132314.h.tsv.gz), subset to the rsIDs in 2026_coloc_with_slopes.csv (lead variants) and 2026_CAD_sQTL_SNPs.csv (LD proxies) |

## longread/ — VSMC long-read transcriptome and predicted proteome

| file | variable | used by | description | provenance |
|---|---|---|---|---|
| `data/longread/260115_sqanti_protein_w_cpm.csv` | `F_SQANTI_PROTEIN_CPM` | 02, 03, 05, 09, 10 | SQANTI Protein classification of the VSMC long-read predicted proteome with transcript category (tx_cat) and CPM per PacBio accession | deposited (this paper) — long-read proteogenomics pipeline output (SQANTI3 / SQANTI Protein on the VSMC Iso-Seq data; SRA BioProject) |
| `data/longread/SMC.demuxed_with_sample_data.tsv` | `F_DEMUX_SAMPLE_COUNTS` | 10 | Per-sample demultiplexed full-length read counts with donor (donor, fl_count, ...) from cDNA Cupcake demux_isoseq_with_genome.py | deposited (this paper) - PacBio long-read RNA-seq, NCBI SRA BioProject |
| `data/longread/SMC.mapped_fl_count.txt` | `F_MAPPED_FL_COUNTS` | 10 | Full-length read counts per collapsed transcript per sample (id, S01-S12; samples S01/S02 = donor 1, ...); Iso-Seq3 collapse output used as the SQANTI3 input | deposited (this paper) - PacBio long-read RNA-seq, NCBI SRA BioProject; Iso-Seq3 output |
| `data/longread/SMC.sqanti_protein_classification.tsv` | `F_SQANTI_PROTEIN_TSV` | 05, 10 | SQANTI Protein classification of the predicted VSMC proteome (pb, pr_splice_cat, pr_splice_subcat, pr_transcripts, tx_cat, ...) | deposited (this paper) — long-read proteogenomics pipeline (SQANTI Protein) |
| `data/longread/SMC_CORRECTED_best_orf.tsv` | `F_BEST_ORF_CPM` | 10 | Best ORF per transcript with gene and CPM (pb_acc, gene, CPM, ...) used for the per-gene isoform composition donuts | deposited (this paper) - long-read proteogenomics pipeline output |
| `data/longread/SMC_CORRECTED_orf_refined.fasta` | `F_SMC_PROTEIN_FASTA` | 02, 03 | Predicted VSMC protein sequences (one per long-read transcript) | deposited (this paper) — long-read proteogenomics pipeline output |
| `data/longread/SMC_CORRECTED_with_cds.gtf` | `F_SMC_GTF` | 02, 03, 04, 05 | Long-read VSMC transcript models with CDS (SQANTI3-corrected; long-read proteogenomics pipeline) | deposited (this paper) — SRA BioProject (processed long-read proteogenomics pipeline output) |
| `data/longread/SMC_CORRECTED_with_cds_exon_coords.txt` | `F_SMC_EXON_COORDS` | 05 | Per-exon table of the long-read models (trans_id = GENE/PB.x.y/CPM) | deposited (this paper) — derived from F_SMC_GTF by the long-read proteogenomics pipeline |
| `data/longread/SMC_CORRECTED_with_cds_intron_coords.txt` | `F_SMC_INTRON_COORDS` | 01, 02, 03 | Intron (junction) coordinates of every VSMC long-read transcript model: chr, intron_start, intron_end, strand, junc_id, trans_id (GENE/PB.x.y/N), gene_id | deposited (this paper) — derived from the long-read transcript GTF (SRA BioProject; long-read proteogenomics pipeline) |
| `data/longread/SMC_classification.txt` | `F_SQANTI3_CLASS` | 10 | SQANTI3 v4.2 transcript classification of the VSMC long-read transcriptome against GENCODE v38 (isoform, associated_gene, associated_transcript, structural_category, length, exons, all_canonical, RTS_stage, ...) | deposited (this paper) - PacBio long-read RNA-seq, NCBI SRA BioProject; SQANTI3 output |
| `data/longread/SMC_corrected.psl` | `F_SMC_PSL` | 02 | Exon block structure per long-read transcript (PSL) | deposited (this paper) — derived from F_SMC_GTF |

## reference/ — public reference annotation

| file | variable | used by | description | provenance |
|---|---|---|---|---|
| `data/reference/BIOGRID-ORGANISM-Homo_sapiens-5.0.252.tab3.txt` | `F_BIOGRID_TAB3` | 07, 09 | BioGRID human interaction records, TAB3 format; physical interactions are used for gene-level coverage and Entrez/symbol/synonym resolution of the bait genes | public — BioGRID release 5.0.252 (https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-5.0.252/) |
| `data/reference/all_SMC_genes_expressed_from_HPA.tsv` | `F_HPA_SMC_GENES` | 07 | Human Protein Atlas genes expressed in smooth muscle cells (columns Gene, Gene description); used as one of three VSMC-expression evidence sources for AP-MS preys | public — Human Protein Atlas, smooth muscle cell expression set (single-cell type atlas; downloaded June 2025) |
| `data/reference/gencode.v38.annotation.gtf` | `F_GENCODE_V38_GTF` | 02, 04, 05, 10 | GENCODE v38 primary annotation GTF (SNP location classification for Fig. 4a; transcripts per gene for Table S11) | public - GENCODE release 38 (https://www.gencodegenes.org) |
| `data/reference/gencode_v38_exon_coords.txt` | `F_GENCODE_V38_EXONS` | 05 | GENCODE v38 per-exon table in the same shape as the long-read exon table | public — derived from the GENCODE v38 GTF |
| `data/reference/gencode_v38_intron_coords_to_isoforms.txt` | `F_GENCODE_V38_INTRONS` | 01 | Intron coordinates of every GENCODE v38 transcript with the list of transcripts (gencode_v38_isoforms) and gene_id per intron | public — derived from the GENCODE v38 annotation GTF |
| `data/reference/uniprot_bed/{unipLocSignal;unipLocExtra;unipLocTransMemb;unipLocCytopl;unipDomain;unipInterest;unipStruct;unipDisulfBond;unipModif;unipRepeat}.bed` | `D_UNIPROT_BED` | 05 | UCSC UniProt track exports (hg38) for the browser feature tracks (10 BED files) | public — UCSC Table Browser UniProt tracks (hg38) |
| `data/reference/uniprot_combined_annot.tsv` | `F_UNIPROT_ANNOT` | 03 | UniProt feature intervals on hg38 (category; name; chrom; chromStart; chromEnd) | public — UniProt (UCSC UniProt genome annotation tracks) |

## protein_features/ — per-isoform protein feature predictions

| file | variable | used by | description | provenance |
|---|---|---|---|---|
| `data/protein_features/AAIndex_flexibility_granthamvolume_kytedoolittlehydropathy_zimmermanpolarity.csv` | `F_AAINDEX` | 03 | Per-residue AAindex scales: flexibility; Grantham volume; Kyte–Doolittle hydropathy; Zimmerman polarity | public — AAindex (Kawashima et al. 2008) |
| `data/protein_features/af3_metrics_per_isoform.csv` | `F_AF3_METRICS_ISO` | 03 | AlphaFold3 per-isoform metrics (pLDDT; radius of gyration; surface charge; disorder and secondary-structure fractions; PTM exposure) | generated by python/compute_af3_metrics.py |
| `data/protein_features/af3_metrics_per_pair.csv` | `F_AF3_METRICS_PAIR` | 03 | AlphaFold3 per-pair metrics (TM-score; alignment identity and gap statistics; junction-region metrics per side) | generated by python/compute_af3_metrics.py |
| `data/protein_features/hmmscan_pfam_SMC_lrRNA_hmmscan_domtblout.txt` | `F_HMMSCAN_DOMTBLOUT` | 03 | hmmscan domain table (Pfam-A v35.0) for the predicted proteome | deposited (this paper) — HMMER v3.3 hmmscan against Pfam-A v35.0 |
| `data/protein_features/iupred_smc_lr_all_residues.tsv` | `F_IUPRED_RESIDUES` | 03 | IUPred3 per-residue scores (pb_id; position; iupred_score) | deposited (this paper) — IUPred3 run on F_SMC_PROTEIN_FASTA |
| `data/protein_features/iupred_smc_lr_summary.csv` | `F_IUPRED_SUMMARY` | 03 | IUPred3 per-isoform summary (percent disordered; mean disorder score) | deposited (this paper) — IUPred3 run on F_SMC_PROTEIN_FASTA |
| `data/protein_features/netsurfp3_smc_master_residue_predictions.csv` | `F_NETSURFP_RESIDUES` | 03 | NetSurfP-3.0 per-residue predictions (pb_id; position; rsa; q3; disorder) | deposited (this paper) — NetSurfP-3.0 run on F_SMC_PROTEIN_FASTA |
| `data/protein_features/netsurfp3_smc_protein_features_slim.csv` | `F_NETSURFP_FEATURES` | 03 | NetSurfP-3.0 per-isoform summary (secondary-structure fractions; mean RSA; mean disorder) | deposited (this paper) — NetSurfP-3.0 run on F_SMC_PROTEIN_FASTA |
| `data/protein_features/smc_deeploc_master.csv` | `F_DEEPLOC_SMC` | 03 | DeepLoc 2.0 per-isoform predictions (Localizations; Signals; Membrane types; compartment probabilities) | deposited (this paper) — DeepLoc 2.0 run on F_SMC_PROTEIN_FASTA |
| `data/protein_features/smc_lr_exons_with_ptms.tsv` | `F_PTM_EXONS` | 03 | PTM-POSE exon table (event_id; Number of PTMs Affected) | deposited (this paper) — PTM-POSE projection of PhosphoSitePlus sites onto the long-read models |
| `data/protein_features/smc_lr_ptms_summary.tsv` | `F_PTM_SUMMARY` | 03 | PTM-POSE per-site table (Region ID; PTM Position in Isoform; Modification Class) | deposited (this paper) — PTM-POSE projection of PhosphoSitePlus sites onto the long-read models |

## proteomics/ — VSMC discovery and IS-PRM proteomics

| file | variable | used by | description | provenance |
|---|---|---|---|---|
| `data/proteomics/231023_SMC_hybrid_gencode_v43_PEPTIDES.csv` | `F_CLEAVER_DIGEST` | 05 | Cleaver in-silico tryptic digest of the hybrid VSMC long-read/GENCODE v43 proteome with isoform membership per peptide | deposited (this paper) — generated with Cleaver from the hybrid protein database |
| `data/proteomics/231030_rts_EDJ_hits_combined_POGO.bed` | `F_POGO_RTS_BED` | 05 | PoGo genomic coordinates of real-time-search (donor B94) peptides | deposited (this paper) — PRIDE (discovery search results mapped with PoGo) |
| `data/proteomics/231031_L101_B94_Z99_all_combined_Pogo.bed` | `F_POGO_DDA_ALL_BED` | 05 | PoGo genomic coordinates of all DDA peptides (donors L101 and Z99 and B94) | deposited (this paper) — PRIDE (discovery search results mapped with PoGo) |
| `data/proteomics/231031_dda_cad_peptides_POGO.bed` | `F_POGO_DDA_CAD_BED` | 05 | PoGo genomic coordinates of DDA peptides in CAD genes | deposited (this paper) — PRIDE (discovery search results mapped with PoGo) |
| `data/proteomics/240710_R_CAD173_rep_TOM_OUTPUT_include_NA_filtered_out.csv` | `F_ISPRM_DETECTED` | 05 | IS-PRM detected peptides after post-acquisition filtering (column Peptide) | deposited (this paper) — PRIDE (IS-PRM results; Alternative-Proteome-Detection-Project) |
| `data/proteomics/240710_R_CAD173_rep_TOM_POGO.bed` | `F_POGO_ISPRM_BED` | 05 | PoGo genomic coordinates of IS-PRM (Tomahto) peptides | deposited (this paper) — PRIDE (IS-PRM processed results; Alternative-Proteome-Detection-Project) |
| `data/proteomics/CAD173_info.csv` | `F_CAD173_TARGETS` | 05 | The 173 IS-PRM target peptides (column AAseq) with gene and specificity | deposited (this paper) — IS-PRM target list (Supplementary Table 12 source) |

## apms/ — AP-MS

| file | variable | used by | description | provenance |
|---|---|---|---|---|
| `data/apms/Bait_Protein_Sequences_CAD.fa` | `F_APMS_BAIT_FASTA` | 10 | Protein sequences of the 26 bait ORFs, headers carrying the PB accession (Table S8) | deposited (this paper) - plasmid construction (Supplementary Table 8) |
| `data/apms/CAD_APMS_Baits_risk_v_nonrisk.csv` | `F_APMS_BAIT_PAIRS` | 02, 06, 09, 10 | The 13 isoform pairs cloned as 3xFLAG AP-MS baits: gene symbol, risk PB accession, nonrisk PB accession, description (header carries a UTF-8 BOM) | experimental source data — bait selection for the AP-MS experiment (Supplementary Table 8) |
| `data/apms/SMC_specificity_table.csv` | `F_SMC_SPECIFICITY_MANUAL` | 07 | Manual VSMC-expression check for the 31 AP-MS preys absent from the HPA smooth-muscle set: PreyGene, SMC Protein (detected in VSMC Z99 discovery proteomics, YES/NO), SMC LR RNAseq (detected in VSMC long-read transcriptome, YES/NO) | experimental source data — manual lookup of each prey in the VSMC discovery proteomics (donor Z99, 28 fractions) and the VSMC long-read transcript catalog |
| `data/apms/comppass_output_10_CAD_set_1.csv` | `F_COMPPASS_SET` | 07 | CompPASS output for AP-MS experimental set 1: Bait, Prey, AvePSM, Z, WD, Entropy | deposited (this paper) — PRIDE (AP-MS search results scored with CompPASS) |
| `data/apms/comppass_output_10_CAD_set_2.csv` | `F_COMPPASS_SET` | 07 | CompPASS output for AP-MS experimental set 2 | deposited (this paper) — PRIDE (AP-MS search results scored with CompPASS) |
| `data/apms/comppass_output_10_CAD_set_3.csv` | `F_COMPPASS_SET` | 07 | CompPASS output for AP-MS experimental set 3 | deposited (this paper) — PRIDE (AP-MS search results scored with CompPASS) |
| `data/apms/intact_isoform_records.csv` | `F_INTACT_RECORDS` | 09 | Stored result of the IntAct PSICQUIC query per bait gene (gene, intact_isoform_record, intact_isoforms_seen); reflects IntAct content at the query date and is regenerated when QUERY_INTACT is TRUE | public — IntAct via PSICQUIC REST (https://www.ebi.ac.uk/intact/); stored query result generated by 09_bait_prior_evidence_biogrid_intact.Rmd |
| `data/apms/saint_express_output_10_CAD_set_1_isoforms.csv` | `F_SAINT_SET` | 07 | SAINTexpress output for AP-MS experimental set 1 (26 isoform baits + RAB11B; GFP/bead pulldowns as controls): Bait, Prey, PreyGene, per-replicate Spec, ctrlCounts, AvgSpec, FoldChange, SaintScore, BFDR | deposited (this paper) — PRIDE (AP-MS search results scored with SAINTexpress) |
| `data/apms/saint_express_output_10_CAD_set_2_isoforms.csv` | `F_SAINT_SET` | 07 | SAINTexpress output for AP-MS experimental set 2 | deposited (this paper) — PRIDE (AP-MS search results scored with SAINTexpress) |
| `data/apms/saint_express_output_10_CAD_set_3_isoforms.csv` | `F_SAINT_SET` | 07 | SAINTexpress output for AP-MS experimental set 3 | deposited (this paper) — PRIDE (AP-MS search results scored with SAINTexpress) |
| `data/apms/srp54_investigation.csv` | `F_SRP54_SPECTRAL` | 10 | Spectral counts of the SRP54 bait and SRP complex members in the nonrisk and risk SRP54 pulldowns and controls (PreyGene, status, Spec*, ctrlCounts*) (Fig. 6b, 6c, Table S19) | deposited (this paper) - AP-MS, ProteomeXchange PRIDE |

## experimental/ — experimental source data

| file | variable | used by | description | provenance |
|---|---|---|---|---|
| `data/experimental/ITGA1_spheroid_compression_raw.xlsx` | `F_ITGA1_SPHEROID_MECH` | 10 | Per-spheroid compression mechanics, first sheet, header on row 3, columns B:K (Type, Date [mddyy batch code], #, Youngs, Shear, Def, Circ) (Fig. 6r-t; Extended Data Fig. 8e; Table S20) | experimental source data - spheroid compression assay |
| `data/experimental/ITGA1_spheroid_raw_data.xlsx` | `F_ITGA1_SPHEROID_IMAGING` | 10 | ITGA1 spheroid imaging: sheets 'Actin per Membrane Scan' (Type, value) and 'SHG' (Type, Value) (Extended Data Fig. 8d, 8g) | experimental source data - immunofluorescence and second-harmonic generation imaging |
| `data/experimental/Isoform_gel_quant_plot_data.xlsx` | `F_GEL_DENSITOMETRY` | 10 | RT-PCR gel densitometry: sheets PARP12, WIPI1, DCLRE1B (Condition, Genotype, Isoform, Average (n=2)) and ITGA1_raw (Condition, Genotype, Raw) (Fig. 2g; Extended Data Fig. 2e, 2g, 2i) | experimental source data - junction-specific RT-PCR band intensities |
| `data/experimental/PARP12_OE_DESeq2_results.csv` | `F_PARP12_DESEQ2` | 10 | DESeq2 results, risk vs nonrisk PARP12 isoform overexpression RNA-seq; required columns gene, log2FoldChange, pvalue, padj (Fig. 2j significance stars) | experimental source data - DESeq2 on the OE RNA-seq (GEO) |
| `data/experimental/PARP12_OE_plots_raw_data.xlsx` | `F_PARP12_OE_DATA` | 10 | PARP12 isoform overexpression: sheets 'Quality control Isoform OE' (qPCR), 'Pathway Enrichment' (HALLMARK Pathway, FDR), 'Gene Diff Exp' (normalized expression per sample), 'NFkB ELISA' (Sample ID, GFP, Canonical OE, Alternate OE) (Fig. 2h-k) | experimental source data - qPCR, GSEA of the OE RNA-seq, p65 phospho-ELISA |

## Software-only dependency

`08_apms_go_enrichment.Rmd` requires the Bioconductor annotation package `org.Hs.eg.db` (version recorded in each script's session information).
