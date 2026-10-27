#!/usr/bin/env python
"""
Compute AlphaFold3 structural metrics for the risk/nonrisk isoform pairs.

Runs between 02_junction_isoform_mapping_and_pairs.Rmd and
03_protein_feature_annotation.Rmd. For every predicted isoform with an AlphaFold3
model it computes whole-isoform metrics (mean pLDDT, radius of gyration, surface
charge, secondary-structure fractions, pLDDT-defined disorder, exposed PTM sites,
AF3 confidence scores); for every pair it computes the structural similarity
(TM-align), a global sequence alignment summary (identity, gaps, gap position
class) and junction-window metrics on the amino-acid windows exported by script 02.

Outputs (read by 03_protein_feature_annotation.Rmd):
    data/protein_features/af3_metrics_per_isoform.csv
    data/protein_features/af3_metrics_per_pair.csv

Requirements: Python >= 3.9 with numpy, pandas, gemmi, biopython; the TM-align and
DSSP (mkdssp) executables on PATH. AlphaFold3 models are expected as
<AF3_OUTPUT_DIR>/<PB accession>/<run>/..._model.cif with the accompanying
*_summary_confidences.json and *_confidences.json.

Usage:
    python python/compute_af3_metrics.py --af3-dir /path/to/af3_outputs \
        --manifest /path/to/isoform_manifest.csv [--repo-root .]

The manifest is a CSV with one column `pb_acc` listing every isoform to process
(all members of the unified pairs).
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import gemmi
from Bio import Align

# ---------------------------------------------------------------------------
# Parameters (as used for the manuscript)
# ---------------------------------------------------------------------------
RSA_SURFACE_CUTOFF = 0.20   # residues with relative solvent accessibility above this are "exposed"
PLDDT_IDR_CUTOFF = 70       # residues with pLDDT below this are counted as disordered
JUNCTION_MIN_LEN = 5        # minimum junction-window length (aa) for junction metrics
CHARGED_RESIDUES = {"ARG": +1, "LYS": +1, "HIS": +1, "ASP": -1, "GLU": -1}
TMALIGN = "TMalign"
DSSP = "mkdssp"

MAX_ASA = {  # Tien et al. 2013 theoretical maxima, used to convert ASA to RSA
    "ALA": 121, "ARG": 265, "ASN": 187, "ASP": 187, "CYS": 148, "GLN": 214, "GLU": 214,
    "GLY": 97, "HIS": 216, "ILE": 195, "LEU": 191, "LYS": 230, "MET": 203, "PHE": 228,
    "PRO": 154, "SER": 143, "THR": 163, "TRP": 264, "TYR": 255, "VAL": 165,
}
ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


# ---------------------------------------------------------------------------
# AlphaFold3 output discovery and parsing
# ---------------------------------------------------------------------------
def find_af3_outputs(af3_dir: Path, pb_acc: str) -> Optional[dict]:
    pb_dir = af3_dir / pb_acc
    if not pb_dir.is_dir():
        return None
    subdirs = [d for d in pb_dir.iterdir() if d.is_dir() and not d.name.startswith("seed-")]
    if not subdirs:
        return None
    inner = subdirs[0]
    cif_files = list(inner.glob("*_model.cif"))
    if not cif_files:
        return None
    return {
        "cif": cif_files[0],
        "summary_json": next(iter(inner.glob("*_summary_confidences.json")), None),
        "confidences_json": next(iter(inner.glob("*_confidences.json")), None),
    }


def parse_cif(cif_path: Path) -> dict:
    """Return the C-alpha trace of the first chain with per-residue pLDDT (B-factor field)."""
    structure = gemmi.read_structure(str(cif_path))
    structure.setup_entities()
    chain = structure[0][0]
    residues, seq_chars = [], []
    for res in chain:
        info = gemmi.find_tabulated_residue(res.name)
        one = info.one_letter_code.upper() if info else "X"
        ca = res.find_atom("CA", "*", gemmi.Element("C"))
        if ca is None:
            continue
        residues.append({"resnum": res.seqid.num, "resname": res.name,
                         "ca_xyz": (ca.pos.x, ca.pos.y, ca.pos.z), "plddt": ca.b_iso, "one": one})
        seq_chars.append(one)
    return {"residues": residues, "sequence": "".join(seq_chars), "n_res": len(residues)}


def run_dssp(cif_path: Path, pb_acc: str, cache_dir: Path) -> Optional[pd.DataFrame]:
    """Secondary structure and relative solvent accessibility per residue (cached)."""
    cache_file = cache_dir / f"{pb_acc}_dssp.csv"
    if cache_file.exists():
        return pd.read_csv(cache_file)
    dssp_out = cache_dir / f"{pb_acc}.dssp"
    try:
        subprocess.run([DSSP, "--output-format", "dssp", str(cif_path), str(dssp_out)],
                       check=True, capture_output=True, timeout=300)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None
    rows = []
    with open(dssp_out) as fh:
        body = False
        for line in fh:
            if line.startswith("  #  RESIDUE"):
                body = True
                continue
            if not body or len(line) < 38:
                continue
            try:
                rows.append({"resnum": int(line[5:10]), "one": line[13],
                             "ss": line[16].strip() or "C", "asa": int(line[34:38])})
            except ValueError:
                continue
    if not rows:
        return None
    df = pd.DataFrame(rows)
    df["resname"] = df["one"].map(ONE_TO_THREE)
    df["rsa"] = df["asa"] / df["resname"].map(MAX_ASA)
    df["ss_simple"] = df["ss"].apply(lambda s: "H" if s in "HGI" else ("E" if s in "EB" else "L"))
    df.to_csv(cache_file, index=False)
    return df


# ---------------------------------------------------------------------------
# Whole-isoform metrics
# ---------------------------------------------------------------------------
def radius_of_gyration(residues) -> float:
    coords = np.array([r["ca_xyz"] for r in residues])
    return float(np.sqrt(np.sum((coords - coords.mean(axis=0)) ** 2, axis=1).mean()))


def surface_charge(residues, dssp) -> float:
    """Net charge of exposed Arg/Lys/His (+1) and Asp/Glu (-1) residues."""
    rsa = dssp.set_index("resnum")["rsa"]
    return float(sum(CHARGED_RESIDUES[r["resname"]] for r in residues
                     if r["resname"] in CHARGED_RESIDUES
                     and r["resnum"] in rsa.index and rsa.loc[r["resnum"]] > RSA_SURFACE_CUTOFF))


def ss_percentages(dssp) -> dict:
    n = len(dssp)
    if n == 0:
        return {"pct_helix": np.nan, "pct_sheet": np.nan, "pct_loop": np.nan}
    counts = dssp["ss_simple"].value_counts()
    return {"pct_helix": 100 * counts.get("H", 0) / n,
            "pct_sheet": 100 * counts.get("E", 0) / n,
            "pct_loop": 100 * counts.get("L", 0) / n}


def ptm_exposure(pb_acc, dssp, ptm_df) -> dict:
    pb_ptms = ptm_df[ptm_df["pb_acc"] == pb_acc]
    if pb_ptms.empty:
        return {"ptm_total": 0, "ptm_buried": 0, "ptm_exposed": 0, "ptm_pct_exposed": np.nan}
    rsa = dssp.set_index("resnum")["rsa"]
    exposed = buried = 0
    for pos in pb_ptms["position"].astype(int):
        if pos not in rsa.index:
            continue
        if rsa.loc[pos] > RSA_SURFACE_CUTOFF:
            exposed += 1
        else:
            buried += 1
    total = exposed + buried
    return {"ptm_total": total, "ptm_buried": buried, "ptm_exposed": exposed,
            "ptm_pct_exposed": 100 * exposed / total if total else np.nan}


def parse_af3_summary(path) -> dict:
    keys = ["ptm", "iptm", "fraction_disordered", "has_clash", "ranking_score"]
    if path is None or not Path(path).exists():
        return {k: np.nan for k in keys}
    with open(path) as fh:
        d = json.load(fh)
    return {k: d.get(k, np.nan) for k in keys}


def parse_af3_confidences(path) -> dict:
    if path is None or not Path(path).exists():
        return {"mean_pae": np.nan}
    with open(path) as fh:
        pae = json.load(fh).get("pae")
    if pae is None:
        return {"mean_pae": np.nan}
    arr = np.array(pae)
    return {"mean_pae": float(arr[~np.eye(arr.shape[0], dtype=bool)].mean())}


# ---------------------------------------------------------------------------
# Pair-level metrics
# ---------------------------------------------------------------------------
def run_tmalign(cif_a: Path, cif_b: Path) -> Optional[dict]:
    """TM-score of the pair, taking the larger of the two length normalisations."""
    try:
        res = subprocess.run([TMALIGN, str(cif_a), str(cif_b)], capture_output=True, text=True,
                             timeout=300, check=True)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None
    tm_a = tm_b = np.nan
    for line in res.stdout.splitlines():
        if "TM-score=" in line and "normalized by length of Structure_1" in line:
            tm_a = float(line.split("TM-score=")[1].split()[0])
        elif "TM-score=" in line and "normalized by length of Structure_2" in line:
            tm_b = float(line.split("TM-score=")[1].split()[0])
    return {"tm_score": max([x for x in (tm_a, tm_b) if not np.isnan(x)], default=np.nan)}


def sequence_alignment_metrics(seq_a: str, seq_b: str) -> dict:
    """Global alignment (match 1, mismatch 0, gap open -1, extend -0.5) summarised as
    identity within the aligned region, aligned length, gap counts, coverage of the shorter
    sequence, length ratio and the position of substantial gaps (> 5 % of the longer
    sequence) relative to the longer sequence."""
    nan_result = {"pct_identity": np.nan, "aligned_length": 0, "n_gap_openings": 0,
                  "total_gap_length": 0, "query_coverage_pct": np.nan, "length_ratio": np.nan,
                  "gap_position_class": "alignment_failed"}
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 1, 0
    aligner.open_gap_score, aligner.extend_gap_score = -1, -0.5
    try:
        aln = aligner.align(seq_a, seq_b)[0]
    except Exception:
        return nan_result
    len_a, len_b = len(seq_a), len(seq_b)
    if len_a == 0 or len_b == 0:
        return nan_result
    blocks_a, blocks_b = aln.aligned
    if len(blocks_a) == 0:
        return nan_result
    n_match = n_aligned = 0
    for (a0, a1), (b0, b1) in zip(blocks_a, blocks_b):
        for i in range(a1 - a0):
            n_aligned += 1
            n_match += seq_a[a0 + i] == seq_b[b0 + i]
    n_gaps = total_gap = 0
    for i in range(len(blocks_a) - 1):
        for gap in (blocks_a[i + 1][0] - blocks_a[i][1], blocks_b[i + 1][0] - blocks_b[i][1]):
            if gap > 0:
                n_gaps += 1
                total_gap += gap
    for lead in (blocks_a[0][0], blocks_b[0][0]):
        if lead > 0:
            n_gaps += 1
            total_gap += lead
    for trail in (len_a - blocks_a[-1][1], len_b - blocks_b[-1][1]):
        if trail > 0:
            n_gaps += 1
            total_gap += trail
    ref_blocks, ref_len = (blocks_a, len_a) if len_a >= len_b else (blocks_b, len_b)
    n_term = ref_blocks[0][0]
    c_term = ref_len - ref_blocks[-1][1]
    internal = sum(ref_blocks[i + 1][0] - ref_blocks[i][1] for i in range(len(ref_blocks) - 1))
    thr = 0.05 * ref_len
    substantial = [n_term > thr, c_term > thr, internal > thr]
    if sum(substantial) == 0:
        gap_class = "minimal"
    elif sum(substantial) > 1:
        gap_class = "scattered"
    elif n_term > thr:
        gap_class = "N-terminal"
    elif c_term > thr:
        gap_class = "C-terminal"
    else:
        gap_class = "internal"
    shorter, longer = min(len_a, len_b), max(len_a, len_b)
    return {"pct_identity": 100 * n_match / max(n_aligned, 1), "aligned_length": n_aligned,
            "n_gap_openings": n_gaps, "total_gap_length": total_gap,
            "query_coverage_pct": 100 * n_aligned / shorter, "length_ratio": shorter / longer,
            "gap_position_class": gap_class}


def junction_window_metrics(residues, dssp, aa_start, aa_end) -> Optional[dict]:
    """Secondary structure, surface charge and disorder within [aa_start, aa_end]."""
    if pd.isna(aa_start) or pd.isna(aa_end):
        return None
    aa_start, aa_end = int(aa_start), int(aa_end)
    length = aa_end - aa_start + 1
    if length < JUNCTION_MIN_LEN:
        return None
    window = set(range(aa_start, aa_end + 1))
    res = [r for r in residues if r["resnum"] in window]
    ds = dssp[dssp["resnum"].isin(window)]
    if not res or ds.empty:
        return None
    ss = ss_percentages(ds)
    plddt = np.array([r["plddt"] for r in res])
    return {"junc_length": length, "junc_pct_helix": ss["pct_helix"], "junc_pct_sheet": ss["pct_sheet"],
            "junc_surface_charge": surface_charge(res, ds),
            "junc_idr_pct": float(100 * (plddt < PLDDT_IDR_CUTOFF).mean()),
            "junc_mean_plddt": float(plddt.mean())}


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def load_aa_mapping(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df[df["status"].str.startswith("ok", na=False)].copy()


def lookup_aa_window(aa_map, pb_acc, junction_key):
    hit = aa_map[(aa_map["pb_acc"] == pb_acc) & (aa_map["junction_key"] == junction_key)]
    if hit.empty:
        return None, None
    return hit.iloc[0]["aa_window_start"], hit.iloc[0]["aa_window_end"]


def load_ptm_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["pb_acc"] = df["Region ID"].str.extract(r"(PB\.[0-9]+\.[0-9]+)")
    df = df.rename(columns={"PTM Position in Isoform": "position"}).dropna(subset=["pb_acc", "position"])
    df["position"] = df["position"].astype(int)
    return df[["pb_acc", "position", "Modification Class"]]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def process_isoform(pb_acc, af3_dir, cache_dir, ptm_df):
    outs = find_af3_outputs(af3_dir, pb_acc)
    if outs is None:
        return None
    parsed = parse_cif(outs["cif"])
    dssp = run_dssp(outs["cif"], pb_acc, cache_dir)
    if dssp is None:
        return None
    record = {"pb_acc": pb_acc, "n_residues": parsed["n_res"], "sequence": parsed["sequence"],
              "cif_path": str(outs["cif"]),
              "mean_plddt": float(np.mean([r["plddt"] for r in parsed["residues"]])),
              "radius_of_gyration": radius_of_gyration(parsed["residues"]),
              "surface_charge": surface_charge(parsed["residues"], dssp),
              "idr_fraction_pct": float(100 * (np.array([r["plddt"] for r in parsed["residues"]])
                                               < PLDDT_IDR_CUTOFF).mean())}
    record.update(ss_percentages(dssp))
    record.update(ptm_exposure(pb_acc, dssp, ptm_df))
    record.update(parse_af3_summary(outs["summary_json"]))
    record.update(parse_af3_confidences(outs["confidences_json"]))
    with open(cache_dir / f"{pb_acc}_residues.json", "w") as fh:
        json.dump({"residues": parsed["residues"], "sequence": parsed["sequence"]}, fh)
    return record


def process_pair(row, iso_lookup, aa_map, cache_dir):
    risk_pb, nonrisk_pb = row["risk_pb"], row["nonrisk_pb"]
    junction, pair_type = row["cad_sqtl_junction"], row["pair_type"]
    risk, nonrisk = iso_lookup.get(risk_pb), iso_lookup.get(nonrisk_pb)
    if risk is None or nonrisk is None:
        return None
    try:
        with open(cache_dir / f"{risk_pb}_residues.json") as fh:
            risk_c = json.load(fh)
        with open(cache_dir / f"{nonrisk_pb}_residues.json") as fh:
            nr_c = json.load(fh)
    except FileNotFoundError:
        return None
    risk_dssp = run_dssp(Path(risk["cif_path"]), risk_pb, cache_dir)
    nr_dssp = run_dssp(Path(nonrisk["cif_path"]), nonrisk_pb, cache_dir)
    tm = run_tmalign(Path(risk["cif_path"]), Path(nonrisk["cif_path"])) or {}
    aln = sequence_alignment_metrics(risk_c["sequence"], nr_c["sequence"])

    # Junction-window metrics use the amino-acid windows from script 02. A salvaged
    # counterpart does not contain the focal junction and therefore has no window;
    # manual AP-MS pairs have no junction context. Both are left as NA.
    junc = {}
    if pair_type != "manual_apms" and junction not in (None, "manual_apms"):
        for side, pb, cache, dssp in (("risk", risk_pb, risk_c, risk_dssp),
                                      ("nonrisk", nonrisk_pb, nr_c, nr_dssp)):
            start, end = lookup_aa_window(aa_map, pb, junction)
            if start is not None:
                m = junction_window_metrics(cache["residues"], dssp, start, end)
                if m:
                    junc.update({f"junc_{side}_{k}": v for k, v in m.items()})

    record = {"gene_name": row["gene_name"], "junction": junction, "pair_type": pair_type,
              "risk_pb": risk_pb, "nonrisk_pb": nonrisk_pb, "tm_score": tm.get("tm_score", np.nan)}
    record.update(aln)
    record.update({
        "delta_length": risk["n_residues"] - nonrisk["n_residues"],
        "delta_rg": risk["radius_of_gyration"] - nonrisk["radius_of_gyration"],
        "delta_surface_charge": risk["surface_charge"] - nonrisk["surface_charge"],
        "delta_pct_helix": risk["pct_helix"] - nonrisk["pct_helix"],
        "delta_pct_sheet": risk["pct_sheet"] - nonrisk["pct_sheet"],
        "delta_idr_pct": risk["idr_fraction_pct"] - nonrisk["idr_fraction_pct"],
        "delta_ptm_pct_exposed": risk.get("ptm_pct_exposed", np.nan) - nonrisk.get("ptm_pct_exposed", np.nan),
        "delta_mean_plddt": risk["mean_plddt"] - nonrisk["mean_plddt"],
    })
    record.update(junc)
    return record


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--af3-dir", required=True, type=Path, help="directory of AlphaFold3 model outputs")
    ap.add_argument("--manifest", required=True, type=Path, help="CSV with a pb_acc column")
    ap.add_argument("--repo-root", default=Path("."), type=Path)
    args = ap.parse_args()

    root = args.repo_root
    pairs_file = root / "results/02_junction_isoform_mapping_and_pairs/leafcutter_isoform_pairs_unified_flagged.csv"
    aa_map_file = root / "results/02_junction_isoform_mapping_and_pairs/junction_aa_mapping.csv"
    ptm_file = root / "data/protein_features/smc_lr_ptms_summary.tsv"
    out_iso = root / "data/protein_features/af3_metrics_per_isoform.csv"
    out_pair = root / "data/protein_features/af3_metrics_per_pair.csv"
    out_fail = root / "data/protein_features/af3_metrics_failures.csv"
    cache_dir = root / "results/af3_metric_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest)
    target_pbs = manifest["pb_acc"].unique().tolist()
    unified = pd.read_csv(pairs_file)
    aa_map = load_aa_mapping(aa_map_file)
    ptm_df = load_ptm_table(ptm_file)
    print(f"{len(target_pbs)} isoforms, {len(unified)} pairs, {len(aa_map)} junction windows")

    iso_results, failures = [], []
    for pb in target_pbs:
        rec = process_isoform(pb, args.af3_dir, cache_dir, ptm_df)
        (iso_results if rec else failures).append(
            rec or {"pb_acc": pb, "stage": "isoform", "reason": "no_model_or_dssp_failure"})
    pd.DataFrame(iso_results).to_csv(out_iso, index=False)
    iso_lookup = {r["pb_acc"]: r for r in iso_results}

    pair_results = []
    for row in unified.to_dict("records"):
        rec = process_pair(row, iso_lookup, aa_map, cache_dir)
        if rec is None:
            failures.append({"pb_acc": f"{row['risk_pb']}|{row['nonrisk_pb']}", "stage": "pair",
                             "reason": f"missing_isoform ({row['pair_type']})"})
        else:
            pair_results.append(rec)
    pd.DataFrame(pair_results).to_csv(out_pair, index=False)
    if failures:
        pd.DataFrame(failures).to_csv(out_fail, index=False)
    print(f"per-isoform: {len(iso_results)} rows -> {out_iso}")
    print(f"per-pair:    {len(pair_results)} rows -> {out_pair}")
    if failures:
        print(f"failures:    {len(failures)} -> {out_fail}")


if __name__ == "__main__":
    main()
