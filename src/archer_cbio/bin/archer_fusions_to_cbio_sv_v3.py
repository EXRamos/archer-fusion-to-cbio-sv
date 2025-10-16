#!/usr/bin/env python3
"""
archer_fusions_to_cbio_sv_v3.py

Enhancements over v2:
- Default filtering: EXCLUDE intronic events and EXCLUDE out-of-frame events
  (override with --keep-intron and/or --keep-out-of-frame).
- Adds a "Percent_Coverage" column sourced from best-available Archer fields
  to help compare AF-like values across datasets.

Usage:
  python archer_fusions_to_cbio_sv_v3.py -o data_sv.tsv sample.full_results.txt
  python archer_fusions_to_cbio_sv_v3.py -o data_sv.tsv --keep-intron sample.txt
  python archer_fusions_to_cbio_sv_v3.py -o data_sv.tsv --keep-out-of-frame sample.txt
  python archer_fusions_to_cbio_sv_v3.py -o data_sv.tsv --keep-intron --keep-out-of-frame sample.txt

Other options:
  --ncbi-build  GRCh37 (default) or GRCh38
  --strong-only include only FCA_*_STRONG_EVIDENCE_ABERRATION=TRUE

Notes:
- "Percent_Coverage" is populated by checking several Archer keys at the FCA or FC level, e.g.:
    FCA_{idx}_PERCENT_FUSION
    FCA_{idx}_PERCENT_GSP2_COVERAGE
    FC_{n}_PERCENT_FUSION
    FC_{n}_PERCENT_RNA_READS
  The first present key (in that priority order) is used; otherwise left blank.
- We still emit one row per supporting FC isoform (or one per FCA with counts=0 if no FC IDs found).
"""

import re
import csv
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ANNO_RE = re.compile(
    r"^(?P<gene>[A-Za-z0-9._-]+)\((?P<strand>[+-])\)"
    r"(?P<tx>(?:[A-Z]{2}_\d+(?:\.\d+)?(?:/[A-Z]{2}_\d+(?:\.\d+)?)*)?)\|"
    r"(?P<region>exon|intron):(?P<region_no>\d+)\|"
    r"(?P<chr1>chr[\w]+):(?P<pos1>\d+)(?:,(?P<chr2>chr[\w]+):(?P<pos2>\d+))?$"
)

GA_RE = re.compile(r"^(chr[\w]+):(\d+):(\d+):([+-])$")

HEADER = [
    "Sample_Id","SV_Status",
    "Site1_Hugo_Symbol","Site1_Ensembl_Transcript_Id","Site1_Entrez_Gene_Id",
    "Site1_Region_Number","Site1_Region","Site1_Chromosome","Site1_Contig","Site1_Position","Site1_Description",
    "Site2_Hugo_Symbol","Site2_Ensembl_Transcript_Id","Site2_Entrez_Gene_Id",
    "Site2_Region_Number","Site2_Region","Site2_Chromosome","Site2_Contig","Site2_Position","Site2_Description",
    "Site2_Effect_On_Frame","NCBI_Build","Class",
    "Tumor_Split_Read_Count","Tumor_Paired_End_Read_Count",
    "Event_Info","Annotation","DNA_Support","RNA_Support","SV_Length",
    "Normal_Read_Count","Tumor_Read_Count","Normal_Variant_Count","Tumor_Variant_Count","Normal_Paired_End_Read_Count",
    "Percent_Coverage"
]

def parse_annotation(s: str) -> Dict[str, Optional[str]]:
    s = (s or "").strip()
    m = ANNO_RE.match(s)
    if not m:
        m2 = re.match(
            r"^(?P<gene>[A-Za-z0-9._-]+)\((?P<strand>[+-])\)\|(?P<region>exon|intron):(?P<region_no>\d+)\|(?P<chr1>chr[\w]+):(?P<pos1>\d+)$",
            s,
        )
        if not m2:
            return {"gene": None, "strand": None, "tx": None, "region": None, "region_no": None, "chr": None, "pos": None, "raw": s}
        gd = m2.groupdict()
        return {"gene": gd["gene"], "strand": gd["strand"], "tx": None, "region": gd["region"],
                "region_no": gd["region_no"], "chr": gd["chr1"], "pos": gd["pos1"], "raw": s}
    gd = m.groupdict()
    chrom = gd.get("chr2") or gd.get("chr1")
    pos = gd.get("pos2") or gd.get("pos1")
    tx = gd.get("tx").lstrip("/") if gd.get("tx") else None
    return {
        "gene": gd.get("gene"),
        "strand": gd.get("strand"),
        "tx": tx,
        "region": gd.get("region"),
        "region_no": gd.get("region_no"),
        "chr": chrom,
        "pos": pos,
        "raw": s,
    }

def score_tx(s: str) -> int:
    if "NM_" in s: return 3
    if "XM_" in s: return 2
    if "NR_" in s: return 1
    return 0

def pick_best_fca_annotation(kv: Dict[str, str], idx: int, side: int) -> Optional[str]:
    prefix = f"FCA_{idx}_ANNOTATION_{side}-"
    anns = [(k, v) for k, v in kv.items() if k.startswith(prefix)]
    if not anns:
        return None
    anns.sort(key=lambda kvp: (-score_tx(kvp[1]), kvp[0]))
    return anns[0][1]

def load_kv(p: Path) -> Dict[str, str]:
    d = {}
    with p.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "\t" not in line: 
                continue
            k, v = line.rstrip("\n").split("\t", 1)
            d[k] = v
    return d

def fc_ga_list(kv: Dict[str, str], fc_id: str):
    out = []
    for k, v in kv.items():
        m = re.match(rf"^{re.escape(fc_id)}_GENOMIC_ANNOTATION_(\d+)$", k)
        if m:
            out.append((v, int(m.group(1)), GA_RE.match(v.strip())))
    out.sort(key=lambda t: t[1])
    return out

def fallback_chr_pos_from_ga(kv: Dict[str, str], fc_id: str, desired_chr: Optional[str]):
    for val, order, mm in fc_ga_list(kv, fc_id):
        if not mm: continue
        chr_, p1, p2, strand = mm.groups()
        if desired_chr and chr_ != desired_chr:
            continue
        return chr_, p2
    for val, order, mm in fc_ga_list(kv, fc_id):
        if mm:
            chr_, p1, p2, strand = mm.groups()
            return chr_, p2
    return "", ""

def get_percent_coverage(kv: Dict[str, str], fca_idx: int, fc_id: Optional[str]) -> str:
    for k in (f"FCA_{fca_idx}_PERCENT_FUSION", f"FCA_{fca_idx}_PERCENT_GSP2_COVERAGE"):
        if k in kv and kv[k] != "":
            return kv[k]
    if fc_id:
        for k in (f"{fc_id}_PERCENT_FUSION", f"{fc_id}_PERCENT_RNA_READS"):
            if k in kv and kv[k] != "":
                return kv[k]
    return ""

def archer_to_cbio_sv(kv: Dict[str, str], ncbi_build="GRCh37", strong_only=False,
                      keep_intron=False, keep_out_of_frame=False) -> List[Dict[str, str]]:
    samp = kv.get("SAMPLE_NAME", "UNKNOWN_SAMPLE")
    rows: List[Dict[str, str]] = []

    fca_idxs = sorted(int(m.group(1)) for k in kv for m in [re.match(r"^FCA_(\d+)_GENES_UNIQUE$", k)] if m)
    for idx in fca_idxs:
        genes = kv.get(f"FCA_{idx}_GENES_UNIQUE", "")
        if ":" not in genes:
            continue
        g1, g2 = genes.split(":", 1)

        strong = kv.get(f"FCA_{idx}_STRONG_EVIDENCE_ABERRATION", "FALSE").upper() == "TRUE"
        if strong_only and not strong:
            continue

        inframe_flag = kv.get(f"FCA_{idx}_HAS_INFRAME_TRANSLATION", "").upper()
        inframe = inframe_flag == "TRUE"

        ann1_str = pick_best_fca_annotation(kv, idx, side=1)
        ann2_str = pick_best_fca_annotation(kv, idx, side=2)
        a1 = parse_annotation(ann1_str or "")
        a2 = parse_annotation(ann2_str or "")

        # Default filtering: drop any intronic OR out-of-frame events unless overridden
        region1 = (a1.get("region") or "").lower()
        region2 = (a2.get("region") or "").lower()
        if not keep_intron and ("intron" in region1 or "intron" in region2):
            continue
        if not keep_out_of_frame and (not inframe):
            continue

        fc_ids = [kv[k] for k in kv if re.match(rf"^FCA_{idx}_SUPPORTING_ISOFORM_ID_\d+$", k) and kv[k].startswith("FC_")]
        if not fc_ids and f"FC_{idx}_EITHER_R1_OR_R2" in kv:
            fc_ids = [f"FC_{idx}"]
        if not fc_ids:
            fc_ids = [None]

        for fc in fc_ids:
            def get_count(name: str) -> int:
                if not fc: return 0
                val = kv.get(f"{fc}_{name}", "0")
                try: return int(float(val))
                except: return 0

            split_reads = get_count("EITHER_R1_OR_R2")
            pair_reads  = get_count("BOTH_R1_AND_R2")

            site1_chr = a1.get("chr") or ""
            site1_pos = a1.get("pos") or ""
            site2_chr = a2.get("chr") or ""
            site2_pos = a2.get("pos") or ""

            if (not site1_chr or not site1_pos) and fc:
                chr_, pos_ = fallback_chr_pos_from_ga(kv, fc, a1.get("chr"))
                site1_chr = site1_chr or chr_
                site1_pos = site1_pos or pos_
            if (not site2_chr or not site2_pos) and fc:
                chr_, pos_ = fallback_chr_pos_from_ga(kv, fc, a2.get("chr"))
                if chr_ == site1_chr and pos_ == site1_pos:
                    chr2, pos2 = fallback_chr_pos_from_ga(kv, fc, None)
                    chr_, pos_ = chr2, pos2
                site2_chr = site2_chr or chr_
                site2_pos = site2_pos or pos_

            sv_class = "Translocation" if g1 != g2 else "Intragenic"
            percent_cov = get_percent_coverage(kv, idx, fc)

            rows.append({
                "Sample_Id": samp,
                "SV_Status": "SOMATIC",
                "Site1_Hugo_Symbol": g1 or (a1.get("gene") or ""),
                "Site1_Ensembl_Transcript_Id": "",
                "Site1_Entrez_Gene_Id": "",
                "Site1_Region_Number": a1.get("region_no") or "",
                "Site1_Region": (a1.get("region") or "").capitalize() if a1.get("region") else "",
                "Site1_Chromosome": site1_chr,
                "Site1_Contig": "",
                "Site1_Position": site1_pos,
                "Site1_Description": a1.get("raw") or "",
                "Site2_Hugo_Symbol": g2 or (a2.get("gene") or ""),
                "Site2_Ensembl_Transcript_Id": "",
                "Site2_Entrez_Gene_Id": "",
                "Site2_Region_Number": a2.get("region_no") or "",
                "Site2_Region": (a2.get("region") or "").capitalize() if a2.get("region") else "",
                "Site2_Chromosome": site2_chr,
                "Site2_Contig": "",
                "Site2_Position": site2_pos,
                "Site2_Description": a2.get("raw") or "",
                "Site2_Effect_On_Frame": "In-frame" if inframe else "Out-of-frame",
                "NCBI_Build": ncbi_build,
                "Class": sv_class,
                "Tumor_Split_Read_Count": str(split_reads),
                "Tumor_Paired_End_Read_Count": str(pair_reads),
                "Event_Info": f"Fusion: {g1}-{g2}",
                "Annotation": "; ".join([a1.get("raw") or "", a2.get("raw") or ""]).strip("; ").strip(),
                "DNA_Support": "No",
                "RNA_Support": "Yes" if split_reads > 0 else "No",
                "SV_Length": "",
                "Normal_Read_Count": "",
                "Tumor_Read_Count": "",
                "Normal_Variant_Count": "",
                "Tumor_Variant_Count": "",
                "Normal_Paired_End_Read_Count": "",
                "Percent_Coverage": percent_cov,
            })
    return rows

def write_tsv(rows: List[Dict[str, str]], out_path: Path):
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)

def main():
    ap = argparse.ArgumentParser(description="Convert Archer full_results to cBioPortal data_sv.tsv (v3 with filters & percent coverage)")
    ap.add_argument("-o", "--out", required=True, help="Output TSV path (e.g., data_sv.tsv)")
    ap.add_argument("--ncbi-build", default="GRCh37", choices=["GRCh37","GRCh38"], help="Reference build label")
    ap.add_argument("--strong-only", action="store_true", help="Include only Archer 'Strong Evidence' fusions")
    ap.add_argument("--keep-intron", action="store_true", help="Include intronic events (default: excluded)")
    ap.add_argument("--keep-out-of-frame", action="store_true", help="Include out-of-frame events (default: excluded)")
    ap.add_argument("inputs", nargs="+", help="Archer *.full_results.txt file(s)")
    args = ap.parse_args()

    all_rows: List[Dict[str, str]] = []
    for path in args.inputs:
        kv = load_kv(Path(path))
        rows = archer_to_cbio_sv(
            kv,
            ncbi_build=args.ncbi_build,
            strong_only=args.strong_only,
            keep_intron=args.keep_intron,
            keep_out_of_frame=args.keep_out_of_frame,
        )
        all_rows.extend(rows)

    write_tsv(all_rows, Path(args.out))
    print(f"Wrote {len(all_rows)} rows to {args.out} (keep_intron={args.keep_intron}, keep_out_of_frame={args.keep_out_of_frame})")

if __name__ == "__main__":
    main()
