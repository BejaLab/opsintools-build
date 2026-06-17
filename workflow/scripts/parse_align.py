import csv
import json
from pathlib import Path

from Bio import AlignIO

GAP_CHARS = {"-", "."}

def load_ref(path):
    with Path(path).open("r", encoding="utf-8") as handle:
        ref = json.load(handle)
    if "id" not in ref or "seq" not in ref or "tms" not in ref:
        raise ValueError(f"Reference file {path} is missing required keys (id, seq, tms)")
    return ref

def parse_alignment(path):
    aln = AlignIO.read(str(path), "clustal")
    if len(aln) < 2:
        raise ValueError(f"Alignment {path} must contain at least two sequences")
    return aln

def normalize_seq(seq):
    return "".join(seq.split()).upper()


def ungap(seq):
    return "".join(ch for ch in seq if ch not in GAP_CHARS)


def unique_substring_start(full_seq, trimmed_seq, ref_id):
    if not trimmed_seq:
        raise ValueError(f"Aligned sequence for {ref_id} has no residues after removing gaps")
    first = full_seq.find(trimmed_seq)
    if first == -1:
        raise ValueError(f"Could not map trimmed aligned sequence to full sequence for {ref_id}")
    if full_seq.find(trimmed_seq, first + 1) != -1:
        raise ValueError(f"Trimmed aligned sequence maps ambiguously in full sequence for {ref_id}")
    return first


def match_record_for_ref(
    ref, records, used_record_ids
):
    ref_id = str(ref["id"])
    full_seq = normalize_seq(str(ref["seq"]))

    candidates = []
    for rec in records:
        if rec.id in used_record_ids:
            continue
        aligned_seq = normalize_seq(str(rec.seq))
        trimmed = ungap(aligned_seq)
        try:
            start = unique_substring_start(full_seq, trimmed, ref_id)
        except ValueError:
            continue

        desc = f"{rec.id} {rec.description}".upper()
        name_match = ref_id.upper() in desc
        candidates.append((rec, aligned_seq, start, name_match))

    if not candidates:
        raise ValueError(f"No aligned sequence could be matched to reference {ref_id}")

    named = [c for c in candidates if c[3]]
    if len(named) == 1:
        rec, aligned_seq, start, _ = named[0]
        return rec, aligned_seq, start
    if len(named) > 1:
        raise ValueError(f"Multiple alignment records match reference name {ref_id}")

    if len(candidates) == 1:
        rec, aligned_seq, start, _ = candidates[0]
        return rec, aligned_seq, start
    raise ValueError(f"Multiple alignment records could map to reference {ref_id}; cannot disambiguate")


def tm_label(tms, pos_1based):
    if pos_1based is None:
        return ""
    for helix, interval in tms.items():
        if not isinstance(interval, list) or len(interval) != 2:
            raise ValueError(f"Invalid TM interval for {helix}: {interval}")
        start, end = int(interval[0]), int(interval[1])
        if start <= pos_1based <= end:
            return str(helix)
    return ""


aln_dir = Path(snakemake.input.aln)
aln_path = aln_dir / "t_coffee.aln"
if not aln_path.exists():
    raise FileNotFoundError(f"Alignment file not found: {aln_path}")

ref_paths = [Path(p) for p in snakemake.input.refs]
if len(ref_paths) != 2:
    raise ValueError(f"Expected exactly 2 reference files, got {len(ref_paths)}")
ref1, ref2 = (load_ref(ref_paths[0]), load_ref(ref_paths[1]))

alignment = parse_alignment(aln_path)
records = list(alignment)

used = set()
rec1, aln_seq1, start1 = match_record_for_ref(ref1, records, used)
used.add(rec1.id)
rec2, aln_seq2, start2 = match_record_for_ref(ref2, records, used)

if len(aln_seq1) != len(aln_seq2):
    raise ValueError(
        f"Matched alignment records have different lengths: {rec1.id}={len(aln_seq1)}, {rec2.id}={len(aln_seq2)}"
    )

out_path = Path(snakemake.output[0])

ref1_name = str(ref1["id"])
ref2_name = str(ref2["id"])
pos_counter1 = 0
pos_counter2 = 0

with out_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(
        [
            "ref1",
            "ref2",
            "ref1_pos",
            "ref2_pos",
            "ref1_res",
            "ref2_res",
            "ref1_transmembrane_helix",
            "ref2_transmembrane_helix",
        ]
    )

    for c1, c2 in zip(aln_seq1, aln_seq2):
        is_gap1 = c1 in GAP_CHARS
        is_gap2 = c2 in GAP_CHARS

        if is_gap1 and is_gap2:
            continue

        ref1_pos = None
        ref2_pos = None
        ref1_res = ""
        ref2_res = ""

        if not is_gap1:
            pos_counter1 += 1
            ref1_pos = start1 + pos_counter1
            ref1_res = c1

        if not is_gap2:
            pos_counter2 += 1
            ref2_pos = start2 + pos_counter2
            ref2_res = c2

        writer.writerow(
            [
                ref1_name,
                ref2_name,
                "" if ref1_pos is None else ref1_pos,
                "" if ref2_pos is None else ref2_pos,
                ref1_res,
                ref2_res,
                tm_label(ref1["tms"], ref1_pos),
                tm_label(ref2["tms"], ref2_pos),
            ]
        )
