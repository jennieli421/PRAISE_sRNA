"""Map human MODOMICS tRNAs onto hq-mature.fa by exact parent-DNA sequence.

Isotype + anticodon only narrow candidates. Match requires identical parent DNA
(U→T). No CCA trim, reverse complement, or length edit.
"""

import pandas as pd
from pathlib import Path

MODOMICS_FASTA = Path("/home/lixiang/data/database/modomics/modified_tRNA_unicode.fasta")
SYMBOL_MAP = Path("/home/lixiang/data/database/modomics/fasta_symbol_map.tsv")
REFERENCE_FASTA = Path("/home/lixiang/data/sports_database/hsa/tRNA/hq-mature.fa")
OUT_TSV = Path("/home/lixiang/data/database/modomics/join_modomics_tRNA-hq-mature.tsv")
QC_TSV = Path("/home/lixiang/data/database/modomics/join_modomics_tRNA-hq-mature.qc.tsv")

PARENT = {
    "Adenosine": "A",
    "Cytidine": "C",
    "Guanosine": "G",
    "Uridine": "U",
}
ISOTYPE_ALIAS = {"Ini": "iMet", "Sec": "SeC"}


def to_dna(text, char2parent):
    """Mod / RNA letters -> DNA. Unknown or multi-parent symbols become N."""
    dna = []
    unknown = []
    for ch in text:
        if ch in "ACGT":
            dna.append(ch)
        elif ch in "ACGU":
            dna.append("T" if ch == "U" else ch)
        else:
            # Convert modification characters to parent bases.
            parent = char2parent.get(ch)
            if parent is None:
                unknown.append(ch)
                dna.append("N")
            else:
                dna.append("T" if parent == "U" else parent)
    return "".join(dna), "".join(dict.fromkeys(unknown))


def main():
    # symbol map covers sequence body and Feature (anticodon) characters.
    symbol_df = pd.read_csv(SYMBOL_MAP, sep="\t", dtype=str)
    nucleoside = symbol_df[symbol_df["moiety type"] == "nucleoside"]
    char2parent = {}
    conflicts = []
    for row in nucleoside.to_dict("records"):
        ch = row["fasta_char"]
        parent = PARENT.get(row["reference moiety"])
        if ch in char2parent and char2parent[ch] != parent:
            conflicts.append(f"{ch!r}: {char2parent[ch]!r} vs {parent!r}")
        char2parent[ch] = parent
    if conflicts:
        raise ValueError("conflicting symbol parents:\n  " + "\n  ".join(conflicts))

    # MODOMICS records (Feature and sequence share char2parent).
    modomics = []
    header = None
    for line in MODOMICS_FASTA.read_text().splitlines():
        if line.startswith(">"):
            header = dict(p.split(":", 1) for p in line[1:].split("|") if ":" in p)
            continue
        if not line.strip() or header is None:
            continue
        seq = line.strip()
        feature = header.get("Feature", "")
        subtype = header.get("Subtype", "")
        anticodon, anticodon_unknown = to_dna(feature, char2parent)
        parent_dna, sequence_unknown = to_dna(seq, char2parent)
        modomics.append({
            "modomics_id": header.get("id", ""),
            "subtype": subtype,
            "isotype": ISOTYPE_ALIAS.get(subtype, subtype),
            "feature": feature,
            "anticodon": anticodon,
            "modomics_sequence": seq,
            "parent_dna": parent_dna,
            "anticodon_unknown": anticodon_unknown,
            "sequence_unknown": sequence_unknown,
        })

    # hq-mature reference.
    reference = []
    header = None
    for line in REFERENCE_FASTA.read_text().splitlines():
        if line.startswith(">"):
            header = line[1:].strip()
            continue
        if not line.strip() or header is None:
            continue
        fields = header.split(",")
        copy_id = fields[3]
        body = copy_id.split("tRNA-", 1)[1]
        if body.startswith("mt_"):
            body = body[3:]
        isotype, anticodon, _ = body.split("-", 2)
        reference.append({
            "copy_id": copy_id,
            "ref_isotype": isotype,
            "ref_anticodon": anticodon,
            "ref_compartment": fields[2],
            "hq_mature_sequence": line.strip(),
        })
        header = None

    # Ambiguous MODOMICS -> QC; the rest are indexed for matching.
    qc_rows = []
    by_trna = {}
    for rec in modomics:
        if rec["sequence_unknown"]:
            qc_rows.append({**rec, "unknown_symbols": rec["sequence_unknown"],
                            "status": "unmatched_ambiguous_sequence"})
            continue
        if rec["anticodon_unknown"]:
            qc_rows.append({**rec, "unknown_symbols": rec["anticodon_unknown"],
                            "status": "unmatched_ambiguous_anticodon"})
            continue
        by_trna.setdefault((rec["isotype"], rec["anticodon"]), []).append(rec)

    empty = {
        "modomics_id": "", "subtype": "", "isotype": "",
        "feature": "", "anticodon": "", "modomics_sequence": "",
    }
    rows = []
    for ref in reference:
        key = (ref["ref_isotype"], ref["ref_anticodon"])
        hits = [m for m in by_trna.get(key, []) if m["parent_dna"] == ref["hq_mature_sequence"]]
        if not hits:
            rows.append({**ref, **empty})
            continue
        for m in hits:
            rows.append({
                **ref,
                "modomics_id": m["modomics_id"],
                "subtype": m["subtype"],
                "isotype": m["isotype"],
                "feature": m["feature"],
                "anticodon": m["anticodon"],
                "modomics_sequence": m["modomics_sequence"],
            })

    join = pd.DataFrame(rows)[[
        "copy_id", "ref_isotype", "ref_anticodon", "ref_compartment",
        "hq_mature_sequence",
        "modomics_id", "subtype", "isotype", "feature", "anticodon",
        "modomics_sequence",
    ]]
    qc = pd.DataFrame(qc_rows)
    if len(qc):
        qc = qc[[
            "modomics_id", "subtype", "isotype", "feature", "anticodon",
            "modomics_sequence", "parent_dna", "unknown_symbols", "status",
        ]]
    else:
        qc = pd.DataFrame(columns=[
            "modomics_id", "subtype", "isotype", "feature", "anticodon",
            "modomics_sequence", "parent_dna", "unknown_symbols", "status",
        ])

    join.to_csv(OUT_TSV, sep="\t", index=False)
    qc.to_csv(QC_TSV, sep="\t", index=False)

    matched = join[join["modomics_id"] != ""]
    print(f"join rows {len(join)}; copy_id {join['copy_id'].nunique()}")
    print(
        f"matched rows {len(matched)} "
        f"(copy_id {matched['copy_id'].nunique()}, "
        f"modomics_id {matched['modomics_id'].nunique()})"
    )
    print(f"unmatched reference {(join['modomics_id'] == '').sum()}; QC {len(qc)}")
    if len(qc):
        print(qc["status"].value_counts().to_string())
    print(OUT_TSV)
    print(QC_TSV)


if __name__ == "__main__":
    main()
