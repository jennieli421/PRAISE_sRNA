## mpileup → summary table for plotting
# seqID columns are a second pass: mpileup has no QNAME, so they cannot
# be counted in the '*' loop. Needs pysam (praise env).

import csv
import sys
from collections import defaultdict
from pathlib import Path
import pandas as pd
import pysam

if len(sys.argv) != 2:
    raise SystemExit(
        "usage: mpileup_summarize.py infile.txt\n"
        "run from the outdir that contains out.mpileup and realigned.bam"
    )
infile_txt = sys.argv[1]

def count_star(bases: str) -> int:
    """Ψ del：只数本位 '*'。跳过 ^Q、$；+/- 后的长度+序列只跳过、不计入 del。"""
    n = i = 0
    while i < len(bases):
        b = bases[i]
        if b == "^":
            i += 2
        elif b == "$":
            i += 1
        elif b in "+-":
            i += 1
            num = 0
            while i < len(bases) and bases[i].isdigit():
                num = num * 10 + int(bases[i])
                i += 1
            i += num
        elif b == "*":
            n += 1
            i += 1
        else:
            i += 1
    return n

inp, outp = "out.mpileup", "out.mod.tsv"
with open(inp) as f, open(outp, "w") as o:
    _ = o.write("chr_name\tchr_index\tref_base\tcoverage\tdel_count\tdelRate\n")
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 5: 
          print("c < 5") 
          break 
        cov = int(c[3])
        del_n = count_star(c[4])
        _ = o.write(f"{c[0]}\t{c[1]}\t{c[2]}\t{cov}\t{del_n}\t{del_n / cov:.6f}\n")
        
# -------------
df = pd.read_csv("out.mod.tsv", sep="\t")
print(df["ref_base"].value_counts())

df_T = df[(df["ref_base"] == "T") & (df["delRate"] >= 0.1) & (df["coverage"] >= 10)]
print(df_T.shape[0], "rows")

def fa_to_df(fa_path):
    rows, name, seq = [], None, []
    for line in Path(fa_path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                rows.append((name, "".join(seq).upper()))
            name, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if name is not None:
        rows.append((name, "".join(seq).upper()))
    return pd.DataFrame(rows, columns=["chr_name", "ref_fa"])


def add_del_seqid(plot_tsv, bam_path, infile_txt):
    """Map bowtie IDs to seqIDs and annotate deletion/non-deletion positions."""


    # 1. bowtieID -> seqID
    with open(infile_txt, newline="") as f:
        reader = csv.DictReader(f)
        bowtie_col = reader.fieldnames[0]
        bowtie_to_seq = {
            row[bowtie_col]: row["id"]
            for row in reader
        }

    deleted = defaultdict(set)
    kept = defaultdict(set)

    # CIGAR operations
    MATCH_OPS = {0, 7, 8}                # M, =, X
    # 记录到 kept
    # M: 比对位置匹配操作，可能有 match 或 mismatch
    # =: 与参考序列完全匹配
    # X: 与参考序列不完全匹配
    REF_CONSUMING_OPS = {0, 2, 3, 7, 8}  # M, D, N, =, X

    # 2. Collect seqIDs at each reference position
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.reference_name is None or not read.cigartuples:
                continue

            seq_id = bowtie_to_seq.get(read.query_name) # query_name is bowtieID
            if seq_id is None:
                continue

            if read.is_secondary:
                status = "second"
            elif read.is_supplementary:
                status = "supp"
            else:
                status = "primary"
            tag = f"{seq_id}:{status}"

            ref_name = read.reference_name
            ref_pos = read.reference_start + 1  # 1-based

            for op, length in read.cigartuples:
                if op == 2:  # D
                    store = deleted
                elif op in MATCH_OPS:  # M, =, X
                    store = kept
                else:
                    store = None

                if store is not None:
                    for pos in range(ref_pos, ref_pos + length):
                        store[(ref_name, pos)].add(tag)

                if op in REF_CONSUMING_OPS:
                    ref_pos += length

    # 3. Annotate the TSV
    df = pd.read_csv(plot_tsv, sep="\t")
    positions = list(zip(df["chr_name"], df["chr_index"].astype(int)))

    def format_tags(store, key):
        return ",".join(sorted(store.get(key, ())))

    df["del_seqid"] = [format_tags(deleted, key) for key in positions]
    df["non_del_seqid"] = [format_tags(kept, key) for key in positions]

    df.to_csv(plot_tsv, sep="\t", index=False)

    n_del = df["del_seqid"].ne("").sum()
    n_keep = df["non_del_seqid"].ne("").sum()
    print(
        f"{plot_tsv}: "
        f"del_seqid {n_del}/{len(df)}, "
        f"non_del_seqid {n_keep}/{len(df)}"
    )

ref = fa_to_df("/home/lixiang/data/PRAISE/sRNA_ref/all_rna.fa")
merged = df_T.merge(ref, on="chr_name", how="left")
merged.to_csv("plot.tsv", sep="\t", index=False)
add_del_seqid("plot.tsv", "realigned.bam", infile_txt)


# Sanity check ------------------------------------------------------------
base_at = merged.apply(
    lambda r: r["ref_fa"][int(r["chr_index"]) - 1]
    if isinstance(r["ref_fa"], str) and 1 <= int(r["chr_index"]) <= len(r["ref_fa"])
    else None,
    axis=1,
)
print("ref mismatch:", (merged["ref_base"] != base_at).sum())
print("ref NA (df有、FA无):", merged["ref_fa"].isna().sum())
