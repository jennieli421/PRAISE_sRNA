import pandas as pd
from pathlib import Path

fasta_path = Path("/home/lixiang/data/database/modomics/modified_tRNA_unicode.fasta")
csv_path = Path("/home/lixiang/data/database/modomics/natural_modifications.csv")
out_path = Path("/home/lixiang/data/database/modomics/fasta_symbol_map.tsv")

# Characters from sequence bodies and from Feature (anticodon) in the header.
chars = []
seen = set()
for line in fasta_path.read_text().splitlines():
    if not line:
        continue
    if line.startswith(">"):
        # Handel anticodon in the header.
        header = dict(part.split(":", 1) for part in line[1:].split("|") if ":" in part)
        text = header.get("Feature", "")
    else:
        # Handle sequence bodies.
        text = line.strip()
    for ch in text:
        if ch not in seen:
            seen.add(ch)
            chars.append(ch)

mods = pd.read_csv(csv_path, dtype=str).fillna("")
keep = [
    "ID", "name", "short name", "abbrev", "new abbrev",
    "reference moiety", "moiety type",
]

rows = []
for ch in chars:
    hit = mods[(mods["abbrev"] == ch) | (mods["new abbrev"] == ch)]
    for rec in hit[keep].to_dict("records"):
        rec["name"] = " ".join(rec["name"].split())
        rows.append({
            "fasta_char": ch,
            "codepoint": f"U+{ord(ch):04X}",
            **rec,
        })

table = pd.DataFrame(rows)
unmatched = [ch for ch in chars if ch not in set(table["fasta_char"])]
table.to_csv(out_path, sep="\t", index=False)

print(f"characters {len(chars)}")
print(f"rows {len(table)}")
print(f"unmatched {unmatched}")
print(out_path)
