# Sports → PRAISE realignment（全 RNA）

srun: 
```bash
/usr/bin/srun --partition=local --gres=shard:1 --cpus-per-task=4 --mem=512G --time=24:00:00 --pty bash 

/usr/bin/squeue --user="$USER"

/usr/bin/scancel 586
```


ID 链：

```
seqname_map.txt : seqID → Illumina QNAME(s)
infile.txt      : bowtieID → seqID, reads
infile.fa       : bowtieID → seq
output_* SAM    : QNAME=bowtieID, RNAME=FA header
hsa/*/*.fa      : header == RNAME（暂假定全局不重复）
```

NH 不参与 realignment 计分；副本仅防 `get_tag("NH")` 崩溃。默认 `-rm off`，不跑 `remove_multi_mapping`。  
丰度：realign 后再按 `infile.txt` 的 reads 展开。


---

## `realignment_forward_sports.py` 相对原版的唯一改动

副本：`/home/lixiang/data/PRAISE/realignment_forward_sports.py`  
（相对原版仅：`NH` 缺省时写 `1`）

```python
# 四处 tags 写出：
("NH", read.get_tag("NH") if read.has_tag("NH") else 1)
```


---

## build sRNA reference 

```bash
cd /home/lixiang/data/PRAISE/sRNA_ref

# SAM stem ↔ FA（与 sports 命名一致；勿含 output_match_genome / genome.fa）
declare -a PAIRS=(
  "rRNA_high-qual:rRNA/high-qual.fa"
  "tRNA_hq-mature:tRNA/hq-mature.fa"
  "tRNA_hq-pre:tRNA/hq-pre.fa"
  "miRNA_hq-mature:miRNA/hq-mature.fa"
  "miRNA_hq-hairpin:miRNA/hq-hairpin.fa"
  "snRNA_high-qual:snRNA/high-qual.fa"
  "snoRNA_high-qual:snoRNA/high-qual.fa"
  "yRNA_high-qual:yRNA/high-qual.fa"
  "piRNA_high-qual:piRNA/high-qual.fa"
  "mRNA_high-qual:mRNA/high-qual.fa"
  "other_combine:other/combine.fa"
  "rRNA_low-qual:rRNA/low-qual.fa"
  "tRNA_low-qual:tRNA/low-qual.fa"
  "miRNA_low-qual:miRNA/low-qual.fa"
  "snRNA_low-qual:snRNA/low-qual.fa"
  "snoRNA_low-qual:snoRNA/low-qual.fa"
  "yRNA_low-qual:yRNA/low-qual.fa"
  "piRNA_low-qual:piRNA/low-qual.fa"
  "mRNA_low-qual:mRNA/low-qual.fa"
)

# ---- A) 合并 FA（一次性；已有则跳过）----
ref=/home/lixiang/data/PRAISE/sRNA_ref/all_rna.fa

mkdir -p "$(dirname "$ref")"
if [[ ! -s "$ref" ]]; then
  : > "$ref"
  for p in "${PAIRS[@]}"; do
    fa=${p#*:}
    cat "$DB/$fa" >> "$ref"
  done
fi

# 检查重复header
grep "^>" "$ref" | sed 's/^>//' | sort | uniq -d


# 去除含有X的序列
ref=/home/lixiang/data/PRAISE/sRNA_ref/all_rna.fa
ref_nx=/home/lixiang/data/PRAISE/sRNA_ref/all_rna.noX.fa
seqkit grep -v -s -i -p X "$ref" -o "$ref_nx"
seqkit grep -s -i -p X "$ref_nx" | grep -c "^>"   # 应为 0

```


## build SAM from SPORTS output (bothAD)
勿把 `output_match_genome` 与未含 `genome.fa` 的 RNA ref 混用。若只要 HQ，从 `PAIRS` 删掉 `*_low-qual` 即可。

### 1. setup 
```bash
# ---- Define sample ----
sample=control_1
root=/home/lixiang/data/20260930_293T_PRAISE/${sample}
bt=${root}/sports2_out_bothAD/.bowtie_result/${sample}_R1_bothAD
prep=${root}/sports2_out_bothAD/.preprocess/${sample}_R1_bothAD

# # ---- Define database ----
# PRAISE=/home/lixiang/data/PRAISE
# DB=/home/lixiang/data/sports_database/hsa

# # ---- A) merge all RNA SAM ----
# outdir=${root}/praise_from_sports/all_rna
# mkdir -p "$outdir"
# cd "$outdir"

# : > ${outdir}/all_rna.sam
# for p in "${PAIRS[@]}"; do
#   stem=${p%%:*}
#   cat "${bt}/output_${stem}" >> ${outdir}/all_rna.sam
# done

# # headerless SAM → name-sorted BAM
# ref_nx=/home/lixiang/data/PRAISE/sRNA_ref/all_rna.noX.fa
# samtools view -b -T "$ref_nx" ${outdir}/all_rna.sam \
#   | samtools sort -n -@ 8 -o ${outdir}/in.name.bam


# ---- B) merge tRNA SAM ----
outdir=${root}/praise_from_sports/tRNA_hq_mature
mkdir -p "$outdir"
cd "$outdir"

# mature tRNA's bowtie SAM
stem=tRNA_hq-mature
cp "${bt}/output_${stem}" ${outdir}/tRNA.sam

# headerless SAM → name-sorted BAM
ref_nx=${DB}/tRNA/hq-mature.fa
samtools view -b -T "$ref_nx" ${outdir}/tRNA.sam \
  | samtools sort -n -@ 8 -o ${outdir}/in.name.bam
```


### 2. realignment
```bash
# 2) realignment（-x 必须与 BAM 同一 FA；-rm 保持默认 off）
# 全 RNA 内存大：建议 -t 1~2；失败后用 --resume 续跑（跳过已有 .realign.bam）
cd "$PRAISE"
python realignment_forward_sports.py --fast -t 4 -ms 4.8 \
  -x "$ref_nx" \
  -i ${outdir}/in.name.bam \
  -o ${outdir}/realigned.bam \
  -f ${outdir}/filtered.bam 
  # --keep-tmp

# 若中断 / 缺 shard，对已有 tmp 目录续跑（只重做缺 .realign 的片）：
python realignment_forward_sports.py --fast -t 1 -ms 4.8 \
  -x $ref_nx \
  -i ${outdir}/in.name.bam \
  -o ${outdir}/realigned.bam \
  -f ${outdir}/filtered.bam \
  --resume /home/lixiang/data/PRAISE/tmp_AmCih2My_in.name.bam \
  --keep-tmp
```

### 3. expand by abundance 
```bash
# 3) 按 infile.txt 展开丰度 → coord-sort
python3 - "${outdir}" "${prep}/infile.txt" <<'PY'
import csv, sys, pysam
from pathlib import Path

outdir, counts_path = Path(sys.argv[1]), Path(sys.argv[2])
bam_in, bam_out = outdir / "realigned.bam", outdir / "expanded.bam"

counts = {}
with open(counts_path) as f:
    r = csv.DictReader(f)
    hcol = r.fieldnames[0]
    for row in r:
        counts[row[hcol]] = int(row["reads"])

inp = pysam.AlignmentFile(bam_in, "rb")
out = pysam.AlignmentFile(bam_out, "wb", header=inp.header)
n_skip = 0
for read in inp.fetch(until_eof=True):
    n = counts.get(read.query_name)
    if n is None:
        n_skip += 1
        continue
    for i in range(n):
        r = pysam.AlignedSegment(out.header)
        r.query_name = f"{read.query_name}_{i}"
        r.flag = read.flag
        r.reference_id = read.reference_id
        r.reference_start = read.reference_start
        r.mapping_quality = read.mapping_quality
        r.cigar = read.cigar
        r.next_reference_id = read.next_reference_id
        r.next_reference_start = read.next_reference_start
        r.template_length = read.template_length
        r.query_sequence = read.query_sequence
        r.query_qualities = read.query_qualities
        r.set_tags(read.get_tags())
        out.write(r)
inp.close(); out.close()
print(f"done expand; skipped_no_count={n_skip}")
PY

samtools sort -@ 8 -o ${outdir}/expanded.sorted.bam ${outdir}/expanded.bam
samtools index ${outdir}/expanded.sorted.bam
```

### 4. mpileup
```bash
samtools mpileup -BQ0 -f "$ref_nx" \
  ${outdir}/expanded.sorted.bam \
  -o ${outdir}/out.mpileup --ff UNMAP,QCFAIL 

# 参数-aa不需要；-d 不需要
# -B 禁用 BAQ（Base Alignment Quality）调整碱基质量值
# -Q0	设置最低碱基质量（base quality）为 0，即不因碱基质量低于默认阈值 13 而过滤
# --ff UNMAP,QCFAIL	排除带有 UNMAP 或 QCFAIL 标志的 reads。默认本来会丢掉 SECONDARY，加上该参数后 secondary 会进入 pileup，本位的 * 被算进 del_count
```

---
## mpileup → 统计表

mpileup 格式
`coverage` = `$4`；`del_count` = 本位 `*` 个数（Ψ：该位点被删）；**不算** substitution；**不把** 前一位的 `-` / `+` 并入 del。  
`delRate = del_count/coverage`。FA join：`/home/lixiang/data/PRAISE/sRNA_ref/all_rna.fa`。

`del_seqid` / `non_del_seqid`：该位点上有 deletion、以及盖住该位点但没有 deletion 的 seqID。每条写成 `seqID:primary`、`seqID:secondary` 或 `seqID:unknown`（supplementary），逗号分隔。来自 realigned.bam 的 bowtieID，经 infile.txt 对应到 seqID。图上不画 seqID。


```bash
# 在 outdir 下运行。praise 环境提供 pysam。
# 写出 plot.tsv 后另读 realigned.bam，填 del_seqid 和 non_del_seqid。
python3 ~/data/20260930_293T_PRAISE/mpileup_summarize.py ${prep}/infile.txt
```


## R plot

```bash
mm activate r-plot 
Rscript plot_praise_mod.R
```

与modomics映射
```bash
python /home/lixiang/data/20260930_293T_PRAISE/join_modomics_trna.py
```




---

## BMAT (optional)
```bash
PARSE=/home/lixiang/data/PRAISE/parse-mpileup.py
python "$PARSE" -p 8 \
  -i ${outdir}/out.mpileup \
  -o ${outdir}/out.bmat
```