# PRAISE 分析笔记（293T）20260930 

工作目录：`/home/lixiang/data/20260811_293T_PRAISE`  
环境：`export PATH=/home/lixiang/envs/praise/bin:/home/lixiang/software/bin:$PATH`  
仓库脚本：`PRAISE/` 

文库结构

```
5AD—NNNNN—Insert—NNNNN—3AD—index—P7
```

读段方向（PE150）：

| Read | 5′→3′ 读段内容 |
|------|----------------|
| **R1** | Insert — NNNNN(3′UMI) — 3AD — index — P7 |
| **R2** | Insert — NNNNN(5′UMI) — 5AD反向互补 — P5反向互补 |

接头序列（建库正向定义，供 cutadapt `-a`/`-A`）：

- **3AD**：`AGATCGGAAGAGCACACGTCT`
- **5AD** 恒定区：`GTTCAGAGTTCTACAGTCCGACGATC`（R2 读段中为反向互补 `GATCGTCGGACTGTAGAACTCTGAAC`）
- Insert **15–200 nt** 
- 两侧各 **5N UMI**（共 10 nt）



### 分析策略：仅保留「双端接头完整」read pairs

因建库经特殊处理，**只用 R1、R2 均切到预期接头的 read pair** 进入下游（非 `psu.md` 默认的「仅 R1 切 3AD」）。

--------

srun: 
```bash
/usr/bin/srun --partition=local --gres=shard:1 --cpus-per-task=8 --mem=512G --time=24:00:00 --pty bash 
/usr/bin/squeue --user="$USER"
/usr/bin/scancel 586

# 修复 mamba 环境调用
bash --noprofile --norc
```


## 数据预处理
在堡垒机上运行


### raw data 

```bash
file                       format  type    num_seqs        sum_len  min_len  avg_len  max_len
control_1.R1.raw.fastq.gz  FASTQ   DNA   14,602,198  2,204,931,898      151      151      151
control_1.R2.raw.fastq.gz  FASTQ   DNA   14,602,198  2,204,931,898      151      151      151
control_2.R1.raw.fastq.gz  FASTQ   DNA   16,021,812  2,419,293,612      151      151      151
control_2.R2.raw.fastq.gz  FASTQ   DNA   16,021,812  2,419,293,612      151      151      151
praise_1.R1.raw.fastq.gz   FASTQ   DNA    9,707,609  1,465,848,959      151      151      151
praise_1.R2.raw.fastq.gz   FASTQ   DNA    9,707,609  1,465,848,959      151      151      151
praise_2.R1.raw.fastq.gz   FASTQ   DNA   10,443,376  1,576,949,776      151      151      151
praise_2.R2.raw.fastq.gz   FASTQ   DNA   10,443,376  1,576,949,776      151      151      151
```

### cutadapt

第一步：质量修剪
-q 20：常规质量修剪。从reads的3'端开始，切除质量值（Phred score）低于20的碱基，直到遇到质量≥20的碱基为止。
--nextseq-trim=20：NextSeq/NovaSeq专用质量修剪。这两个平台采用双色化学，当测序反应“脱落”时会产生长串高质量但错误的G碱基（poly-G尾巴）。此参数会忽略G碱基的质量分数，将其视为低质量并切除。阈值20代表切除质量低于20的区域。


第二步：接头匹配与修剪
-O（最小重叠长度）指接头序列与read必须至少有10个碱基完全匹配（在计算错配之前），Cutadapt才会启动修剪。
-e 允许接头匹配的错配率。计算公式为：允许错配数（2.1） = e（0.1）× 接头长度（21 or 26）


第三步：长度筛选
--minimum-length 25:25：修完后，如果R1的长度小于25 或 R2的长度小于25，则该read对被标记为“太短”。格式 25:25 分别对应R1和R2的阈值。
--maximum-length 210:210：同理，如果R1长度大于210 或 R2长度大于210，则该read对被标记为“太长”。这通常用于排除接头二聚体或未正确断裂的嵌合体。

第四步：配对保留
--discard-untrimmed：如果某条read未被修剪，则丢弃该read。对于双端数据，强制双端同时有接头才保留。
--pair-filter=both：如果R1和R2都不合格才丢弃

参考：[FAQ: How should my NEBNext Small RNA library be trimmed?](https://www.neb.com/en-us/faqs/how-should-my-nebnext-small-rna-library-be-trimmed)


```bash
sample=control_2
mkdir -p $sample/cutadapt

# R1切3AD 
cutadapt -j 32 \
  -a 3AD=ATCACGAGATCGGAAGAGCACACGTCT \
  -a 3AD=ATCACGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -O 5 -e 0.1 \
  -q 20 --nextseq-trim=20 \
  -m 10 \
  --info-file ${sample}/cutadapt/info_R1_cut3AD.tsv \
  -o ${sample}/cutadapt/${sample}_R1.fastq.gz \
  rawdata/${sample}.R1.raw.fastq.gz > ${sample}/cutadapt/log-3AD.log 2>&1


# R2切5AD
cutadapt -j 32 \
  -a 5AD=GATCGTCGGACTGTAGAACTCTGAAC \
  -O 5 -e 0.1 \
  -q 20 --nextseq-trim=20 \
  -m 10 \
  --info-file ${sample}/cutadapt/info_R2_cut5AD.tsv \
  -o ${sample}/cutadapt/${sample}_R2.fastq.gz \
  rawdata/${sample}.R2.raw.fastq.gz > ${sample}/cutadapt/log-5AD.log 2>&1

cd $sample/cutadapt 
运行 process_cutadapt_info_file.py 

```




### 修复配对
```bash
mkdir ../categorize && cd $_
# umi extract 要去严格配对关系
repair.sh \
    in1=../cutadapt/${sample}_R1.fastq.gz \
    in2=../cutadapt/${sample}_R2.fastq.gz \
    out1=${sample}_R1_paired.fastq.gz \
    out2=${sample}_R2_paired.fastq.gz \
    outsingle=${sample}_unpaired.fastq.gz

```


按 ID 提取 bothAD：

```bash
# 安装 BBMap
# mm activate cutadapt
# mm install -c bioconda bbmap
# include=t：提取名单中的 reads（f 则排除）

cd ${sample}/categorize

in1=${sample}_R1_paired.fastq.gz
in2=${sample}_R2_paired.fastq.gz

filter=../cutadapt/valid_ids_bothAD.txt
out1=${sample}_R1_bothAD.fq.gz
out2=${sample}_R2_bothAD.fq.gz

seqkit grep -f "$filter" "$in1" -o "$out1"
seqkit grep -f "$filter" "$in2" -o "$out2"
```


### cutadapt round2 
接头若缺了自己的 5′ 端，文档写明默认找不到，例如 APTERmysequence。要这种匹配，需要 -a "ADAPTER;anywhere"。

错误 -a 的逻辑：找到 adapter 后，把 adapter 本身以及 adapter 后面的所有序列全部剪掉。
<!-- ```bash
# R1继续切
cd cutadapt-round2
cutadapt -j 32 \
  -a 5AD=GTTCTACAGTCCGACGATC \
  -O 5 -e 0.1 \
  -q 20 --nextseq-trim=20 \
  -m 10 \
  -o ${sample}_R1_bothAD_cut5AD.fastq.gz \
  ../categorize/${sample}_R1_bothAD.fq.gz > log-R1-bothAD-5AD.log 2>&1
GTTCAGAGTTCTACAGTCCGACGATC

cutadapt -j 32 \
  -a 5AD=GTTCTACAGTCCGACGATC \
  -O 5 -e 0.1 \
  -q 20 --nextseq-trim=20 \
  -m 10 \
  -o ${sample}_R1_bothAD_cut5AD_cut5AD.fastq.gz \
  ${sample}_R1_bothAD_cut5AD.fastq.gz > log-R1-bothAD-5AD-5AD.log 2>&1
``` -->

try 2  正确
-g 逻辑：找到 adapter 后，删除 adapter 以及 adapter 前面的序列，保留 adapter 后面的序列。
```bash
cutadapt -j 32 \
  -g "5AD=GTTCTACAGTCCGACGATC" \
  -O 5 -e 0.1 \
  -q 20 --nextseq-trim=20 \
  -m 10 \
  -o ${sample}_R1_bothAD_cut5AD_g.fastq.gz \
  ../categorize/${sample}_R1_bothAD.fq.gz > log-R1-bothAD-5AD-g.log 2>&1


### compare two ways of cutadapt 
cd /home/lixiang/data/20260930_293T_PRAISE/control_1/cutadapt-round2
A=control_1_R1_bothAD_cut5AD_cut5AD.fastq.gz
B=control_1_R1_bothAD_cut5AD_g.fastq.gz
RAW=/home/lixiang/data/20260930_293T_PRAISE/rawdata/control_1.R1.raw.fastq.gz
OUT=compare_cut5AD_cut5AD_vs_g.tsv

# -i: 只要 ID（空格前）；默认输出 id \t seq \t qual
seqkit fx2tab -i "$A" | cut -f1,2 > a.tsv
seqkit fx2tab -i "$B" | cut -f1,2 > b.tsv
seqkit fx2tab -i "$RAW" | cut -f1,2 \
| gawk -F'\t' -v OFS='\t' '
    ARGIND==1 { a[$1]=$2; next }
    ARGIND==2 { b[$1]=$2; next }
    ARGIND==3 {
      if (($1 in a) || ($1 in b)) raw[$1]=$2
    }
    END {
      print "read_id","in_A","in_B","in_raw","len_A","len_B","len_raw","same_AB","seq_A","seq_B","seq_raw"
      for (id in a) ids[id]
      for (id in b) ids[id]
      for (id in ids) {
        sa=(id in a)?a[id]:""; sb=(id in b)?b[id]:""; sr=(id in raw)?raw[id]:""
        same=(sa!="" && sb!="" && sa==sb)?1:0
        print id, (sa!=""), (sb!=""), (sr!=""), length(sa), length(sb), length(sr), same, sa, sb, sr
      }
    }
  ' a.tsv b.tsv - > compare_cut5AD_cut5AD_vs_g.tsv

# 
cd /home/lixiang/data/20260930_293T_PRAISE/control_1/cutadapt-round2
IN=compare_cut5AD_cut5AD_vs_g.tsv
OUT=compare_split
mkdir -p "$OUT"

awk -F'\t' -v d="$OUT" '
NR==1 {
  print > d"/only_A.tsv"
  print > d"/only_B.tsv"
  print > d"/both_same.tsv"
  print > d"/both_diff.tsv"
  next
}
# 1) A=0, B≠0
$2==0 && $3!=0 {
  n_onlyB++
  print > d"/only_B.tsv"
  print $1 > d"/only_B.ids"
  next
}
# 2) A≠0, B=0
$2!=0 && $3==0 {
  n_onlyA++
  print > d"/only_A.tsv"
  print $1 > d"/only_A.ids"
  next
}
# 3) 两边都有且序列相同
$2!=0 && $3!=0 && $8==1 {
  n_same++
  print > d"/both_same.tsv"
  print $1 > d"/both_same.ids"
  next
}
# 4) 两边都有但序列不同
$2!=0 && $3!=0 && $8==0 {
  n_diff++
  print > d"/both_diff.tsv"
  print $1 > d"/both_diff.ids"
  next
}
{ n_other++ }
END {
  tot = n_onlyA + n_onlyB + n_same + n_diff + n_other
  printf "category\tn\tpct\n"
  printf "only_A\t%d\t%.4f\n", n_onlyA+0, 100*(n_onlyA+0)/tot
  printf "only_B\t%d\t%.4f\n", n_onlyB+0, 100*(n_onlyB+0)/tot
  printf "both_same\t%d\t%.4f\n", n_same+0, 100*(n_same+0)/tot
  printf "both_diff\t%d\t%.4f\n", n_diff+0, 100*(n_diff+0)/tot
  printf "other\t%d\t%.4f\n", n_other+0, 100*(n_other+0)/tot
  printf "total\t%d\t100\n", tot
}
' "$IN" | tee "$OUT/counts.tsv"
```




UMI extract (cutadapt-round2)
``` bash
# control_1
umi_tools extract \
  --extract-method=regex \
  --bc-pattern="^(?P<umi_1>.{5})(?P<seq>.+)(?P<umi_2>.{5})$" \
  --stdin=cutadapt-round2/${sample}_R1_bothAD_cut5AD_g.fastq.gz \
  --stdout=cut_umi/${sample}_R1_bothAD_cut5AD_g.fq.gz \
  --log=cut_umi/log/umi-from-cutadapt-g.log

# sports only takes unzipped file 
cd cut_umi
gunzip *.gz


mm activate sports 

# umi-from-cut-adapt 单端比对
sample=control_1
output=/home/lixiang/data/20260930_293T_PRAISE/${sample}/sports2_out_cutumi_from_cut5ad_g
fq1=/home/lixiang/data/20260930_293T_PRAISE/${sample}/cut_umi/${sample}_R1_bothAD_cut5AD_g.fq # control_1

sports -i $fq1 \
       -ry ~/data/sports_database/hsa/ref.yaml \
       -l 10 -L 200 \
       -V -z -d -p 32 -o $output \
       --minscore L,-0.8,-0.8 \
       --seedmismatch 1  \
       --seedlen 8 


```


#### umi extract

```env setup 
mm create -n praise -c conda-forge -c bioconda python=3.10 cutadapt=4.8 seqkit=0.13.2 umi_tools=1.1.6 -y
mm activate praise 
```

``` bash
mkdir -p cut_umi/log/

# control_1
# bothAD paired 提取
umi_tools extract \
  --extract-method=regex \
  --bc-pattern="^(?P<umi_1>.{5})(?P<seq>.+)(?P<umi_2>.{5})$" \
  --stdin=categorize/${sample}_R1_bothAD.fq.gz \
  --stdout=cut_umi/${sample}_R1_bothAD.fq.gz \
  --log=cut_umi/log/umi-from-bothAD-paired.log


# praise_1 
# bothAD paired：R1 两端各提取 5nt UMI
umi_tools extract \
  --extract-method=regex \
  --bc-pattern="^(?P<umi_1>.{5})(?P<seq>.+)(?P<umi_2>.{5})$" \
  --stdin=categorize/${sample}_R1_bothAD.fq.gz \
  --stdout=cut_umi/${sample}_R1_bothAD.fq.gz \
  --log=cut_umi/log/umi.log

# sports only takes unzipped file 
cd cut_umi
gunzip *.gz
```

<!-- 
# 同时cut ad umi
```bash
mkdir cutad_umi
cutadapt \
    -j 0 \
    -a 3ad_umi=NNNNNATCACGAGATCGGAAGAGCACACGTCT \
    -u 5 \
    -q 20 \
    --nextseq-trim=20 \
    -m 15 \
    -M 200 \
    --info-file ${sample}/cutad_umi/info_3AD.tsv \
    -o ${sample}/cutad_umi/${sample}_R1.fastq.gz \
    rawdata/${sample}.R1.raw.fastq.gz > ${sample}/cutad_umi/log-3ad-umi.log 2>&1
``` -->


## sports 比对

bothAD noAD 分开比对
```bash
mm activate sports 

# umi-from-cut-adapt 单端比对
sample=control_1
output=/home/lixiang/data/20260930_293T_PRAISE/${sample}/sports2_out_cutumi_from_cutadapt
fq1=/home/lixiang/data/20260930_293T_PRAISE/${sample}/cut_umi/${sample}_R1.fastq # control_1

# umi-from-bothAD-paired.log 单端比对
sample=control_1
output=/home/lixiang/data/20260930_293T_PRAISE/${sample}/sports2_out_bothAD
fq1=/home/lixiang/data/20260930_293T_PRAISE/${sample}/cut_umi/${sample}_R1_bothAD.fq # control_1


sample=praise_1
output=/home/lixiang/data/20260930_293T_PRAISE/${sample}/sports2_out_bothAD
fq1=/home/lixiang/data/20260930_293T_PRAISE/${sample}/cut_umi/${sample}_R1_bothAD.fq # praise_1

sports -i $fq1 \
       -ry ~/data/sports_database/hsa/ref.yaml \
       -l 10 -L 200 \
       -V -z -d -p 32 -o $output \
       --minscore L,-0.8,-0.8 \
       --seedmismatch 1  \
       --seedlen 8 

```


```bash
# --- SE（bothAD 等）---
cd control_1/sports2_out_bothAD 

awk -F'\t' 'NR>1 { a[$6]+=$3; t+=$3 }
END {
  for (c in a) {
    n=c; if(n=="NO_Annotation") n="unannotated"
    printf "%s\t%.0f\t%.2f\n", n, a[c], a[c]*100/t
  }
}' */output.tsv | sort -k2,2nr > annotation_summary_all.tsv


# --- PE（noAD 等）---
cd /home/lixiang/data/20260830_293T_15_200nt_PRAISE/control_2/sports2_out_noAD

awk -F'\t' 'NR>1 {
  a[$7]+=$4; t+=$4
}
END {
  for (c in a) {
    n=c; if(n=="NO_Annotation") n="unannotated"
    printf "%s\t%.0f\t%.2f\n", n, a[c], a[c]*100/t
  }
}' control_2_noAD/pe/output.tsv | sort -k2,2nr > annotation_summary_all_pe.tsv
```



## sports 比对结果可视化

```bash
mm activate my_r_env 

sample=praise_1
sample=praise_2
sample=control_1

sample=control_2

Rscript plot_sports_pe_no_unannotated.R \
  ${sample}/sports2_out_bothAD/${sample}_R1_bothAD/output.tsv \
  ${sample}/sports2_out_bothAD/plots_annotated_only \
  ${sample}

Rscript plot_sports_pe_no_unannotated.R \
  ${sample}/sports2_out_combined_2/combined_R1/pe/output.tsv \
  ${sample}/sports2_out_combined_2/plots_annotated_only_pe



samples=("praise_1" "praise_2" "control_1" "control_2")

for sample in "${samples[@]}"; do
    # 第一个 Rscript：noAD 目录，没有第三个参数
    Rscript plot_sports_pe_no_unannotated.R \
        "${sample}/sports2_out_noAD/${sample}_R1_noAD/pe/output.tsv" \
        "${sample}/sports2_out_noAD/plots_annotated_only"

    # 第二个 Rscript：bothAD 目录，有第三个参数（样本名 + _R1_bothAD）
    Rscript plot_sports_pe_no_unannotated.R \
        "${sample}/sports2_out_bothAD/${sample}_R1_bothAD/output.tsv" \
        "${sample}/sports2_out_bothAD/plots_annotated_only" \
        "${sample}_R1_bothAD"
done

```


我现在的上游处理流程与hisat2不同。我用sports比对。我们的比对基于小RNA数据库，而非带坐标的基因组。因此我们的输出有所不同，但我希望能与praise的处理流程兼容，或者我们选择修改realignment脚本。以下是我认为有用的信息，如果你有其他问题随时问我。
在sports流程中，我们将完全一样的序列合并再进行比对。
.preprocess/seqname_map中有seqID -> illumina ID（Qname）
.preprocess/infile.txt: bowtieID -> seqID（唯一对应，基于序列经过hash function得到bowtieID）
.preprocess/infile.fa: bowtieID -> seq 

SAM: 
.bowtie_results/output_RNA: bowtieID -> 比对结果，第三列即RNAME，第一列可用bowtieID当做QNAME，后续映射回Illumina ID
REF: 
/home/lixiang/data/sports_database/hsa/*/*.fa: header与 bowtie结果的第三列一样。将数据库中所有fa合并即可。
