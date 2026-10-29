

# 堡垒机跑PRAISE 20260930 

## First round trim adapter 
```bash
sample=praise_1
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
运行 process_cutadapt_info_parquet.py 
```

## second round trim adapter 
```bash
# times 参数(best)
mkdir -p $sample/cutadapt-round2
cutadapt -j 32 \
  -g "5AD=ACAGTCCGACGATC" \
  -g "5AD=GTTCTACAGTCCGACGATC" \
  -O 5 -e 0.1 --times 3 \
  -q 20 --nextseq-trim=20 \
  -m 10 \
  -o ${sample}_R1_bothAD_cut5AD_g_3times.fastq.gz \
  ../categorize/${sample}_R1_bothAD.fq.gz > log-R1-bothAD-5AD-g-3times.log 2>&1
  ```


### compare two ways of cutadapt 
```bash
cd /home/lixiang/data/20260930_293T_PRAISE/control_1/cutadapt-round2
A=control_1_R1_bothAD_cut5AD_g_g.fastq.gz
B=control_1_R1_bothAD_cut5AD_g.fastq.gz
RAW=/home/lixiang/data/20260930_293T_PRAISE/rawdata/control_1.R1.raw.fastq.gz
OUT=compare_cut5AD_g.tsv

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


## 提umi
```bash
mkdir -p cut_umi/log

umi_tools extract \
  --extract-method=regex \
  --bc-pattern="^(?P<umi_1>.{5})(?P<seq>.+)(?P<umi_2>.{5})$" \
  --stdin=cutadapt-round2/${sample}_R1_bothAD_cut5AD_g_3times.fastq.gz \
  --stdout=cut_umi/${sample}_R1_bothAD_cut5AD_g_3times.fq.gz \
  --log=cut_umi/log/umi-from-cutadapt-g-3times.log

cd cut_umi
gunzip *.gz
```


# SPORTS 单端比对

## 自定义版本

sports2 创建新环境：
```bash
cd 文件夹
mm create -n sports python==3.10
mm install -c conda-forge -c bioconda bowtie bowtie2 seqkit taxonkit trimmomatic r-base r-ggplot2 r-stringr r-data.table gmap 
pip install .
```
sports2 运行
```bash
mm activate sports
sample=praise_1
output=/home/lixiang/data/20260930_293T_15_200nt_PRAISE/${sample}/sports2_out_cutumi_from_cut5ad_g_3times
fq1=/home/lixiang/data/20260930_293T_15_200nt_PRAISE/${sample}/cut_umi/${sample}_R1_bothAD_cut5AD_g_3times.fq 

sports -i $fq1        -ry /mnt/genome/sports2_reference/hsa/ref.yaml       -l 15 -L 200        -V -z -d -p 32 -o $output        --minscore L,-0.8,-0.8        --seedmismatch 1         --seedlen 8 
```
注意：需要修改单双引号的问题


## 公共版本
```bash
sample=control_1
output=/home/lixiang/data/20260930_293T_15_200nt_PRAISE/${sample}/sports2_out_cutumi_from_cut5ad_g_3times
fq1=/home/lixiang/data/20260930_293T_15_200nt_PRAISE/${sample}/cut_umi/${sample}_R1_bothAD_cut5AD_g_3times.fq 

mm activate /mnt/usrData/xielefanData/miniforge3/envs/ncrna

sports -i $fq1        -ry /mnt/genome/sports2_reference/hsa/ref.yaml       -l 15 -L 200        -V -z -d -p 32 -o $output        --minscore L,-0.8,-0.8        --seedmismatch 1         --seedlen 8 
```


# 总结 SPORTS2输出
less output.tsv | grep NO_Annotation | less
less output.tsv | grep NO_Annotation | awk '{print ">" $1; print $2}' > to_blast

# --- SE（bothAD 等）---
cd control_1/sports2_out_bothAD 

awk -F'\t' 'NR>1 { a[$6]+=$3; t+=$3 }
END {
  for (c in a) {
    n=c; if(n=="NO_Annotation") n="unannotated"
    printf "%s\t%.0f\t%.2f\n", n, a[c], a[c]*100/t
  }
}' */output.tsv | sort -k2,2nr > annotation_summary_all.tsv


control-1
rsRNA   380607  37.38
tsRNA   325593  31.98
unannotated     167648  16.47
msRNA   60838   5.98
snsRNA  31197   3.06
snosRNA 20850   2.05
lncsRNA 9554    0.94
piRNA   9015    0.89
miRNA   5762    0.57
ysRNA   2803    0.28
SRP_RNA,lncsRNA 2645    0.26
RNase_MRP_RNA,lncsRNA,ribozyme  388     0.04
lncsRNA,vault_RNA       273     0.03
RNase_P_RNA,lncsRNA,ribozyme    249     0.02
RNase_MRP_RNA,lncsRNA   187     0.02
RNase_P_RNA,lncsRNA     177     0.02
lncsRNA,pseudogene      126     0.01
sncRNA  110     0.01
pseudogene      66      0.01
ncRNA   44      0.00
SRP_RNA 31      0.00
ncRNA,sncRNA    13      0.00
tmRNA   6       0.00
lncsRNA,telomerase_RNA  4       0.00
RNase_MRP_RNA   3       0.00
RNase_MRP_RNA,ribozyme  3       0.00
RNase_P_RNA     1       0.00
