#!/usr/bin/env Rscript
# 基于 sports2 output.tsv，去掉 NO_Annotation 后重画
# 兼容 PE（Seq_5end / Seq_3end）与 SE（Sequence）
#
# 用法:
#   Rscript plot_sports_pe_no_unannotated.R <input.tsv> [outdir] [title]
#
# 示例:
#   Rscript plot_sports_pe_no_unannotated.R \
#     /home/lixiang/data/20260830_293T_15_200nt_PRAISE/control_2/sports2_out_combined/combined_R1/pe/output.tsv
#
#   Rscript plot_sports_pe_no_unannotated.R \
#     /path/to/pe/output.tsv \
#     /path/to/outdir \
#     "control_2 combined PE"
#
# 若不指定 outdir，默认写到 <input 所在目录>/plots_annotated_only
#
# 输出:
#   annotation_bar.png
#   length_stack.png      — RPM（分母 = annotated 总 reads）
#   length_facet.png      — raw Reads
#   annotation_summary.tsv

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})
source("/home/lixiang/scripts/theme_jennie.R")

# 固定配色 dictionary（name = color）；新 class 只在此追加一对即可
# 顺序即 stack/legend 默认层级
RNA_COLORS <- c(
  rsRNA         = "#377EB8",  # 蓝色 (Blue)
  rRNA          = "#377EB8",  # 蓝色 (Blue)
  msRNA         = "#4DAF4A",  # 绿
  mRNA          = "#4DAF4A",  # 绿
  tsRNA         = "#BF1534",  # 红色 (Red)
  tRNA          = "#BF1534",  # 红色 (Red)
  snsRNA        = "#984EA3",  # 紫色 (Purple)
  snRNA         = "#984EA3",  # 紫色 (Purple)
  snosRNA       = "#B15928",  # 橙褐色 (Orange-brown)
  snoRNA        = "#B15928",  # 橙褐色 (Orange-brown)
  lncsRNA       = "#E6AB02",  # 金黄色 (Gold)
  lncRNA        = "#E6AB02",  # 金黄色 (Gold)
  piRNA         = "#FBB4C4",  # 浅粉色 (Light pink)
  miRNA         = "#03A62C",  # 棕
  ysRNA         = "#66C2A5",  # 蓝绿色 (Teal)
  pseudogene    = "#666666",  # 灰色 (Gray)
  SRP_RNA       = "#1B9E77",  # 深青色 (Dark cyan)
  sncRNA        = "#7570B3",  # 紫灰色 (Purple-gray)
  RNase_MRP_RNA = "#D95F02",  # 橙色 (Orange)
  ncRNA         = "#E7298A",  # 粉红色 (Pink)
  ribozyme      = "#66A61E",  # 黄绿色 (Yellow-green)
  RNase_P_RNA   = "#A6761D",  # 土黄色 (Ochre)
  circRNA       = "#1F78B4",  # 亮蓝色 (Bright blue)
  vault_RNA     = "#B2DF8A",   # 浅绿色 (Light green)
  telomerase_RNA = "#EE5D05",  # 橘色
  tmRNA         = "#FF7F00"  # 橙色
)
RNA_CLASS_LEVELS <- names(RNA_COLORS)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || args[1] %in% c("-h", "--help")) {
  message("Usage: Rscript plot_sports_pe_no_unannotated.R <input.tsv> [outdir] [title]")
  quit(status = if (length(args) < 1) 1 else 0)
}

input <- normalizePath(args[1], mustWork = TRUE)
outdir <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  file.path(dirname(input), "plots_annotated_only")
}
plot_title <- if (length(args) >= 3 && nzchar(args[3])) {
  args[3]
} else {
  paste(rev(rev(strsplit(dirname(input), "/")[[1]])[1:2]), collapse = "/")
}

message("Input : ", input)
message("Outdir: ", outdir)
message("Title : ", plot_title)

# ---------- 读入并过滤 ----------
df <- read.delim(input, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("Reads", "Length", "Annotation")
miss <- setdiff(need, names(df))
if (length(miss)) stop("缺少列: ", paste(miss, collapse = ", "))

# PE: Seq_5end; SE: Sequence；均无则无法兜底 Length=unknown
has_pe <- "Seq_5end" %in% names(df)
has_se <- "Sequence" %in% names(df)

df <- df %>%
  filter(!Annotation %in% c("NO_Annotation", "unannotated", "Unannotated")) %>%
  mutate(
    Length_num = suppressWarnings(as.numeric(Length)),
    Length_num = ifelse(
      is.na(Length_num),
      if (has_pe) {
        nchar(Seq_5end)
      } else if (has_se) {
        nchar(Sequence)
      } else {
        NA_real_
      },
      Length_num
    )
  ) %>%
  filter(!is.na(Length_num), Length_num > 0, Reads > 0)

if (nrow(df) == 0) stop("过滤 NO_Annotation 后无剩余序列")

df_exp <- df %>%
  mutate(Annotation = gsub("\\s+", "", Annotation)) %>%
  separate_rows(Annotation, sep = ",") %>%
  group_by(ID) %>%
  mutate(n_anno = n(), Reads_frac = Reads / n_anno) %>%
  ungroup() %>%
  mutate(Class = Annotation)

# 画图前只检查一次：未知 class 直接终止（请手动追加到 RNA_COLORS）
unknown <- setdiff(unique(df_exp$Class), names(RNA_COLORS))
if (length(unknown) > 0) {
  stop(
    "未知 RNA class，已终止（请手动添加到 RNA_COLORS）。\n",
    "未知: ", paste(sort(unknown), collapse = ", "), "\n",
    "已知: ", paste(names(RNA_COLORS), collapse = ", ")
  )
}
cls_levels <- intersect(names(RNA_COLORS), unique(df_exp$Class))
message("RNA classes: ", paste(cls_levels, collapse = ", "))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

total_anno <- sum(df$Reads)
message(sprintf("Annotated unique seqs: %s", comma(nrow(df))))
message(sprintf("Annotated total reads: %s", comma(total_anno)))

# ---------- 1) annotation 汇总 + 柱状图 ----------
sum_cls <- df_exp %>%
  group_by(Class) %>%
  summarise(Reads = round(sum(Reads_frac),2), .groups = "drop") %>%
  mutate(
    Pct = round(100 * Reads / sum(Reads),2),
    Class = reorder(Class, -Reads)
  ) %>%
  arrange(desc(Reads))

write.table(
  sum_cls,
  file.path(outdir, "annotation_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 柱状图横轴：仍按丰度；颜色按固定色板（与种类绑定）
sum_cls_bar <- sum_cls %>%
  arrange(desc(Reads)) %>%
  mutate(Class = factor(as.character(Class), levels = as.character(Class)))

p_bar <- ggplot(sum_cls_bar, aes(x = Class, y = Reads, fill = Class)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Pct)), vjust = -0.3, size = 3.2) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = RNA_COLORS, breaks = cls_levels, drop = FALSE) +
  labs(
    title = paste0("Annotated RNA classes — ", plot_title),
    subtitle = sprintf("NO_Annotation removed; total annotated reads = %s", comma(total_anno)),
    x = NULL, y = "Reads"
  ) +
  theme_jennie() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.subtitle = element_text(hjust = 0.5))


save_png(p_bar, file.path(outdir, "annotation_bar"), width = 8, height = 5)

# ---------- 2) 长度分布 ----------
len_cls <- df_exp %>%
  mutate(Length = as.integer(Length_num)) %>%
  group_by(Class, Length) %>%
  summarise(Reads = sum(Reads_frac), .groups = "drop") %>%
  mutate(RPM = Reads * 1e6 / total_anno)

# stack/facet：固定色板顺序，跨样本堆积层级一致
len_cls$Class <- factor(as.character(len_cls$Class), levels = cls_levels)

len_range <- seq(min(len_cls$Length), max(len_cls$Length))
len_full <- tidyr::expand_grid(
  Class = factor(cls_levels, levels = cls_levels),
  Length = len_range
) %>%
  left_join(len_cls, by = c("Class", "Length")) %>%
  mutate(
    RPM = ifelse(is.na(RPM), 0, RPM),
    Reads = ifelse(is.na(Reads), 0, Reads)
  )

p_stack <- ggplot(len_full, aes(x = Length, y = RPM, fill = Class)) +
  geom_col(width = 1, position = "stack") +
  scale_fill_manual(values = RNA_COLORS, breaks = cls_levels, drop = FALSE) +
  labs(
    title = paste0("Length distribution (annotated only) — ", plot_title),
    subtitle = "Y = RPM among annotated reads; unannotated excluded",
    x = "Length (nt)", y = "RPM"
  ) +
  theme_jennie() +
  theme(legend.title = element_blank(), legend.position = "top") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(plot.subtitle = element_text(hjust = 0.5))

p_facet <- ggplot(len_full, aes(x = Length, y = Reads, fill = Class)) +
  geom_col(width = 1) +
  facet_wrap(~Class, scales = "free", ncol = 4) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = RNA_COLORS, breaks = cls_levels, drop = FALSE) +
  labs(
    title = paste0("Length distribution by class — ", plot_title),
    x = "Length (nt)", y = "Reads"
  ) +
  theme_jennie() +
  theme(legend.position = "none")


save_png(p_stack, file.path(outdir, "length_stack"), width = 12, height = 9)
save_png(p_facet, file.path(outdir, "length_facet"), width = 15, height = 10)

message("Done. Figures written to: ", outdir)
