suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
  library(grid)
  library(stringr)
  library(UpSetR)
})
source("/home/lixiang/scripts/theme_jennie.R")

out_dir <- "/home/lixiang/data/20260930_293T_PRAISE/plot_mod"
dir.create(out_dir, showWarnings = FALSE)

## 1. Data loading ----------------
praise_rep1 <- read.table(
  "/home/lixiang/data/20260930_293T_PRAISE/praise_1/praise_from_sports/tRNA_hq_mature/plot.tsv",
  header = TRUE,
  colClasses = c(ref_base = "character")
)
control_rep1 <- read.table(
  "/home/lixiang/data/20260930_293T_PRAISE/control_1/praise_from_sports/tRNA_hq_mature/plot.tsv",
  header = TRUE,
  colClasses = c(ref_base = "character")
)

samples <- c("praise_rep1", "control_rep1")

parse_plot_tsv <- function(raw, sample) {
  raw %>%
    tidyr::separate(
      chr_name,
      into = c("rna_type", "isotype", "source", "copy_id",
               "symbol", "high_qual", "tRNA_len", "src_db"),
      sep = ",", convert = TRUE
    ) %>%
    dplyr::rename(pos = chr_index) %>%
    mutate(
      sample = sample,
      isotype_class = if_else(str_detect(isotype, "^mature-mt_"), "mt", "cytosolic"),
      isotype = sub("^mature-(mt_)?", "", isotype),
      mod_id = str_c(copy_id, pos, sep = ":"),
      type = "Psi",
      unmod_count = coverage - del_count
    )
}

df <- parse_plot_tsv(praise_rep1, "praise_rep1")
df2 <- parse_plot_tsv(control_rep1, "control_rep1")
df_all <- bind_rows(df, df2)
processed_list <- list(praise_rep1 = df, control_rep1 = df2)


## basic stats ----------------
mod_stats <- imap_dfr(processed_list, ~{
  tibble(
    Sample = .y,
    n_tsRNA = n_distinct(.x$copy_id),
    n_psi_mod = nrow(.x),
    mod_per_tsRNA = round(n_psi_mod / n_tsRNA, 2)
  )
})
print(mod_stats)

mod_by_isotype <- df_all %>%
  group_by(sample, isotype_class, isotype) %>%
  summarise(
    n_copies = n_distinct(copy_id),
    mod_sites = sum(del_count > 0),
    .groups = "drop"
  ) %>%
  arrange(sample, desc(mod_sites))
print(mod_by_isotype)
write.table(
  mod_by_isotype,
  file.path(out_dir, "mod_sites_by_isotype.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


## Plot UpSet ----------------
plot_upset <- function(data_list,
                       sets_order = names(data_list),
                       outfile = file.path(out_dir, "upset_plot.pdf")) {
  upset_matrix <- fromList(data_list)
  pdf(outfile, width = 8, height = 5)
  # upset() only returns an object. Drawing happens in print(),
  # which RStudio does automatically and a function call does not.
  print(UpSetR::upset(
    upset_matrix,
    sets = sets_order,
    order.by = "freq",
    keep.order = TRUE,
    sets.bar.color = "#f2a93b",
    main.bar.color = "#2c7fb8",
    text.scale = 1.2
  ))
  # print() leaves a layout viewport pushed; go back to the page
  # before adding the title, otherwise the text replaces the plot.
  popViewport(0)
  grid.text(
    expression("tRNA " * psi * " modification"),
    x = 0.5, y = 0.98,
    gp = gpar(fontsize = 14, fontface = "bold")
  )
  dev.off()
}

upset_data <- map(processed_list, ~ .x$mod_id)
plot_upset(upset_data, sets_order = samples)


## GtRNAdb mature reference, one sequence per copy ----------------
# 名字形如 Homo_sapiens_tRNA-Ala-AGC-1-1，对应 copy_id mature-tRNA-Ala-AGC-1-1。
# 文件里是 RNA（U）。横轴改成 T，和测序碱基一致。
read_gtrna_mature <- function(path) {
  lines <- readLines(path)
  header_at <- which(startsWith(lines, ">"))
  gtrna_id <- sub("^>Homo_sapiens_", "", lines[header_at])
  gtrna_id <- sub(" .*", "", gtrna_id)

  ref_seq <- character(length(header_at))
  for (i in seq_along(header_at)) {
    start <- header_at[i] + 1
    end <- if (i < length(header_at)) header_at[i + 1] - 1 else length(lines)
    ref_seq[i] <- paste(lines[start:end], collapse = "")
  }

  tibble(gtrna_id = gtrna_id, ref_seq = chartr("U", "T", ref_seq))
}

gtrna_ref <- read_gtrna_mature("/home/lixiang/database/GtRNAdb/hg38-mature-tRNAs.fa")

# SPORTS 的成熟序列 = GtRNAdb 序列 + 末尾 CCA。His 还在 5' 多一个 G。
# offset 把 plot.tsv 的 pos 换成 GtRNAdb 上的位置。线粒体不在这个 fasta 里，会被丢掉。
df_plot <- df_all %>%
  mutate(gtrna_id = sub("^mature-", "", copy_id)) %>%
  inner_join(gtrna_ref, by = "gtrna_id") %>%
  mutate(
    matches_at_1 = str_starts(ref_fa, ref_seq),
    matches_at_2 = str_starts(str_sub(ref_fa, 2), ref_seq),
    offset = case_when(
      matches_at_1 ~ 0L,
      matches_at_2 ~ 1L,
      TRUE ~ NA_integer_
    ),
    gtrna_pos = pos - offset
  ) %>%
  filter(!is.na(offset), gtrna_pos >= 1, gtrna_pos <= nchar(ref_seq))


## position plots ----------------
# 一张图一种 tRNA（如 tRNA-Glu-TTC）。
# 里面每个 copy 一块；praise 和 control 各占一行，避免 y 轴挤在一起。
focus_aa <- c("Ala", "Cys", "Glu", "Gln")
types <- df_plot %>%
  filter(str_detect(isotype, paste(focus_aa, collapse = "|"))) %>%
  distinct(isotype) %>%
  arrange(isotype) %>%
  pull(isotype)

for (iso in types) {
  copy_ids <- df_plot %>%
    filter(isotype == iso) %>%
    distinct(gtrna_id) %>%
    arrange(gtrna_id) %>%
    pull(gtrna_id)

  plots <- list()
  n_rows <- numeric(length(copy_ids))

  for (i in seq_along(copy_ids)) {
    one <- df_plot %>% filter(gtrna_id == copy_ids[i])
    bases <- strsplit(one$ref_seq[1], "", fixed = TRUE)[[1]]
    seq_len <- length(bases)
    n_rows[i] <- n_distinct(one$sample)

    stacked <- one %>%
      pivot_longer(c(del_count, unmod_count), names_to = "group", values_to = "count") %>%
      mutate(
        group = if_else(group == "del_count", "Psi", "Unmodified"),
        group = factor(group, levels = c("Unmodified", "Psi")),
        sample = factor(sample, levels = samples)
      )

    plots[[i]] <- ggplot(stacked, aes(gtrna_pos, count, fill = group)) +
      geom_col(width = 0.9, position = "stack") +
      geom_text(
        data = one,
        aes(gtrna_pos, coverage, label = del_count),
        inherit.aes = FALSE,
        vjust = -0.4,
        size = 5
      ) +
      facet_wrap(~ sample, ncol = 1, scales = "free_y") +
      scale_fill_manual(values = c("Unmodified" = "grey80", "Psi" = "#1f78b4")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
      scale_x_continuous(
        breaks = 1:seq_len,
        labels = bases,
        limits = c(0.5, seq_len + 0.5),
        expand = c(0, 0)
      ) +
      labs(title = copy_ids[i], x = NULL, y = "Coverage") +
      theme_jennie() +
      theme(
        legend.position = if (i == length(copy_ids)) "bottom" else "none",
        legend.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(margin = margin(r = 10)),
        panel.spacing.y = unit(16, "pt"),
        plot.margin = margin(10, 12, 10, 18),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
      )
  }

  png(
    file.path(out_dir, paste0("plot_", iso, ".png")),
    width = 16,
    height = 2.6 * sum(n_rows),
    units = "in",
    res = 150
  )
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(length(plots), 1, heights = n_rows)))
  for (i in seq_along(plots)) {
    print(
      plots[[i]],
      vp = viewport(layout.pos.row = i, layout.pos.col = 1),
      newpage = FALSE
    )
  }
  dev.off()
}
