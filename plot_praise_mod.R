# Plotting libraries and the shared theme.
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
  header = TRUE, sep = "\t", quote = "",
  colClasses = c(ref_base = "character", del_seqid = "character", non_del_seqid = "character")
)
control_rep1 <- read.table(
  "/home/lixiang/data/20260930_293T_PRAISE/control_1/praise_from_sports/tRNA_hq_mature/plot.tsv",
  header = TRUE, sep = "\t", quote = "",
  colClasses = c(ref_base = "character", del_seqid = "character", non_del_seqid = "character")
)

samples <- c("praise_rep1", "control_rep1")

# Split chr_name, tag the sample, and split coverage into deletion vs unmodified.
parse_plot_tsv <- function(raw, sample) {
  raw %>%
    tidyr::separate(
      chr_name,
      into = c("rna_type", "isotype", "source", "copy_id",
               "symbol", "high_qual", "tRNA_len", "src_db"),
      sep = ",", convert = TRUE, remove = FALSE
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

## 2. Basic statistics ----------------
# How many modified tRNAs and sites each sample has.
mod_stats <- imap_dfr(processed_list, ~{
  tibble(
    Sample = .y,
    n_tsRNA = n_distinct(.x$copy_id),
    n_psi_mod = nrow(.x),
    mod_per_tsRNA = round(n_psi_mod / n_tsRNA, 2)
  )
})
# print(mod_stats)

mod_by_isotype <- df_all %>%
  group_by(sample, isotype_class, isotype) %>%
  summarise(
    n_copies = n_distinct(copy_id),
    mod_sites = sum(del_count > 0),
    .groups = "drop"
  ) %>%
  arrange(sample, desc(mod_sites))
# print(mod_by_isotype)
write.table(
  mod_by_isotype,
  file.path(out_dir, "mod_sites_by_isotype.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## 3. Upset plot ----------------
# Which modification sites are shared by praise and control.
plot_upset <- function(data_list,
                       sets_order = names(data_list),
                       outfile = file.path(out_dir, "upset_plot.pdf")) {
  upset_matrix <- fromList(data_list)
  pdf(outfile, width = 8, height = 5)
  print(UpSetR::upset(
    upset_matrix,
    sets = sets_order,
    order.by = "freq",
    keep.order = TRUE,
    sets.bar.color = "#f2a93b",
    main.bar.color = "#2c7fb8",
    text.scale = 1.2
  )) 
  popViewport(0)  # leave UpSetR's viewport before adding the title
  grid.text(
    expression("tRNA " * psi * " modification"), x = 0.5, y = 0.98,
    gp = gpar(fontsize = 14, fontface = "bold")
  )
  dev.off()
}

upset_data <- map(processed_list, ~ .x$mod_id)
plot_upset(upset_data, sets_order = samples)


# X axis is ref_fa. del_seqid / non_del_seqid are not drawn.
df_plot <- df_all

# Exact-sequence MODOMICS matches (empty modomics_sequence rows dropped).
mod_join <- read.table(
  "/home/lixiang/data/database/modomics/join_modomics_tRNA-hq-mature.tsv",
  header = TRUE, sep = "\t", quote = "",
  colClasses = "character", comment.char = ""
) %>%
  filter(nzchar(modomics_sequence)) %>%
  distinct(copy_id, modomics_id, modomics_sequence)

## 4. Plotting ----------------
# One PNG per isotype. Each copy is its own plot: samples stacked.

focus_aa <- c("Ala", "Cys", "Glu", "Gln")
types <- df_plot %>%
  filter(str_detect(isotype, paste(focus_aa, collapse = "|"))) %>%
  distinct(isotype) %>%
  arrange(isotype) %>%
  pull(isotype)

for (iso in types) {
  copy_ids <- df_plot %>%
    filter(isotype == iso) %>%
    distinct(copy_id) %>%
    arrange(copy_id) %>%
    pull(copy_id)

  plots <- list()
  n_rows <- numeric(length(copy_ids))

  for (i in seq_along(copy_ids)) {
    one <- df_plot %>% filter(copy_id == copy_ids[i])
    bases <- strsplit(one$ref_fa[1], "", fixed = TRUE)[[1]]
    seq_len <- length(bases)
    n_rows[i] <- n_distinct(one$sample)

    # X axis is ref_fa, add modomics sequences
    x_labels <- paste0(bases, "\n", seq_along(bases))
    plot_title <- one$copy_id[1]
    mod_hit <- mod_join %>% filter(copy_id == copy_ids[i])
    # One equal-length match: add MODOMICS row on the x axis and id in the title.
    if (nrow(mod_hit) == 1L && nchar(mod_hit$modomics_sequence[1]) == seq_len) {
      mod_bases <- strsplit(mod_hit$modomics_sequence[1], "", fixed = TRUE)[[1]]
      x_labels <- paste0(bases, "\n", seq_along(bases), "\n", mod_bases)
      plot_title <- paste0(one$copy_id[1], " / MODOMICS ", mod_hit$modomics_id[1])
    }
    # Sites with deletion: paint the whole x-tick (ref / pos / modomics) red.
    del_pos <- sort(unique(one$pos[one$del_count > 0]))
    axis_cols <- ifelse(seq_along(bases) %in% del_pos, "red", "black")

    stacked <- one %>%
      pivot_longer(c(del_count, unmod_count), names_to = "group", values_to = "count") %>%
      mutate(
        group = if_else(group == "del_count", "Psi", "Unmodified"),
        group = factor(group, levels = c("Unmodified", "Psi")),
        sample = factor(sample, levels = samples)
      )

    plots[[i]] <- ggplot(stacked, aes(pos, count, fill = group)) +
      geom_col(width = 0.9, position = "stack") +
      # bar label texts: deletion/coverage
      geom_text(
        data = one, aes(pos, coverage, label = paste0(del_count, "/", coverage)),
        inherit.aes = FALSE, vjust = -0.4, size = 4
      ) +
      # Sample name inside the panel, top-left. vjust > 1 moves it down from the top edge.
      geom_text(
        data = one %>% distinct(sample), aes(x = 0.7, y = Inf, label = sample),
        inherit.aes = FALSE, hjust = -0.1, vjust = 1.9, size = 5
      ) +
      facet_wrap(~ sample, ncol = 1, scales = "free_y") +
      coord_cartesian(clip = "off") +
      scale_fill_manual(values = c("Unmodified" = "grey80", "Psi" = "#1f78b4")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
      scale_x_continuous(
        breaks = 1:seq_len,
        labels = x_labels,
        limits = c(0.5, seq_len + 0.5),
        expand = c(0, 0)
      ) +
      labs(title = plot_title, x = NULL, y = "Coverage") +
      theme_jennie() +
      theme(
        legend.position = if (i == length(copy_ids)) "bottom" else "none",
        legend.title = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size = 10, lineheight = 0.85, colour = axis_cols),
        axis.title.y = element_text(size = 15,margin = margin(r = 10)),
        panel.spacing.y = unit(10, "pt"),   # space between the praise and control panels
        plot.title = element_text(size = 15, margin = margin(b = 2)),
        plot.margin = margin(10, 12, 14, 18),  # top, right, bottom, left
        panel.border = element_rect(colour = "gray50", fill = NA, linewidth = 0.6)
      )
  }

  png(
    file.path(out_dir, paste0("plot_", iso, ".png")),
    width = 16, height = 3.2 * sum(n_rows), units = "in", res = 150
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
