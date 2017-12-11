library(tidyverse)

load_csv_results <- function(path_to_csv_files) {
  if (!grepl(".+/$", path_to_csv_files)) {
    path_to_csv_files <- paste0(path_to_csv_files, "/")
  }
  csv_results <- list.files(path_to_csv_files)

  # keep csv files only
  csv_results <- csv_results[grepl("^.+\\.csv$", csv_results)]
  csv_results <- paste0(path_to_csv_files, csv_results)

  results_df <- NULL
  n_incomplete <- 0
  n_read <- 0
  for (i in 1:length(csv_results)) {
    filename <- csv_results[i]
    print(paste(i, "- reading:", filename))
    tmp_df <- read_csv(filename)

    # check (dataframe size hard coded...)
    if (nrow(tmp_df) != 4 * 5) {
      print(paste(filename, "is an incomplete file."))
      n_incomplete <- n_incomplete + 1
      next
    }

    if (is.null(results_df)) {
      results_df <- tmp_df
    } else {
      results_df <- bind_rows(results_df, tmp_df)
    }

    n_read <- n_read + 1
    # only need to read in 100 dataframes
    if (n_read == 100) break
  }

  print(paste("Done reading", n_read, "of", length(csv_results), "files."))
  print(paste(n_incomplete, "files are incomplete."))

  return(results_df)
}

simulation_results_df <- load_csv_results("./results/high_dim/")
n_row <- unique(simulation_results_df$n_row)

# proportion of TP identified
simulation_results_df %>%
  filter(n_signif_blocks %in% c(1, 2, 3, 15)) %>%
  mutate(sensitivity = TP / n_signif_Y) %>% # aka recall or TPR
  arrange(as.numeric(n_signif_blocks)) %>%
  mutate(n_signif_blocks = factor(n_signif_blocks, levels = unique(n_signif_blocks),
         labels = paste(unique(n_signif_blocks), "cross-correlated blocks"),
         ordered = TRUE)) %>%
  ggplot(aes(x = sensitivity, y = (..count..)/(sum(..count..)/4))) +
    geom_histogram(bins = 15) +
    facet_wrap(~n_signif_blocks, scales = "free_y") +
    theme_bw() +
    labs(x = "TPP(v)", y = "Proportion")

ggsave("./figures/sensitivity_1_2_3_15_signif_blocks.png")
ggsave("./figures/sensitivity_1_2_3_15_signif_blocks.pdf")
