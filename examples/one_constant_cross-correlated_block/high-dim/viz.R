library(tidyverse)

load_csv_results <- function(subdir) {
  results_dir <- "./results/"
  if (!grepl(".+/$", subdir)) { subdir <- paste0(subdir, "/") }
  subdir_full_path <- paste0(results_dir, subdir)
  csv_results <- list.files(subdir_full_path)

  # keep csv files only
  csv_results <- csv_results[grepl("^.+\\.csv$", csv_results)]
  csv_results <- paste0(subdir_full_path, csv_results)

  results_df <- NULL
  for (filename in csv_results) {
    tmp_df <- read_csv(filename, guess_max = 1e6)
    # (setting guess_max to avoid readr from guessing that a column is integer when it is double)
    if ( !("iter" %in% colnames(tmp_df)) ) {
      print(paste("'iter' not recorded for file", filename))
      tmp_df$iter <- 1:nrow(tmp_df)
    }
    if (is.null(results_df)) {
      results_df <- tmp_df
    } else {
      results_df <- bind_rows(results_df, tmp_df)
    }
  }

  return(results_df)
}


FDRcorrectedSCCA_df <- load_csv_results("constant_high_dimensional_FDRcorrectedSCCA")
SCCA_5foldCV_df <- load_csv_results("constant_high_dimensional_5foldCV")
SCCA_PMA_df <- load_csv_results("constant_high_dimensional_PMA")
ordinary_SCCA_df <- load_csv_results("constant_high_dimensional_ordinary_SCCA")

dim(FDRcorrectedSCCA_df)
# [1] 77000    10
dim(SCCA_5foldCV_df)
# [1] 31000    11
dim(SCCA_PMA_df)
# [1] 26000    11
dim(ordinary_SCCA_df)
# [1] 51000    11

colnames(FDRcorrectedSCCA_df)
colnames(SCCA_5foldCV_df)
colnames(SCCA_PMA_df)
colnames(ordinary_SCCA_df)

FDRcorrectedSCCA_df <- mutate(FDRcorrectedSCCA_df, method = "FDR-corrected SCCA") %>%
  mutate(method = paste0(method, " (q = ", target_FDR, ")"))
SCCA_5foldCV_df <- mutate(SCCA_5foldCV_df, method = "L1 SCCA: 5-fold CV")
SCCA_PMA_df <- mutate(SCCA_PMA_df, method = "L1 SCCA: permutation-based")
ordinary_SCCA_df <- mutate(ordinary_SCCA_df, method = "L1 SCCA: fixed penalty") %>%
  mutate(method = paste0(method, " (", lambda_X, ")"))

simulation_results_df <- bind_rows(FDRcorrectedSCCA_df,
                                   SCCA_5foldCV_df,
                                   SCCA_PMA_df,
                                   ordinary_SCCA_df)

dim(simulation_results_df)
# [1] 185000     13
table((simulation_results_df$lambda_X == simulation_results_df$lambda_Y) |
      (is.na(simulation_results_df$lambda_X) & is.na(simulation_results_df$lambda_Y)))

# check if number of iterations is consistently the same (500)

simulation_results_df %>%
  select(-n_row, -n_col_X, -n_col_Y) %>%
  group_by(n_signif_X, n_signif_Y, method) %>%
  summarize(n_iter = n(), max_iter = max(iter)) %>%
  filter(n_iter != 500)

# remove entries that were generated twice within two overlapping simulation runs

simulation_results_df <- simulation_results_df %>%
  distinct(n_signif_X, n_signif_Y, method, iter, .keep_all = TRUE)

# check again if number of iterations is consistently the same

simulation_results_df %>%
  select(-n_row, -n_col_X, -n_col_Y) %>%
  group_by(n_signif_X, n_signif_Y, method) %>%
  summarize(n_iter = n(), max_iter = max(iter)) %>%
  filter(n_iter != 500)

dim(simulation_results_df)
# [1] 171500     13
table((simulation_results_df$lambda_X == simulation_results_df$lambda_Y) |
      (is.na(simulation_results_df$lambda_X) & is.na(simulation_results_df$lambda_Y)))

#--- Estimated FDR

custom_cols = c("FDR-corrected SCCA (q = 0.2)" = "#a50026",
                "FDR-corrected SCCA (q = 0.1)" = "#d73027",
                "FDR-corrected SCCA (q = 0.05)" = "#fdae61",
                "L1 SCCA: 5-fold CV" = "#4575b4",
                "L1 SCCA: permutation-based" = "#c51b7d",
                "L1 SCCA: fixed penalty (0.3)" = "#91cf60",
                "L1 SCCA: fixed penalty (0.01)" = "#1a9850")
custom_line_size = c("FDR-corrected SCCA (q = 0.2)" = 1,
                     "FDR-corrected SCCA (q = 0.1)" = 1,
                     "FDR-corrected SCCA (q = 0.05)" = 1,
                     "L1 SCCA: 5-fold CV" = 1/2,
                     "L1 SCCA: permutation-based" = 1/2,
                     "L1 SCCA: fixed penalty (0.3)" = 1/2,
                     "L1 SCCA: fixed penalty (0.01)" = 1/2)
custom_shape = c("FDR-corrected SCCA (q = 0.2)" = 1,
                 "FDR-corrected SCCA (q = 0.1)" = 1,
                 "FDR-corrected SCCA (q = 0.05)" = 1,
                 "L1 SCCA: 5-fold CV" = 2,
                 "L1 SCCA: permutation-based" = 3,
                 "L1 SCCA: fixed penalty (0.3)" = 4,
                 "L1 SCCA: fixed penalty (0.01)" = 5)

simulation_results_df %>%
  select(method, n_signif_X, n_signif_Y, FDP) %>%
  mutate(method = factor(method, levels = names(custom_cols))) %>%
  # (put the method names in desired order)
  group_by(n_signif_X, n_signif_Y, method) %>%
  summarize(FDR = mean(FDP), SE = sd(FDP) / sqrt(n())) %>%
  tbl_df() %>%
  mutate(n_signif_X = factor(n_signif_X,
                             levels = sort(unique(n_signif_X)),
                             labels = paste0("s[X]==", sort(unique(n_signif_X))),
                             ordered = TRUE)) %>%
  select(n_signif_X, n_signif_Y, FDR, SE, method) %>%
  mutate(lwr = FDR - 2*SE, upr = FDR + 2*SE) %>%
  mutate(lwr = ifelse(lwr < 0, 0, lwr)) %>%
  ggplot(aes(x = n_signif_Y, y = FDR)) +
    geom_segment(aes(x = 0, y = 0.05, xend = max(n_signif_Y), yend = 0.05), lty = 2, color = "grey", lwd = 0.2) +
    geom_segment(aes(x = 0, y = 0.1, xend = max(n_signif_Y), yend = 0.1), lty = 2, color = "grey", lwd = 0.2) +
    geom_segment(aes(x = 0, y = 0.2, xend = max(n_signif_Y), yend = 0.2), lty = 2, color = "grey", lwd = 0.2) +
    geom_line(aes(color = method, linetype = method, size = method)) +
    geom_point(aes(color = method, shape = method)) +
#    geom_ribbon(aes(ymin = FDR - 2*SE, ymax = FDR + 2*SE), alpha = 0.3) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, color = method),
                  width = 2, size = 0.4, show.legend = FALSE) +
    scale_y_sqrt(breaks = c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.00)) +
    scale_color_manual(values = custom_cols) +
    scale_size_manual(values = custom_line_size) +
    scale_shape_manual(values = custom_shape) +
    facet_wrap(~n_signif_X, nrow = 2,
               labeller = labeller(n_signif_X = label_parsed, n_col = label_parsed)) +
    xlab(expression(s[Y])) + ylab("Estimated FDR(v) +/- 2SE") +
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.position = c(0.88, 0.25),
          legend.key.width = unit(3, "lines"))
#          legend.background = element_rect(color = "black", size = 0.2))

ggsave("./viz/high_dim_FDR_combined.pdf")

#--- Mean sensitivity (TPP)

simulation_results_df %>%
  select(method, n_signif_X, n_signif_Y, TP) %>%
  mutate(method = factor(method, levels = names(custom_cols))) %>%
  # (put the method names in desired order)
  mutate(TPP = TP / n_signif_Y) %>%
  group_by(n_signif_X, n_signif_Y, method) %>%
  summarize(TPR = mean(TPP), SE = sd(TPP) / sqrt(n())) %>%
  tbl_df() %>%
  mutate(n_signif_X = factor(n_signif_X,
                             levels = sort(unique(n_signif_X)),
                             labels = paste0("s[X]==", sort(unique(n_signif_X))),
                             ordered = TRUE)) %>%
  select(n_signif_X, n_signif_Y, TPR, SE, method) %>%
  ggplot(aes(x = n_signif_Y, y = TPR)) +
    geom_line(aes(color = method, linetype = method, size = method)) +
    geom_point(aes(color = method, shape = method)) +
    geom_errorbar(aes(ymin = TPR - 2*SE, ymax = TPR + 2*SE, color = method),
                  width = 2, size = 0.4, show.legend = FALSE) +
    scale_color_manual(values = custom_cols) +
    scale_size_manual(values = custom_line_size) +
    scale_shape_manual(values = custom_shape) +
    facet_wrap(~n_signif_X, nrow = 2,
               labeller = labeller(n_signif_X = label_parsed, n_col = label_parsed)) +
    xlab(expression(s[Y])) + ylab("Estimated TPR(v) +/- 2SE") +
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.position = c(0.88, 0.25),
          legend.key.width = unit(3, "lines"))
#          legend.background = element_rect(color = "black", size = 0.2))

ggsave("./viz/high_dim_TPR_combined.pdf")
