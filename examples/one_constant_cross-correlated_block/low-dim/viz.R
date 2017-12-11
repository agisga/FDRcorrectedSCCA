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
    if ( !("iter" %in% colnames(tmp_df)) ) {
      warning(paste("'iter' not recorded for file", filename))
      tmp_df$iter <- 1:nrow(tmp_df)
    }
    # (setting guess_max to avoid readr from guessing that a column is integer when it is double)
    if (is.null(results_df)) {
      results_df <- tmp_df
    } else {
      results_df <- bind_rows(results_df, tmp_df)
    }
  }

  return(results_df)
}


FDRcorrectedSCCA_df <- load_csv_results("constant_low_dimensional_FDRcorrectedSCCA")
SCCA_5foldCV_df <- load_csv_results("constant_low_dimensional_5foldCV")
SCCA_PMA_df <- load_csv_results("constant_low_dimensional_PMA")
ordinary_SCCA_df <- load_csv_results("constant_low_dimensional_ordinary_SCCA")

dim(FDRcorrectedSCCA_df)
# [1] 24000    13
dim(SCCA_5foldCV_df)
# [1] 8000   14
dim(SCCA_PMA_df)
# [1] 8000   14
dim(ordinary_SCCA_df)
# [1] 16000    14

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
# [1] 56000    16
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

df_for_FDR_viz <- simulation_results_df %>%
  select(method, n_signif_X, n_signif_Y, FDP_u, FDP_v) %>%
  mutate(method = factor(method, levels = names(custom_cols))) %>%
  # (put the method names in desired order)
  gather(key = "u_or_v", value = "FDP", -starts_with("n_"), -method) %>%
  mutate(u_or_v = gsub("FDP_(.)", "Estimated FDR(\\1)", u_or_v)) %>%
  group_by(n_signif_X, n_signif_Y, method, u_or_v) %>%
  summarize(FDR = mean(FDP), SE = sd(FDP) / sqrt(n()), n_repl = n()) %>%
  tbl_df() %>%
  select(n_signif_X, n_signif_Y, FDR, SE, method, u_or_v) %>%
  mutate(lwr = FDR - 2*SE, upr = FDR + 2*SE)

# plot of FDR(u)

df_for_FDR_viz %>%
  filter(u_or_v == "Estimated FDR(u)") %>%
  mutate(n_signif_Y = factor(n_signif_Y,
                             levels = sort(unique(n_signif_Y)),
                             labels = paste0("s[Y]==", sort(unique(n_signif_Y))),
                             ordered = TRUE)) %>%
  mutate(lwr = ifelse(lwr < 0, 0, lwr)) %>%
  ggplot(aes(x = n_signif_X, y = FDR)) +
    geom_segment(aes(x = 0, y = 0.05, xend = max(n_signif_X), yend = 0.05), lty = 2, color = "grey", lwd = 0.2) +
    geom_segment(aes(x = 0, y = 0.1, xend = max(n_signif_X), yend = 0.1), lty = 2, color = "grey", lwd = 0.2) +
    geom_segment(aes(x = 0, y = 0.2, xend = max(n_signif_X), yend = 0.2), lty = 2, color = "grey", lwd = 0.2) +
    geom_line(aes(color = method, linetype = method, size = method)) +
    geom_point(aes(color = method, shape = method)) +
#    geom_ribbon(aes(ymin = FDR - 2*SE, ymax = FDR + 2*SE), alpha = 0.3) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, color = method),
                  width = 2, size = 0.4, show.legend = FALSE) +
    scale_y_sqrt(breaks = c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.00)) +
    scale_color_manual(values = custom_cols) +
    scale_size_manual(values = custom_line_size) +
    scale_shape_manual(values = custom_shape) +
    facet_wrap(~n_signif_Y, nrow = 1,
               labeller = labeller(n_signif_Y = label_parsed, n_col = label_parsed)) +
    xlab(expression(s[Y])) + ylab("Estimated FDR(u) +/- 2SE") +
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.position = "none")

ggsave("./viz/low_dim_FDR(u)_combined_all_methods_nolegend.pdf")

# plot of FDR(v)

df_for_FDR_viz %>%
  filter(u_or_v == "Estimated FDR(v)") %>%
  mutate(n_signif_X = factor(n_signif_X,
                             levels = sort(unique(n_signif_X)),
                             labels = paste0("s[X]==", sort(unique(n_signif_X))),
                             ordered = TRUE)) %>%
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
    facet_wrap(~n_signif_X, nrow = 1,
               labeller = labeller(n_signif_X = label_parsed, n_col = label_parsed)) +
    xlab(expression(s[Y])) + ylab("Estimated FDR(v) +/- 2SE") +
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.position = "none")

ggsave("./viz/low_dim_FDR(v)_combined_all_methods_nolegend.pdf")

#--- Mean sensitivity

df_for_TPR_viz <- simulation_results_df %>%
  select(method, n_signif_X, n_signif_Y, TP_u, TP_v) %>%
  mutate(method = factor(method, levels = names(custom_cols))) %>%
  # (put the method names in desired order)
  gather(key = "u_or_v", value = "TP", -starts_with("n_"), -method) %>%
  mutate(u_or_v = gsub("TP_(.)", "Estimated TPR(\\1)", u_or_v)) %>%
  mutate(TPP = ifelse(u_or_v == "Estimated TPR(u)",
                      TP / n_signif_X,
                      TP / n_signif_Y)) %>%
  group_by(n_signif_X, n_signif_Y, method, u_or_v) %>%
  summarize(TPR = mean(TPP), SE = sd(TPP) / sqrt(n()), count = n()) %>%
  tbl_df() %>%
  select(n_signif_X, n_signif_Y, TPR, SE, method, u_or_v)

# plot of TPR(v)

df_for_TPR_viz %>%
  filter(u_or_v == "Estimated TPR(v)") %>%
  mutate(n_signif_X = factor(n_signif_X,
                             levels = sort(unique(n_signif_X)),
                             labels = paste0("s[X]==", sort(unique(n_signif_X))),
                             ordered = TRUE)) %>%
  ggplot(aes(x = n_signif_Y, y = TPR)) +
    geom_line(aes(color = method, linetype = method, size = method)) +
    geom_point(aes(color = method, shape = method)) +
    geom_errorbar(aes(ymin = TPR - 2*SE, ymax = TPR + 2*SE, color = method),
                  width = 2, size = 0.4, show.legend = FALSE) +
    scale_color_manual(values = custom_cols) +
    scale_size_manual(values = custom_line_size) +
    scale_shape_manual(values = custom_shape) +
    facet_wrap(~n_signif_X, nrow = 1,
               labeller = labeller(n_signif_X = label_parsed, n_col = label_parsed)) +
    xlab(expression(s[Y])) + ylab("Estimated TPR(v) +/- 2SE") +
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.position = "none")

ggsave("./viz/low_dim_TPR(v)_combined_all_methods_nolegend.pdf")

# plot of TPR(u)

df_for_TPR_viz %>%
  filter(u_or_v == "Estimated TPR(u)") %>%
  mutate(n_signif_Y = factor(n_signif_Y,
                             levels = sort(unique(n_signif_Y)),
                             labels = paste0("s[Y]==", sort(unique(n_signif_Y))),
                             ordered = TRUE)) %>%
  ggplot(aes(x = n_signif_X, y = TPR)) +
    geom_line(aes(color = method, linetype = method, size = method)) +
    geom_point(aes(color = method, shape = method)) +
    geom_errorbar(aes(ymin = TPR - 2*SE, ymax = TPR + 2*SE, color = method),
                  width = 2, size = 0.4, show.legend = FALSE) +
    scale_color_manual(values = custom_cols, guide = guide_legend(nrow = 3)) +
    scale_size_manual(values = custom_line_size, guide = guide_legend(nrow = 3)) +
    scale_shape_manual(values = custom_shape, guide = guide_legend(nrow = 3)) +
    facet_wrap(~n_signif_Y, nrow = 1,
               labeller = labeller(n_signif_Y = label_parsed, n_col = label_parsed)) +
    xlab(expression(s[Y])) + ylab("Estimated TPR(u) +/- 2SE") +
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.key.width = unit(3, "lines"))
#          legend.background = element_rect(color = "black", size = 0.2))

ggsave("./viz/low_dim_TPR(u)_combined_all_methods_legend.pdf")
