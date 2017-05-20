# FDRcorrectedSCCA  Copyright (C) 2017  Alexej Gossmann

library(dplyr)
library(tidyr)
library(ggplot2)

load(file = "./results/with_piece-wise_constant_Sigma.RData")
low_dimensional <- simulation_results_df
low_dimensional <- mutate(low_dimensional, n_col_X = n_col_X, n_col_Y = n_col_Y)

load(file = "./results/with_piece-wise_constant_Sigma_and_high-dimensional_FDR_5.RData")
high_dimensional_1 <- simulation_results_df
high_dimensional_1 <- mutate(high_dimensional_1, n_col_X = n_col_X, n_col_Y = n_col_Y)

load(file = "./results/with_piece-wise_constant_Sigma_and_high-dimensional_FDR_10.RData")
high_dimensional_2 <- simulation_results_df
high_dimensional_2 <- mutate(high_dimensional_2, n_col_X = n_col_X, n_col_Y = n_col_Y)

load(file = "./results/with_piece-wise_constant_Sigma_and_high-dimensional_FDR_20.RData")
high_dimensional_3 <- simulation_results_df
high_dimensional_3 <- mutate(high_dimensional_3, n_col_X = n_col_X, n_col_Y = n_col_Y)

simulation_results_df <- bind_rows(low_dimensional,
                                   high_dimensional_1,
                                   high_dimensional_2,
                                   high_dimensional_3)

# Estimated FDR

simulation_results_df %>%
  filter(n_row == 500) %>%
  mutate(target_FDR = factor(target_FDR, levels = c(0.05, 0.1, 0.2)),
         n_col = paste("list(n==", n_row, ", p[X]==", n_col_X, ", p[Y]==", n_col_Y, ")")) %>%
  mutate(n_col_vs_n_signif_X_vs_n_signif_Y_vs_target_FDR =
           paste(n_col, n_signif_X, n_signif_Y, target_FDR, sep = "_")) %>%
  mutate(n_col_vs_n_signif_X_vs_n_signif_Y =
           factor(n_col_vs_n_signif_X_vs_n_signif_Y_vs_target_FDR,
                  levels = rev(sort(unique(n_col_vs_n_signif_X_vs_n_signif_Y_vs_target_FDR))),
                  ordered = TRUE)) %>%
  group_by(n_col_vs_n_signif_X_vs_n_signif_Y_vs_target_FDR) %>%
  summarize(FDR = mean(FDP), SE = sd(FDP) / sqrt(length(FDP)),
            n_col = unique(n_col), target_FDR = unique(target_FDR),
            n_signif_X = unique(n_signif_X), n_signif_Y = unique(n_signif_Y)) %>%
  mutate(n_col = factor(n_col, rev(sort(unique(n_col))), ordered = TRUE)) %>%
  mutate(n_signif_X = factor(paste0("s[X]==", n_signif_X))) %>%
  select(n_signif_X, n_signif_Y, FDR, target_FDR, SE, n_col) %>%
  ggplot(aes(x = n_signif_Y, y = FDR)) +
    geom_line(aes(color = target_FDR, linetype = target_FDR)) +
    scale_color_manual(values = 1:3, breaks = c(0.05, 0.1, 0.2),
                       labels = list(bquote(q==0.05), bquote(q==0.1), bquote(q==0.2))) +
    scale_linetype_manual(values = c(3, 1, 2),
                          labels = list(bquote(q==0.05), bquote(q==0.1), bquote(q==0.2))) +
    geom_errorbar(aes(ymin = FDR - 2*SE, ymax = FDR + 2*SE, color = target_FDR),
                  width = 2, size = 0.4, show.legend = FALSE) +
    facet_grid(n_col ~ n_signif_X, scales = "free_y",
               labeller = labeller(n_signif_X = label_parsed, n_col = label_parsed)) +
    xlab(expression(s[Y])) + ylab("Estimated FDR(v) +/- 2SE") +
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.position = "bottom")

ggsave("./results/with_piece-wise_constant_Sigma_FDR_graph.pdf")

# bar graph: mean number of true positives
simulation_results_df %>%
  filter(n_row == 500 & target_FDR == 0.1) %>%
  mutate(n_col = paste("list(n==", n_row, ", p[X]==", n_col_X, ", p[Y]==", n_col_Y, ")")) %>%
  mutate(n_col_vs_n_signif_X_vs_n_signif_Y = paste(n_col, n_signif_X,
                                                   n_signif_Y, sep = "_")) %>%
  mutate(n_col_vs_n_signif_X_vs_n_signif_Y = factor(n_col_vs_n_signif_X_vs_n_signif_Y,
               levels = rev(sort(unique(n_col_vs_n_signif_X_vs_n_signif_Y))),
               ordered = TRUE)) %>%
  group_by(n_col_vs_n_signif_X_vs_n_signif_Y) %>%
  summarize(mean_TP = mean(TP), n_col = unique(n_col),
            se_TP = sd(TP) /  sqrt(length(TP)),
            n_signif_X = unique(n_signif_X), n_signif_Y = unique(n_signif_Y)) %>%
  mutate(n_signif_X = factor(paste0("s[X]==", n_signif_X))) %>%
  mutate(n_col = factor(n_col, rev(sort(unique(n_col))), ordered = TRUE)) %>%
  ggplot(aes(x = factor(n_signif_Y), y = mean_TP)) +
    geom_bar(stat = "identity", fill = "lightgrey", color = "black") +
    geom_errorbar(aes(ymin = mean_TP - 2 * se_TP, ymax = mean_TP + 2 * se_TP),
                  width = 0.2, size = 0.4, color = "black") +
    facet_grid(n_col ~ n_signif_X,
               labeller = labeller(n_signif_X = label_parsed, n_col = label_parsed)) +
    xlab(expression(s[Y])) + ylab("TP(v) +/- 2SE") +
    theme_bw()

ggsave("./results/with_piece-wise_constant_Sigma_TP_graph.pdf")

# histograms of proportions of true positives and false positives (combined for all n_signif conditions across the board)
simulation_results_df %>%
  filter(n_row == 500 & target_FDR == 0.1) %>%
  mutate(TPP = TP / n_signif_Y) %>% # aka recall or sensitivity
  select(TPP, FDP) %>%
  gather() %>%
  ggplot(aes(value)) +
    geom_histogram(bins = 20, fill = "lightgrey", color = "black") +
    facet_wrap(~key) +
    theme_bw() +
    theme(axis.title = element_blank())

ggsave("./results/with_piece-wise_constant_Sigma_FDP_and_TPP_hist.pdf")
