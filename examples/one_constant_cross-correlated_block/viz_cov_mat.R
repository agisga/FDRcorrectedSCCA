library(tidyverse)

n_col_X <- 5
n_col_Y <- 5
n_signif_X <- 2
n_signif_Y <- 2
n_col <- n_col_X + n_col_Y
cor_btwn <- 0.4
cor_wthn <- 0.1

# generate the correlation matrix of matrix [X | Y]
Sigma <- matrix(0, n_col, n_col)
Sigma[1:n_col_X, 1:n_col_X] <- cor_wthn
Sigma[(n_col_X+1):n_col, (n_col_X+1):n_col] <- cor_wthn
Sigma[1:n_signif_X, (n_col_X+1):(n_col_X+n_signif_Y)] <- cor_btwn
Sigma[(n_col_X+1):(n_col_X+n_signif_Y), 1:n_signif_X] <- cor_btwn
Sigma[1:n_signif_X, 1:n_signif_X] <- Sigma[1:n_signif_X, 1:n_signif_X] + cor_btwn
Sigma[(n_col_X+1):(n_col_X+n_signif_Y), (n_col_X+1):(n_col_X+n_signif_Y)] <- Sigma[(n_col_X+1):(n_col_X+n_signif_Y), (n_col_X+1):(n_col_X+n_signif_Y)] + cor_btwn
diag(Sigma) <- 1

Sigma_df <- as_data_frame(Sigma) %>%
  mutate(y = 1:n_col) %>%
  gather(x, value, -y) %>%
  mutate(x = as.numeric(gsub("V", "", x)))

Sigma_df %>% ggplot(aes(x, y)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value)) +
  scale_y_reverse() +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme_void() +
  theme(legend.title = element_blank())

ggsave("Sigma.pdf")
