library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)

df <- read_csv("analysis/trophic_data.csv")
df$type <- as.factor(df$type)
df$type <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", df$type, perl=TRUE)

vol_simp <- function(S_ind, K=3, m=rep(1.0, 3)) {
  (K*S_ind)^K / (factorial(K) * prod(m))
}

df <- df %>%
  mutate(norm_vol = vol / vol_simp(S_ind)) %>%
  group_by(S_ind, type) %>%
  mutate(rep_id = row_number()) %>%
  ungroup()

median_df <- df %>%
  group_by(S_ind, type) %>%
  summarise(median_norm = median(norm_vol), .groups = "drop")

p <- ggplot() +
  geom_line(data = df,
    aes(x = S_ind, y = norm_vol, group = interaction(type, rep_id), color = type),
    alpha = 0.1, linewidth = 0.5
  ) +
  geom_line(data = median_df,
    aes(x = S_ind, y = median_norm, color = type),
    alpha = 0.8, linewidth = 1
  ) +
  scale_color_manual(values = c("blue", "red", "green", "purple")) +
  coord_cartesian(xlim = c(0,10.0), ylim = c(0,0.8)) +
  labs(x = "Ind. Power Cap", y = "Norm. Vol.") +
  facet_wrap(~ type) +
  theme_pubr(base_size=12) +
  theme(legend.position = "none")

ggsave("figures/compare_0.pdf", p, width = 6, height = 4.5, dpi = 450)



df_plain <- read_csv("analysis/plain_data.csv")
df_plain$type <- as.factor(df_plain$type)

df_plain <- df_plain %>%
  mutate(norm_vol = vol / vol_simp(S_ind)) %>%
  group_by(S_ind, type) %>%
  mutate(rep_id = row_number()) %>%
  ungroup()

median_df <- df_plain %>%
  group_by(S_ind, type) %>%
  summarise(median_norm = median(norm_vol), .groups = "drop")

p1 <- ggplot() +
  geom_line(data = df_plain,
    aes(x = S_ind, y = norm_vol, group = interaction(type, rep_id)),
    alpha = 0.1, linewidth = 0.5, color = "black"
  ) +
  geom_line(data = median_df,
    aes(x = S_ind, y = median_norm),
    alpha = 0.8, linewidth = 1, color = "black"
  ) +
  # scale_color_manual(values = c("blue", "red", "green", "purple")) +
  coord_cartesian(xlim = c(0,10.0), ylim = c(0,0.8)) +
  labs(x = "Ind. Power Cap", y = "Norm. Vol.") +
  # facet_wrap(~ type) +
  theme_pubr(base_size=12) +
  theme(legend.position = "none")

ggsave("figures/plain.pdf", p1, width = 6, height = 4.5, dpi = 450)
