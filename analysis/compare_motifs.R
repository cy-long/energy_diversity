library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)
library(purrr)

df <- read_csv("analysis/trophic_data_quadratic.csv")
df$type <- df$tr_type %>% gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", ., perl=TRUE) %>% as.factor()
df$norm_vol[is.nan(df$norm_vol)] <- 0
df <- df %>% rename(S_tot = S) # temporary

median_df <- df %>%
  group_by(S_tot, type) %>%
  summarise(median_norm = median(norm_vol), .groups = "drop")

p <- ggplot() +
  geom_line(data = df,
    aes(x = S_tot, y = vol_d, group = interaction(type, rep_id), color = type),
    alpha = 0.1, linewidth = 0.5
  ) +
  # geom_line(data = median_df,
  #   aes(x = S_tot, y = 10*median_norm, color = type),
  #   alpha = 0.8, linewidth = 1
  # ) +
  scale_color_manual(values = c("blue", "red", "green", "purple")) +
  coord_cartesian(xlim = c(0,5.0), ylim = c(0,5.0)) +
  labs(x = "Total Power Cap", y = "Rescaled Vol.") +
  facet_wrap(~ type) +
  theme_pubr(base_size=12) +
  theme(legend.position = "none")

print(p)

ggsave("figures/compare_2.pdf", p, width = 6, height = 4.5, dpi = 450)


df_plain <- read_csv("analysis/plain_data.csv")
df_plain$type <- as.factor(df_plain$type)

df_plain <- df_plain %>%
  mutate(norm_vol = vol / vol_simp(S_tot)) %>%
  group_by(S_tot, type) %>%
  mutate(rep_id = row_number()) %>%
  ungroup()

median_df <- df_plain %>%
  group_by(S_tot, type) %>%
  summarise(median_norm = median(norm_vol), .groups = "drop")

p1 <- ggplot() +
  geom_line(data = df_plain,
    aes(x = S_tot, y = norm_vol, group = interaction(type, rep_id)),
    alpha = 0.1, linewidth = 0.5, color = "black"
  ) +
  geom_line(data = median_df,
    aes(x = S_tot, y = median_norm),
    alpha = 0.8, linewidth = 1, color = "black"
  ) +
  # scale_color_manual(values = c("blue", "red", "green", "purple")) +
  coord_cartesian(xlim = c(0,10.0), ylim = c(0,0.8)) +
  labs(x = "Ind. Power Cap", y = "Norm. Vol.") +
  # facet_wrap(~ type) +
  theme_pubr(base_size=12) +
  theme(legend.position = "none")

ggsave("figures/plain.pdf", p1, width = 6, height = 4.5, dpi = 450)
