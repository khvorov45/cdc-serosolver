library(tidyverse)

sim_data <- read_csv("sim/sim-data.csv", col_types = cols())

titre_plot <- sim_data %>%
  ggplot(aes(sample_year, titre_rescaled)) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_continuous(
    "Sample time",
    breaks = seq(2014, 2014.75, 0.25),
    labels = glue::glue("2014 Q{1:4}")
  ) +
  scale_y_log10("Titre", breaks = 5 * 2^(0:15)) +
  facet_wrap(~virus_year) +
  geom_point(alpha = 0.5) +
  geom_line(aes(group = individual), alpha = 0.5)

ggsave("sim/titres.pdf", titre_plot, width = 20, height = 20, units = "cm")
