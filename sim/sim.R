library(tidyverse)

# SECTION Antigenic map

virus_key <- c(
  "HK68" = 1968, "EN72" = 1972, "VI75" = 1975, "TX77" = 1977,
  "BK79" = 1979, "SI87" = 1987, "BE89" = 1989, "BJ89" = 1989,
  "BE92" = 1992, "WU95" = 1995, "SY97" = 1997, "FU02" = 2002,
  "CA04" = 2004, "WI05" = 2005, "PE06" = 2006
)

antigenic_map <- read_csv(
  system.file("extdata", "fonville_map_approx.csv", package = "serosolver"),
  col_types = cols()
) %>%
  select(strain_name = Strain, x_coord = X, y_coord = Y) %>%
  mutate(
    strain_year = virus_key[strain_name],
    strain_quarter = strain_year * 4,
  )

# NOTE(sen) The map does not have all the desired quaters to they need to be
# predicted somehow
strain_year_min <- 2000
strain_year_max <- 2002
strain_quarters_desired <- seq(strain_year_min * 4, strain_year_max * 4 - 1)

# NOTE(sen) This dodgy manipulation is from
# serosolver::generate_antigenic_map_flexible
antigenic_map_xy_fit <-
  smooth.spline(antigenic_map$x_coord, antigenic_map$y_coord, spar = 0.3)
antigenic_map_x_quarter_fit <- lm(x_coord ~ strain_quarter, antigenic_map)
antigenic_map_predicted <- tibble(
  strain_quarter = strain_quarters_desired,
  x_coord = predict(antigenic_map_x_quarter_fit, tibble(strain_quarter)),
  y_coord = predict(antigenic_map_xy_fit, tibble(x_coord = x_coord))$y$x_coord
)

# NOTE(sen) Verify that it's the same as serosolver's
antigenic_map_predicted_serosolver <-
  serosolver::generate_antigenic_map_flexible(
    antigenic_map %>%
      select(Strain = strain_year, X = x_coord, Y = y_coord),
    buckets = 4,
    year_min = strain_year_min,
    year_max = strain_year_max
  )
all(antigenic_map_predicted_serosolver$x_coord == antigenic_map_predicted$x_coord)
all(antigenic_map_predicted_serosolver$y_coord == antigenic_map_predicted$y_coord)

# NOTE(sen) The manipulation above happens to work well enough for _this_ map
# but it's not a robust way of predicting titres for _any_ map.

write_csv(antigenic_map_predicted, "sim/sim-agmap.csv")

antigenic_map_coords_plot <- antigenic_map %>%
  ggplot(aes(x_coord, y_coord)) +
  theme_bw() +
  geom_point() +
  geom_text(aes(label = strain_name), hjust = 0.5, vjust = 1, nudge_y = -0.2) +
  geom_line(data = antigenic_map_predicted)

antigenic_map_x_time_plot <- antigenic_map %>%
  select(strain_quarter, x_coord) %>%
  mutate(source = "real") %>%
  bind_rows(
    antigenic_map_predicted %>%
      select(strain_quarter, x_coord) %>%
      mutate(source = "predicted-serosolver")
  ) %>%
  ggplot(aes(strain_quarter, x_coord, col = source)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null")
  ) +
  geom_point()

ggdark::ggsave_dark(
  "sim/antigenic-map.pdf",
  ggpubr::ggarrange(
    antigenic_map_coords_plot, antigenic_map_x_time_plot,
    ncol = 1
  ),
  width = 10, height = 15, units = "cm"
)

# SECTION Simulation

n_individuals <- 50

sim_data_ages <- tibble(
  pid = 1:n_individuals,
  dob = runif(
    n_individuals,
    lubridate::ymd("1960-01-01"),
    lubridate::ymd("1999-12-31")
  ) %>%
    lubridate::as_date(),
  dob_quarter = lubridate::year(dob) * 4 + ceiling(lubridate::month(dob) / 3)
)

infection_histories <- antigenic_map_predicted %>%
  mutate(infection_prob = runif(nrow(antigenic_map_predicted), 0, 1)) %>%
  select(strain_quarter, infection_prob) %>%
  slice(rep(1:n(), each = n_individuals)) %>%
  bind_cols(
    sim_data_ages %>%
      select(pid, dob_quarter) %>%
      slice(rep(1:n(), times = nrow(antigenic_map_predicted)))
  ) %>%
  mutate(
    infection_prob = if_else(dob_quarter <= strain_quarter, infection_prob, 0),
    infected = rbinom(n(), 1, infection_prob)
  )

# NOTE(sen) Nobody should be infected with strains that circulated before their
# birth
infection_histories %>%
  filter(dob_quarter > strain_quarter) %>%
  pull(infected) %>%
  `==`(0) %>%
  all()

# NOTE(sen) Serosolver's `melt_antigenic_coords` actually just finds the
# distance matrix but in its own very special way

antigenic_coords_matrix <- antigenic_map_predicted %>%
  select(x_coord, y_coord) %>%
  as.matrix()
rownames(antigenic_coords_matrix) <- antigenic_map_predicted$strain_quarter

antigenic_distances <- dist(antigenic_coords_matrix) %>% as.matrix()

# NOTE(sen) Serosolver thought that a c++ function was warranted for this
# amazingly complex job
calc_boosting <- function(distances, parameter) {
  pmax(1 - distances * parameter, 0)
}

parameter_long_term_boosting <- 0.1 # NOTE(sen) sigma1
parameter_short_term_boosting <- 0.03 # NOTE(sen) sigma2

long_term_boosting <-
  calc_boosting(antigenic_distances, parameter_long_term_boosting)
short_term_boosting <-
  calc_boosting(antigenic_distances, parameter_short_term_boosting)

# NOTE(sen) So the "boosting" is as follows. If the antigenic distance is more
# than 1 / paramenter then it's 0. Otherwise it's 1 minus distance times the
# appropriate parameter So when the parameter is 0 then the boosting is actually
# 1 for every virus. If the parameter is large (infinite) then the boosting is 0
# for every virus

param_titre_contribution_long <- 1.8 # NOTE(sen) mu
param_titre_contribution_short <- 2.7 # NOTE(sen) mu_short

# NOTE(sen) Each column represents the titre contributions for that strain from
# all the other strains (and itself) if the individual was infected with all the
# strains and the sample was taken after all the infections and waning isn't a
# factor
titre_contribution_long <- long_term_boosting * param_titre_contribution_long
titre_contribution_short <- short_term_boosting * param_titre_contribution_short

sampling_quarters <- strain_quarters_desired[strain_quarters_desired %% 2 == 0]

parameter_wane_per_quarter <- 0.2 # NOTE(sen) wane
parameter_seniority <- 0.05 # NOTE(sen) tau

clamp01 <- function(x) {
  pmax(pmin(x, 1), 0)
}

sim_titres <- tibble(
  pid = unique(infection_histories$pid) %>% rep(each = length(sampling_quarters)),
  sampling_quarter = rep(sampling_quarters, length(unique(pid))),
) %>%
  slice(rep(1:n(), each = length(strain_quarters_desired))) %>%
  mutate(
    measured_strain_quarter = rep(strain_quarters_desired, length.out = n()),
  ) %>%
  inner_join(infection_histories, c("pid", "measured_strain_quarter" = "strain_quarter")) %>%
  mutate(infected_at_sampling_time = if_else(measured_strain_quarter <= sampling_quarter, infected, 0L)) %>%
  group_by(pid, sampling_quarter) %>%
  mutate(
    quarters_from_infection_to_measurement = sampling_quarter - measured_strain_quarter,
    wane = clamp01(1 - parameter_wane_per_quarter * quarters_from_infection_to_measurement),
    # NOTE(sen) titre_contribution matrices need to correspond to
    # measured_strain_quarter in terms of where each strain is positioned
    max_titre_contribution_long = map_dbl(
      measured_strain_quarter,
      ~ sum(titre_contribution_long[as.character(.x), ] * infected_at_sampling_time)
    ),
    max_titre_contribution_short = map_dbl(
      measured_strain_quarter,
      ~ sum(titre_contribution_short[as.character(.x), ] * infected_at_sampling_time * wane)
    ),
    titre_index =
      (max_titre_contribution_long + max_titre_contribution_short) *
        clamp01(1 - parameter_seniority * sum(infected_at_sampling_time))
  ) %>%
  ungroup() %>%
  mutate(
    titre = 5 * 2^titre_index
  )

titre_plot <- sim_titres %>%
  ggplot(aes(sampling_quarter, titre)) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_continuous(
    "Sample time",
    # breaks = seq(2014, 2014.75, 0.25),
    # labels = glue::glue("2014 Q{1:4}")
  ) +
  scale_y_log10("Titre", breaks = 5 * 2^(0:15)) +
  facet_wrap(~measured_strain_quarter) +
  geom_point(alpha = 0.5) +
  geom_line(aes(group = pid), alpha = 0.5)

ggsave("sim/titres.pdf", titre_plot, width = 20, height = 20, units = "cm")
