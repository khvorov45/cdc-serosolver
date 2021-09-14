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

# NOTE(sen) So the "boosting" is this as follows. If the antigenic distance is
# more than 1 / paramenter then it's 0. Otherwise it's 1 minus distance times
# the appropriate parameter So when the parameter is 0 then the boosting is
# actually 1 for every virus. If the parameter is large (infinite) then the
# boosting is 0 for every virus

param_titre_contribution_long <- 1.8 # NOTE(sen) mu
param_titre_contribution_short <- 2.7 # NOTE(sen) mu_short

# NOTE(sen) Each column represents the titre contributions for that strain from
# all the other strains (including itself) if the individual was infected with
# all the strains and the sample was taken after all the infections and waning
# isn't a factor
titre_contribution_short <- short_term_boosting * param_titre_contribution_short
titre_contribution_long <- long_term_boosting * param_titre_contribution_long

sampling_quarters <- strain_quarters_desired[strain_quarters_desired %% 2 == 0]

parameter_wane_per_quarter <- 0.2 # NOTE(sen) wane
parameter_seniority <- 0.05 # NOTE(sen) tau

Rcpp::sourceCpp("sim/simtitre.cpp")

simulate_individual_titre_cpp(
  strain_quarters_desired,
  titre_contribution_long["8000", ],
  titre_contribution_short["8000", ],
  infection_histories %>% filter(pid == 1) %>% pull(infected),
  8001,
  parameter_wane_per_quarter,
  parameter_seniority
)

simulate_individual_titre_multiple_timepoints_cpp(
  strain_quarters_desired,
  titre_contribution_long["8000", ],
  titre_contribution_short["8000", ],
  infection_histories %>% filter(pid == 1) %>% pull(infected),
  c(8000, 8001),
  parameter_wane_per_quarter,
  parameter_seniority
)

titre_contribution_long["8000", "8000"]

sim_titres <- map_dfr(sampling_quarters, function(sampling_quarter) {
  infection_histories_before_sampling <- infection_histories %>%
    mutate(infected = if_else(strain_quarter <= sampling_quarter, infected, 0L))
  infection_histories_before_sampling %>%
    group_by(pid) %>%
    group_map(function(data, key) {
      strain_quarters <- data$strain_quarter
      infected <- data$infected

      infected_matrix <- as.matrix(infected)
      rownames(infected_matrix) <- strain_quarters

      horizontal_matrix_with_ones <-
        matrix(rep(1, length(strain_quarters)), nrow = 1)
      colnames(horizontal_matrix_with_ones) <- strain_quarters

      # NOTE(sen) Remove the contribution of the strains that the individual
      # wasn't infected by
      titre_contribution_mask <- infected_matrix %*% horizontal_matrix_with_ones

      titre_contribution_long_masked <-
        titre_contribution_long * titre_contribution_mask
      titre_contribution_short_masked <-
        titre_contribution_short * titre_contribution_mask

      sampling_quarters_from_strain <-
        sampling_quarter - matrix(strain_quarters) %*% horizontal_matrix_with_ones

      wane_amount <-
        pmin(pmax(1 - sampling_quarters_from_strain * parameter_wane_per_quarter, 0), 1)

      titre_contribution_short_with_waning <-
        titre_contribution_short_masked * wane_amount

      titre_contribution_total_masked <-
        titre_contribution_long_masked + titre_contribution_short_masked

      seniority <- max(1 - parameter_seniority * (sum(infected) - 1), 0)

      titre_contribution_total_per_strain <-
        colSums(titre_contribution_total_masked) * seniority

      tibble(
        pid = key$pid,
        sampling_quarter,
        strain_quarter = names(titre_contribution_total_per_strain),
        titre_unit = titre_contribution_total_per_strain
      )
    })
})


sim_titres


infection_histories %>%
  filter(pid == 34)

sim_titres %>%
  filter(pid == 34, sampling_quarter == 8000)

matrix(-5:4, ncol = 2) %>% pmax(0)

sim_titres















temp_one_infection_history <- infection_histories %>%
  filter(pid == first(pid))

temp_one_infection_history_matrix <- temp_one_infection_history %>%
  pull(infected) %>%
  as.matrix()

rownames(temp_one_infection_history_matrix) <-
  temp_one_infection_history$strain_quarter

temp_one_infection_history_matrix %*%
  matrix(rep(1, nrow(temp_one_infection_history_matrix)), nrow = 1)

param_mu <- 1.8
param_mu_short <- 2.7
short_term_boosting %>%
  as.matrix() %>%
  `*`(param_mu_short) %>%
  colSums()


param_mu





c(temp)
str(antigenic_distances)

c(antigenic_distances)
c(temp)
antigenic_distances[lower.tri(antigenic_distances, diag = TRUE)]

data(example_par_tab, package = "serosolver")




# NOTE(sen) strain isolation times are indices of Q1 of virus years, assume that
# these sampled viruses were sampled in Q2 of each year.
sampled_viruses <- 1999 * 4
# seq(min(strain_isolation_times) + 1, max(strain_isolation_times), by = 4)

serum_samples_per_year <- 2
resolution_quarterly <- 4 # NOTE(sen) 4 per year, i.e. quarterly resolution_quarterly

# NOTE(sen) One year's worth of sampling times
sampling_times <- strain_isolation_times[c(1, 3)] # seq(2014 * resolution_quarterly, 2015 * resolution_quarterly - 1, by = 1)


length(attack_rates)
length(strain_isolation_times)
serosolver::check_attack_rates(attack_rates, strain_isolation_times)

all_simulated_data <- serosolver::simulate_data(
  par_tab = example_par_tab,
  n_indiv = 30,
  buckets = resolution_quarterly,
  strain_isolation_times = strain_isolation_times,
  measured_strains = sampled_viruses,
  sampling_times = sampling_times,
  nsamps = serum_samples_per_year,
  antigenic_map = example_antigenic_map,
  age_min = 20 * resolution_quarterly,
  age_max = 65 * resolution_quarterly,
  attack_rates = attack_rates,
  repeats = 1
)

# NOTE(sen) Unique values from each column
iwalk(
  all_simulated_data$data,
  function(vec, name) {
    cat(glue::glue("{name}: {paste(unique(vec), collapse = ' ')}\n\n"))
  }
)

sim_data <- all_simulated_data$data %>%
  inner_join(all_simulated_data$ages, "individual") %>%
  as_tibble() %>%
  mutate(
    sample_year = samples / 4,
    virus_year = (virus - 1) / 4,
    titre_rescaled = 5 * 2^(titre), # NOTE(sen) I guess
    year_of_birth = DOB / 4
  )

write_csv(sim_data, "sim/sim-data.csv")
