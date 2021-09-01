library(tidyverse)

fit_par_tab <- read_csv("sim/fit-par-tab.csv", col_types = cols())
sim_data <- read_csv("sim/sim-data.csv", col_types = cols())
sim_agmap <- read_csv("sim/sim-agmap.csv", col_types = cols())

# NOTE(sen) Amateur hour programming continues
library(coda)
library(data.table)
library(serosolver)
all_chains <- load_mcmc_chains(
  "sim",
  thin = 1,
  burnin = 0,
  par_tab = fit_par_tab,
  unfixed = FALSE,
  convert_mcmc = TRUE
)
detach("package:coda")
detach("package:data.table")

parameters <- all_chains$theta_chain %>% as.data.frame()
infection_histories <- all_chains$inf_chain %>% as.data.frame()

library(plyr)
library(data.table)
titre_predictions <- serosolver::get_titre_predictions(
  chain = parameters,
  infection_histories = infection_histories,
  titre_dat = sim_data %>% as.data.frame(),
  individuals = unique(sim_data$individual),
  antigenic_map = sim_agmap %>% as.data.frame(),
  par_tab = fit_par_tab %>% as.data.frame(),
  for_res_plot = FALSE,
  titre_before_infection = FALSE
)
detach("package:plyr")
detach("package:data.table")

titre_predictions$predictions %>%
  as_tibble() %>%
  print(n = 100) # %>%
filter(individual == 1, virus == first(virus))

get_titre_predictions <- function(chain = parameters,
                                  infection_histories,
                                  titre_dat = sim_data %>% as.data.frame(),
                                  individuals = unique(sim_data$individual),
                                  antigenic_map = sim_agmap %>% as.data.frame(),
                                  strain_isolation_times = NULL,
                                  par_tab = fit_par_tab %>% as.data.frame(),
                                  nsamp = 10,
                                  add_residuals = FALSE,
                                  mu_indices = NULL,
                                  measurement_indices_by_time = NULL,
                                  for_res_plot = FALSE,
                                  expand_titredat = FALSE,
                                  titre_before_infection = FALSE, titres_for_regression = FALSE) {
  ## Need to align the iterations of the two MCMC chains
  ## and choose some random samples
  samps <- intersect(unique(infection_histories$sampno), unique(chain$sampno))
  chain <- chain[chain$sampno %in% samps, ]
  infection_histories <- infection_histories[infection_histories$sampno %in% samps, ]

  ## Take subset of individuals
  titre_dat <- titre_dat[titre_dat$individual %in% individuals, ]
  infection_histories <- infection_histories[infection_histories$i %in% individuals, ]

  titre_dat$individual <- match(titre_dat$individual, individuals)
  infection_histories$i <- match(infection_histories$i, individuals)

  ## Format the antigenic map to solve the model
  if (!is.null(antigenic_map)) {
    strain_isolation_times <- unique(antigenic_map$inf_times) # How many strains are we testing against and what time did they circulate
  } else {
    antigenic_map <- data.frame("x_coord" = 1, "y_coord" = 1, "inf_times" = strain_isolation_times)
  }
  nstrain <- length(strain_isolation_times)
  n_indiv <- length(individuals)

  ## Empty data structures to save output to
  infection_history_dens <- NULL
  tmp_samp <- sample(samps, nsamp)

  ## See the function in posteriors.R
  titre_dat1 <- titre_dat

  if (expand_titredat) {
    titre_dat1 <- expand.grid(
      individual = unique(titre_dat$individual),
      samples = unique(titre_dat$samples),
      titre = 0, run = 1
    )
    titre_dat2 <- unique(titre_dat[, c("individual", "virus", "group", "DOB")])
    titre_dat1 <- merge(titre_dat1, titre_dat2)
    titre_dat1 <- titre_dat1[
      order(titre_dat1$group, titre_dat1$individual, titre_dat1$samples, titre_dat1$virus),
      c("individual", "samples", "virus", "titre", "run", "group", "DOB")
    ]
  }
  model_func <- serosolver::create_posterior_func(par_tab, titre_dat1, antigenic_map, 100,
    mu_indices = mu_indices, version = 2,
    measurement_indices_by_time = measurement_indices_by_time, function_type = 4,
    titre_before_infection = titre_before_infection
  )

  predicted_titres <- residuals <- residuals_floor <-
    observed_predicted_titres <- matrix(nrow = nrow(titre_dat1), ncol = nsamp)
  samp_record <- numeric(nsamp)


  ## For each sample, take values for theta and infection histories and simulate titres
  inf_hist_all <- list(nsamp)
  for (i in 1:nsamp) {
    index <- tmp_samp[i]
    pars <- get_index_pars(chain, index)
    pars <- pars[!(names(pars) %in% c(
      "lnlike", "likelihood", "prior_prob",
      "sampno", "total_infections", "chain_no"
    ))]
    ## pars <- pars[names(pars) %in% par_tab$names]
    tmp_inf_hist <- infection_histories[infection_histories$sampno == index, ]
    tmp_inf_hist <- as.matrix(Matrix::sparseMatrix(i = tmp_inf_hist$i, j = tmp_inf_hist$j, x = tmp_inf_hist$x, dims = c(n_indiv, nstrain)))
    predicted_titres[, i] <- model_func(pars, tmp_inf_hist)
    observed_predicted_titres[, i] <- add_noise(predicted_titres[, i], pars, NULL, NULL)
    inf_hist_all[[i]] <- tmp_inf_hist
    ## Get residuals between observations and predictions
    residuals[, i] <- titre_dat1$titre - floor(predicted_titres[, i])
    residuals_floor[, i] <- titre_dat1$titre - observed_predicted_titres[, i]
    samp_record[i] <- index
  }
  colnames(predicted_titres) <- tmp_samp

  ## If generating for residual plot, can return now
  if (for_res_plot) {
    return(list(
      residuals, samp_record, titre_dat1,
      predicted_titres,
      observed_predicted_titres,
      residuals_floor
    ))
  }

  # residuals <- cbind(titre_dat1, residuals)

  ## Get 95% credible interval and means
  dat2 <- t(apply(predicted_titres, 1, function(x) quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975))))

  ## Get 95% credible interval and means of observations
  obs_dat <- t(apply(observed_predicted_titres, 1, function(x) quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975))))

  residuals <- t(apply(residuals, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))
  residuals <- cbind(titre_dat1, residuals)

  ## Find multivariate posterior mode estimate from the chain
  best_pars <- get_best_pars(chain)
  best_pars <- best_pars[!(names(best_pars) %in% c(
    "lnlike", "likelihood", "prior_prob",
    "sampno", "total_infections", "chain_no"
  ))]
  # best_pars <- best_pars[names(best_pars) %in% par_tab$names]
  best_I <- chain$sampno[which.max(chain$lnlike)]
  best_inf <- infection_histories[infection_histories$sampno == best_I, ]
  best_inf <- as.matrix(Matrix::sparseMatrix(i = best_inf$i, j = best_inf$j, x = best_inf$x, dims = c(n_indiv, nstrain)))

  ## Generate trajectory for best parameters
  best_traj <- model_func(best_pars, best_inf)
  best_residuals <- titre_dat1$titre - floor(best_traj)
  best_residuals <- cbind(titre_dat1, best_residuals, "sampno" = best_I)
  dat2 <- as.data.frame(dat2)
  obs_dat <- as.data.frame(obs_dat)

  colnames(dat2) <- colnames(obs_dat) <- c("lower", "lower_50", "median", "upper_50", "upper")
  dat2$max <- best_traj
  dat2 <- cbind(titre_dat1, dat2)
  obs_dat <- cbind(titre_dat1, obs_dat)
  tmp_inf_chain <- data.table(subset(infection_histories, sampno %in% tmp_samp))

  ## Get infection history density for each individual and each epoch
  data.table::setkey(tmp_inf_chain, "i", "j")
  infection_history_dens <- tmp_inf_chain[, list(V1 = sum(x) / length(tmp_samp)), by = key(tmp_inf_chain)]
  infection_history_dens$j <- strain_isolation_times[infection_history_dens$j]
  colnames(infection_history_dens) <- c("individual", "variable", "value")
  infection_history_final <- infection_history_dens
  best_inf <- data.frame(best_inf)
  best_inf$individual <- 1:nrow(best_inf)
  best_inf$individual <- individuals[best_inf$individual]

  dat2$individual <- individuals[dat2$individual]
  infection_history_final$individual <- individuals[infection_history_final$individual]
  if (titres_for_regression) {
    return(list(
      "all_predictions" = predicted_titres, "all_inf_hist" = inf_hist_all,
      "summary_titres" = dat2, "best_inf_hist" = best_inf, "predicted_observations" = obs_dat
    ))
  }

  if (add_residuals) {
    result <- list(
      "predictions" = dat2, "histories" = infection_history_final,
      "residuals" = residuals, "bestRes" = best_residuals, "best_infhist" = best_inf,
      "predicted_observations" = obs_dat
    )
  } else {
    result <- list(
      "predictions" = dat2, "histories" = infection_history_final,
      "best_infhist" = best_inf, "predicted_observations" = obs_dat
    )
  }
  return(result)
}
