#include "Rcpp.h"

using namespace Rcpp;

// [[Rcpp::export]]
double simulate_individual_titre_cpp(
    const IntegerVector strain_quarters,
    const NumericVector titre_contribution_long,
    const NumericVector titre_contribution_short,
    const LogicalVector infections_in_quarters,
    const int measurement_quarter,
    const double parameter_wane_per_quarter,
    const double parameter_seniority
) {
    double titre = 0;
    int n_prior_infections = 0;
    for (int strain_index = 0;
        strain_index < strain_quarters.length();
        ++strain_index) {
        if (infections_in_quarters[strain_index]) {
            const int infection_quarter = strain_quarters[strain_index];
            const int quarters_since_infection =
                measurement_quarter - infection_quarter;
            if (quarters_since_infection >= 0) {
                double infection_titre_contribution = 0;

                infection_titre_contribution +=
                    titre_contribution_long[strain_index];

                double wane =
                    1.0 - parameter_wane_per_quarter * quarters_since_infection;
                if (wane < 0) {
                    wane = 0;
                }

                infection_titre_contribution +=
                    titre_contribution_short[strain_index] * wane;

                double seniority =
                    1.0 - parameter_seniority * n_prior_infections;
                if (seniority < 0) {
                    seniority = 0;
                }
                infection_titre_contribution *= seniority;
                ++n_prior_infections;

                titre += infection_titre_contribution;
            }
        }
    }
    return titre;
}

// [[Rcpp::export]]
NumericVector simulate_individual_titre_multiple_timepoints_cpp(
    const IntegerVector strain_quarters,
    const NumericVector titre_contribution_long,
    const NumericVector titre_contribution_short,
    const LogicalVector infections_in_quarters,
    const IntegerVector measurement_quarters,
    const double parameter_wane_per_quarter,
    const double parameter_seniority
) {
    NumericVector titres(measurement_quarters.length());
    for (int index = 0; index < measurement_quarters.length(); ++index) {
        titres[index] = simulate_individual_titre_cpp(
            strain_quarters,
            titre_contribution_long,
            titre_contribution_short,
            infections_in_quarters,
            measurement_quarters[index],
            parameter_wane_per_quarter,
            parameter_seniority
        );
    }
    return titres;
}
