#include "InverseCompton.h"
#include <omp.h>
#include "Eigen/Core"
#include "Utils.h"
#include "constants.h"

double CRad::InverseCompton::emissivity_(double epsilon, double electron_energy,
                                         double epsilon_prime) const {
    double lorentz = electron_energy / constants::m_e;
    double e1 = epsilon_prime / electron_energy;  // E_1 in Eq<2.47>
    double Gamma = 4 * epsilon * lorentz / constants::m_e;
    double q = e1 / (Gamma * (1. - e1));
    double f = 1. / (4 * lorentz * lorentz);
    if (f > 0.1)
        return 0.;
    if (q > 1 || q < f) {
        return 0;
    }
    double bracket = 2. * q * log(q) + (1. + 2. * q) * (1 - q) +
                     0.5 * (1. - q) * Gamma * q * Gamma * q / (1. + Gamma * q);
    double integrand = 2 * constants::pi * constants::eRadius *
                       constants::eRadius * constants::m_e *
                       constants::c_speed / lorentz / epsilon * bracket;
    integrand =
        integrand / lorentz / constants::m_e;  // Convert dE1 to d\epsilon_prim
    integrand = integrand * constants::ln10 * epsilon;  // d10^e = ln1010^e de
    /*
     * Be careful! here we don't consider the photon density n(\epsilon) 
     */
    return integrand;
}
double CRad::InverseCompton::InverseCompton::emissivity_photon_integrated_sum(
    double electron_energy, double epsilon_prime) const {
    double lorentz = electron_energy / constants::m_e;
    double E_1 = epsilon_prime / electron_energy;
    double up_bound = epsilon_prime / (1 - E_1);
    double low_bound = up_bound / (4 * lorentz * lorentz);

    auto energy_range = dat.GetMaxMinTargetEnergy();
    up_bound = up_bound > energy_range.first ? energy_range.first : up_bound;
    low_bound =
        low_bound < energy_range.second ? energy_range.second : low_bound;
    if (lorentz < 1) {
        return 0;
    }
    if (low_bound >= up_bound) {
        return 0;
    }
    /*
     * Integration on the log10(energy)
     */
    up_bound = log10(up_bound);
    low_bound = log10(low_bound);
    auto f = [&](double x) -> double {
        return dat.GetSumPhotonDensity(pow(10, x)) *
               emissivity_(pow(10, x), electron_energy, epsilon_prime);
    };
    double res = Utils::integrate_1d(f, low_bound, up_bound);
    return res;
}
double CRad::InverseCompton::differential_emission(double epsilon_prime) const {
    auto electron_energy_range = dat.GetMinMaxElectronEnergy();
    double min_energy = electron_energy_range.first;
    double max_energy = electron_energy_range.second;
    if (epsilon_prime > max_energy) {
        return 0;
    }
    if (epsilon_prime > min_energy) {
        min_energy = epsilon_prime;
    }
    auto f = [&](double x) -> double {
        return emissivity_photon_integrated_sum(pow(10, x), epsilon_prime) *
               dat.GetElectronDensity(pow(10, x)) * constants::ln10 *
               pow(10, x);
    };
    double res = Utils::integrate_1d(f, log10(min_energy), log10(max_energy));
    return res;
}

Eigen::VectorXd CRad::InverseCompton::CalculateDifferentialSpectrum(
    const Eigen::Ref<const Eigen::VectorXd>& epsilon_prime) const {
    Eigen::VectorXd differential_spectrum(epsilon_prime.size());
#pragma omp parallel for
    for (int i = 0; i < epsilon_prime.size(); i++) {
        differential_spectrum[i] = differential_emission(epsilon_prime[i]);
    }
    return differential_spectrum;
}
