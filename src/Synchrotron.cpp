#include "Synchrotron.h"
#include "Eigen/Core"
#include "Utils.h"
#include "constants.h"
#include <gsl/gsl_sf_bessel.h>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace  CRad
{
    double Synchrontron:: K_(double nu, double x)
    {
        if( x <= 0. || x > 700)
        {   
            return 0;
        }
        else {
            return gsl_sf_bessel_Knu(nu, x);
        }

    }
    double Synchrontron::Nu_b_(double B_field)
    {
        return  constants::el_charge * B_field * constants::c_speed/ (2 * constants::pi * constants::m_e);
    }
    double Synchrontron::Nu_c_(double B_field, double gamma)
    {
        return 3.0/2 * Nu_b_(B_field) * pow(gamma, 2);
    }
    double Synchrontron::R_(double x)
    {
        double term_1_num = 1.808 * std::pow(x, 1.0/3);
        double term_1_denom = std::sqrt(1 + 3.4 * std::pow(x, 2.0/3));
        double term_2_num = 1 + 2.21 * std::pow(x, 2.0/3) + 0.347 * std::pow(x, 4.0/3);
        double term_2_denom = 1 + 1.353 * std::pow(x, 2.0/3) + 0.217 * std::pow(x, 4.0/3);
        return term_1_num / term_1_denom * term_2_num / term_2_denom * std::exp(-1 * x);
        
    }
    double Synchrontron::emissivity_(double B_field, double electron_energy, double epsilon)
    {
        double nu = epsilon / constants::hp;
        double gamma = electron_energy/constants::m_e;
        double nu_b  = Nu_b_(B_field);
        double x  = nu/Nu_c_(B_field, gamma) / 2;
       if( nu < nu_b)
        {
            return 0;
        }
        double K_1 = K_(1./3., x);
        double K_4 = K_(4./3., x);
        double first_term = 4 * std::sqrt(3) * constants::pi * constants::e_radius * constants::m_e_g * constants::c_speed * nu_b * pow(x, 2);
        double second_term = K_4  * K_1 - 0.6 * x * (K_4 + K_1) * (K_4 - K_1);
        double res = first_term * second_term/ constants::hp / (epsilon);
        return res;
    }

    double Synchrontron::emissivity_approx_(double B_field, double electron_energy, double epsilon)
    {
        double nu = epsilon / constants::hp;
        double gamma = electron_energy/constants::m_e;
        double x  = nu/Nu_c_(B_field, gamma);
        double prefactor = std::sqrt(3) * std::pow(constants::el_charge, 3) * B_field / (constants::m_e);
        return prefactor * R_(x) / constants::hp / epsilon;

    }
    double Synchrontron::differential_emission(double epsilon)const
    {
        auto electron_energy_range = dat.GetMinMaxElectronEnergy();
        double min_energy = electron_energy_range.first;
        double max_energy = electron_energy_range.second;
if( epsilon > max_energy)
        {
            return 0;
        }
        if( epsilon > min_energy)
        {
            min_energy = epsilon;
        }
        auto f = [&](double x)->double{
            return emissivity(dat.GetB_Field(), pow(10, x), epsilon) * constants::ln10 * pow(10, x) * dat.GetElectronDensity(pow(10, x));
        };
        double res = Utils::integrate_1d(f, std::log10(min_energy), std::log10(max_energy));
        return res;
    }
    Eigen::VectorXd Synchrontron::CalculateDifferentialSpectrum(const Eigen::Ref<const Eigen::VectorXd>& epsilon) const
    {
        Eigen::VectorXd res(epsilon.size());
        size_t index = 0;
//#pragma omp parallel for
        for(auto i = 0; i < epsilon.size(); i++)
        {
            res[i]  = differential_emission(epsilon[i]);
        }
        return res;
    }

    double ProtonSynchrotron::x_(double nu, double B_field, double gamma)
    {
        double nu_c = Synchrontron::Nu_c_(B_field, gamma) ;
        double x = nu/ nu_c * (constants::m_p)/constants::m_e;
        return x;
    }
    double ProtonSynchrotron::Nu_b_(double B_field)
    {
        return Synchrontron::Nu_b_(B_field) * constants::m_e / constants::m_p;
    }
    double ProtonSynchrotron::emissivity_(double B_field, double proton_energy, double epsilon)
    {
        double nu   = epsilon /  constants::hp;
        double gamma = proton_energy / constants::m_p;
        double nu_b = Nu_b_(B_field);
        double x    = x_(nu, B_field, gamma) * 0.5;
        if( nu < nu_b)
        {
            return 0;
        }
        double K_1 = Synchrontron::K_(1./3., x);
        double K_4 = Synchrontron::K_(4./3., x);
        double first_term = 4 * std::sqrt(3) * constants::pi * constants::e_radius * constants::m_e_g * constants::c_speed * nu_b * pow(x, 2);
        double second_term = K_4 * K_1 - 0.6 * x * (K_4 + K_1) * (K_4 - K_1);
        // Convert the dP/dnu to dN/dE
        return first_term * second_term / constants::hp / epsilon;
    }
    double ProtonSynchrotron::emissivity_approx_(double B_field, double proton_energy, double epsilon)
    {
        double nu = epsilon / constants::hp;
        double gamma = proton_energy / constants::m_p;
        double x = x_(nu,B_field, gamma);
        double prefactor = std::sqrt(3) * std::pow(constants::el_charge, 3) / (constants::m_p);
        return prefactor * Synchrontron::R_(x) / constants::hp / epsilon;
    }
    double ProtonSynchrotron::differential_emission(double epsilon) const
    {
        auto proton_energy_range = dat.GetMinMaxProtonEnergy();
        double min_energy = proton_energy_range.first;
        double max_energy = proton_energy_range.second;
        if( epsilon > max_energy)
        {
            return 0;
        }
        if( epsilon > min_energy)
        {
            min_energy = epsilon;
        }

        auto f = [&](double x)->double{
            return emissivity(dat.GetB_Field(), pow(10, x), epsilon) * constants::ln10 * pow(10, x)
                * dat.GetProtonDensity(pow(10, x));
        };
        double res = Utils::integrate_1d(f, std::log10(min_energy), std::log10(max_energy));
        return res;
    }
    double Synchrontron::energy_loss_(double energy, double B_field)
    {
        double gamma = energy/constants::m_e;
        double beta2 = 1 - 1/std::pow(gamma, 2);
        double U_b   = B_field * B_field / 8 * constants::pi;
        return 4.0/3 * constants::c_speed * constants::sigma_T * U_b * beta2 * std::pow(gamma, 2);
    }
    Eigen::VectorXd ProtonSynchrotron::CalculateDifferentialSpectrum(const Eigen::Ref<const Eigen::VectorXd>& epsilon) const
    {
        Eigen::VectorXd res(epsilon.size());
#pragma  omp parallel for
        for(auto i = 0; i < epsilon.size(); i++)
        {
            res[i] = differential_emission(epsilon[i]);
        }
        return res;

    }

}
