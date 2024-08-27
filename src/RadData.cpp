#include "RadData.h"
#include <gsl/gsl_spline.h>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Eigen/Core"
#include "Utils.h"
#include "constants.h"
#include <iostream>

CRad::RadData::RadData() {
    electron_distribution = Eigen::MatrixX2d::Zero(0, 2);
    proton_distribution = Eigen::MatrixX2d::Zero(0, 2);
    B_field = 0.0;
}

CRad::RadData::~RadData() {
    target_photons_density_log_.clear();
    target_photons_energy_log_.clear();
    target_photons.clear();
}
void CRad::RadData::SetElectronDis(const Eigen::VectorXd& energy,
                                   const Eigen::VectorXd& density) {
    if (energy.size() != density.size()) {
        throw std::runtime_error(
            "The dimension of electron energy and electron distribution is not "
            "the Same!");
    }
    if (verbose) {
        // print some info;
    }
    // For negative energy, program exit
    if( (energy.array() <= 0).any() )
    {
        throw std::runtime_error("Have nagetive Energy! Program Exit");
    }
    if((density.array() <= 0).any())
    {
        std::cout << "Warning: the density have zero or negative numver !" << std::endl;
    }
    Eigen::VectorXd positive_density = density.array().max(1e-240);
    // In order to get rid of the situation, density have 0.0
    electron_distribution.resize(energy.size(), 2);
    electron_distribution.col(0) = energy;
    electron_distribution.col(1) = positive_density;

    electron_energy_log_ = energy.array().log10();
    electron_density_log_ = positive_density.array().log10();

    electron_interpolate = std::make_unique<GSLInterpData>(
        electron_energy_log_, electron_density_log_);
}
double CRad::RadData::GetElectronDensity(double energy) {
    double log_energy = std::log10(energy);
    if (log_energy > electron_energy_log_.maxCoeff() ||
        log_energy < electron_energy_log_.minCoeff()) {
        return 0;
    }
    double log10_res = gsl_spline_eval(electron_interpolate->spline, log_energy,
                                       electron_interpolate->acc);
    return pow(10, log10_res);
}
std::pair<double, double> CRad::RadData::GetMinMaxElectronEnergy() {
    double min_energy = electron_distribution.col(0).minCoeff();
    double max_energy = electron_distribution.col(0).maxCoeff();
    return std::make_pair(min_energy, max_energy);
}

void CRad::RadData::SetProtonDis(const Eigen::VectorXd& energy,
                                 const Eigen::VectorXd& density) {
    if (energy.size() != density.size()) {
        throw std::runtime_error("");
    }
    if (verbose) {
        //
    }
    if( (energy.array() <= 0).any() )
    {
        throw std::runtime_error("Have nagetive Energy! Program Exit");
    }
    if((density.array() <= 0).any())
    {
        std::cout << "Warning: the density have zero or negative numver !" << std::endl;
    }
    Eigen::VectorXd positive_density = density.array().max(1e-240);
    proton_distribution.resize(energy.size(), 2);
    proton_distribution.col(0) = energy;
    proton_distribution.col(1) = positive_density;
    proton_energy_log_ = energy.array().log10();
    proton_density_log_ = positive_density.array().log10();
    proton_interpolate  = std::make_unique<GSLInterpData>(proton_energy_log_, proton_density_log_);
}

//TODO: Exists a lof of repeated code!!! need fix
double CRad::RadData::GetProtonDensity(double energy)
{
    double log_energy = std::log10(energy);
    if( log_energy > proton_energy_log_.maxCoeff() ||
            log_energy < proton_energy_log_.minCoeff())
    {
        return 0;
    }
    double log10_res = proton_interpolate->Interpolate(log_energy);
    return pow(10, log10_res);
}

std::pair<double, double> CRad::RadData::GetMinMaxProtonEnergy()
{
    double min_energy = proton_distribution.col(0).minCoeff();
    double max_energy = proton_distribution.col(0).maxCoeff();
    return std::make_pair(min_energy, max_energy);
}

void CRad::RadData::AddTargetPhotons(const Eigen::VectorXd& energy,
                                     const Eigen::VectorXd& density) {
    Eigen::MatrixX2d photon_distribution(energy.size(), 2);
    Eigen::MatrixX2d photon_distribution_log(energy.size(), 2);
    photon_distribution.col(0) = energy;
    photon_distribution.col(1) = density;
    target_photons.push_back(std::move(photon_distribution));
}
void CRad::RadData::AddThermalTargetPhotons(const double temperture, int bins,
                                            double energy_density) {
    double low_boundary = constants::kb * temperture * 1e-12;
    double high_boundary = constants::kb * temperture * 1e4;

    auto energy = Utils::GetLogspaceVec(low_boundary, high_boundary, bins);

    Eigen::VectorXd density;
    if (energy_density > 0) {
        density = Utils::GreyBody(energy, temperture, energy_density);
    } else {
        density = Utils::BlackBody(energy, temperture);
    }
    Eigen::MatrixX2d thermal_photons(energy.size(), 2);
    thermal_photons.col(0) = energy;
    thermal_photons.col(1) = density;
    /* precompute the log10 version */
    Eigen::VectorXd energy_log10 = energy.array().log10();
    Eigen::VectorXd density_log10 = density.array().log10();

    target_photons_energy_log_.push_back(std::move(energy_log10));
    target_photons_density_log_.push_back(std::move(density_log10));
    target_photons.push_back(std::move(thermal_photons));

    target_photons_interpolate.push_back(std::make_unique<GSLInterpData>(
        target_photons_energy_log_.back(), target_photons_density_log_.back()));
}

Eigen::VectorXd CRad::RadData::GetPhotonDensity(double energy) {
    Eigen::VectorXd res(target_photons.size());
    double log_energy = std::log10(energy);
    for (auto index = 0; index < target_photons.size(); ++index) {
        if (log_energy > target_photons_energy_log_[index].maxCoeff() ||
            log_energy < target_photons_energy_log_[index].minCoeff()) {
            res[index] = 0;
        } else {
            res[index] = pow(
                10, target_photons_interpolate[index]->Interpolate(log_energy));
        }
    }
    return res;
}
double CRad::RadData::GetSumPhotonDensity(double energy) {
    auto res = GetPhotonDensity(energy);
    return res.sum();
}
std::pair<double, double> CRad::RadData::GetMaxMinTargetEnergy() {
    double max_v = 0;
    double min_v = 1e8;

    for (const auto& photon_dis : target_photons) {
        if (photon_dis.row(0).maxCoeff() > max_v) {
            max_v = photon_dis.row(0).maxCoeff();
        }
        if (photon_dis.row(0).minCoeff() < min_v) {
            min_v = photon_dis.row(0).minCoeff();
        }
    }
    return std::make_pair(max_v, min_v);
}
