#include "RadData.h"
#include <stdexcept>
#include <utility>
#include <vector>
#include "Eigen/Core"
#include "Utils.h"
#include "constants.h"

CRad::RadData::RadData() {
    electron_distribution = Eigen::MatrixX2d::Zero(0, 2);
    electron_distribution_log = Eigen::MatrixX2d::Zero(0, 2);
    proton_distribution = Eigen::MatrixX2d::Zero(0, 2);
    B_field = 0.0;
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
    electron_distribution.resize(energy.size(), 2);
    electron_distribution.col(0) = energy;
    electron_distribution.col(1) = density;

    electron_distribution_log.resize(energy.size(), 2);
    electron_distribution_log.col(0) = energy.array().log10();
    electron_distribution_log.col(1) = density.array().log10();
}
double CRad::RadData::GetElectronDensity(double energy)
{
    return pow(10,Utils::interpolate(electron_distribution_log.col(0), electron_distribution_log.col(1), std::log10(energy)));
}
    
void CRad::RadData::SetProtonDis(const Eigen::VectorXd& energy,
                                 const Eigen::VectorXd& density) {
    if (energy.size() != density.size()) {
        throw std::runtime_error("");
    }
    if (verbose) {
        //
    }
    proton_distribution.resize(energy.size(), 2);
    proton_distribution.col(1) = energy;
    proton_distribution.col(2) = density;
}
void CRad::RadData::AddThermalTargetPhotons(const double temperture, int bins,
                                            double energy_density) {
    double low_boundary = constants::kb * temperture * 1e-12;
    double high_boundary = constants::kb * temperture * 1e6;

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

    target_photons.push_back(std::move(thermal_photons));
}

Eigen::VectorXd CRad::RadData::GetPhotonDensity(double energy) {
    Eigen::VectorXd res(target_photons.size());
    int index = 0;
    for (const auto& photon_dis : target_photons) {
        Eigen::VectorXd x = photon_dis.col(0).array().log10();
        Eigen::VectorXd y = photon_dis.col(1).array().log10();
        res[index++] = pow(10, Utils::interpolate(x, y, std::log10(energy)));
    }
    return res;
}
double CRad::RadData::GetSumPhotonDensity(double energy)
{
    auto res = GetPhotonDensity(energy);
    return res.sum();
}
std::pair<double, double> CRad::RadData::GetMinMaxTargetEnergy() {
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
