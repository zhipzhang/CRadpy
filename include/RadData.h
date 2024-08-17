/*
 *
 * Store the Initial Data from the radiation computation
 *
 */
#pragma once
#include <memory>
#include <vector>
#include "Eigen/Core"
#include "Utils.h"

namespace CRad {
class RadData {

   public:
    RadData();
    ~RadData();
    bool HaveTargetPhotons() {
        if (target_photons.size() > 0) {
            return true;
        }
        return false;
    }

    bool HaveElectron() {
        if (electron_distribution.rows() > 0) {
            return true;
        }
        return false;
    }

    /* Set the spectra of electron:  array of energy[erg], array of number density[erg^-1]*/
    void SetElectronDis(const Eigen::VectorXd& energy,
                        const Eigen::VectorXd& density);

    /* return the max/min energy of electron spectra */
    std::pair<double, double> GetMinMaxElectronEnergy();

    double GetElectronDensity(double energy);

    /* Set the spectra of proton:  array of energy[erg], array of number density[erg^-1]*/
    void SetProtonDis(const Eigen::VectorXd& energy,
                      const Eigen::VectorXd& density);

    /* Add arbitary isotropic target photon distribution*: array of energy[erg] array of density */
    void AddTargetPhotons(const Eigen::VectorXd& energy,
                          const Eigen::VectorXd& density);

    /* Add blackbody/Graybody target photons */
    void AddThermalTargetPhotons(const double temerature, int bins = 500,
                                 double energy_density = 0);

    std::pair<double, double> GetMaxMinTargetEnergy();

    /* return a VectorXd of multiple components photon density  energy[erg]*/
    Eigen::VectorXd GetPhotonDensity(double energy);

    /* target photons for InverseCompton. return the sum of multiple components unit:[erg, erg^-1]*/
    double GetSumPhotonDensity(double energy);

    std::vector<Eigen::MatrixX2d> target_photons;

    /* magnetic field  unit:[gauss]*/
    double B_field;

    /* electron distribution unit:[erg, erg^-1]*/
    Eigen::MatrixX2d electron_distribution;

    /* proton distribution unit: [erg: erg^-1]*/
    Eigen::MatrixX2d proton_distribution;

    bool verbose = 0;

   private:
    Eigen::VectorXd electron_energy_log_;
    Eigen::VectorXd electron_density_log_;
    std::unique_ptr<GSLInterpData> electron_interpolate;
    std::vector<std::unique_ptr<GSLInterpData>> target_photons_interpolate;
    std::vector<Eigen::VectorXd> target_photons_energy_log_;
    std::vector<Eigen::VectorXd> target_photons_density_log_;
};

}  // namespace CRad
