/*
 *
 * Store the Initial Data from the radiation computation
 *
 */
#pragma once
#include <vector>
#include "Eigen/Core"

namespace CRad {
class RadData {

   public:
    RadData();
    /* Set the spectra of electron:  array of energy[erg], array of number density[erg^-1]*/
    void SetElectronDis(const Eigen::VectorXd& energy,
                        const Eigen::VectorXd& density);
    std::pair<double, double> GetMinMaxElectronEnerg();
    double GetElectronDensity(double energy);

    /* Set the spectra of proton:  array of energy[erg], array of number density[erg^-1]*/
    void SetProtonDis(const Eigen::VectorXd& energy,
                      const Eigen::VectorXd& density);
    /* Add arbitary isotropic target photon distribution*: array of energy[erg] array of density */
    void AddTargetPhotons(const Eigen::VectorXd& energy,
                          const Eigen::VectorXd& density);
    /* Add blackbody/Graybody target photons */
    void AddThermalTargetPhotons(const double temerature, int bins = 1000,
                                 double energy_density = 0);
    std::pair<double, double> GetMinMaxTargetEnergy();
    /* return a VectorXd of multiple components photon density  energy[erg]*/
    Eigen::VectorXd GetPhotonDensity(double energy);
    /* target photons for InverseCompton. May have multiple components unit:[erg, erg^-1]*/
    double GetSumPhotonDensity(double energy);
    std::vector<Eigen::Matrix2Xd> target_photons;
    /* magnetic field  unit:[gauss]*/
    double B_field;
    /* electron distribution unit:[erg, erg^-1]*/
    Eigen::MatrixX2d electron_distribution;
    /* Convert the N(E)dE to ln10EN(x)dx*/
    Eigen::MatrixX2d electron_distribution_log;
    /* proton distribution unit: [erg: erg^-1]*/
    Eigen::MatrixX2d proton_distribution;
    bool verbose = 0;
};

}  // namespace CRad
