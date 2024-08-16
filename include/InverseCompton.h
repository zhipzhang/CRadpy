#pragma once
#include <functional>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "RadData.h"

namespace CRad {
class InverseCompton {
   public:
    InverseCompton(RadData& data) : dat{data} {
        emissivity =
            std::bind(&InverseCompton::emissivity_, this, std::placeholders::_1,
                      std::placeholders::_2, std::placeholders::_3);
    };
    void initialize_kernel();
    Eigen::VectorXd CalculateDifferentialSpectrum(
        const Eigen::Ref<const Eigen::VectorXd>& epsilon_prime);

   private:
    //------------------------------------------------------------------------
    //! Compute the emissivity of InverseCompton
    //! This is come from the Eq.48 of the Gould(1970)
    //! Assuming a isotropic photon field and electron distribution in Lab Frame
    //!
    //! @param epsilon     target photon energy [erg]
    //! @param electron_energy     electron energy [erg]
    //! @param epsilon_prime     scattered energy [erg]
    //! @return        emissivity: dN(E_e, E_photon)/dt/depsilon_prime
    //------------------------------------------------------------------------
    double emissivity_(double epsilon, double electron_energy,
                       double epsilon_prime);

    //------------------------------------------------------------------------
    //! Integrated emissivity over target photon distribution.
    //!
    //! @return        return_description
    //------------------------------------------------------------------------
    double emissivity_photon_integrated_sum(double electron_energy,
                                            double epsilon_prime);
    double differential_emission(double epsilon_prime);

    std::vector<Eigen::MatrixXd> ic_spec_kernel;
    /* Kernel to compute the energy loss term  i: electron_energy j: photon_energy*/
    Eigen::MatrixXd ic_loss_kernel;
    RadData& dat;
    std::function<double(double, double, double)> emissivity;
};
}  // namespace CRad
