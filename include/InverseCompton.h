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

    Eigen::VectorXd CalculateDifferentialSpectrum(
        const Eigen::Ref<const Eigen::VectorXd>& epsilon_prime) const;

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
                       double epsilon_prime) const;

    //------------------------------------------------------------------------
    //! Integrated emissivity over target photon distribution.
    //!
    //! @return        integration on the target photon distribution
    //------------------------------------------------------------------------
    double emissivity_photon_integrated_sum(double electron_energy,
                                            double epsilon_prime) const;
    //------------------------------------------------------------------------
    //! Integrated over electron distribution
    //!
    //! @param epsilon_prime     scattered photon energy [erg]
    //! @return        dN/dt/d(epsilon_prime)  [s^-1 erg^-1]
    //------------------------------------------------------------------------
    double differential_emission(double epsilon_prime) const;

    RadData& dat;
    std::function<double(double, double, double)> emissivity;
};
}  // namespace CRad
