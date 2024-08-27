#include "Eigen/Core"
#include "RadData.h"
#include <functional>

namespace CRad
{   
    class Synchrontron
    {
        friend class ProtonSynchrotron;
        public:
        Synchrontron(RadData& data):dat{data}{
            emissivity = 
                std::bind(&Synchrontron::emissivity_, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        };
        void SetApproxEmissivity()
        {
            emissivity = std::bind(&Synchrontron::emissivity_approx_, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        }
        Eigen::VectorXd CalculateDifferentialSpectrum(const Eigen::Ref<const Eigen::VectorXd>& epsilon) const;
        ~Synchrontron(){};

        protected:

            // modified Bessel functions
            inline static double K_(double nu, double x);
            // compute the gyrofrequency
            inline static double Nu_b_(double B_field);
            // compute the nu_c = 2/3 \nu_b \gamma^2
            inline static  double Nu_c_(double B_field, double gamma);
            //------------------------------------------------------------------------
            //! Compute the emissivity of Synchrontron: exact version
            //! Euquation 7.43 in Derner2009.
            //! @param B_field     Magnetic Field [G]
            //! @param electron_energy     electron energy [erg]
            //! @param epsilon     Photon Energy [erg]
            //! @return        emissivity_ dN(E_e)/dt/depsilon_prime [erg/s/erg]
            //------------------------------------------------------------------------
            double emissivity_(double B_field, double electron_energy, double epsilon);

            double emissivity_approx_(double B_field, double electron_energy, double epsilon);
            
            // using numerical approximation of R(x): EqD.7 of [Aharonian2010]_, is used.
            static double R_(double x);
            //------------------------------------------------------------------------
            //! Integration over emissivity_ over electron distribution
            //!
            //! @param epsilon     output photon energy [erg]
            //! @return        dN/dt/depsilon_prime [erg/s/erg]
            //------------------------------------------------------------------------
            double differential_emission(double epsilon)  const;
            
            static double energy_loss_(double energy, double B_field);
            RadData& dat;
            std::function<double(double, double, double)> emissivity;

    };
    
    // Similar to Electron Synchrontron, we will check two ways of computation.
    // TODO : find the way to incorporate this two Synchrontron into one uniform way.
    class ProtonSynchrotron
    {
        public:
            ProtonSynchrotron(RadData& data):dat{data} {
                emissivity = std::bind(&ProtonSynchrotron::emissivity_, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
            };
            ~ProtonSynchrotron(){};
        void SetApproxEmissivity()
        {
            emissivity = std::bind(&ProtonSynchrotron::emissivity_approx_, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        }
        Eigen::VectorXd CalculateDifferentialSpectrum(const Eigen::Ref<const Eigen::VectorXd>& epsilon) const;
        private:
            double emissivity_(double B_field, double proton_energy, double epsilon);
            double emissivity_approx_(double B_field, double proton_energy, double epsilon);
            double Nu_b_(double B_field);
            RadData& dat;
            std::function<double(double, double, double)> emissivity;
            double x_(double nu, double B_field, double gamma);

            double differential_emission(double epsilon)  const;
            static double energy_loss_(double energy);



    };
}
