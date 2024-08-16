#pragma  once
#include <limits>
#include <variant>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "InverseCompton.h"
namespace CRad {
    class Radiation
    {
        public:
            void SetProtonDistribution(const Eigen::Ref<const Eigen::VectorXd>& energy, const Eigen::Ref<const Eigen::VectorXd>& density);
            void SetElectronDistribution(const Eigen::Ref<const Eigen::VectorXd>& energy, const Eigen::Ref<const Eigen::VectorXd>& density);
            void AddBlackBodyPhotons(double temperature);
            void AddArbitaryPhotons(const Eigen::Ref<const Eigen::VectorXd>& energy, const Eigen::Ref<const Eigen::VectorXd>& density);


        private:
            std::vector<std::variant<std::monostate, InverseCompton>> radiation_mechanism_;
            std::vector<Eigen::VectorXd> differential_spectrum;

            template<typename T>
                Eigen::VectorXd GetSpectrum() const
                {
                    for(size_t i = 0 ; i < radiation_mechanism_.size(); i++)
                    {
                        if(std::holds_alternative<T>(radiation_mechanism_[i]))
                        {
                            return differential_spectrum[i];
                        }
                    }
                    return Eigen::VectorXd::Constant(differential_spectrum.empty()? 0: differential_spectrum[0].size(), std::numeric_limits<double>::quiet_NaN());
                }

    };

}
