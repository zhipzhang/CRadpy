#pragma once
#include <limits>
#include <memory>
#include <variant>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "InverseCompton.h"
#include "RadData.h"
#include "constants.h"
namespace CRad {
class Radiation {
   public:
    Radiation() { radiation_data_ = std::make_unique<RadData>(); };
    void SetProtonDistribution(
        const Eigen::Ref<const Eigen::VectorXd>& energy,
        const Eigen::Ref<const Eigen::VectorXd>& density);
    void SetElectronDistribution(
        const Eigen::Ref<const Eigen::VectorXd>& energy,
        const Eigen::Ref<const Eigen::VectorXd>& density);
    void SetDistance(double distance) { distance_ = distance; }
    void AddBlackBodyPhotons(double temperature);
    void AddArbitaryPhotons(const Eigen::Ref<const Eigen::VectorXd>& energy,
                            const Eigen::Ref<const Eigen::VectorXd>& density);
    void CalculateDifferentialSpectrum(
        const Eigen::Ref<const Eigen::VectorXd>& epsilon_prime);
    Eigen::VectorXd GetTotalSpectrum();
    Eigen::VectorXd GetICSpectrum() {
        Eigen::VectorXd luminosity = GetSpectrum<InverseCompton>();
        if (distance_ == 0)
            return luminosity;
        else
            return luminosity.array() / 4 / constants::pi / distance_ /
                   distance_;
    }

   private:
    std::vector<std::variant<InverseCompton>> radiation_mechanism_;
    std::vector<Eigen::VectorXd> differential_spectrum;

    template <typename T>
    Eigen::VectorXd GetSpectrum() const {
        for (size_t i = 0; i < radiation_mechanism_.size(); i++) {
            if (std::holds_alternative<T>(radiation_mechanism_[i])) {
                return differential_spectrum[i];
            }
        }
        return Eigen::VectorXd::Constant(
            differential_spectrum.empty() ? 0 : differential_spectrum[0].size(),
            std::numeric_limits<double>::quiet_NaN());
    }
    std::unique_ptr<RadData> radiation_data_;
    Eigen::VectorXd PhotonOutEnergy_;
    double distance_ = 0;
};

}  // namespace CRad
