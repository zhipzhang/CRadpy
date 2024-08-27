#pragma once
#include <limits>
#include <memory>
#include <variant>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "InverseCompton.h"
#include "RadData.h"
#include "constants.h"
#include "Synchrotron.h"
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
    void SetB_Field(double B)
    {
        radiation_data_->SetB_Field(B);
    }
    void CalculateDifferentialSpectrum(
        const Eigen::Ref<const Eigen::VectorXd>& epsilon_prime);
    Eigen::VectorXd GetTotalSpectrum();
    Eigen::VectorXd GetICSpectrum() {
        Eigen::VectorXd luminosity = GetSpectrum<InverseCompton>();
        return luminosity;
    }
    Eigen::VectorXd GetSynSpectrum(){
        Eigen::VectorXd luminosity = GetSpectrum<Synchrontron>();
            return luminosity;
    }
    Eigen::VectorXd GetProtonSynSpectrum()
    {
        Eigen::VectorXd luminosity = GetSpectrum<ProtonSynchrotron>();
        return luminosity;
    }
    void SetSynApprox()
    {
        Syn_approx_ = true;
    };
    //  This is only used for python test!!
    const RadData* GetRadiationData()
    {
        return radiation_data_.get();
    }

   private:
    std::vector<std::variant<InverseCompton, Synchrontron, ProtonSynchrotron>> radiation_mechanism_;
    std::vector<Eigen::VectorXd> differential_spectrum;

    template <typename T>
    Eigen::VectorXd GetSpectrum() const {
        for (size_t i = 0; i < radiation_mechanism_.size(); i++) {
            if (std::holds_alternative<T>(radiation_mechanism_[i])) {
                if (distance_ == 0)
                return differential_spectrum[i];
                else
                 return differential_spectrum[i].array() / 4 / constants::pi / distance_ / distance_;
            }
        }
        return Eigen::VectorXd::Constant(
            differential_spectrum.empty() ? 0 : differential_spectrum[0].size(),
            std::numeric_limits<double>::quiet_NaN());
    }
    std::unique_ptr<RadData> radiation_data_;
    Eigen::VectorXd PhotonOutEnergy_;
    double distance_ = 0;
    bool Syn_approx_ = false;
};

}  // namespace CRad
