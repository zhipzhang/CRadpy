#include "Radiation.h"
#include <utility>
#include "Eigen/Core"
#include "RadData.h"

namespace CRad {
void Radiation::SetProtonDistribution(
    const Eigen::Ref<const Eigen::VectorXd>& energy,
    const Eigen::Ref<const Eigen::VectorXd>& density) {
    radiation_data_->SetProtonDis(energy, density);
}
void Radiation::SetElectronDistribution(
    const Eigen::Ref<const Eigen::VectorXd>& energy,
    const Eigen::Ref<const Eigen::VectorXd>& density) {
    radiation_data_->SetElectronDis(energy, density);
}
void Radiation::AddBlackBodyPhotons(double temperature) {
    radiation_data_->AddThermalTargetPhotons(temperature);
}
void Radiation::AddArbitaryPhotons(
    const Eigen::Ref<const Eigen::VectorXd>& energy,
    const Eigen::Ref<const Eigen::VectorXd>& density) {
    radiation_data_->AddTargetPhotons(energy, density);
}
void Radiation::CalculateDifferentialSpectrum(
    const Eigen::Ref<const Eigen::VectorXd>& epsilon_prime) {
    const int size = epsilon_prime.size();
    PhotonOutEnergy_ = epsilon_prime;
    if (radiation_data_->HaveElectron()) {
        if (radiation_data_->HaveTargetPhotons()) {
             radiation_mechanism_.emplace_back(
                std::in_place_type<InverseCompton>, *radiation_data_);
            differential_spectrum.emplace_back(Eigen::VectorXd::Zero(size));
        }
        if( radiation_data_->GetB_Field() > 0)
        {
            auto& synchrontron = radiation_mechanism_.emplace_back(std::in_place_type<Synchrontron>, *radiation_data_);
            if(Syn_approx_)
            {
                std::get<Synchrontron>(synchrontron).SetApproxEmissivity();
            }
            differential_spectrum.emplace_back(Eigen::VectorXd::Zero(size));
        }
    }
    if(radiation_data_->HaveProton())
    {
        if(radiation_data_->GetB_Field() > 0)
        {
            auto& proton_syn = radiation_mechanism_.emplace_back(std::in_place_type<ProtonSynchrotron>, *radiation_data_);
            if( Syn_approx_)
                std::get<ProtonSynchrotron>(proton_syn).SetApproxEmissivity();
            differential_spectrum.emplace_back(Eigen::VectorXd::Zero(size));
        }
    }
    for (size_t i = 0; i < radiation_mechanism_.size(); i++) {
        std::visit(
            [this, i, epsilon_prime](const auto& mech) {
                differential_spectrum[i] =
                    mech.CalculateDifferentialSpectrum(epsilon_prime);
            },
            radiation_mechanism_[i]);
    }
}
Eigen::VectorXd Radiation::GetTotalSpectrum() {
    Eigen::VectorXd total_spectrum =
        Eigen::VectorXd::Zero(PhotonOutEnergy_.size());
    for (const auto& ispectrum : differential_spectrum) {
        total_spectrum += ispectrum;
    }
    return total_spectrum;
}

}  // namespace CRad
