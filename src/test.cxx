#include <Utils.h>
#include "Eigen/Dense"
#include "Radiation.h"

using namespace CRad;

int main(int argc, char** argv) {
    auto energy = Utils::GetLogspaceVec(1e-3, 1e3, 400);
    Eigen::VectorXd density = 1e4 * pow((energy.array()), -2);
    auto rad = new Radiation();
    rad->SetElectronDistribution(energy, density);
    rad->AddBlackBodyPhotons(2.7);
    auto spectrum_energy = Utils::GetLogspaceVec(1e-1, 1e14, 500);
    rad->CalculateDifferentialSpectrum(spectrum_energy);
}
