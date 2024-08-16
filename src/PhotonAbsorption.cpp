#include "PhotonAbsorption.h"
#include "Eigen/Core"
#include "constants.h"

Eigen::ArrayXd CRad::PhotonAbsorption::sigma_(const Eigen::ArrayXd& s) {
    Eigen::ArrayXd values = Eigen::ArrayXd::Zero(s.size());

    Eigen::ArrayXd valid_s = s >= 1;
    if (valid_s.any()) {
        auto s_filtered = s(valid_s);
        auto v_cm = (1 - 1.0 / s_filtered).sqrt();
        auto prefactor = 3.0 / 16.0 * constants::sigma_T * (1 - v_cm.pow(2));
        auto term1 = (3 - v_cm.pow(4)) * ((1 + v_cm) / (1 - v_cm)).log();
        auto term2 = 2 * v_cm * (v_cm.pow(2) - 2);
        values(valid_s) = prefactor * (term1 + term2);
    }
    return values;
}
Eigen::ArrayXd CRad::PhotonAbsorption::AverageSigma_(
    const double epsilon, const Eigen::ArrayXd& bgphotone) {
    auto s = (epsilon * bgphotone) / (constants::m_e * constants::m_e);
    Eigen::ArrayXd res = Eigen::ArrayXd::Zero(s.size());
    Eigen::ArrayXd valid_s = s >= 1;
    if (valid_s.any()) {
        auto term1 = 3.0 / (2 * s.pow(2)) * constants::sigma_T;
        auto term2 = (s + 0.5 * s.log() - 1.0 / 6 + 0.5 / s);
        auto term3 = (s.sqrt() + (1 - 1 / s).sqrt()).log();
        auto term4 = (s + 4.0 / 9.0 - 1.0 / (9.0 * s)) * (1 - 1.0 / s).sqrt();
        res(valid_s) = term1 * (term2 * term3 - term4);
    }
    return res;
}
