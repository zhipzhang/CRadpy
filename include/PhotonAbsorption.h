#include <Eigen/Dense>
#include <memory>
#include "Eigen/Core"
namespace CRad {

class PhotonDistribution;
class PhotonAbsorption {
   public:
    PhotonAbsorption(std::shared_ptr<PhotonDistribution> bgPhotons);

   private:
    //------------------------------------------------------------------------
    //! compute the photon photon cross-section
    //!
    //! @param     s      s = 0.5 * (epsilon * epsilon_1) * (1 - cos\theta) [constant]
    //! @return           cross-section [cm^2]
    //------------------------------------------------------------------------
    static Eigen::ArrayXd sigma_(const Eigen::ArrayXd& s);

    //------------------------------------------------------------------------
    //! photon-photon cross-section for solid-angle averaged
    //!
    //! @param epsilon     incident photon energy [erg]
    //! @param bgphotone     background photon energy [erg]
    //! @return        cross-section [cm^2]
    //------------------------------------------------------------------------
    static Eigen::ArrayXd AverageSigma_(const double epsilon,
                                        const Eigen::ArrayXd& bgphotone);
    std::shared_ptr<PhotonDistribution> bgPhotons_;
};
}  // namespace CRad
