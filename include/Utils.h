#include "Eigen/Core"

namespace CRad {
class Utils {
    using fp = double (*)(double);

   public:
    static Eigen::VectorXd GetLogspaceVec(double start, double end, int num);
    static Eigen::VectorXd BlackBody(
        const Eigen::Ref<const Eigen::MatrixXd>& energy, double temperture);
    static Eigen::VectorXd GreyBody(
        const Eigen::Ref<const Eigen::MatrixXd>& energy, double temperture,
        double energy_density);
    static double interpolate(const Eigen::Ref<const Eigen::VectorXd>& x,
                              const Eigen::Ref<const Eigen::VectorXd>& y, double x0);
    static double integrate_1d(const std::function<double(double)>& f,
                               double low_bound, double up_bound,
                               double eps_rel = 1e-2, double eps_abs = 0);

   private:
    static double gsl_function_wrapper_(double x, void* params);
};
}  // namespace CRad
