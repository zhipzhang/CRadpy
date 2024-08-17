#pragma once
#include <gsl/gsl_spline.h>
#include "Eigen/Core"

namespace CRad {
struct GSLInterpData {
    gsl_interp_accel* acc;
    gsl_spline* spline;
    GSLInterpData() : acc(nullptr), spline(nullptr) {}
    GSLInterpData(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
        acc = gsl_interp_accel_alloc();
        spline = gsl_spline_alloc(gsl_interp_cspline, x.size());
        gsl_spline_init(spline, x.data(), y.data(), x.size());
    }

    double Interpolate(double x0) { return gsl_spline_eval(spline, x0, acc); }
    ~GSLInterpData() {
        if (spline)
            gsl_spline_free(spline);
        if (acc)
            gsl_interp_accel_free(acc);
    }

   private:
};
class Utils {
   public:
    static Eigen::VectorXd GetLogspaceVec(double start, double end, int num);
    static Eigen::VectorXd BlackBody(
        const Eigen::Ref<const Eigen::MatrixXd>& energy, double temperture);
    static Eigen::VectorXd GreyBody(
        const Eigen::Ref<const Eigen::MatrixXd>& energy, double temperture,
        double energy_density);
    static double interpolate(const Eigen::Ref<const Eigen::VectorXd>& x,
                              const Eigen::Ref<const Eigen::VectorXd>& y,
                              double x0);
    static double integrate_1d(const std::function<double(double)>& f,
                               double low_bound, double up_bound,
                               double eps_rel = 1e-2, double eps_abs = 0);

   private:
    static double gsl_function_wrapper_(double x, void* params);
    static size_t hashEigenVectors(const Eigen::VectorXd& x,
                                   const Eigen::VectorXd& y) {
        std::hash<double> hasher;
        size_t seed = 0;
        for (int i = 0; i < x.size(); ++i) {
            seed ^= hasher(x[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        for (int i = 0; i < y.size(); ++i) {
            seed ^= hasher(y[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};
}  // namespace CRad
