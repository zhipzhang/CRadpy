#include "Utils.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <cmath>
#include <iostream>
#include <memory>
#include "Eigen/Core"
#include "constants.h"
Eigen::VectorXd CRad::Utils::GetLogspaceVec(double start, double end, int num) {
    double log_start = std::log10(start);
    double log_end = std::log10(end);
    auto log_vec = Eigen::VectorXd::LinSpaced(num, log_start, log_end);
    return pow(10, log_vec.array()).matrix();
}

Eigen::VectorXd CRad::Utils::BlackBody(
    const Eigen::Ref<const Eigen::MatrixXd>& energy, double temperture) {
    auto nv = energy.array() / constants::hp;
    Eigen::ArrayXd exp_part = exp(energy.array()/(constants::kb * temperture)) - 1;
    Eigen::ArrayXd deno =
        pow(constants::c_speed, 3) * exp_part;
    Eigen::ArrayXd numerator = 8 * constants::pi / constants::hp * pow(nv, 2);
    /* Assure all element is greater than 0 */
    Eigen::ArrayXd result = (deno.inverse() * numerator).max(1e-243);
    return result.matrix();
}

Eigen::VectorXd CRad::Utils::GreyBody(
    const Eigen::Ref<const Eigen::MatrixXd>& energy, double temperture,
    double energy_density) {
    Eigen::ArrayXd exp_part =
        exp(energy.array() / (constants::kb * temperture)) - 1;
    auto density = energy_density * 15 / pow(constants::pi, 4) /
                   pow(constants::kb * temperture, 4) * pow(energy.array(), 2);
    Eigen::ArrayXd result = (exp_part.inverse() * density).max(1e-150);
    return result;
}

namespace CRad {

std::unordered_map<size_t, std::unique_ptr<GSLInterpData>> cache;

double Utils::interpolate(const Eigen::Ref<const Eigen::VectorXd>& x,
                          const Eigen::Ref<const Eigen::VectorXd>& y,
                          double x0) {
    if (x0 > x.maxCoeff() || x0 < x.minCoeff()) {
        return 0;
    }
    size_t hash_key = hashEigenVectors(x, y);

    // 检查是否已有缓存
    auto it = cache.find(hash_key);
    if (it != cache.end()) {
        return gsl_spline_eval(it->second->spline, x0, it->second->acc);
    }

    // 创建并插入 GSLInterpData 对象
    auto new_data = std::make_unique<GSLInterpData>(x, y);
    double result = gsl_spline_eval(new_data->spline, x0, new_data->acc);
    cache.emplace(hash_key, std::move(new_data));

    return result;
}
double Utils::integrate_1d(const std::function<double(double)>& f,
                           double low_bound, double up_bound, double eps_rel,
                           double eps_abs) {
    gsl_function gsl_func;
    gsl_func.function = &gsl_function_wrapper_;
    gsl_func.params = const_cast<void*>(reinterpret_cast<const void*>(&f));

    gsl_integration_workspace* workspace =
        gsl_integration_workspace_alloc(10000);
    double result, abserr;
    gsl_integration_qag(&gsl_func, low_bound, up_bound, eps_abs, eps_rel, 10000,
                        6, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}
double Utils::gsl_function_wrapper_(double x, void* params) {
    // 将参数转换为 std::function<double(double)>* 类型
    std::function<double(double)>* func =
        static_cast<std::function<double(double)>*>(params);
    // 调用 std::function 并返回结果
    return (*func)(x);
}

}  // namespace CRad
