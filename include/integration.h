#include <omp.h>
#include <cmath>
#include <functional>

/*
 * Transform the variable x --> u : x = u/1-u
 */
double transform_infinity(double u);

double transfrom_derivative(double u);

double simpson_integration(const std::function<double(double)>& f, double a,
                           double b, int n, bool use_log_space = false);

double parallel_simpson_integration(const std::function<double(double)>& f,
                                    double a, double b, int n,
                                    bool use_log_space = false);

double simpson_double_integration(
    const std::function<double(double, double)>& f, double ax, double bx,
    double ay, double by, int nx, int ny, bool use_log_spacex = false,
    bool use_log_spacey = false);
double simpson_integration_infinity(const std::function<double(double)>& f,
                                    double a);
