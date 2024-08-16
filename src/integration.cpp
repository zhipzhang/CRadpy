#include "integration.h"
#include <cmath>
#include <functional>
#include <iostream>
double transform_infinity(double u) {
    return u / (1.0 - u);  // Map [0, 1) -> [0, ∞)
}

double transform_variable(double x) {
    return x / (x + 1);
}
double transform_derivative(double u) {
    double denom = (1.0 - u) * (1.0 - u);
    return 1.0 / denom;
}

double simpson_integration(const std::function<double(double)>& f, double a,
                           double b, int n, bool use_log_space) {
    double h, s;

    // Handle logarithmic spacing
    if (use_log_space) {
        double log_a = std::log(a);
        double log_b = std::log(b);
        h = (log_b - log_a) / n;
        s = f(a) * a + f(b) * b;

        for (int i = 1; i < n; i += 2) {
            double u = log_a + i * h;
            s += 4 * f(std::exp(u)) * std::exp(u);
        }

        for (int i = 2; i < n - 1; i += 2) {
            double u = log_a + i * h;
            s += 2 * f(std::exp(u)) * std::exp(u);
        }

        return s * h / 3.0;
    } else {
        h = (b - a) / n;
        s = f(a) + f(b);

        for (int i = 1; i < n; i += 2) {
            s += 4 * f(a + i * h);
        }

        for (int i = 2; i < n - 1; i += 2) {
            s += 2 * f(a + i * h);
        }

        return s * h / 3.0;
    }
}

double parallel_simpson_integration(const std::function<double(double)>& f,
                                    double a, double b, int n,
                                    bool use_log_space) {
    double sum = 0.0;
    int num_threads = omp_get_max_threads();
    int local_n = n / num_threads < 50 ? 50 : n / num_threads;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        double local_sum = 0.0;
        double local_a = a + tid * (b - a) / num_threads;
        double local_b = local_a + (b - a) / num_threads;

        local_sum =
            simpson_integration(f, local_a, local_b, local_n, use_log_space);

#pragma omp critical
        { sum += local_sum; }
    }

    return sum;
}

double parallel_double_integral(const std::function<double(double, double)>& f,
                                double ax, double bx, double ay, double by,
                                int nx, int ny, bool use_log_space_x,
                                bool use_log_space_y) {
    double result = 0.0;

#pragma omp parallel for reduction(+ : result)
    for (int i = 0; i < nx; ++i) {
        double x_start = ax + i * (bx - ax) / nx;
        double x_end = x_start + (bx - ax) / nx;

        auto integrand = [&](double x) -> double {
            auto inner_f = [&](double y) -> double {
                return f(x, y);
            };
            return parallel_simpson_integration(inner_f, ay, by, ny,
                                                use_log_space_y);
        };

        double local_result =
            simpson_integration(integrand, x_start, x_end, 1,
                                use_log_space_x);  // 1 step for outer integral
        result += local_result;
    }

    return result;
}

double simpson_integration_infinity(const std::function<double(double)>& f,
                                    double a) {
    std::function<double(double)> f_transformed = [&](double u) -> double {
        double x = transform_infinity(u);
        return f(x) * transform_derivative(u);
    };
    return parallel_simpson_integration(f_transformed, transform_variable(a),
                                        1 - 1e-8, 500);
};
