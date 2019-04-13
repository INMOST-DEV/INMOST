//
// Created by Dmitri Bagaev on 2019-04-13.
//

#include "dynamic_r2.h"

DynamicR2::DynamicR2() : distribution(-0.04, 0.04) {
    unsigned int seed = static_cast<unsigned int>(time(NULL));
    generator.seed(seed);
}

double DynamicR2::invoke(double x, double y, int iteration) const noexcept {
    std::pair<double, double> min_point = GetMinimumPoint(iteration);

    double mx = min_point.first;
    double my = min_point.second;

    double tmp1 = std::log10(x / mx);
    double tmp2 = (17.5 * (y - my)) / (7.5 + y - my);

    double f = ((16.0 / 25.0) * (tmp1 * tmp1) + 1) * ((1.0 / 25.0) * (tmp2 * tmp2) + 1);

    return f * (1 + distribution(generator));
}

double DynamicR2::GetMinimumValue(int iteration) const noexcept {
    return 1.0;
}

std::pair<double, double> DynamicR2::GetMinimumPoint(int iteration) const noexcept {
    double kx = 200;
    double ky = 800;

    // double x  = std::pow(10, -2 - std::cos((2 * M_PI * iteration) / kx));
    double x = 1.0 / (1.0 / (std::pow(10, -3) - std::pow(10, -1)) - ((double) iteration / 20.0)) + std::pow(10, -1);
    double y = 2 + std::cos((2 * M_PI * iteration) / ky);
    return std::make_pair(x, y);
}

std::pair<double, double> DynamicR2::GetXRange() const noexcept {
    return std::make_pair(-3.5, -0.5);
}

std::pair<double, double> DynamicR2::GetYRange() const noexcept {
    return std::make_pair(0, 4);
}

TTSP::OptimizationParameterType DynamicR2::GetXParameterType() const noexcept {
    return TTSP::OptimizationParameterType::PARAMETER_TYPE_EXPONENT;
}

