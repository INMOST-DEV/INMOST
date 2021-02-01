//
// Created by Dmitri Bagaev on 2019-03-29.
//

#include "static_sin_r2.h"

#include <cmath>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

double StaticSinR2::invoke(double x, double y, int iteration) const noexcept {
    return -1.0 * std::sin(x) * std::cos(y) + 1.0;
}

double StaticSinR2::GetMinimumValue(int iteration) const noexcept {
    return 0;
}

std::pair<double, double> StaticSinR2::GetMinimumPoint(int iteration) const noexcept {
    return std::make_pair(M_PI / 2.0, 0);
}

std::pair<double, double> StaticSinR2::GetXRange() const noexcept {
    return std::make_pair(0, M_PI);
}

std::pair<double, double> StaticSinR2::GetYRange() const noexcept {
    return std::make_pair(-M_PI / 2.0, M_PI / 2.0);
}
