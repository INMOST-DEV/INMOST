//
// Created by Dmitri Bagaev on 2019-03-25.
//

#include "static_sin.h"

#include <cmath>

double StaticSin::invoke(double x, double y, int iteration) const noexcept {
    return 1.0 - std::sin(x);
}

std::pair<double, double> StaticSin::GetXRange() const noexcept {
    return std::make_pair(0, M_PI);
}

std::pair<double, double> StaticSin::GetYRange() const noexcept {
    return std::make_pair(0, M_PI);
}

double StaticSin::GetMinimumValue(int iteration) const noexcept {
    return 0;
}

std::pair<double, double> StaticSin::GetMinimumPoint(int iteration) const noexcept {
    return std::make_pair(M_PI / 2, 0.0);
}
