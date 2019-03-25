//
// Created by Dmitri Bagaev on 2019-03-25.
//

#include "dynamic_x2.h"

#include <cmath>

double DynamicX2::invoke(double x, double y, int iteration) const noexcept {

    double xm = GetMinimumPoint(iteration).first;

    if (x <= xm) {
        double a = 10.0 / (xm * xm);
        double b = -20.0 / xm;
        double c = 10.0;

        return a * x * x + b * x + c;
    } else {
        double a = 10.0 / (xm * xm - 8.0 * xm + 16);
        double b = -2.0 * a * xm;
        double c = b * b / (4.0 * a);

        return a * x * x + b * x + c;
    }
}

double DynamicX2::GetMinimumValue(int iteration) const noexcept {
    return 0.0;
}

std::pair<double, double> DynamicX2::GetMinimumPoint(int iteration) const noexcept {
    return std::make_pair(2 + std::sin(iteration / 50.0), 0.0);
}

std::pair<double, double> DynamicX2::GetXRange() const noexcept {
    return std::make_pair(0, 4);
}

std::pair<double, double> DynamicX2::GetYRange() const noexcept {
    return std::make_pair(0, 1); // unused
}
