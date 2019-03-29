//
// Created by Dmitri Bagaev on 2019-03-25.
//

#ifndef INMOST_DYNAMIC_X2_H
#define INMOST_DYNAMIC_X2_H

#include "dynamic_function.h"
#include <random>

class DynamicX2 : public DynamicFunction {
private:
    mutable std::mt19937                     generator;
    mutable std::uniform_real_distribution<> distribution;
public:
    DynamicX2();

    double invoke(double x, double y, int iteration) const noexcept override;

    double GetMinimumValue(int iteration) const noexcept override;

    std::pair<double, double> GetMinimumPoint(int iteration) const noexcept override;

    std::pair<double, double> GetXRange() const noexcept override;

    std::pair<double, double> GetYRange() const noexcept override;
};


#endif //INMOST_DYNAMIC_X2_H
