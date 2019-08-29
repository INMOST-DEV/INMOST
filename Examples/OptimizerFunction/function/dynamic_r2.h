//
// Created by Dmitri Bagaev on 2019-04-13.
//

#ifndef INMOST_DYNAMIC_R2_H
#define INMOST_DYNAMIC_R2_H

#include "dynamic_function.h"
#include <random>

class DynamicR2 : public DynamicFunction {
private:
    mutable std::mt19937                     generator;
    mutable std::uniform_real_distribution<> distribution;
public:
    DynamicR2();

    double invoke(double x, double y, int iteration) const noexcept override;

    double GetMinimumValue(int iteration) const noexcept override;

    std::pair<double, double> GetMinimumPoint(int iteration) const noexcept override;

    std::pair<double, double> GetXRange() const noexcept override;

    std::pair<double, double> GetYRange() const noexcept override;

    INMOST::OptimizationParameterType GetXParameterType() const noexcept override;
};

#endif //INMOST_DYNAMIC_R2_H
