//
// Created by Dmitri Bagaev on 2019-03-29.
//

#ifndef INMOST_STATIC_SIN_R2_H
#define INMOST_STATIC_SIN_R2_H


#include "dynamic_function.h"

class StaticSinR2 : public DynamicFunction {
public:
    double invoke(double x, double y, int iteration) const noexcept override;

    double GetMinimumValue(int iteration) const noexcept override;

    std::pair<double, double> GetMinimumPoint(int iteration) const noexcept override;

    std::pair<double, double> GetXRange() const noexcept override;

    std::pair<double, double> GetYRange() const noexcept override;
};


#endif //INMOST_STATIC_SIN_R2_H
