//
// Created by Dmitri Bagaev on 2019-03-25.
//

#ifndef INMOST_DYNAMIC_SIN_H
#define INMOST_DYNAMIC_SIN_H

#include "dynamic_function.h"

class StaticSin : public DynamicFunction {
public:
    double GetMinimumValue(int iteration) const noexcept override;

    std::pair<double, double> GetMinimumPoint(int iteration) const noexcept override;

    std::pair<double, double> GetXRange() const noexcept override;

    std::pair<double, double> GetYRange() const noexcept override;

    double invoke(double x, double y, int iteration) const noexcept override;
};


#endif //INMOST_DYNAMIC_SIN_H
