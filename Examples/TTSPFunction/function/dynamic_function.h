//
// Created by Dmitri Bagaev on 2019-03-25.
//

#ifndef INMOST_DYNAMIC_FUNCTION_H
#define INMOST_DYNAMIC_FUNCTION_H

#include <utility>

class DynamicFunction {
public:
    virtual double invoke(double x, double y, int iteration) const noexcept = 0;

    virtual double GetMinimumValue(int iteration) const noexcept = 0;

    virtual std::pair<double, double> GetMinimumPoint(int iteration) const noexcept = 0;

    virtual std::pair<double, double> GetXRange() const noexcept = 0;

    virtual std::pair<double, double> GetYRange() const noexcept = 0;
};


#endif //INMOST_DYNAMIC_FUNCTION_H
