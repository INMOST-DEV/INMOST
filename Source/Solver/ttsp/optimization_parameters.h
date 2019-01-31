//
// Created by bvdmitri on 31.01.19.
//

#ifndef INMOST_OPTIMIZATION_PARAMETERS_H
#define INMOST_OPTIMIZATION_PARAMETERS_H

#include <string>
#include <vector>

namespace TTSP {

    typedef std::pair<double, double> OptimizationParameterRange;
    typedef std::pair<std::string, double> OptimizationParameterPoint;
    typedef std::vector<OptimizationParameterPoint> OptimizationParameterPoints;

    class OptimizationParameter {
    private:
        std::string name;
        std::vector<double> values;
        double default_value;

        static void swap(OptimizationParameter &left, OptimizationParameter &right);
    public:
        OptimizationParameter(const std::string &name, const OptimizationParameterRange &range,
                              double step, double default_value);
        OptimizationParameter(const std::string &name, const std::vector<double> &values, double default_value);

        OptimizationParameter(const OptimizationParameter &other);
        OptimizationParameter(OptimizationParameter &&other) noexcept;

        OptimizationParameter &operator=(const OptimizationParameter &other);

        const std::string &GetName() const;
        const std::vector<double> &GetValues() const;
        const double &GetDefaultValue() const;
    };

    typedef std::pair<OptimizationParameter, double> OptimizationParametersEntry;
    typedef std::vector<OptimizationParametersEntry> OptimizationParameters;

    class OptimizationParametersSpace {
    private:
        std::string solver_name;
        std::string solver_prefix;
        OptimizationParameters parameters;

        static void swap(OptimizationParametersSpace &left, OptimizationParametersSpace &right);
    public:
        OptimizationParametersSpace(const std::string &solver_name, const std::string &solver_prefix,
                                    const OptimizationParameters &parameters);

        OptimizationParametersSpace(const OptimizationParametersSpace &other);
        OptimizationParametersSpace(OptimizationParametersSpace &&other) noexcept;

        OptimizationParametersSpace &operator=(const OptimizationParametersSpace &other);

        bool isSolverNameMatch(const std::string &solver_name) const;
        bool isSolverPrefixMatch(const std::string &solver_prefix) const;

        const std::string &GetSolverName() const;
        const std::string &GetSolverPrefix() const;
        const OptimizationParameters &GetParameters() const;
        const OptimizationParameterPoints GetPoints() const;

        void Update(const OptimizationParameterPoints &update);
    };

}


#endif //INMOST_OPTIMIZATION_PARAMETERS_H
