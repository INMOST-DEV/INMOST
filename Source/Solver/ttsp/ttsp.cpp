//
// Created by bvdmitri on 04.02.19.
//

#include "ttsp.h"

namespace TTSP {

    void OptimizationParameter::swap(OptimizationParameter &left, OptimizationParameter &right) {
        std::swap(left.name, right.name);
        std::swap(left.values, right.values);
        std::swap(left.default_value, right.default_value);
    }

    OptimizationParameter::OptimizationParameter(const std::string &name,
                                                 const OptimizationParameterRange &range,
                                                 double step, double default_value) :
            name(name),
            values(static_cast<unsigned long>((range.second - range.first) / step)),
            default_value(default_value) {
        int index = 0;
        double value = range.first;
        while (value < range.second) {
            values[index++] = value;
            value += step;
        }
    }

    OptimizationParameter::OptimizationParameter(const std::string &name, const std::vector<double> &values,
                                                 double default_value) :
            name(name), values(values), default_value(default_value) {}

    OptimizationParameter::OptimizationParameter(const OptimizationParameter &other) :
            name(other.name), values(other.values), default_value(other.default_value) {}

    OptimizationParameter::OptimizationParameter(OptimizationParameter &&other) noexcept {
        OptimizationParameter::swap(*this, other);
    }

    OptimizationParameter &TTSP::OptimizationParameter::operator=(const OptimizationParameter &other) {
        OptimizationParameter tmp(other);
        OptimizationParameter::swap(*this, tmp);
        return *this;

    }

    const std::string &OptimizationParameter::GetName() const {
        return name;
    }

    const std::vector<double> &OptimizationParameter::GetValues() const {
        return values;
    }

    const double &OptimizationParameter::GetDefaultValue() const {
        return default_value;
    }

    void OptimizationParametersSpace::swap(OptimizationParametersSpace &left, OptimizationParametersSpace &right) {
        std::swap(left.solver_name, right.solver_name);
        std::swap(left.solver_prefix, right.solver_prefix);
        std::swap(left.parameters, right.parameters);
    }

    OptimizationParametersSpace::OptimizationParametersSpace(const std::string &solver_name,
                                                             const std::string &solver_prefix,
                                                             const OptimizationParameters &parameters) :
            solver_name(solver_name), solver_prefix(solver_prefix), parameters(parameters) {}

    OptimizationParametersSpace::OptimizationParametersSpace(const OptimizationParametersSpace &other) :
            solver_name(other.solver_name), solver_prefix(other.solver_prefix), parameters(other.parameters) {}

    OptimizationParametersSpace::OptimizationParametersSpace(OptimizationParametersSpace &&other) noexcept {
        OptimizationParametersSpace::swap(*this, other);
    }

    OptimizationParametersSpace &OptimizationParametersSpace::operator=(const OptimizationParametersSpace &other) {
        OptimizationParametersSpace tmp(other);
        OptimizationParametersSpace::swap(*this, tmp);
        return *this;
    }

    bool OptimizationParametersSpace::isSolverNameMatch(const std::string &solver_name) const {
        return this->solver_name == solver_name;
    }

    bool OptimizationParametersSpace::isSolverPrefixMatch(const std::string &solver_prefix) const {
        return this->solver_prefix == solver_prefix;
    }

    const std::string &OptimizationParametersSpace::GetSolverName() const {
        return solver_name;
    }

    const std::string &OptimizationParametersSpace::GetSolverPrefix() const {
        return solver_prefix;
    }

    const OptimizationParameters &OptimizationParametersSpace::GetParameters() const {
        return parameters;
    }

    const OptimizationParameterPoints OptimizationParametersSpace::GetPoints() const {
        OptimizationParameterPoints points(parameters.size());
        std::transform(parameters.begin(), parameters.end(), points.begin(), [](const OptimizationParametersEntry &p) {
            return std::make_pair(p.first.GetName(), p.second);
        });
        return points;
    }

    void OptimizationParametersSpace::Update(const OptimizationParameterPoints &update) {
        for (int i = 0; i < update.size(); ++i) {
            parameters[i].second = update[i].second;
        }
    }

    void OptimizationParameterResult::swap(OptimizationParameterResult &left, OptimizationParameterResult &right) {
        std::swap(left.points, right.points);
        std::swap(left.preconditioner_time, right.preconditioner_time);
        std::swap(left.solve_time, right.solve_time);
        std::swap(left.is_solved, right.is_solved);
    }

    OptimizationParameterResult::OptimizationParameterResult(const OptimizationParameterPoints &points,
                                                             double preconditioner_time, double solve_time, bool is_solved) :
            points(points), preconditioner_time(preconditioner_time), solve_time(solve_time), is_solved(is_solved) {}

    OptimizationParameterResult::OptimizationParameterResult(const OptimizationParameterResult &other) :
            points(other.points), preconditioner_time(other.preconditioner_time), solve_time(other.solve_time), is_solved(other.is_solved) {}

    OptimizationParameterResult::OptimizationParameterResult(OptimizationParameterResult &&other) noexcept {
        OptimizationParameterResult::swap(*this, other);
    }

    OptimizationParameterResult &OptimizationParameterResult::operator=(const OptimizationParameterResult &other) {
        OptimizationParameterResult tmp(other);
        OptimizationParameterResult::swap(*this, tmp);
        return *this;
    }

    const OptimizationParameterPoints &OptimizationParameterResult::GetPoints() const {
        return points;
    }

    const double &OptimizationParameterResult::GetPreconditionerTime() const {
        return preconditioner_time;
    }

    const double &OptimizationParameterResult::GetSolveTime() const {
        return solve_time;
    }

    double OptimizationParameterResult::GetTime() const {
        return preconditioner_time + solve_time;
    }

    bool OptimizationParameterResult::IsSolved() const {
        return is_solved;
    }

    void OptimizationParameterResultsBuffer::swap(OptimizationParameterResultsBuffer &left,
                                                  OptimizationParameterResultsBuffer &right) {
        std::swap(left.buffer, right.buffer);
        std::swap(left.capacity, right.capacity);
    }

    OptimizationParameterResultsBuffer::OptimizationParameterResultsBuffer(std::size_t capacity) : capacity(capacity) {}

    OptimizationParameterResultsBuffer::OptimizationParameterResultsBuffer(const OptimizationParameterResultsBuffer &other) :
            buffer(other.buffer), capacity(other.capacity) {}

    OptimizationParameterResultsBuffer::OptimizationParameterResultsBuffer(OptimizationParameterResultsBuffer &&other) noexcept {
        OptimizationParameterResultsBuffer::swap(*this, other);
    }

    OptimizationParameterResultsBuffer &OptimizationParameterResultsBuffer::operator=(const OptimizationParameterResultsBuffer &other) {
        OptimizationParameterResultsBuffer tmp(other);
        OptimizationParameterResultsBuffer::swap(*this, tmp);
        return *this;
    }

    std::deque<OptimizationParameterResult>::const_reverse_iterator OptimizationParameterResultsBuffer::begin() const {
        return buffer.crbegin();
    }

    std::deque<OptimizationParameterResult>::const_reverse_iterator OptimizationParameterResultsBuffer::end() const {
        return buffer.crend();
    }

    void OptimizationParameterResultsBuffer::push(const OptimizationParameterResult &result) {
        if (buffer.size() == capacity) {
            buffer.pop_front();
        }
        buffer.push_back(result);
    }

    std::size_t OptimizationParameterResultsBuffer::size() const {
        return buffer.size();
    }

    bool OptimizerInterface::Solve(INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix, INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL,
                                   GetPreconditionerTimeFromSolverLambda preconditioner_time, GetSolveTimeFromSolverLambda solve_time) {
        const OptimizationParameterPoints &current = space.GetPoints();

        std::for_each(current.begin(), current.end(), [&solver](const OptimizationParameterPoint &p) {
            solver.SetParameter(p.first, INMOST::to_string(p.second));
        });

        solver.SetMatrix(matrix);
        bool solved = solver.Solve(RHS, SOL);

        SaveResult(current, solver, solved, preconditioner_time, solve_time);

        const OptimizationParameterPoints &update = this->MakeOptimizationIteration(solver, matrix, RHS);

        space.Update(update);

        return solved;
    }

    void OptimizerInterface::SaveResult(const OptimizationParameterPoints &points, const INMOST::Solver &solver, bool is_solved,
                                        GetPreconditionerTimeFromSolverLambda preconditioner_time, GetSolveTimeFromSolverLambda solve_time) {
        results.push(OptimizationParameterResult(points, preconditioner_time(solver), solve_time(solver), is_solved));
    }

    const OptimizationParametersSpace &OptimizerInterface::GetSpace() const {
        return space;
    }
    const OptimizationParameterResultsBuffer &OptimizerInterface::GetResults() const {
        return results;
    };

    double OptimizerInterface::DefaultGetPreconditionerTime(const INMOST::Solver &solver) {
        return atof(solver.GetParameter("time_prec").c_str());
    }

    double OptimizerInterface::DefaultGetSolveTime(const INMOST::Solver &solver) {
        return atof(solver.GetParameter("time_iter").c_str());
    }
}
