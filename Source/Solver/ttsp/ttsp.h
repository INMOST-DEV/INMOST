//
// Created by bvdmitri on 04.02.19.
//

#ifndef INMOST_TTSP_H
#define INMOST_TTSP_H

#include <vector>
#include <algorithm>
#include <string>
#include <deque>
#include <inmost_solver.h>
#include <Source/Misc/utils.h>

namespace TTSP {

    /// This class is used to define a range of possible values for some parameter
    /// Usage: OptimizationParameterRange range = std::make_pair(0.0, 1.0)
    typedef std::pair<double, double> OptimizationParameterRange;

    /// This class is used to define a solver parameter for optimization algortihms
    /// Usage:
    ///    1.  OptimizationParameter tau("tau", range, 0.1, 1e-3);
    ///    2.  OptimizationParameter kovl("kovl", { 0, 1, 2, 3 }, 0);
    class OptimizationParameter {
    private:
        std::string name;          /// A name of a parameter
        std::vector<double> values;        /// List of possible values for a paramater
        double default_value; /// Default value for this parameter

        static void swap(OptimizationParameter &left, OptimizationParameter &right);
    public:
        /// Default constructor to define an OptimizationParameter with range of possible values and associated step
        /// @param name          - Name of a parameter
        /// @param range         - Range of possible values for this parameter
        /// @param step          - Step associated with the range (used to split up values in a list)
        /// @param default_value - Default value for this parameter
        OptimizationParameter(const std::string &name, const OptimizationParameterRange &range,
                              double step, double default_value);

        /// Default constructor to define an OptimizationParameter with list of possible values
        /// @param name          - Name of a parameter
        /// @param values        - List of possible values for this parameter
        /// @param default_value - Default value for this parameter
        OptimizationParameter(const std::string &name, const std::vector<double> &values, double default_value);

        /// Copy constructor
        /// @param other - OptimizationParameter to make copy of
        OptimizationParameter(const OptimizationParameter &other);

        /// Move constructor
        /// @param other - OptimizationParameter to swap with
        OptimizationParameter(OptimizationParameter &&other) noexcept;

        /// Assignment operator
        /// @param other - OptimizationParameter to assign
        OptimizationParameter &operator=(const OptimizationParameter &other);

        /// Getter for name of a parameter
        const std::string &GetName() const;

        /// Getter for values of a parameter
        const std::vector<double> &GetValues() const;

        /// Getter for default_value of a parameter
        const double &GetDefaultValue() const;
    };

    /// This class is used to store an OptimizationParameter and associated current value
    /// @see OptimizationParameter
    typedef std::pair<OptimizationParameter, double> OptimizationParametersEntry;

    /// This class is used to store a list of OptimizationParameter and associated current value
    /// @see OptimizationParameter
    /// @see OptimizationParametersEntry
    typedef std::vector<OptimizationParametersEntry> OptimizationParameters;

    /// This class is used to store a slice of OptimizationParametersSpace for some optimization parameter
    /// @see OptimizationParameter
    /// @see OptimizationParametersSpace
    typedef std::pair<std::string, double> OptimizationParameterPoint;

    /// This class is used to store a slice of OptimizationParametersSpace for all optimization parameters
    /// @see OptimizationParameter
    /// @see OptimizationParameterPoint
    /// @see OptimizationParametersSpace
    typedef std::vector<OptimizationParameterPoint> OptimizationParameterPoints;

    /// This class is used to define an Optimization parameters space
    /// Usage: OptimizationParametersSpace space("fcbiilu2", "prefix", parameters);
    class OptimizationParametersSpace {
    private:
        std::string solver_name;           /// A name of a solver to which parameters are belong to
        std::string solver_prefix;         /// Solver prefix
        OptimizationParameters parameters; /// List of optimization parameters in this space

        static void swap(OptimizationParametersSpace &left, OptimizationParametersSpace &right);
    public:
        /// Default constructor to define an OptimizationParametersSpace
        /// @param solver_name   - A name of a solver to which parameters are belong to
        /// @param solver_prefix - Solver prefix
        /// @param parameters    - List of optimization parameters in this space
        OptimizationParametersSpace(const std::string &solver_name, const std::string &solver_prefix,
                                    const OptimizationParameters &parameters);

        /// Copy constructor
        /// @param other - OptimizationParametersSpace to make copy of
        OptimizationParametersSpace(const OptimizationParametersSpace &other);

        /// Move constructor
        /// @param other - OptimizationParametersSpace to swap with
        OptimizationParametersSpace(OptimizationParametersSpace &&other) noexcept;

        /// Assignment operator
        /// @param other - OptimizationParametersSpace to assign
        OptimizationParametersSpace &operator=(const OptimizationParametersSpace &other);

        /// Utility method to check against solver name
        /// @param solver_name - A name of a solver to check against
        bool isSolverNameMatch(const std::string &solver_name) const;

        /// Utility method to check against solver prefix
        /// @param solver_prefix - Solver prefix to check against
        bool isSolverPrefixMatch(const std::string &solver_prefix) const;

        /// Getter for solver name
        const std::string &GetSolverName() const;

        /// Getter for solver prefix
        const std::string &GetSolverPrefix() const;

        /// Getter for parameters of this space
        const OptimizationParameters &GetParameters() const;

        /// Getter for slice of parameters of this space
        const OptimizationParameterPoints GetPoints() const;

        void Update(const OptimizationParameterPoints &update);
    };

    /// This class is used to store timings result for some optimization parameter points
    class OptimizationParameterResult {
    private:
        OptimizationParameterPoints points;              /// Optimization parameter points to store result for
        double preconditioner_time; /// Preconditioner timings
        double solve_time;          /// Iterations timings
        bool   is_solved;           /// Is problem solved or not

        static void swap(OptimizationParameterResult &left, OptimizationParameterResult &right);
    public:
        /// Default constructor to define an OptimizationParameterResult
        /// @param points              - Optimization parameter points to store result for
        /// @param preconditioner_time - Preconditioner timings
        /// @param solve_time          - Solve timings
        /// @param is_solved           - Is problem solved or not
        OptimizationParameterResult(const OptimizationParameterPoints &points,
                                    double preconditioner_time, double solve_time, bool is_solved);

        /// Copy constructor
        /// @param other - OptimizationParameterResult to make copy of
        OptimizationParameterResult(const OptimizationParameterResult &other);

        /// Move constructor
        /// @param other - OptimizationParameterResult to swap with
        OptimizationParameterResult(OptimizationParameterResult &&other) noexcept;

        /// Assignment operator
        /// @param other - OptimizationParameterResult to assign
        OptimizationParameterResult &operator=(const OptimizationParameterResult &other);

        /// Getter for points of this result
        const OptimizationParameterPoints &GetPoints() const;

        /// Getter for preconditioner time
        const double &GetPreconditionerTime() const;

        /// Getter for solve time
        const double &GetSolveTime() const;

        /// Getter for is solved status
        bool IsSolved() const;

        /// Getter for cumulative time
        double GetTime() const;
    };

    /// This class is used to store last N solve results and used parameters on it
    /// Implemented as a circular (ring) container
    class OptimizationParameterResultsBuffer {
    private:
        std::size_t capacity;                            /// Max capacity of the results buffer
        std::deque<OptimizationParameterResult> buffer;  /// Buffer to store results in

        static void swap(OptimizationParameterResultsBuffer &left, OptimizationParameterResultsBuffer &right);
    public:
        /// Default constructor to define an OptimizationParameterResultsBuffer
        /// @param capacity - Max capacity of the results buffer
        OptimizationParameterResultsBuffer(std::size_t capacity);

        /// Copy constructor
        /// @param other - OptimizationParameterResultsBuffer to make copy of
        OptimizationParameterResultsBuffer(const OptimizationParameterResultsBuffer &other);

        /// Move constructor
        /// @param other - OptimizationParameterResultsBuffer to swap with
        OptimizationParameterResultsBuffer(OptimizationParameterResultsBuffer &&other) noexcept;

        /// Assignment operator
        /// @param other - OptimizationParameterResultsBuffer to assign
        OptimizationParameterResultsBuffer &operator=(const OptimizationParameterResultsBuffer &other);

        /// Begin iterator of the buffer
        std::deque<OptimizationParameterResult>::const_reverse_iterator begin() const;

        /// End iterator of the buffer
        std::deque<OptimizationParameterResult>::const_reverse_iterator end() const;

        /// Utility method to push result to the start of ring buffer
        /// @param result - Result to store
        /// @see OptimizationParameterResult
        void push(const OptimizationParameterResult &result);

        /// Current size of the buffer (<= capacity)
        std::size_t size() const;
    };

    typedef double (*GetPreconditionerTimeFromSolverLambda)(const INMOST::Solver &solver);
    typedef double (*GetSolveTimeFromSolverLambda)(const INMOST::Solver &solver);

    class OptimizerInterface {
    protected:
        OptimizationParameterResultsBuffer results;
        OptimizationParametersSpace space;
    public:
        OptimizerInterface(const OptimizationParametersSpace &space, std::size_t buffer_capacity) : space(space), results(buffer_capacity) {};

        bool Solve(INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix, INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL,
                   GetPreconditionerTimeFromSolverLambda preconditioner_time = OptimizerInterface::DefaultGetPreconditionerTime,
                   GetSolveTimeFromSolverLambda solve_time = OptimizerInterface::DefaultGetSolveTime);

        void SaveResult(const OptimizationParameterPoints &points, const INMOST::Solver &solver, bool is_solved,
                        GetPreconditionerTimeFromSolverLambda preconditioner_time = OptimizerInterface::DefaultGetPreconditionerTime,
                        GetSolveTimeFromSolverLambda solve_time = OptimizerInterface::DefaultGetSolveTime);

        virtual OptimizationParameterPoints MakeOptimizationIteration(INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix,
                                                                      INMOST::Sparse::Vector &RHS) const = 0;

        const OptimizationParametersSpace &GetSpace() const;

        const OptimizationParameterResultsBuffer &GetResults() const;

        virtual ~OptimizerInterface() {};

        static double DefaultGetPreconditionerTime(const INMOST::Solver &solver);

        static double DefaultGetSolveTime(const INMOST::Solver &solver);
    };

};


#endif //INMOST_TTSP_H
