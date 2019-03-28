#ifndef INMOST_INMOST_TTSP_H
#define INMOST_INMOST_TTSP_H

#include "inmost.h"

#if defined(USE_TTSP)

#include <vector>
#include <algorithm>
#include <string>
#include <deque>
#include <unordered_map>
#include <functional>

namespace TTSP {

    /// This class is used to defined a type of a parameter
    /// - DEFAULT  - parameter values should be used as is
    /// - EXPONENT - parameter values should be converted to 10^v
    enum class OptimizationParameterType {
        PARAMETER_TYPE_DEFAULT,
        PARAMETER_TYPE_EXPONENT
    };

    /// This class is used to store a slice of OptimizationParameters for some optimization parameter
    /// @see OptimizationParameter
    /// @see OptimizationParametersSpace
    class OptimizationParameterPoint {
    private:
        std::string name;
        double      value;

        OptimizationParameterPoint(const std::string &name, double value, OptimizationParameterType type);

        static void swap(OptimizationParameterPoint &left, OptimizationParameterPoint &right);
    public:
        OptimizationParameterPoint(const OptimizationParameterPoint &other);

        OptimizationParameterPoint(OptimizationParameterPoint &&other) noexcept;

        OptimizationParameterPoint &operator=(const OptimizationParameterPoint &other);

        const std::string &GetName() const noexcept;

        double GetValue() const noexcept;

        ~OptimizationParameterPoint();

        static double convert(double value, OptimizationParameterType type);

        friend class OptimizationParameters;
    };

    /// This class is used to store a slice of OptimizationParameters for all optimization parameters
    /// @see OptimizationParameter
    /// @see OptimizationParameterPoint
    /// @see OptimizationParametersSpace
    typedef std::vector<OptimizationParameterPoint> OptimizationParameterPoints;

    /// This class is used to define a range of possible values for some parameter
    /// Usage: OptimizationParameterRange range = std::make_pair(0.0, 1.0)
    typedef std::pair<double, double> OptimizationParameterRange;

    /// This class is used to define a solver parameter for optimization algortihms
    /// Usage:
    ///    1.  OptimizationParameter tau("tau", range, 0.1, 1e-3);
    ///    2.  OptimizationParameter kovl("kovl", { 0, 1, 2, 3 }, 0);
    class OptimizationParameter {
    private:
        std::string               name;          /// A name of a parameter
        std::vector<double>       values;        /// List of possible values for a parameter
        double                    default_value; /// Default value for this parameter
        OptimizationParameterType type;          /// Parameter type

        static void swap(OptimizationParameter &left, OptimizationParameter &right);

    public:
        /// Default constructor to define an OptimizationParameter with range of possible values and associated step
        /// @param name          - Name of a parameter
        /// @param range         - Range of possible values for this parameter
        /// @param step          - Step associated with the range (used to split up values in a list)
        /// @param default_value - Default value for this parameter
        /// @param type          - Parameter type
        OptimizationParameter(const std::string &name, const OptimizationParameterRange &range, double step,
                              double default_value, OptimizationParameterType type = OptimizationParameterType::PARAMETER_TYPE_DEFAULT);

        /// Default constructor to define an OptimizationParameter with list of possible values
        /// @param name          - Name of a parameter
        /// @param values        - List of possible values for this parameter
        /// @param default_value - Default value for this parameter
        /// @param type          - Parameter type
        OptimizationParameter(const std::string &name, const std::vector<double> &values,
                              double default_value, OptimizationParameterType type = OptimizationParameterType::PARAMETER_TYPE_DEFAULT);

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
        const std::string &GetName() const noexcept;

        /// Getter for values of a parameter
        const std::vector<double> &GetValues() const noexcept;

        /// Getter for minimal value of a parameter
        double GetMinimalValue() const noexcept;

        /// Getter for maximum value of a parameter
        double GetMaximumValue() const noexcept;

        double GetMinStrictBound(double v, double modifier = 1.0) const noexcept;

        double GetMaxStrictBound(double v, double modifier = 1.0) const noexcept;

        /// Getter for closest value of a parameter to specified values in creation stage
        double GetClosestTo(double to) const noexcept;

        /// Getter for closest value index of a parameter to specified values in creation stage
        std::size_t GetClosestIndexTo(double to) const noexcept;

        /// Getter for number of values of a parameter
        std::size_t GetValuesCount() const noexcept;

        /// Getter for default_value of a parameter
        double GetDefaultValue() const noexcept;

        /// Getter for type of a parameter
        OptimizationParameterType GetType() const noexcept;

        double ExtractValueFromPoint(const OptimizationParameterPoint &point) const noexcept;
    };

    /// This class is used to store an OptimizationParameter and associated current value
    /// @see OptimizationParameter
    typedef std::pair<OptimizationParameter, double> OptimizationParametersEntry;

    /// This class is used to store a list of OptimizationParameterEntry
    /// @see OptimizationParameter
    /// @see OptimizationParametersEntry
    typedef std::vector<OptimizationParametersEntry> OptimizationParameterEntries;

    typedef std::pair<std::size_t, double> OptimizationAlgorithmSuggestion;

    /// This class is used to store a result of optimization iterations with suggested points and changed parameter
    /// @see OptimizationParameter
    /// @see OptimizationParameterPoint
    /// @see OptimizationParameterPoints
    /// @see OptimizationParameters
    class OptimizationParametersSuggestion {
    private:
        std::size_t                 changed_index;
        double                      changed_value;
        OptimizationParameterPoints before;
        double                      metrics_before;
        OptimizationParameterPoints after;
    public:
        OptimizationParametersSuggestion(std::size_t changed_index, double changed_value,
                                         const OptimizationParameterPoints &before, double metrics_before,
                                         const OptimizationParameterPoints &after);

        std::size_t GetChangedParameterIndex() const noexcept;

        double GetChangedValue() const noexcept;

        const OptimizationParameterPoints &GetPointsBefore() const noexcept;

        double GetMetricsBefore() const noexcept;

        const OptimizationParameterPoints &GetPointsAfter() const noexcept;
    };

    /// This class is used to define an Optimization parameters entries and associated metrics with it
    /// Usage: OptimizationParametersSpace space("fcbiilu2", "prefix", parameters);
    class OptimizationParameters {
    private:
        OptimizationParameterEntries entries;            /// List of optimization parameters entries
        double                       metrics;            /// Metrics value for current optimization parameters

        static void swap(OptimizationParameters &left, OptimizationParameters &right);

    public:
        /// Default constructor to define an OptimizationParametersSpace
        /// @param entries       - List of optimization parameter entries
        /// @param metrics       - Default metrics value
        OptimizationParameters(const OptimizationParameterEntries &entries, double metrics);

        /// Copy constructor
        /// @param other - OptimizationParametersSpace to make copy of
        OptimizationParameters(const OptimizationParameters &other);

        /// Move constructor
        /// @param other - OptimizationParametersSpace to swap with
        OptimizationParameters(OptimizationParameters &&other) noexcept;

        /// Assignment operator
        /// @param other - OptimizationParametersSpace to assign
        OptimizationParameters &operator=(const OptimizationParameters &other);

        /// Getter for parameters of this space
        const OptimizationParameterEntries &GetParameterEntries() const noexcept;

        /// Getter for parameter entry of this space by index
        const OptimizationParametersEntry &GetParameterEntry(std::size_t index) const;

        /// Getter for parameter reference of this space by index
        const OptimizationParameter &GetParameter(std::size_t index) const;

        /// Getter for parameters count
        std::size_t Size() const;

        /// Getter for slice of parameters of this space
        const OptimizationParameterPoints GetPoints() const noexcept;

        /// Getter for metrics for current parameters of this space
        double GetMetrics() const noexcept;

        /// Getter for slice of parameters of this space with changed parameter
        const OptimizationParameterPoints GetPointsWithChangedParameter(const OptimizationParameter &parameter, double value) const noexcept;

        void Update(const std::vector<double> &update, double metrics);

        void Update(std::size_t index, double value, double metrics);
    };

    /// This class is used to store timings result for some optimization parameter points
    class OptimizationParameterResult {
    private:
        std::size_t                 changed_index;       /// Changed optimization parameter index
        double                      changed_value;       /// New changed value for 'changed' parameter
        OptimizationParameterPoints before;              /// Optimization parameter points before changing
        double                      metrics_before;      /// Optimization metrics before changing
        OptimizationParameterPoints after;               /// Optimization parameter points after changing
        double                      metrics_after;       /// Optimization metrics after changing
        bool                        is_good;             /// Is problem solved or not with new parameters

        static void swap(OptimizationParameterResult &left, OptimizationParameterResult &right);

    public:
        /// Default constructor to define an OptimizationParameterResult
        OptimizationParameterResult(std::size_t changed_index, double changed_value,
                                    const OptimizationParameterPoints &before, const OptimizationParameterPoints &after,
                                    double metrics_before, double metrics_after, bool is_good);

        /// Copy constructor
        /// @param other - OptimizationParameterResult to make copy of
        OptimizationParameterResult(const OptimizationParameterResult &other);

        /// Move constructor
        /// @param other - OptimizationParameterResult to swap with
        OptimizationParameterResult(OptimizationParameterResult &&other) noexcept;

        /// Assignment operator
        /// @param other - OptimizationParameterResult to assign
        OptimizationParameterResult &operator=(const OptimizationParameterResult &other);

        /// Getter for points of this result before changing
        const OptimizationParameterPoints &GetPointsBefore() const noexcept;

        /// Getter for points of this result after changing
        const OptimizationParameterPoints &GetPointsAfter() const noexcept;

        /// Getter for metrics before changing
        double GetMetricsBefore() const noexcept;

        /// Getter for metrics after changing
        double GetMetricsAfter() const noexcept;

        /// Getter for is solved status
        bool IsGood() const noexcept;

        /// Getter for changed parameter
        std::size_t GetChangedParameterIndex() const noexcept;

        /// Getter for new changed value
        double GetChangedValue() const noexcept;
    };

    /// This class is used to store last N solve results and used parameters on it
    /// Implemented as a circular (ring) container
    class OptimizationParameterResultsBuffer {
    private:
        std::size_t                             capacity;                            /// Max capacity of the results buffer
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
        std::deque<OptimizationParameterResult>::const_reverse_iterator cbegin() const noexcept;

        const OptimizationParameterResult &at(std::size_t index) const;

        /// End iterator of the buffer
        std::deque<OptimizationParameterResult>::const_reverse_iterator cend() const noexcept;

        /// Utility method to push result to the start of ring buffer
        /// @param result - Result to store
        /// @see OptimizationParameterResult
        void push(const OptimizationParameterResult &result);

        /// Current size of the buffer (<= capacity)
        std::size_t size() const noexcept;

        /// Getter for emptiness status
        bool IsEmpty() const noexcept;

        /// Getter for last result successfulness status
        bool IsLastResultSuccessful() const noexcept;

        bool IsSuccessfulResultExist() const noexcept;

        /// Getter for last successful result
        const OptimizationParameterResult &GetLastSuccessfulResult() const;

        /// Getter for last unique point entries
        const std::vector<OptimizationParameterResult> GetLastUniqueEntries(std::size_t max) const;
    };

    typedef std::map<std::string, std::string> OptimizerProperties;

    enum class OptimizerVerbosityLevel {
        Level0,
        Level1,
        Level2,
        Level3
    };

    enum class OptimizerRestartStrategy {
        RESTART_STRATEGY_NO_RESTART,
        RESTART_STRATEGY_WITH_BEST
    };

    typedef std::pair<bool, double> OptimizationFunctionInvokeResult;

    class OptimizerInterface {
    private:
        int mpi_rank;

        std::string             name;
        OptimizerVerbosityLevel verbosity = OptimizerVerbosityLevel::Level0;

        OptimizerRestartStrategy restart_strategy = OptimizerRestartStrategy::RESTART_STRATEGY_NO_RESTART;
        std::size_t              max_fails        = std::numeric_limits<std::size_t>::max();
        std::size_t              fails_count      = 0;

        std::deque<OptimizationParametersSuggestion> suggestions;

        void RestartWithBestStrategy();

    protected:
        OptimizationParameterResultsBuffer results;
        OptimizationParameters             parameters;
        const OptimizerProperties          properties;

        virtual bool UpdateSpaceWithLatestResults();

        virtual OptimizationAlgorithmSuggestion AlgorithmMakeSuggestion(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                             const OptimizationParameterPoints &,
                                                                                                                             void *)> &invoke, void *data) const = 0;

    public:
        OptimizerInterface(const std::string &name, const OptimizationParameters &parameters, const OptimizerProperties &properties, std::size_t buffer_capacity) :
                name(name), parameters(parameters), properties(properties), results(buffer_capacity), mpi_rank(0) {
#ifdef USE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
        };

        OptimizationParametersSuggestion Suggest(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &, const OptimizationParameterPoints &,
                                                                                                      void *)> &invoke = [](const OptimizationParameterPoints &,
                                                                                                                            const OptimizationParameterPoints &,
                                                                                                                            void *) -> OptimizationFunctionInvokeResult {
            return std::make_pair(false, 0);
        }, void *data = nullptr);

        void SaveResult(std::size_t changed_index, double changed_value,
                        const OptimizationParameterPoints &before, double metrics_before,
                        const OptimizationParameterPoints &after, double metrics_after, bool is_good);

        void SaveResult(std::size_t changed_index, double changed_value, const OptimizationParameterPoints &after, double metrics_after, bool is_good);

        void SaveResult(const OptimizationParametersSuggestion &suggestion, double metrics, bool is_good);

        const OptimizationParameterResultsBuffer &GetResults() const noexcept;

        const OptimizationParameterPoints GetPoints() const noexcept;

        void SetVerbosityLevel(OptimizerVerbosityLevel level) noexcept;

        OptimizerVerbosityLevel GetVerbosityLevel() const noexcept;

        bool HasProperty(const std::string &name) const noexcept;

        const std::string &GetProperty(const std::string &name) const;

        const OptimizationParametersSuggestion &GetLastSuggestion() const;

        void SetRestartStrategy(OptimizerRestartStrategy strategy, std::size_t max_fails) noexcept;

        const std::string &GetName() const noexcept;

        virtual ~OptimizerInterface() {};
    };

    typedef std::unordered_map<std::string, OptimizerInterface *> SavedOptimizersMap;

    class Optimizers {
    private:
        static SavedOptimizersMap optimizers;
    public:
        static bool IsOptimizerAvailable(const std::string &type);

        static std::vector<std::string> GetAvailableOptimizers();

        static OptimizerInterface *GetOptimizer(const std::string &name, const std::string &type,
                                                const OptimizationParameters &parameters, const OptimizerProperties &properties, std::size_t buffer_capacity);

        static void SaveOptimizerOrReplace(const std::string &name, const std::string &type,
                                           const OptimizationParameters &parameters, const OptimizerProperties &properties, std::size_t buffer_capacity);

        static OptimizerInterface *GetSavedOptimizer(const std::string &name);
    };

};

#endif

#endif //INMOST_INMOST_TTSP_H
