//
// Created by Dmitri Bagaev on 2019-07-09.
//

#ifndef INMOST_TTSP_CONFIGURATION_H
#define INMOST_TTSP_CONFIGURATION_H

#include <vector>
#include <string>
#include <inmost_optimizer.h>

namespace INMOST {

    namespace TTSP {

        class TTSPConfigurationParameterEntry {
        private:
            std::string               ref;
            std::string               name;
            double                    initial;
            std::vector<double>       values;
            OptimizationParameterType type;

            static void swap(TTSPConfigurationParameterEntry &left, TTSPConfigurationParameterEntry &right);

        public:
            TTSPConfigurationParameterEntry(const TTSPConfigurationParameterEntry &other);

            TTSPConfigurationParameterEntry(TTSPConfigurationParameterEntry &&other) noexcept;

            TTSPConfigurationParameterEntry &operator=(const TTSPConfigurationParameterEntry &other);

            TTSPConfigurationParameterEntry(const XMLReader::XMLTree &entry);

            const std::string &GetRef() const noexcept;

            const std::string &GetName() const noexcept;

            double GetInitial() const noexcept;

            const std::vector<double> &GetValues() const noexcept;

            OptimizationParameterType GetParameterType() const noexcept;

            ~TTSPConfigurationParameterEntry();

            friend std::ostream &operator<<(std::ostream &out, const TTSPConfigurationParameterEntry &c);
        };

        class TTSPConfigurationSolverPrefixEntry {
        private:
            bool                                         enabled;
            std::string                                  prefix;
            std::string                                  optimizer;
            INMOST::OptimizerVerbosityLevel              verbosity;
            std::size_t                                  buffer_capacity;
            std::vector<TTSPConfigurationParameterEntry> parameters;

            static void swap(TTSPConfigurationSolverPrefixEntry &left, TTSPConfigurationSolverPrefixEntry &right);

        public:
            TTSPConfigurationSolverPrefixEntry(const TTSPConfigurationSolverPrefixEntry &other);

            TTSPConfigurationSolverPrefixEntry(TTSPConfigurationSolverPrefixEntry &&other) noexcept;

            TTSPConfigurationSolverPrefixEntry &operator=(const TTSPConfigurationSolverPrefixEntry &other);

            TTSPConfigurationSolverPrefixEntry(const XMLReader::XMLTree &entry, const std::vector<TTSPConfigurationParameterEntry> &global);

            bool IsEnabled() const noexcept;

            const std::string &GetPrefix() const noexcept;

            const std::string &GetOptimizer() const noexcept;

            INMOST::OptimizerVerbosityLevel GetVerbosityLevel() const noexcept;

            std::size_t GetBufferCapacity() const noexcept;

            const std::vector<TTSPConfigurationParameterEntry> &GetParameters() const noexcept;

            ~TTSPConfigurationSolverPrefixEntry();

            friend std::ostream &operator<<(std::ostream &out, const TTSPConfigurationSolverPrefixEntry &c);
        };

        class TTSPConfigurationSolverEntry {
        private:
            bool                                            enabled;
            std::string                                     solver;
            std::vector<TTSPConfigurationSolverPrefixEntry> prefixes;

            static void swap(TTSPConfigurationSolverEntry &left, TTSPConfigurationSolverEntry &right);

        public:
            TTSPConfigurationSolverEntry(const TTSPConfigurationSolverEntry &other);

            TTSPConfigurationSolverEntry(TTSPConfigurationSolverEntry &&other) noexcept;

            TTSPConfigurationSolverEntry &operator=(const TTSPConfigurationSolverEntry &other);

            TTSPConfigurationSolverEntry(const XMLReader::XMLTree &entry, const std::vector<TTSPConfigurationParameterEntry> &global);

            bool IsEnabled() const noexcept;

            const std::string &GetSolver() const noexcept;

            const std::vector<TTSPConfigurationSolverPrefixEntry> &GetPrefixes() const noexcept;

            ~TTSPConfigurationSolverEntry();

            friend std::ostream &operator<<(std::ostream &out, const TTSPConfigurationSolverEntry &c);
        };

        class TTSPConfiguration {
        private:
            std::vector<TTSPConfigurationParameterEntry> global_parameter_refs;
            std::vector<TTSPConfigurationSolverEntry>    solvers;

            void parse(std::ifstream &stream, bool verbose = false);

            TTSPConfiguration(const TTSPConfiguration &other) = delete;

            TTSPConfiguration(TTSPConfiguration &&other) = delete;

            TTSPConfiguration &operator=(const TTSPConfiguration &other) = delete;

        public:
            TTSPConfiguration(std::ifstream &stream);

            TTSPConfiguration(const std::string &path);

            const std::vector<TTSPConfigurationSolverEntry> &GetSolvers() const noexcept;

            const TTSPConfigurationSolverPrefixEntry &FindForSolverAndPrefix(const std::string &solver_name, const std::string &solver_prefix) const;
        };

    }

}


#endif //INMOST_TTSP_CONFIGURATION_H
