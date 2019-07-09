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
            std::string                         ref;
            std::string                         name;
            double                              initial;
            std::vector<double>                 values;
            OptimizationParameterType           type;

            static void swap(TTSPConfigurationParameterEntry &left, TTSPConfigurationParameterEntry &right);

        public:
            TTSPConfigurationParameterEntry(const TTSPConfigurationParameterEntry &other);

            TTSPConfigurationParameterEntry(TTSPConfigurationParameterEntry &&other) noexcept;

            TTSPConfigurationParameterEntry &operator=(const TTSPConfigurationParameterEntry &other);

            TTSPConfigurationParameterEntry(const XMLReader::XMLTree &entry);

            const std::string &getRef() const noexcept;

            const std::string &getName() const noexcept;

            double getInitial() const noexcept;

            const std::vector<double> &getValues() const noexcept;

            OptimizationParameterType getParameterType() const noexcept;

            ~TTSPConfigurationParameterEntry();

            friend std::ostream & operator << (std::ostream &out, const TTSPConfigurationParameterEntry &c);
        };

        class TTSPConfiguration {
        private:
            std::vector<TTSPConfigurationParameterEntry> global_parameter_refs;

            static void swap(TTSPConfiguration &left, TTSPConfiguration &right);

            void parse(std::ifstream &stream, bool verbose = false);

            TTSPConfiguration(const TTSPConfiguration &other) = delete;
            TTSPConfiguration(TTSPConfiguration &&other) = delete;
            TTSPConfiguration &operator=(const TTSPConfiguration &other) = delete;
        public:
            TTSPConfiguration(std::ifstream &stream);
            TTSPConfiguration(const std::string &path);
        };

    }

}


#endif //INMOST_TTSP_CONFIGURATION_H
