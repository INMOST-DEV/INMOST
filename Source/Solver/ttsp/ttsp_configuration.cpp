//
// Created by Dmitri Bagaev on 2019-07-09.
//

#include <Source/Misc/utils.h>
#include "ttsp_configuration.h"

namespace INMOST {

    namespace TTSP {

        void TTSPConfigurationParameterEntry::swap(INMOST::TTSP::TTSPConfigurationParameterEntry &left, INMOST::TTSP::TTSPConfigurationParameterEntry &right) {
            std::swap(left.ref, right.ref);
            std::swap(left.name, right.name);
            std::swap(left.initial, right.initial);
            std::swap(left.values, right.values);
            std::swap(left.type, right.type);
        }

        TTSPConfigurationParameterEntry::TTSPConfigurationParameterEntry(const TTSPConfigurationParameterEntry &other) :
                ref(other.ref), name(other.name), initial(other.initial), values(other.values),
                type(other.type) {
        }

        TTSPConfigurationParameterEntry::TTSPConfigurationParameterEntry(TTSPConfigurationParameterEntry &&other) noexcept {
            TTSPConfigurationParameterEntry::swap(*this, other);
        }

        TTSPConfigurationParameterEntry &TTSPConfigurationParameterEntry::operator=(const TTSPConfigurationParameterEntry &other) {
            TTSPConfigurationParameterEntry tmp(other);
            TTSPConfigurationParameterEntry::swap(*this, tmp);
            return *this;
        }

        TTSPConfigurationParameterEntry::TTSPConfigurationParameterEntry(const XMLReader::XMLTree &entry) {

            int num_attrib = entry.NumAttrib();
            int namei      = entry.FindAttrib("name");
            if (namei == num_attrib) {
                throw "Bad parameter in XML solver optimization configuration file: Missing 'name' property!";
            }

            int initiali = entry.FindAttrib("initial");
            if (initiali == num_attrib) {
                throw "Bad parameter in XML solver optimization configuration file: Missing 'initial' property!";
            }

            this->name    = entry.GetAttrib(namei).value;
            this->initial = atof(entry.GetAttrib(initiali).value.c_str());

            int refi = entry.FindAttrib("ref");
            if (refi != num_attrib) {
                this->ref = entry.GetAttrib(refi).value;
            }

            int typei = entry.FindAttrib("type");
            if (typei != num_attrib) {
                this->type = (entry.GetAttrib(typei).value == "log10" ? OptimizationParameterType::PARAMETER_TYPE_EXPONENT
                                                                      : OptimizationParameterType::PARAMETER_TYPE_DEFAULT);
            } else {
                this->type = OptimizationParameterType::PARAMETER_TYPE_DEFAULT;
            }


            if (entry.children.size() != 0) {
                for (std::vector<const XMLReader::XMLTree>::iterator value = entry.children.begin(); value != entry.children.end(); ++value) {
                    if ((*value).tag.name == "value") {
                        this->values.push_back(atof((*value).contents.c_str()));
                    } else {
                        throw "Bad parameter in XML solver optimization configuration file: Missing 'value' tag!";
                    }
                }
            } else {
                double from = 0.0;
                double to   = 0.0;
                double step = 0.0;

                int fromi = entry.FindAttrib("from");
                if (fromi != num_attrib) {
                    from = atof(entry.GetAttrib(fromi).value.c_str());
                }

                int toi = entry.FindAttrib("to");
                if (toi != num_attrib) {
                    to = atof(entry.GetAttrib(toi).value.c_str());
                }

                int stepi = entry.FindAttrib("step");
                if (stepi != num_attrib) {
                    step = atof(entry.GetAttrib(stepi).value.c_str());
                }

                if (from >= to) {
                    throw "Bad parameter in XML solver optimization configuration file: Invalid 'from' and 'to' properties!";
                } else if (step <= 0) {
                    throw "Bad parameter in XML solver optimization configuration file: Invalid 'step' property!";
                } else {
                    double v = from;
                    while (v <= to) {
                        this->values.push_back(v);
                        v += step;
                    }
                }

            }
        }

        TTSPConfigurationParameterEntry::~TTSPConfigurationParameterEntry() {}

        std::ostream &operator<<(std::ostream &out, const TTSPConfigurationParameterEntry &c) {
            out << "name:\t\t" << c.name << std::endl;
            if (!c.ref.empty()) out << "ref:\t\t" << c.ref << std::endl;
            out << "type:\t\t" << (c.type == OptimizationParameterType::PARAMETER_TYPE_EXPONENT ? "log10" : "default") << std::endl;
            out << "init:\t\t" << c.initial << std::endl;
            out << "values:\t\t[";
            std::for_each(c.values.cbegin(), c.values.cend(), [&out](double v) {
                out << v << ",";
            });
            out << "]" << std::endl;
            return out;
        }

        const std::string &TTSPConfigurationParameterEntry::getRef() const noexcept {
            return ref;
        }

        const std::string &TTSPConfigurationParameterEntry::getName() const noexcept {
            return name;
        }

        double TTSPConfigurationParameterEntry::getInitial() const noexcept {
            return initial;
        }

        const std::vector<double> &TTSPConfigurationParameterEntry::getValues() const noexcept {
            return values;
        }

        OptimizationParameterType TTSPConfigurationParameterEntry::getParameterType() const noexcept {
            return type;
        }

        void TTSPConfiguration::swap(TTSPConfiguration &left, TTSPConfiguration &right) {
            std::swap(left.global_parameter_refs, right.global_parameter_refs);
        }

        void TTSPConfiguration::parse(std::ifstream &stream, bool verbose) {
            try {
                INMOST::XMLReader          reader("", stream);
                INMOST::XMLReader::XMLTree root = reader.ReadXML();

                if (root.tag.name != "SolverOptimization") {
                    if (verbose) {
                        std::cout << __FILE__ << ":" << __LINE__ << ": Bad XML solver optimization configuration file: Missing open SolverOptimization tag!" << std::endl;
                    }
                    return;
                }

                global_parameter_refs.clear();

                for (INMOST::xml_reader_tree_iterator_t entry = root.children.begin(); entry != root.children.end(); ++entry) {
                    std::string entry_type = INMOST::string_to_lower((*entry).tag.name);

                    if (entry_type == "parameter") {
                        TTSPConfigurationParameterEntry e = TTSPConfigurationParameterEntry(*entry);
                        global_parameter_refs.emplace_back(e);
                        if (verbose) {
                            std::cout << "Parsed optimization parameter: " << std::endl << e << std::endl;
                        }
                    }
                }


            } catch (const char *e) {
                std::cout << "Error while parsing solver optimization configuration file: " << e << std::endl;
                std::exit(-1);
            } catch (const std::string &e) {
                std::cout << "Error while parsing solver optimization configuration file: " << e << std::endl;
                std::exit(-1);
            } catch (...) {
                std::cout << "Error while parsing solver optimization configuration file" << std::endl;
                std::exit(-1);
            }

        }

        TTSPConfiguration::TTSPConfiguration(std::ifstream &stream) {
            int rank = MPIGetRank();
            parse(stream, rank == 0);
        }

        TTSPConfiguration::TTSPConfiguration(const std::string &path) {
            int rank = MPIGetRank();

            if (path.empty()) {
                if (rank == 0) {
                    std::cout << __FILE__ << ":" << __LINE__ << ": XML solver optimization configuration file path is empty." << std::endl;
                }
                return;
            }

            std::ifstream stream;
            stream.open(path.c_str());

            if (stream.fail()) {
                if (rank == 0) {
                    std::cout << __FILE__ << ":" << __LINE__ << ": XML solver optimization configuration file path " << path << " not found." << std::endl;
                }
                return;
            }

            parse(stream, rank == 0);
        }

    }
}