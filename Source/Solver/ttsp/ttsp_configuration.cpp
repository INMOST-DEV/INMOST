//
// Created by Dmitri Bagaev on 2019-07-09.
//

#include "../../Misc/utils.h"
#include "ttsp_configuration.h"

#if defined(USE_OPTIMIZER)

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
                for (std::vector<XMLReader::XMLTree>::const_iterator value = entry.children.begin(); value != entry.children.end(); ++value) {
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
            out << "      -------------------------------" << std::endl;
            out << "      name:     " << c.name << std::endl;
            out << "      ref:      " << c.ref << std::endl;
            out << "      type:     " << (c.type == OptimizationParameterType::PARAMETER_TYPE_EXPONENT ? "log10" : "default") << std::endl;
            out << "      initial:  " << c.initial << std::endl;
            out << "      values:   " << "[";
            std::for_each(c.values.cbegin(), c.values.cend(), [&out](double v) {
                out << v << ",";
            });
            out << "]" << std::endl;
            return out;
        }

        const std::string &TTSPConfigurationParameterEntry::GetRef() const noexcept {
            return ref;
        }

        const std::string &TTSPConfigurationParameterEntry::GetName() const noexcept {
            return name;
        }

        double TTSPConfigurationParameterEntry::GetInitial() const noexcept {
            return initial;
        }

        const std::vector<double> &TTSPConfigurationParameterEntry::GetValues() const noexcept {
            return values;
        }

        OptimizationParameterType TTSPConfigurationParameterEntry::GetParameterType() const noexcept {
            return type;
        }

        void TTSPConfigurationSolverPrefixEntry::swap(TTSPConfigurationSolverPrefixEntry &left, TTSPConfigurationSolverPrefixEntry &right) {
            std::swap(left.enabled, right.enabled);
            std::swap(left.prefix, right.prefix);
            std::swap(left.optimizer, right.optimizer);
            std::swap(left.verbosity, right.verbosity);
            std::swap(left.buffer_capacity, right.buffer_capacity);
            std::swap(left.parameters, right.parameters);
        }

        TTSPConfigurationSolverPrefixEntry::TTSPConfigurationSolverPrefixEntry(const TTSPConfigurationSolverPrefixEntry &other) :
                enabled(other.enabled), prefix(other.prefix), optimizer(other.optimizer), verbosity(other.verbosity), buffer_capacity(other.buffer_capacity),
                parameters(other.parameters) {}

        TTSPConfigurationSolverPrefixEntry::TTSPConfigurationSolverPrefixEntry(TTSPConfigurationSolverPrefixEntry &&other) noexcept {
            TTSPConfigurationSolverPrefixEntry::swap(*this, other);
        }

        TTSPConfigurationSolverPrefixEntry &TTSPConfigurationSolverPrefixEntry::operator=(const TTSPConfigurationSolverPrefixEntry &other) {
            TTSPConfigurationSolverPrefixEntry tmp(other);
            TTSPConfigurationSolverPrefixEntry::swap(*this, tmp);
            return *this;
        }

        TTSPConfigurationSolverPrefixEntry::TTSPConfigurationSolverPrefixEntry(const XMLReader::XMLTree &entry, const std::vector<TTSPConfigurationParameterEntry> &global) {
            this->prefix = entry.tag.name;

            int num_attr = entry.NumAttrib();

            int optimizeri = entry.FindAttrib("optimizer");
            if (optimizeri != num_attr) {
                this->optimizer = entry.GetAttrib(optimizeri).value;
            } else {
                throw "Bad parameter in XML solver optimization configuration file: Missing 'optimizer' property in prefix configuration!";
            }

            int enabledi = entry.FindAttrib("enabled");
            if (enabledi != num_attr) {
                this->enabled = entry.GetAttrib(enabledi).value == "true";
            } else {
                this->enabled = true;
            }

            int verbosityi = entry.FindAttrib("verbosity");
            if (verbosityi != num_attr) {
                std::string v = entry.GetAttrib(verbosityi).value;
                if (v == "0") {
                    this->verbosity = OptimizerVerbosityLevel::Level0;
                } else if (v == "1") {
                    this->verbosity = OptimizerVerbosityLevel::Level1;
                } else if (v == "2") {
                    this->verbosity = OptimizerVerbosityLevel::Level2;
                } else if (v == "3") {
                    this->verbosity = OptimizerVerbosityLevel::Level3;
                } else {
                    throw "Bad parameter in XML solver optimization configuration file: Invalid 'verbosity' property in prefix configuration!";
                }
            } else {
                this->verbosity = OptimizerVerbosityLevel::Level1;
            }

            std::size_t buffer_capacityi = entry.FindAttrib("buffer");
            if (buffer_capacityi != num_attr) {
                this->buffer_capacity = static_cast<std::size_t>(atol(entry.GetAttrib(buffer_capacityi).value.c_str()));
            } else {
                this->buffer_capacity = 15;
            }

            for (std::vector<XMLReader::XMLTree>::const_iterator param = entry.children.begin(); param != entry.children.end(); ++param) {

                if ((*param).tag.name == "parameter") {
                    if ((*param).FindAttrib("ref") != (*param).NumAttrib()) {
                        std::string ref  = (*param).GetAttrib("ref");
                        auto        fref = std::find_if(global.cbegin(), global.cend(), [&ref](const TTSPConfigurationParameterEntry &e) {
                            return e.GetRef() == ref;
                        });
                        if (fref != global.cend()) {
                            TTSPConfigurationParameterEntry e = TTSPConfigurationParameterEntry(*fref);
                            this->parameters.emplace_back(e);
                        } else {
                            throw "Bad parameter in XML solver optimization configuration file: Cannot find global parameter by ref: " + ref;
                        }
                    } else {
                        TTSPConfigurationParameterEntry e = TTSPConfigurationParameterEntry(*param);
                        this->parameters.emplace_back(e);
                    }
                } else {
                    throw "Bad parameter in XML solver optimization configuration file: Invalid 'parameter' child tag in prefix configuration!";
                }

            }

        }

        bool TTSPConfigurationSolverPrefixEntry::IsEnabled() const noexcept {
            return enabled;
        }

        const std::string &TTSPConfigurationSolverPrefixEntry::GetPrefix() const noexcept {
            return prefix;
        }

        const std::string &TTSPConfigurationSolverPrefixEntry::GetOptimizer() const noexcept {
            return optimizer;
        }

        INMOST::OptimizerVerbosityLevel TTSPConfigurationSolverPrefixEntry::GetVerbosityLevel() const noexcept {
            return verbosity;
        }

        std::size_t TTSPConfigurationSolverPrefixEntry::GetBufferCapacity() const noexcept {
            return buffer_capacity;
        }

        const std::vector<TTSPConfigurationParameterEntry> &TTSPConfigurationSolverPrefixEntry::GetParameters() const noexcept {
            return parameters;
        }

        TTSPConfigurationSolverPrefixEntry::~TTSPConfigurationSolverPrefixEntry() {}

        std::ostream &operator<<(std::ostream &out, const TTSPConfigurationSolverPrefixEntry &c) {
            out << "    =================================" << std::endl;
            out << "    enabled:    " << c.enabled << std::endl;
            out << "    prefix:     " << c.prefix << std::endl;
            out << "    optimizer:  " << c.optimizer << std::endl;
            out << "    parameters: " << std::endl;
            std::for_each(c.parameters.cbegin(), c.parameters.cend(), [&out](const TTSPConfigurationParameterEntry &e) {
                out << e;
            });
            out << std::endl;
            return out;
        }

        void TTSPConfigurationSolverEntry::swap(TTSPConfigurationSolverEntry &left, TTSPConfigurationSolverEntry &right) {
            std::swap(left.enabled, right.enabled);
            std::swap(left.solver, right.solver);
            std::swap(left.prefixes, right.prefixes);
        }

        TTSPConfigurationSolverEntry::TTSPConfigurationSolverEntry(const TTSPConfigurationSolverEntry &other) :
                enabled(other.enabled), solver(other.solver), prefixes(other.prefixes) {}

        TTSPConfigurationSolverEntry::TTSPConfigurationSolverEntry(TTSPConfigurationSolverEntry &&other) noexcept {
            TTSPConfigurationSolverEntry::swap(*this, other);
        }

        TTSPConfigurationSolverEntry &TTSPConfigurationSolverEntry::operator=(const TTSPConfigurationSolverEntry &other) {
            TTSPConfigurationSolverEntry tmp(other);
            TTSPConfigurationSolverEntry::swap(*this, tmp);
            return *this;
        }

        TTSPConfigurationSolverEntry::TTSPConfigurationSolverEntry(const XMLReader::XMLTree &entry, const std::vector<TTSPConfigurationParameterEntry> &global) {
            this->solver = entry.tag.name;

            int num_attr = entry.NumAttrib();

            int enabledi = entry.FindAttrib("enabled");
            if (enabledi != num_attr) {
                this->enabled = entry.GetAttrib(enabledi).value == "true";
            } else {
                this->enabled = true;
            }

            for (std::vector<XMLReader::XMLTree>::const_iterator prefix = entry.children.begin(); prefix != entry.children.end(); ++prefix) {
                TTSPConfigurationSolverPrefixEntry p = TTSPConfigurationSolverPrefixEntry(*prefix, global);
                this->prefixes.emplace_back(p);
            }

        }

        bool TTSPConfigurationSolverEntry::IsEnabled() const noexcept {
            return enabled;
        }

        const std::string &TTSPConfigurationSolverEntry::GetSolver() const noexcept {
            return solver;
        }

        const std::vector<TTSPConfigurationSolverPrefixEntry> &TTSPConfigurationSolverEntry::GetPrefixes() const noexcept {
            return prefixes;
        }

        TTSPConfigurationSolverEntry::~TTSPConfigurationSolverEntry() {}

        std::ostream &operator<<(std::ostream &out, const TTSPConfigurationSolverEntry &c) {
            out << "  enabled:    " << c.enabled << std::endl;
            out << "  solver:     " << c.solver << std::endl;
            out << "  prefixes:   " << std::endl;
            std::for_each(c.prefixes.cbegin(), c.prefixes.cend(), [&out](const TTSPConfigurationSolverPrefixEntry &e) {
                out << e;
            });
            out << std::endl;
            return out;
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
                    } else {
                        TTSPConfigurationSolverEntry s = TTSPConfigurationSolverEntry(*entry, global_parameter_refs);
                        this->solvers.emplace_back(s);
                        if (verbose) {
                            std::cout << "Parsed optimization configuration for solver: " << std::endl;
                            std::cout << s << std::endl;
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

        const std::vector<TTSPConfigurationSolverEntry> &TTSPConfiguration::GetSolvers() const noexcept {
            return solvers;
        }

        const TTSPConfigurationSolverPrefixEntry &TTSPConfiguration::FindForSolverAndPrefix(const std::string &solver_name, const std::string &solver_prefix) const {
            std::vector<TTSPConfigurationSolverEntry>::const_iterator its = std::find_if(this->solvers.cbegin(), this->solvers.cend(),
                                                                                         [&solver_name](const TTSPConfigurationSolverEntry &s) {
                                                                                             return s.GetSolver() == solver_name;
                                                                                         });
            if (its == this->solvers.cend()) {
                std::cout << "Unable to find optimization parameters entry for solver: " << solver_name << std::endl;
                throw "Unable to find optimization parameters entry for solver: " + solver_name;
            }

            std::vector<TTSPConfigurationSolverPrefixEntry>::const_iterator ip = std::find_if((*its).GetPrefixes().cbegin(), (*its).GetPrefixes().cend(),
                                                                                              [&solver_prefix](const TTSPConfigurationSolverPrefixEntry &p) {
                                                                                                  return p.GetPrefix() == solver_prefix;
                                                                                              });

            if (ip == (*its).GetPrefixes().cend()) {
                std::cout << "Unable to find optimization parameters entry for prefix: " << solver_prefix << std::endl;
                throw "Unable to find optimization parameters entry for prefix: " + solver_prefix;
            }

            return (*ip);
        }
    }
}

#endif