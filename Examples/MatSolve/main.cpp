#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "inmost.h"
using namespace INMOST;

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

void removeSpaces(char* source) {
    char* i = source;
    char* j = source;
    while(*j != 0)
    {
        *i = *j++;
        if(*i != ' ')
            i++;
    }
    *i = 0;
}

enum OptionType {
    REAL,
    ENUM
};

class InnerOption {
public:

    InnerOption(const std::string &name, const std::string &value, OptionType type) : name(name), value(value),
                                                                                      type(type) { }
    std::string name;
    std::string value;
    OptionType type;
};

class InnerOptions {
public:
    std::vector<InnerOption *> options;
};

InnerOptions *parseInnerDatabaseOptions(char *optionsFile) {
    FILE *databaseFile = fopen(optionsFile, "r");
    if (!databaseFile) {
        std::cout << "Inner options file not found" << std::endl;
        return nullptr;
    }
    InnerOptions *options = new InnerOptions();
    char *tmp = (char *) calloc(256, sizeof(char));
    char *parameterName = (char *) calloc(128, sizeof(char));
    char *parameterValue = (char *) calloc(128, sizeof(char));
    char *type = (char *) calloc(36, sizeof(char));
    while (!feof(databaseFile) && fgets(tmp, 256, databaseFile)) {
        //removeSpaces(tmp);
        char *line = tmp;
        //Comment str
        if (line[0] == '#') continue;
        //First 4 chars is 'real' or 'enum'
        bool isReal = false, isEnum = false;
        sscanf(line, "%s %s %s", type, parameterName, parameterValue);
        if (strncmp(type, "real", 4) == 0) {
            //fprintf(stdout, "Real parameter:");
            isReal = true;
        } else if (strncmp(type, "enum", 4) == 0) {
            //fprintf(stdout, "Enum parameter:");
            isEnum = true;
        } else {
            fprintf(stderr, "Skipping bad line: %s", line);
            continue;
        }
        options->options.push_back(new InnerOption(std::string(parameterName), std::string(parameterValue), isReal ? REAL : ENUM));

    }
    free(parameterValue);
    free(parameterName);
    free(tmp);
    free(type);
    return options;
}

char *findInnerOptions(const char *databaseFilePath) {
    FILE *databaseFile = fopen(databaseFilePath, "r");
    if (databaseFile == NULL) {
        fprintf(stderr, "Database file not found\n");
        #if defined(USE_MPI)
        MPI_Finalize();
        #endif
        exit(1);
    }
    char *tmp = (char *) calloc(256, sizeof(char));
    char *parameterName = (char *) calloc(64, sizeof(char));
    char *parameterValue = (char *) calloc(64, sizeof(char));
    char *fileName = NULL;
    while (!feof(databaseFile) && fgets(tmp, 256, databaseFile)) {
        removeSpaces(tmp);
        char *line = tmp;
        if (strncmp(line, "inner", 5) != 0) continue;
        line += 6;
        size_t fileNameLength = 0;
        while (!isspace(*(line + fileNameLength))) fileNameLength += 1;
        fileName = (char *) calloc(fileNameLength + 1, sizeof(char));
        memcpy(fileName, line, fileNameLength);
        fileName[fileNameLength] = 0;
    }
    free(tmp);
    return fileName;
}

int main(int argc, char ** argv) {
    int processRank = 0, processorsCount = 1;

    #if defined(USE_MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);  // Get the rank of the current process
    MPI_Comm_size(MPI_COMM_WORLD, &processorsCount ); // Get the total number of processors used
    #endif

    if (processRank == 0) {
        std::cout << "Starting MatSolve2" << std::endl;
    }

    {
        std::string matrixFileName = "";
        std::string vectorBFileName = "";
        std::string vectorXFileName = "";
        std::string parametersFileName = "";
        int solverType = 0;

        bool matrixFound = false;
        bool vectorBFound = false;
        bool vectorXFound = false;
        bool parametersFound = false;
        bool typeFound = false;

        //Parse argv parameters
        if (argc == 1) goto helpMessage;
        int i;
        for (i = 1; i < argc; i++) {
            //Print help message and exit
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                helpMessage:
                if (processRank == 0) {
                    std::cout << "Help message: " << std::endl;
                    std::cout << "Command line options: " << std::endl;
                    std::cout << "Required: " << std::endl;
                    std::cout << "-m, --matrix <Matrix file name>" << std::endl;
                    std::cout << "Optional: " << std::endl;
                    std::cout << "-b, --bvector <RHS vector file name>" << std::endl;
                    std::cout << "-x, --xvector <X vector file name>" << std::endl;
                    std::cout << "-p, --parameters <Solver parameters file name>" << std::endl;
                    std::cout << "-t, --type <Solver type index>" << std::endl;
                    std::cout << "    0: INNER_ILU2 " << std::endl;
                    std::cout << "    1: INNER_DDPQILUC " << std::endl;
                    std::cout << "    2: INNER_MPTILUC " << std::endl;
                    std::cout << "    3: INNER_MPTILU2 " << std::endl;
                    std::cout << "    4: Trilinos_Aztec " << std::endl;
                    std::cout << "    5: Trilinos_Belos " << std::endl;
                    std::cout << "    6: Trilinos_ML " << std::endl;
                    std::cout << "    7: Trilinos_Ifpack " << std::endl;
                    std::cout << "    8: PETSc " << std::endl;
                    std::cout << "    9: ANI " << std::endl;
                    std::cout << "    10: FCBIILU2 " << std::endl;
                    std::cout << "    11: K3BIILU2 " << std::endl << std::endl;
                    std::cout << "-h, --help - print this message." << std::endl;
                }
                #if defined(USE_MPI)
                MPI_Finalize();
                #endif
                return 0;
            }
            //Matrix file name found with -m or --matrix options
            if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--matrix") == 0) {
                if (processRank == 0) {
                    std::cout << "Matrix file found: " << argv[i + 1] << std::endl;
                }
                matrixFound = true;
                matrixFileName = std::string(argv[i + 1]);
                i++;
                continue;
            }
            //B vector file name found with -b or --bvector options
            if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--bvector") == 0) {
                if (processRank == 0) {
                    std::cout << "B vector file found: " << argv[i + 1] << std::endl;
                }
                vectorBFound = true;
                vectorBFileName = std::string(argv[i + 1]);
                i++;
                continue;
            }
            //X vector file name found with -x or --xvector options
            if (strcmp(argv[i], "-x") == 0 || strcmp(argv[i], "--xvector") == 0) {
                if (processRank == 0) {
                    std::cout << "X vector file found: " << argv[i + 1] << std::endl;
                }
                vectorXFound = true;
                vectorXFileName = std::string(argv[i + 1]);
                i++;
                continue;
            }
            //Parameters file name found with -p or --parameters options
            if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--parameters") == 0) {
                if (processRank == 0) {
                    std::cout << "Solver parameters file found: " << argv[i + 1] << std::endl;
                }
                parametersFound = true;
                parametersFileName = std::string(argv[i + 1]);
                i++;
                continue;
            }
            //Solver type found with -t ot --type options
            if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--type") == 0) {
                if (processRank == 0) {
                    std::cout << "Solver type index found: " << argv[i + 1] << std::endl;
                }
                typeFound = true;
                solverType = atoi(argv[i + 1]);
                i++;
                continue;
            }
        }

        if (!matrixFound) {
            if (processRank == 0) {
                std::cout <<
                "Matrix not found, you can specify matrix file name using -m or --matrix options, otherwise specify -h option to see all options, exiting...";
            }
            return -1;
        }

        if (!typeFound) {
            if (processRank == 0) {
                std::cout <<
                "Solver type not found in command line, you can specify solver type with -t or --type option, using INNER_ILU2 solver by default." <<
                std::endl;
            }
        }

        if (!vectorBFound) {
            if (processRank == 0) {
                std::cout <<
                "B vector not found, you can specify b vector file name with -b or --bvector option, using identity vector by default." <<
                std::endl;
            }
        }


        Solver::Type type;
        switch (solverType) {
            case 0:
                type = Solver::INNER_ILU2;
                break;
            case 1:
                type = Solver::INNER_DDPQILUC;
                break;
            case 2:
                type = Solver::INNER_MPTILUC;
                break;
            case 3:
                type = Solver::INNER_MPTILU2;
                break;
            case 4:
                type = Solver::Trilinos_Aztec;
                break;
            case 5:
                type = Solver::Trilinos_Belos;
                break;
            case 6:
                type = Solver::Trilinos_ML;
                break;
            case 7:
                type = Solver::Trilinos_Ifpack;
                break;
            case 8:
                type = Solver::PETSc;
                break;
            case 9:
                type = Solver::ANI;
                break;
            case 10:
                type = Solver::FCBIILU2;
                break;
            case 11:
                type = Solver::K3BIILU2;
                break;
                
            default:
                if (processRank == 0) {
                    std::cout << "Invalid solver type index: " << solverType << " , using INNER_ILU2 instead." <<
                    std::endl;
                }
                type = Solver::INNER_ILU2;
                break;
        }

        Solver::Initialize(&argc, &argv, parametersFound ? parametersFileName.c_str()
                                                         : NULL); // Initialize the linear solver in accordance with args

        if (processRank == 0) {
            std::cout << "Solving with " << Solver::TypeName(type) << std::endl;
        }

        Sparse::Matrix mat("A"); // Declare the matrix of the linear system to be solved
        Sparse::Vector b("rhs"); // Declare the right-hand side vector
        Sparse::Vector x("sol"); // Declare the solution vector
        double tempTimer = Timer(), solvingTimer;
        mat.Load(matrixFileName); //if interval parameters not set, matrix will be divided automatically
        BARRIER
        if (processRank == 0) {
            std::cout << "load matrix:    " << Timer() - tempTimer << std::endl;
        }
        tempTimer = Timer();
        if (vectorBFound) {
            b.Load(vectorBFileName); // Load RHS vector
            //b.Save(vectorBFileName + "_saved_test.rhs");
        } else { // Set local RHS to 1 if it was not specified
            INMOST_DATA_ENUM_TYPE mbeg, mend, k;
            mat.GetInterval(mbeg, mend);
            b.SetInterval(mbeg, mend);
            for (k = mbeg; k < mend; ++k) b[k] = 1.0;
        }
        BARRIER
        if (processRank == 0) {
            std::cout << "load vector:    " << Timer() - tempTimer << std::endl;
        }
        solvingTimer = Timer();
        int iters;
        bool success;
        double resid, realresid = 0;
        std::string reason;
        Solver s(type); // Declare the linear solver by specified type

        if (parametersFound) {
            char *fileName = findInnerOptions(parametersFileName.c_str());
            if (fileName != NULL) {
                InnerOptions *options = parseInnerDatabaseOptions(fileName);
                if (options != nullptr) {
                    for (int ii = 0; ii < options->options.size(); ii++) {
                        InnerOption *option = options->options[ii];
                        if (option->type == ENUM) {
                            s.SetParameterEnum(option->name, (unsigned int) atoi(option->value.c_str()));
                        } else {
                            s.SetParameterReal(option->name, atof(option->value.c_str()));
                        };
                    }
                    delete options;
                }
                free(fileName);
            }
        }
        BARRIER

        //s.SetParameterEnum("maximum_iterations", 1000);
        //s.SetParameterEnum("gmres_substeps", 4);
        //s.SetParameterReal("relative_tolerance", 1.0e-6);
        //s.SetParameterReal("absolute_tolerance", 1.0e-16);
        //s.SetParameterReal("divergence_tolerance", 1e+200);

        //s.SetParameterEnum("reorder_nonzeros", 0);
        //s.SetParameterEnum("rescale_iterations", 8);
        //s.SetParameterEnum("adapt_ddpq_tolerance", 0);

        //s.SetParameterReal("drop_tolerance", 3.0e-3);
        //s.SetParameterReal("reuse_tolerance", 1.0e-5);
        //s.SetParameterReal("ddpq_tolerance", 0.0);

        //s.SetParameterEnum("condition_estimation", 1);
        //s.SetParameterEnum("schwartz_overlap", 3);


        tempTimer = Timer();
        s.SetMatrix(mat); // Compute the preconditioner for the original matrix
        BARRIER
        if (processRank == 0) std::cout << "preconditioner time: " << Timer() - tempTimer << std::endl;
        tempTimer = Timer();
        success = s.Solve(b, x); // Solve the linear system with the previously computted preconditioner
        BARRIER
        solvingTimer = Timer() - solvingTimer;
        if (processRank == 0) std::cout << "iterations time:     " << Timer() - tempTimer << std::endl;
        iters = s.Iterations(); // Get the number of iterations performed
        resid = s.Residual();   // Get the final residual achieved
        reason = s.GetReason(); // Get the convergence reason
        //x.Save("output.sol");  // Save the solution if required

        // Compute the true residual
        double aresid = 0, bresid = 0;
        Sparse::Vector test;
        tempTimer = Timer();
        Solver::OrderInfo info;
        info.PrepareMatrix(mat, 0);
        info.PrepareVector(x);
        info.Update(x);

        mat.MatVec(1.0, x, 0.0, test); // Multiply the original matrix by a vector
        {
            INMOST_DATA_ENUM_TYPE mbeg, mend, k;
            info.GetLocalRegion(info.GetRank(), mbeg, mend);
            for (k = mbeg; k < mend; ++k) {
                aresid += (test[k] - b[k]) * (test[k] - b[k]);
                bresid += b[k] * b[k];
            }
        }
        double temp[2] = {aresid, bresid}, recv[2] = {aresid, bresid};
#if defined(USE_MPI)
        MPI_Reduce(temp, recv, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
        realresid = sqrt(recv[0] / (recv[1] + 1.0e-100));
        info.RestoreVector(x);
        if (processRank == 0) {
            std::cout << "||Ax-b||=" << sqrt(recv[0]) << " ||b||=" << sqrt(recv[1]) << " ||Ax-b||/||b||=" <<
            sqrt(recv[0] / (recv[1] + 1.0e-100)) << std::endl;
            std::cout << "norms: " << Timer() - tempTimer << std::endl;
            std::cout << processorsCount << " processors for Solver::type=" << type;
            if (success) {
                std::cout << " solved in " << solvingTimer << " secs";
                std::cout << " with " << iters << " iterations to " << resid << " norm";
            } else std::cout << " failed to solve with " << iters << "iterations";
            std::cout << " matrix \"" << matrixFileName << "\"";
            if (vectorBFound) std::cout << " vector \"" << vectorBFileName << "\"";
            std::cout << " true residual ||Ax-b||/||b||=" << realresid;
            std::cout << std::endl;
            std::cout << "reason: " << reason << std::endl;
        }

        if (vectorXFound) {
            Sparse::Vector ex("exact");  // Declare the exact solution vector
            Sparse::Vector err("error"); // Declare the solution error vector
            INMOST_DATA_ENUM_TYPE mbeg, mend, k;
            mat.GetInterval(mbeg, mend);
            err.SetInterval(mbeg, mend);
            ex.Load(vectorXFileName);
            BARRIER
            double dif1 = 0, dif2 = 0, dif8 = 0, norm = 0;
            for (k = mbeg; k < mend; ++k) {
                double dloc = err[k] = fabs(x[k] - ex[k]);
                dif1 += dloc;
                dif2 += dloc * dloc;
                dif8 = (dif8 > dloc) ? dif8 : dloc;
                norm += fabs(ex[k]);
            }
#if defined(USE_MPI)
            if (processorsCount > 1) {
                double temp1 = dif1;
                MPI_Reduce(&temp1, &dif1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                temp1 = dif2;
                MPI_Reduce(&temp1, &dif2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                temp1 = dif8;
                MPI_Reduce(&temp1, &dif8, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                temp1 = norm;
                MPI_Reduce(&temp1, &norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
#endif
            dif2 = sqrt(dif2);
            norm += 1.0e-100;
            if (processRank == 0) {
                std::cout << "Difference with exact solution \"" << vectorXFileName << "\": " << std::scientific <<
                std::setprecision(6) << std::endl;
                std::cout << "dif1 = " << dif1 << "  dif2 = " << dif2 << "  dif8 = " << dif8 << "  ||ex||_1 = " <<
                norm << std::endl;
                std::cout << "rel1 = " << dif1 / norm << "  rel2 = " << dif2 / norm << "  rel8 = " << dif8 / norm <<
                std::endl;
            }
        }
    }
    BARRIER
	Solver::Finalize(); // Finalize solver and close MPI activity
	return 0;
}
