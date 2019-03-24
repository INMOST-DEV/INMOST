//
// Created by bvdmitri on 10.02.19.
//

#include "series.h"

#include <algorithm>
#include <fstream>
#include <iostream>

// trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}

static inline bool is_file_exist(const std::string &file) {
    std::ifstream infile(file);
    return infile.good();
}

MatrixSeries::MatrixSeries(const std::string &file, const std::string &directory_prefix) : current(0) {
    std::ifstream input(file);
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }

        line = trim(line);
        std::string path = directory_prefix + line.substr(2);

        if (!is_file_exist(path)) {
            std::cerr << "[WARN] File " << path << " does not exist. Skipping..." << std::endl;
            continue;
        }

        switch (line.at(0)) {
            case 'A':
                matrices.push_back(path);
                break;
            case 'b':
                rhss.push_back(path);
                break;
            default:
                std::cerr << "[WARN] Invalid line in matrix series configuration file: " << line << std::endl;
        }
    }
}

bool MatrixSeries::end() const {
    return current == matrices.size();
}

void MatrixSeries::restart() {
    current = 0;
}

std::pair<const char *, const char *> MatrixSeries::next() {
    const char *matrix = current < matrices.size() ? matrices.at(current).c_str() : nullptr;
    const char *rhs = current < rhss.size() ? rhss.at(current).c_str() : nullptr;
    current += 1;
    return std::make_pair(matrix, rhs);
}

int MatrixSeries::size() const {
    return matrices.size();
}
