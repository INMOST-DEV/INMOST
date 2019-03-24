//
// Created by bvdmitri on 10.02.19.
//

#ifndef INMOST_SERIES_H
#define INMOST_SERIES_H

#include <vector>
#include <string>

class MatrixSeries {
private:
    std::size_t current;
    std::vector<std::string> matrices;
    std::vector<std::string> rhss;

    MatrixSeries(const MatrixSeries &other) = delete;
    MatrixSeries(MatrixSeries &&other) = delete;
    MatrixSeries &operator=(const MatrixSeries &other) = delete;
public:
    MatrixSeries(const std::string &file, const std::string &directory_prefix = "");

    bool end() const;

    void restart();

    int size() const;

    std::pair<const char *, const char *> next();
};


#endif //INMOST_SERIES_H
