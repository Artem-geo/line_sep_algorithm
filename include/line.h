#pragma once

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <numbers>

namespace line {
    class Line {
    public:
        int32_t symb;
        Eigen::MatrixXd xy;        
    public:
        Line(int32_t lsymb,
             std::vector<double>& x,
             std::vector<double>& y);
        friend std::ostream& operator<<(std::ostream& os, const Line& line);
    private:
        Line() = default;
    };
}