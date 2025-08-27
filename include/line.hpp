#pragma once

#include <Eigen/Dense>
#include <set>

namespace line {
    struct LR_hood { // left-right line neighbourhood
        std::set<int32_t> left;
        std::set<int32_t> right;
    };
    class Line {
    public:
        int num;
        int32_t symb;
        Eigen::MatrixXd xy; 
        double azimuth {0};
        std::pair<double, double> xlim;
        std::pair<double, double> ylim;
        double proj {0.0};
        std::vector<double> dist_left;
        std::vector<double> dist_right;
    public:
        Line() = default;
        Line(int lnum, int32_t lsymb,
             std::vector<double>& x, std::vector<double>& y);
        friend std::ostream& operator<<(std::ostream& os, const Line& line);
    };
}