#include "line.h"
#include "misc.h"
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/statistics/linear_regression.hpp>
#include <cmath>
#include <format>
#include <fstream>
#include <sstream>

using namespace misc;

namespace line {
    Line::Line(int32_t lsymb,
               std::vector<double>& x,
               std::vector<double>& y)
        : symb(lsymb)
    {
        if (x.size() != y.size())
            throw std::invalid_argument("X and Y columns should have the same size in order to be written to a line");

        xy = Eigen::MatrixXd::Constant(x.size(), 2, rDUMMY);
        Eigen::Map<Eigen::VectorXd> xv(x.data(), x.size());
        Eigen::Map<Eigen::VectorXd> yv(y.data(), y.size());
        xy.col(0) = xv;
        xy.col(1) = yv;
    }
    std::ostream& operator<<(std::ostream& os, const Line& line)
    {
        os << "Line symbol: " << line.symb << "; x/y size: " << line.xy.rows() << std::endl;
        return os;
    }
}

