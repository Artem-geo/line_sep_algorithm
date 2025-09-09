#pragma once

#include "line_sep/line_sep.hpp"
#include <filesystem>
#include <limits>
#include <map>
#include <numbers>
#include <vector>

namespace line_sep::misc {
    constexpr double DOUBLE_MAX {std::numeric_limits<double>::max()};
    constexpr double DOUBLE_MIN {std::numeric_limits<double>::min()};
    const double PI = std::numbers::pi;
    constexpr double rDUMMY = -1e+32;

    namespace io {
        void load_lines(const std::filesystem::path& filepath, std::map<int32_t, line::Line>&lines);
        void safe_lines(const std::filesystem::path& filepath, const std::map<int32_t, line::Line>& lines);
        void safe_lines_with_dist(const std::filesystem::path& filepath, const std::map<int32_t, line::Line>& lines);
        void safe_grid(const std::filesystem::path& filepath, const std::map<int32_t, line::Line>& lines, const std::map<int32_t, std::shared_ptr<line::GCol>> grid);
        std::vector<std::string> tokenise(const std::string& line);
    }
    double calc_circular_mean(const std::vector<double>& azimuths);
}