#pragma once

#include "line.h"
#include <filesystem>
#include <map>
#include <numbers>
#include <set>
#include <vector>

using VLines = std::vector<line::Line>;

namespace misc {

    const double PI = std::numbers::pi;
    const double rDUMMY = -1e+32;

    namespace io {
        void load_lines(const std::filesystem::path& filepath, std::map<int, VLines>& lines, std::set<int>& line_numbers);
        void safe_lines(const std::filesystem::path& filepath, const std::map<int, VLines>& lines);
        std::vector<std::string> tokenise(const std::string& line);
    }
    double calc_azimuth(std::map<int, VLines>& lines);
    double calc_circular_mean(const std::vector<double>& azimuths);
    void rotate_lines(std::map<int, VLines>& lines, double angle);
}