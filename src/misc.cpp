#include "line_sep/line.hpp"
#include "line_sep/misc.hpp"
#include <algorithm>
#include <execution>
#include <fstream>
#include <ranges>

namespace line_sep::misc {
    namespace io {
        void load_lines(const std::filesystem::path& filepath, 
                        std::map<int32_t, line::Line>& lines)
        {
            std::vector<int> lsymbs;
            std::vector<double> xs;
            std::vector<double> ys;
            int lnum {0};
            double x {0.0};
            double y {0.0};
            int lsymb {0};
        
            std::ifstream ifile(filepath);
            if (!ifile.is_open())
                throw std::ios::failure("Can't open input file");
        
            std::string sline;
            for (size_t i {0}; std::getline(ifile, sline); ++i) {
                if (i == 0)
                    continue;
        
                std::vector<std::string> tokens = tokenise(sline);
                lnum = std::stoi(tokens[0]);
                x = std::stod(tokens[1]);
                y = std::stod(tokens[2]);
                lsymb = std::stoi(tokens[3]);
                
                if ((i > 1) && (lsymb != lsymbs.back())) {
                    lines[lsymbs.back()] = line::Line(lnum, lsymbs.back(), xs, ys);
                    lsymbs.clear();
                    xs.clear();
                    ys.clear();
                }
                lsymbs.push_back(lsymb);
                xs.push_back(x);
                ys.push_back(y);
            }
            lines[lsymbs.back()] = line::Line(lnum, lsymbs.back(), xs, ys);
            ifile.close();
        }

        void safe_lines(const std::filesystem::path& filepath, const std::map<int32_t, line::Line>& lines)
        {
            std::ofstream ofile(filepath);
            if (!ofile.is_open())
                throw std::ios::failure("Can't open input file");

            ofile << "LNUM,X,Y,LSYMB\n";
            for (const auto& [symb, ln] : lines) {
                for (int row {0}; row < ln.xy.rows(); ++row)
                    ofile << ln.num << "," << ln.xy.row(row)(0) << "," << ln.xy.row(row)(1) << "," << ln.symb << "\n";
            }
            ofile.close();
        }

        void safe_lines_with_dist(const std::filesystem::path& filepath, const std::map<int32_t, line::Line>& lines)
        {
            std::ofstream ofile(filepath);
            if (!ofile.is_open())
                throw std::ios::failure("Can't open input file");

            ofile << "LNUM,X,Y,LSYMB,DIST_LEFT,DIST_RIGHT\n";
            for (const auto& [symb, ln] : lines) {
                for (int row {0}; row < ln.xy.rows(); ++row) {
                    ofile << ln.num << "," << ln.xy.row(row)(0) << "," << ln.xy.row(row)(1) << "," << ln.symb;
                    ofile << "," << ln.dist_left[row] << "," << ln.dist_right[row] << "\n";
                }
            }
            ofile.close();
        }

        void safe_grid(const std::filesystem::path& filepath, const std::map<int32_t, line::Line>& lines,
                        const std::map<int32_t, std::shared_ptr<line::GCol>> grid)
        {
            std::ofstream ofile(filepath);
            if (!ofile.is_open())
                throw std::ios::failure("Can't open input file");

            ofile << "LNUM_MAIN,LNUM_SEC,X,Y,CELL\n";
            for (const auto& [symb, grid] : grid) {
                for (int row {0}; row < (*grid).size(); ++row) {
                    if ((*grid)[row].size() == 0)
                        continue;
                    for (auto point : (*grid)[row]) {
                        ofile << symb << ","
                                << point.first << ',' 
                                << lines.at(point.first).xy.row(point.second)(0) << ',' 
                                << lines.at(point.first).xy.row(point.second)(1) << ',' << row << '\n';
                    }
                }
            }
            ofile.close();
        }

        std::vector<std::string> tokenise(const std::string& line)
        {
            std::vector<std::string> tokens;
            std::stringstream ss(line);
            std::string token;
            while (std::getline(ss, token, ','))
                tokens.push_back(token);
            return tokens;
        }
    }  
    double calc_circular_mean(const std::vector<double>& azimuths)
    {
        double cos_alpha = 0.0;
        double sin_alpha = 0.0;

        for (const auto& azimuth : azimuths) {
            cos_alpha += std::cos(azimuth);
            sin_alpha += std::sin(azimuth);
        }
        return std::atan2(sin_alpha, cos_alpha);
    }
}