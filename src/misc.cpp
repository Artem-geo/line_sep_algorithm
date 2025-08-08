#include "line.h"
#include "misc.h"
#include <execution>
#include <iostream>
#include <ranges>

using namespace line;

namespace misc {
    namespace io {
        void load_lines(const std::filesystem::path& filepath, 
                        std::map<int, VLines>& lines, 
                        std::set<int>& line_numbers)
        {
            std::vector<int> lnums;
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
                x = std::stof(tokens[1]);
                y = std::stof(tokens[2]);
                lsymb = std::stoi(tokens[3]);
        
                if ((i > 1) && (lsymb != lsymbs.back())) {
                    lines[lnums.back()].emplace_back(lsymbs.back(), xs, ys);
                    line_numbers.emplace(lnums.back());
                    lsymbs.clear();
                    lnums.clear();
                    xs.clear();
                    ys.clear();
                }
                lnums.push_back(lnum);
                xs.push_back(x);
                ys.push_back(y);
                lsymbs.push_back(lsymb);
            }
            lines[lnums.back()].emplace_back(lsymbs.back(), xs, ys);
            line_numbers.emplace(lnums.back());
            ifile.close();
        }

        void safe_lines(const std::filesystem::path& filepath, const std::map<int, VLines>& lines)
        {
            std::ofstream ofile(filepath);
            if (!ofile.is_open())
                throw std::ios::failure("Can't open input file");

            ofile << "LNUM,X,Y,LSYMB\n";
            for (auto& [lnum, plines] : lines) {
                for (const Line& line : plines) {
                    for (int row {0}; row < line.xy.rows(); ++row) {
                        ofile << lnum << "," << line.xy.row(row)(0) << "," << line.xy.row(row)(1) << "," << line.symb << "\n";
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

    double calc_azimuth(std::map<int, VLines>& lines)
    {
        double azimuth {0};
        std::vector<double> azimuths;
        double sx {0};
        double sy {0};
        double sxx {0};
        double syy {0};
        double sxy {0};
        double sxx_prime {0};
        double syy_prime {0};
        double sxy_prime {0};
        for (auto& [lnum, plines] : lines) { // lnum - line number, plines - line partials (if any)
            for (Line& line : plines) {
                // running sums
                sx = line.xy.col(0).mean();
                sy = line.xy.col(1).mean();
                sxx = line.xy.col(0).array().square().mean();
                syy = line.xy.col(1).array().square().mean();
                sxy = (line.xy.col(0).array() * line.xy.col(1).array()).mean();

                // accumulators
                sxx_prime = sxx - std::pow(sx, 2);
                syy_prime = syy - std::pow(sy, 2);
                sxy_prime = sxy - sx * sy;

                // azimuth in rad
                azimuth = 0.5 * std::atan2(2*sxy_prime, sxx_prime - syy_prime);
                azimuth = std::numbers::pi/2 - azimuth; // to get azimuth from [0, pi), where 0 = North, pi = South
                azimuths.push_back(azimuth);
            }
        }

        return calc_circular_mean(azimuths);
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

    void rotate_lines(std::map<int, VLines>& lines, double angle)
    {
        Eigen::Rotation2D<double> r(-angle);
        Eigen::MatrixXd rmatrix = r.toRotationMatrix();
        for (auto& [lnum, plines] : lines) {
            for (auto& line : plines) {
                line.xy = line.xy * rmatrix;
            }
        }
    }
    
}