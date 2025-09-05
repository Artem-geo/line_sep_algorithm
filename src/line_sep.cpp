#include "line.hpp"
#include "line_sep.hpp"
#include "misc.hpp"
#include <algorithm>
#include <execution>
#include <iostream>

using misc::PI;

namespace line_sep {
    namespace geometry {
        double calc_azimuth(std::map<int32_t, line::Line>& lines)
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
            for (auto& [symb, ln] : lines) { // lnum - line number, plines - line partials (if any)
                // running sums
                sx = ln.xy.col(0).mean();
                sy = ln.xy.col(1).mean();
                sxx = ln.xy.col(0).array().square().mean();
                syy = ln.xy.col(1).array().square().mean();
                sxy = (ln.xy.col(0).array() * ln.xy.col(1).array()).mean();

                // accumulators
                sxx_prime = sxx - std::pow(sx, 2);
                syy_prime = syy - std::pow(sy, 2);
                sxy_prime = sxy - sx * sy;

                // azimuth in rad
                azimuth = 0.5 * std::atan2(2 * sxy_prime, sxx_prime - syy_prime);
                azimuth = misc::PI / 2 - azimuth; // to get azimuth from [0, pi), where 0 = North, pi = South
                ln.azimuth = azimuth;
                azimuths.push_back(azimuth);
            }
            return misc::calc_circular_mean(azimuths);
        }

        Eigen::RowVector2d get_pivot(std::map<int32_t, line::Line>& lines)
        {
            Eigen::RowVector2d pivot = {0.0, 0.0};
            for (auto& [symb, ln] : lines)
                pivot += Eigen::RowVector2d({(ln.xlim.first + ln.xlim.second) / 2, (ln.ylim.first + ln.ylim.second) / 2});
            pivot(0) /= lines.size();
            pivot(1) /= lines.size();
            return pivot;
        }

        void rotate_lines(std::map<int32_t, line::Line>& lines, double angle, const Eigen::RowVector2d& pivot)
        {
            Eigen::Rotation2D<double> r(angle);
            Eigen::MatrixXd rmatrix = r.toRotationMatrix();
            for (auto& [symb, ln] : lines) {
                ln.xy = (ln.xy.rowwise() - pivot) * rmatrix.transpose();
                ln.xy = ln.xy.rowwise() + pivot;
                ln.azimuth -= angle;
                ln.azimuth = (ln.azimuth < 0) ? ln.azimuth + PI : ln.azimuth; // so that resultant azimuth is [0, pi) again (without negative values)
                ln.xlim = {ln.xy.col(0).minCoeff(), ln.xy.col(0).maxCoeff()};
                ln.ylim = {ln.xy.col(1).minCoeff(), ln.xy.col(1).maxCoeff()};
            }
        }

        void project_lines(std::map<int32_t, line::Line>& lines, const Eigen::RowVector2d& pivot) // fine tunes (following the main rotation) the lines and projects them onto the X axis
        {
            Eigen::MatrixXd rotated_matrix;
            double angle {0.0};
            for (auto& [symb, ln] : lines) {
                angle = (ln.azimuth > PI / 2) ? ln.azimuth - PI : ln.azimuth; // to choose whether the line should be rotated CW (angle>pi/2) or CCW (angle<pi/2 or = pi/2)
                Eigen::Rotation2D<double> r(angle);
                Eigen::MatrixXd rmatrix = r.toRotationMatrix();
                rotated_matrix = (ln.xy.rowwise() - pivot) * rmatrix.transpose(); // fine-tune line orientation
                rotated_matrix = rotated_matrix.rowwise() + pivot;
                ln.proj = rotated_matrix.col(0).mean(); // projection equals to the mean x value of the fine-tuned matrix
            }
        }

        void set_left_adj(std::map<int32_t, line::LR_adj>& line_adjs, std::vector<std::pair<int32_t, double>> line_signs, double ndst, double tol)
        {
            double tdst = ndst * tol;
            for (auto litr = line_signs.rbegin(); litr != line_signs.rend(); ++litr) { // litr for "line iterator"
                if (litr + 1 == line_signs.rend()) // if next line is the last one
                    break;
                for (auto nlitr = litr + 1; nlitr != line_signs.rend(); ++nlitr) { // nlitr for "next line iterator"
                    if (!is_adjacent(litr, nlitr, ndst)) // if nlitr is further than 1.5 line separation it is considered as not adjacent
                        break;
                    if (!is_partial(litr, nlitr, tdst))
                        line_adjs[(*litr).first].left.insert((*nlitr).first);
                }
            }
        }

        void set_right_adj(std::map<int32_t, line::LR_adj>& line_adjs, std::vector<std::pair<int32_t, double>> line_signs, double ndst, double tol)
        {
            double tdst = ndst * tol;
            for (auto litr = line_signs.begin(); litr != line_signs.end(); ++litr) { // litr for "line iterator"
                if (litr + 1 == line_signs.end()) // if next line is the last one
                    break;
                for (auto nlitr = litr + 1; nlitr != line_signs.end(); ++nlitr) { // nlitr for "next line iterator"
                    if (!is_adjacent(litr, nlitr, ndst)) // if nlitr is further than 1.5 line separation it is considered as not adjacent
                        break;
                    if (!is_partial(litr, nlitr, tdst))
                        line_adjs[(*litr).first].right.insert((*nlitr).first);
                }
            }
        }

        void sort_lines(std::map<int32_t, line::Line>& lines,
                        std::map<int32_t, line::LR_adj>& line_adjs,
                        double ndst, double tol)
        {
            double tdst = ndst * tol; // distance tolerance
            std::vector<std::pair<int32_t, double>> line_signs; // line signatures = {symbol, projection}

            for (auto& [symb, ln] : lines) {
                line_signs.push_back({symb, ln.proj}); // create signatures
                line_adjs[symb] = line::LR_adj(); // init neighbours
            }

            std::sort(line_signs.begin(), line_signs.end(), // line signatures are sorted in ascending order depending on proj
                      [](const std::pair<int32_t, double>& line_sign1, const std::pair<int32_t, double>& line_sign2)
                      {
                          return line_sign1.second < line_sign2.second;
                      });

            set_left_adj(line_adjs, line_signs, ndst, tol);
            set_right_adj(line_adjs, line_signs, ndst, tol);
        }

        std::pair<double, double> get_ylim(const std::map<int32_t, line::Line>& lines)
        {
            double ymin {misc::DOUBLE_MAX};
            double ymax {misc::DOUBLE_MIN};

            for (const auto& [symb, ln] : lines) {
                ymin = (ln.ylim.first < ymin) ? ln.ylim.first : ymin;
                ymax = (ln.ylim.second > ymax) ? ln.ylim.second : ymax;
            }
            return {ymin, ymax};
        }
    }

    namespace gridding {
        void init_grid(const std::map<int32_t, line::LR_adj>& line_adjs, std::map<int32_t, std::shared_ptr<line_sep::GCol>>& grid, int ncells)
        {
            std::set<std::set<int32_t>> lines_groups;
            for (auto [symb, lr_lines] : line_adjs) {
                if (lr_lines.left.size() != 0)
                    lines_groups.insert(lr_lines.left);
                if (lr_lines.right.size() != 0)
                    lines_groups.insert(lr_lines.right);
            }

            for (auto& lines_group : lines_groups) {
                auto ptr_temp = std::make_shared<GCol>(ncells, Cell());
                for (auto& symb : lines_group)
                    grid[symb] = ptr_temp;
                ptr_temp = nullptr;
            }
        }

        void grid_lines(const std::map<int32_t, line::Line>& lines, const std::map<int32_t, line::LR_adj>& line_adjs, 
                        std::map<int32_t, std::shared_ptr<GCol>>& grid, 
                        double ymin_grid, double dy)
        {
            size_t idx {0};
            for (auto& [symb, line] : lines) {
                for (size_t i {0}; i < (size_t)line.xy.rows(); ++i) {
                    idx = (int)std::ceil(((line.xy.row(i)(1) - ymin_grid) / dy) - 1); // index of the cell to put the point into
                    (*(grid.at(symb)))[idx].insert({symb, i});
                }
                idx = 0;
            }
        }

        std::tuple<int, double, double, double> get_number_cells(const std::map<int32_t, line::Line>& lines, double nomdst, double grid_tol)
        {
            const std::pair<double, double>& ylim = geometry::get_ylim(lines);
            double dy = nomdst * grid_tol;
            double ymin_grid = (ylim.first - dy / 2);
            double ymax_grid = (ylim.second + dy / 2);
            int ncells = (int)std::ceil((ymax_grid - ymin_grid) / dy); // number of cells = ((ymax+dy/2) - (ymin-dy/2))/dy rounded to the nearest int up
            return {ncells, ymin_grid, ymax_grid, dy}; 
        }
    }

    namespace distance {
        void calc_line_dists(std::map<int32_t, line::Line>& lines, const std::map<int32_t, line::LR_adj>& line_adjs, 
                             const std::map<int32_t, std::shared_ptr<GCol>>& grid, double ymin_grid, double dy)
        {
            for (auto [symb, line] : lines) {
                auto& left_adj = line_adjs.at(symb).left;
                auto& right_adj = line_adjs.at(symb).right;
                std::cout << "Symbol: " << symb;
                if (!left_adj.empty()) {
                    auto left_symb = left_adj.begin();
                    std::cout << " left: " << *left_symb;
                    calc_dist(lines, LR::left, grid, symb, *left_symb, ymin_grid, dy);
                }
                if (!right_adj.empty()) {
                    auto right_symb = right_adj.begin();
                    std::cout << " right: " << *right_symb;
                    calc_dist(lines, LR::right, grid, symb, *right_symb, ymin_grid, dy);
                }
                std::cout << '\n';
            }
        }

        void calc_dist(std::map<int32_t, line::Line>& lines, LR lr, 
                       const std::map<int32_t, std::shared_ptr<GCol>>& grid, 
                       int32_t lsymb_target, int32_t lsymb_adj, double ymin_grid, double dy)
        {
            auto& xy_target = lines.at(lsymb_target).xy;
            auto& dists_target = (lr == LR::left) ? lines.at(lsymb_target).dist_left : lines.at(lsymb_target).dist_right;
            size_t idx {0};

            auto calc_dist_from_cell = [&](size_t i, Cell& cell_target) // i - index of current line point, idx - adjacent line cell index
                {
                    double dist {misc::DOUBLE_MAX};
                    double temp_dist {0.0};
                    for (auto& point : cell_target) {
                        temp_dist = std::sqrt((xy_target.row(i) - lines.at(point.first).xy.row(point.second)).squaredNorm());
                        dist = (temp_dist <  dist) ? temp_dist : dist;
                    }
                    return dist;
                };

            for (size_t i {0}; i < (size_t)xy_target.rows(); ++i) {
                idx = (int)std::ceil(((xy_target.row(i)(1) - ymin_grid) / dy) - 1);
                auto& cell_target = (*(grid.at(lsymb_adj)))[idx];
                if (cell_target.size() != 0)
                    dists_target[i] = calc_dist_from_cell(i, cell_target);
            }
        }
    }
}