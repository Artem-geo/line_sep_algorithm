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
        void init_grid(const std::map<int32_t, line::LR_adj>& line_adjs, std::map<int32_t, std::shared_ptr<line_sep::GCol>> grid, double nomdst)
        {
            std::set<std::set<int32_t>> lines_groups;
            for (auto [symb, lr_lines] : line_adjs) {
                if (lr_lines.left.size() != 0)
                    lines_groups.insert(lr_lines.left);
                if (lr_lines.right.size() != 0)
                    lines_groups.insert(lr_lines.right);
            }

            for (auto lines_group : lines_groups) {
                auto ptr_temp 
                for ()
            }
        }

        std::tuple<double, double, double> grid_lines(const std::map<int32_t, line::Line>& lines, const std::map<int32_t, line::LR_adj>& line_adjs, 
                                                      std::map<int32_t, std::shared_ptr<GCol>>& grid, double nomdst)
        {
            const std::pair<double, double>& ylim = geometry::get_ylim(lines);
            double dy = nomdst * 0.1; // 10% of line spacing should be fine enought for the analysis
            double ymin_grid = (ylim.first - dy / 2);
            double ymax_grid = (ylim.second + dy / 2);
            int num_cells = (int) std::ceil((ymax_grid - ymin_grid) / dy); // number of cells = ((ymax+dy/2) - (ymin-dy/2))/dy rounded to the nearest int up

            for (auto& [symb, lh] : line_adjs) {
                if (lh.left.size() != 0)
                    create_grid(lines, lh.left, grid, num_cells, dy, ymin_grid);
                if (lh.right.size() != 0)
                    create_grid(lines, lh.right, grid, num_cells, dy, ymin_grid);
            }
            return {dy, ymin_grid, ymax_grid};
        }
        void create_grid(const std::map<int32_t, line::Line>& lines, const std::set<int32_t>& line_adj, 
                         std::map<int32_t, std::shared_ptr<GCol>>& grid, 
                         int num_cells, double dy, double ymin_grid)
        {
            auto grid_column = std::make_shared<GCol>(num_cells, Cell()); // temporary ptr to a grid column of num_cells size
            for (auto symb : line_adj) {
                if (!grid.contains(symb)) { // if the grid item has not been created yet
                    grid[symb] = grid_column;
                    auto& xy = lines.at(symb).xy;
                    int idx {-1};
                    for (size_t nr {0}; nr < (size_t) xy.rows(); ++nr) { // nr is a row number
                        Point point = {symb, nr};
                        idx = (int) std::ceil((xy.row(nr)(1) - ymin_grid) / dy - 1); // index of the cell to put the point into
                        grid_column->at(idx).insert(point);
                    }
                }
            }
            grid_column = nullptr; // disconnect the temporary link from the pointer
        }
    }

    namespace distance {
        void calc_line_dists(std::map<int32_t, line::Line>& lines, const std::map<int32_t, line::LR_adj>& line_adjs,
                             std::map<int32_t, std::shared_ptr<GCol>>& grid, const std::tuple<double, double, double>& grid_params)
        {
            for (auto [symb, line] : lines) {
                auto& left_adj = line_adjs.at(symb).left;
                auto& right_adj = line_adjs.at(symb).right;
                if (!left_adj.empty()) {
                    auto left_symb = left_adj.begin();
                    calc_dist(line.xy, line.dist_left, *grid[*left_symb], lines, grid_params);
                }
                if (!right_adj.empty()) {
                    if (symb == 110) {
                        int temp {0};
                    }
                    auto right_symb = right_adj.begin();
                    calc_dist(line.xy, line.dist_right, *grid[*right_symb], lines, grid_params);
                }
            }
        }

        void calc_dist(const Eigen::MatrixXd& xy, std::vector<double>& distances, const GCol& grid_column,
                       const std::map<int32_t, line::Line>& lines,
                       const std::tuple<double, double, double>& grid_params)
        {
            auto [dy, ymin_grid, ymax_grid] = grid_params;
            size_t idx {0};

            auto calc_dist_from_cell = [&](size_t i, size_t idx) // i - index of current line point, idx - adjacent line cell index
                {
                    double dist {-1.0};
                    double temp_dist {0.0};
                    for (auto& point : grid_column[idx]) {
                        temp_dist = std::sqrt((xy.row(i) - lines.at(point.first).xy.row(point.second)).squaredNorm());
                        dist = (temp_dist > dist) ? temp_dist : dist;
                    }
                    return dist;
                };

            auto calc_dist_from_cells_adjacent = [&](size_t i, size_t idx)
                {
                    // cases to consider: first cell, last cell, empty previous, empty next
                    double dist_cell_1 {-1.0};
                    if (idx != 0) {
                        if (grid_column[idx - 1].size() != 0)
                            dist_cell_1 = calc_dist_from_cell(i, idx - 1);
                    }

                    double dist_cell_2 {-1.0};
                    if (idx != distances.size() - 1) {
                        if (grid_column[idx + 1].size() != 0)
                            dist_cell_2 = calc_dist_from_cell(i, idx + 1);
                    }

                    if ((dist_cell_1 >= 0) && (dist_cell_2 >= 0))
                        return (dist_cell_1 + dist_cell_2) / 2;
                    else
                        return (dist_cell_1 >= 0) ? dist_cell_1 : ((dist_cell_2 >= 0) ? dist_cell_2 : misc::rDUMMY);
                };

            for (size_t i {0}; i < distances.size(); ++i) {
                idx = (int)std::ceil(((xy.row(i)(1) - ymin_grid) / dy) - 1); // index of the cell
                if (grid_column[idx].size() != 0) {
                    distances[i] = calc_dist_from_cell(i, idx);
                }
                else {
                    distances[i] = calc_dist_from_cells_adjacent(i, idx);
                }
            }
        }
    }
}