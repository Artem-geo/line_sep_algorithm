#include "line.hpp"
#include "line_sep.hpp"
#include "misc.hpp"
#include <execution>
#include <algorithm>s

using misc::PI;

namespace line_sep {
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

    void set_left_hood(std::map<int32_t, line::LR_hood>& lhood, std::vector<std::pair<int32_t, double>> line_signs, double ndst, double tol)
    {
        double tdst = ndst * tol;
        for (auto litr = line_signs.rbegin(); litr != line_signs.rend(); ++litr) { // litr for "line iterator"
            if (litr + 1 == line_signs.rend()) // if next line is the last one
                break;
            for (auto nlitr = litr + 1; nlitr != line_signs.rend(); ++nlitr) { // nlitr for "next line iterator"
                if (!is_adjacent(litr, nlitr, ndst)) // if nlitr is further than 1.5 line separation it is considered as not adjacent
                    break;
                if (!is_partial(litr, nlitr, tdst))
                    lhood[(*litr).first].left.insert((*nlitr).first);
            }
        }
    }

    void set_right_hood(std::map<int32_t, line::LR_hood>& lhood, std::vector<std::pair<int32_t, double>> line_signs, double ndst, double tol)
    {
        double tdst = ndst * tol;
        for (auto litr = line_signs.begin(); litr != line_signs.end(); ++litr) { // litr for "line iterator"
            if (litr + 1 == line_signs.end()) // if next line is the last one
                break;
            for (auto nlitr = litr + 1; nlitr != line_signs.end(); ++nlitr) { // nlitr for "next line iterator"
                if (!is_adjacent(litr, nlitr, ndst)) // if nlitr is further than 1.5 line separation it is considered as not adjacent
                    break;
                if (!is_partial(litr, nlitr, tdst))
                    lhood[(*litr).first].right.insert((*nlitr).first);
            }
        }
    }

    void sort_lines(std::map<int32_t, line::Line>& lines,
                    std::map<int32_t, line::LR_hood>& lhood,
                    double ndst, double tol)
    {
        double tdst = ndst * tol; // distance tolerance
        std::vector<std::pair<int32_t, double>> line_signs; // line signatures = {symbol, projection}

        for (auto& [symb, ln] : lines) {
            line_signs.push_back({symb, ln.proj}); // create signatures
            lhood[symb] = line::LR_hood(); // init neighbours
        }

        std::sort(line_signs.begin(), line_signs.end(), // line signatures are sorted in ascending order depending on proj
                  [](const std::pair<int32_t, double>& line_sign1, const std::pair<int32_t, double>& line_sign2)
                  {
                      return line_sign1.second < line_sign2.second;
                  });

        set_left_hood(lhood, line_signs, ndst, tol);
        set_right_hood(lhood, line_signs, ndst, tol);
    }

    void grid_lines(const std::map<int32_t, line::Line>& lines, const std::map<int32_t, line::LR_hood>& lhood, std::map<int32_t, std::shared_ptr<GCol>>& grid, double nomdst)
    {
        const std::pair<double, double>& ylim = get_ylim(lines);
        double dy = nomdst * 0.05; // 5% of line spacing should be fine enought for the analysis
        double ymin_grid = (ylim.first - dy / 2);
        double ymax_grid = (ylim.second + dy / 2);
        int num_cells = (int) std::ceil((ymax_grid - ymin_grid) / dy); // number of cells = ((ymax+dy/2) - (ymin-dy/2))/dy rounded to the nearest int up
        
        for (auto& [symb, lh] : lhood) {
            if (lh.left.size() != 0) {
                create_grid(lines, lh.left, grid, num_cells, dy, ymin_grid);
            }
            if (lh.right.size() != 0) {
                create_grid(lines, lh.right, grid, num_cells, dy, ymin_grid);
            }
        }
    }

    void create_grid(const std::map<int32_t, line::Line>& lines, const std::set<int32_t>& lh, std::map<int32_t, std::shared_ptr<GCol>>& grid, int num_cells, double dy, double ymin_grid)
    {
        auto grid_column = std::make_shared<GCol>(num_cells, Cell()); // temporary ptr to a grid column of num_cells size
        for (auto symb : lh) {
            if (!grid.contains(symb)) { // if the grid item has not been created yet
                grid[symb] = grid_column;
                auto& xy = lines.at(symb).xy;
                auto stride = xy.stride(); // stride is specified to get y coordinate in the column align matrix (default of Eigen)
                int index {-1};

                for (size_t nr {0}; nr < (size_t) xy.rows(); ++nr) { // nr is a row number
                    const double* xc = xy.row(nr).data();
                    const double* yc = xc + stride;
                    Point point = {xc, yc};
                    index = (int) std::ceil((*yc - ymin_grid)/dy - 1); // index of the cell to put the point into
                    grid_column->at(index).insert(point);
                }
            }
        }
        grid_column = nullptr; // disconnect the temporary link from the pointer
    }

    void calc_dist(std::map<int32_t, line::Line>& lines, const std::set<int32_t>& lh, std::map<int32_t, std::shared_ptr<GCol>>& grid)
    {
        
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