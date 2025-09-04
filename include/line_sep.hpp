#pragma once

#include <Eigen/Dense>
#include <map>
#include <memory>

namespace line_sep {
    class line::Line;
    struct line::LR_adj;

    using Point = std::pair<int32_t, size_t>;
    using Cell = std::set<Point>;
    using GCol = std::vector<Cell>;

    enum class LR {left, right};

    namespace geometry {
        double calc_azimuth(std::map<int32_t, line::Line>&lines);
        Eigen::RowVector2d get_pivot(std::map<int32_t, line::Line>&lines);
        void rotate_lines(std::map<int32_t, line::Line>&lines, double angle, const Eigen::RowVector2d & pivot);
        void project_lines(std::map<int32_t, line::Line>&lines, const Eigen::RowVector2d & pivot);

        template<typename IT>
        bool is_partial(IT litr1, IT litr2, double tdst)
        {
            return std::abs((*litr1).second - (*litr2).second) < tdst;
        }

        template<typename IT>
        bool is_adjacent(IT litr1, IT litr2, double ndst)
        {
            return std::abs((*litr1).second - (*litr2).second) < 1.5 * ndst; // lines further than 1.5 line spacings are not adjecent
        }

        void set_left_adj(std::map<int32_t, line::LR_adj>&line_adjs, std::vector<std::pair<int32_t, double>> line_signs, double ndst, double tol);
        void set_right_adj(std::map<int32_t, line::LR_adj>&line_adjs, std::vector<std::pair<int32_t, double>> line_signs, double ndst, double tol);
        void sort_lines(std::map<int32_t, line::Line>&lines, std::map<int32_t, line::LR_adj>&line_adjs, double ndst, double tol); // tolerance = ndst percentage to separate partials

        std::pair<double, double> get_ylim(const std::map<int32_t, line::Line>& lines);
    }
    namespace gridding {
        void init_grid(const std::map<int32_t, line::LR_adj>& line_adjs, std::map<int32_t, std::shared_ptr<line_sep::GCol>>& grid, int ncells);
        void grid_lines(const std::map<int32_t, line::Line>& lines, const std::map<int32_t, line::LR_adj>& line_adjs, std::map<int32_t, std::shared_ptr<GCol>>& grid, double ymin_grid, double dy);
        std::tuple<int, double, double, double> get_number_cells(const std::map<int32_t, line::Line>& lines, double nomdst, double grid_tol);
    }
    namespace distance {
        void calc_line_dists(std::map<int32_t, line::Line>& lines, const std::map<int32_t, line::LR_adj>& line_adjs, const std::map<int32_t, std::shared_ptr<GCol>>& grid, double ymin_grid, double dy);
        void calc_dist(std::map<int32_t, line::Line>& lines, LR lr, const std::map<int32_t, std::shared_ptr<GCol>>& grid, int32_t lsymb_target, int32_t lsymb_adj, double ymin_grid, double dy);
    }  
}