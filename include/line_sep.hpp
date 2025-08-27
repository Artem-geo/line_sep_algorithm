#pragma once

#include <Eigen/Dense>
#include <map>
#include <memory>

namespace line_sep {
    class line::Line;
    struct line::LR_hood;

    using Point = std::pair<const double*, const double*>;
    using Cell = std::set<Point>;
    using GCol = std::vector<Cell>;

    double calc_azimuth(std::map<int32_t, line::Line>& lines);
    Eigen::RowVector2d get_pivot(std::map<int32_t, line::Line>& lines);
    void rotate_lines(std::map<int32_t, line::Line>& lines, double angle, const Eigen::RowVector2d& pivot);
    void project_lines(std::map<int32_t, line::Line>& lines, const Eigen::RowVector2d& pivot);

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

    void set_left_hood(std::map<int32_t, line::LR_hood>& lhood, std::vector<std::pair<int32_t, double>> line_signs, double ndst, double tol);
    void set_right_hood(std::map<int32_t, line::LR_hood>& lhood, std::vector<std::pair<int32_t, double>> line_signs, double ndst, double tol);
    void sort_lines(std::map<int32_t, line::Line>& lines, std::map<int32_t, line::LR_hood>& lhood, double ndst, double tol); // tolerance = ndst percentage to separate partials
    void grid_lines(const std::map<int32_t, line::Line>& lines, const std::map<int32_t, line::LR_hood>& lhood, std::map<int32_t, std::shared_ptr<GCol>>& grid, double nomdst);
    void create_grid(const std::map<int32_t, line::Line>& lines, const std::set<int32_t>& lh, std::map<int32_t, std::shared_ptr<GCol>>& grid, int num_cells, double dy, double ymin_grid);
    void calc_dist(std::map<int32_t, line::Line>& lines, const std::set<int32_t>& lh, std::map<int32_t, std::shared_ptr<GCol>>& grid);
    
    std::pair<double, double> get_ylim(const std::map<int32_t, line::Line>& lines);
}

// grid lines function
// calculate ymin and ymax to create a grid
// calculate dy of the grid as 5% of the nominal line spacing
// calculate number of cells in a grid column
// create a grid map where key = lsymb, value = shared_ptr to the grid
// for every line in the hood
    // if left hood set is not empty
        // create grid
    // if right hood set is not empty
        // create grid

// create grid function
// create a temp shared_ptr
// 