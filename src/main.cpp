#include "line.hpp"
#include "line_sep.hpp"
#include "misc.hpp"
#include <chrono>
#include <filesystem>
#include <iostream>

int main()
{
    try {
        std::string fpath = "C:\\Dev\\line_sep_algorithm\\qgis\\points_315_approx_symb.csv";
        std::string appdx = "_rotated";
        std::filesystem::path filepath(fpath);
        std::filesystem::path ofilepath = filepath;
        std::filesystem::path ofilename = filepath.stem().concat(appdx);
        std::filesystem::path oextension = filepath.extension();

        ofilepath.replace_filename(ofilename);
        ofilepath.replace_extension(oextension);

        std::map<int32_t, line::Line> lines;
        misc::io::load_lines(filepath, lines);

        double azimuth = line_sep::geometry::calc_azimuth(lines);
        Eigen::RowVector2d pivot = line_sep::geometry::get_pivot(lines);

        line_sep::geometry::rotate_lines(lines, azimuth, pivot);
        line_sep::geometry::project_lines(lines, pivot); // projection is written to proj attribute of a line

        std::map<int32_t, line::LR_adj> line_adjs;
        double ndst = 500.0; // nominal line separation
        line_sep::geometry::sort_lines(lines, line_adjs, ndst, 0.15);


        //std::map<int32_t, std::shared_ptr<line_sep::GCol>> grid;
        //line_sep::gridding::init_grid(line_adjs, grid, 300);
        //auto grid_params = line_sep::gridding::grid_lines(lines, line_adjs, grid, 300);


        //line_sep::calc_line_dists(lines, lhood, grid, grid_params);
        //misc::io::safe_grid(ofilepath, lines, grid);
        //misc::io::safe_lines(ofilepath, lines);

        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        //std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;

    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }
}