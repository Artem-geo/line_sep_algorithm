#include "line.h"
#include "misc.h"
#include <chrono>
#include <filesystem>
#include <iostream>

int main()
{
    using namespace line;
    using namespace misc;
    using namespace misc::io;

    try {
        std::string fpath = "C:\\Dev\\line_sep_algorithm\\qgis\\points_315_approx_symb.csv";
        std::string appdx = "_rotated";
        std::filesystem::path filepath(fpath);
        std::filesystem::path ofilepath = filepath;
        std::filesystem::path ofilename = filepath.stem().concat(appdx);
        std::filesystem::path oextension = filepath.extension();

        ofilepath.replace_filename(ofilename);
        ofilepath.replace_extension(oextension);
        
        std::set<int> lnums;
        std::map<int, VLines> lines;
        
        load_lines(filepath, lines, lnums);

        double azimuth = calc_azimuth(lines);
        rotate_lines(lines, azimuth);

        //safe_lines(ofilepath, lines);
    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }
}