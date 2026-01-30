/**
 * @file convert_to_csv.cc
 * @author Kevin Scheuer
 * @brief
 * @date 2026-01-30
 *
 * TODO:
 *  this is basically replacement for python script, so it will need to loop over args
 *  but instead provide a vector of files to the individual scripts
 */

#include <string>
#include <vector>

int main(std::vector<std::string_view> args(argv, argv + argc))
{

    // parse args
    std::vector<std::string> input_files;
    std::string output_file;
    for (const auto &arg : args)
    {

        // parse input files
        if (arg.find("--input") != std::string_view::npos || arg.find("-i") != std::string_view::npos)
        {
            // Skip the flag itself and collect all following arguments until next flag
            auto it = std::find(args.begin(), args.end(), arg);
            if (it != args.end()) {
                ++it;
                while (it != args.end() && !it->starts_with("-")) {
                    input_files.push_back(*it);
                    ++it;
                }
            }
        }
        else if (arg.find("--output") != std::string_view::npos || arg.find("-o") != std::string_view::npos)
        {
            // parse output file
        }
    }

    // TEMP: print out input files
    for (const auto &file : input_files)
    {
        std::cout << "Input file: " << file << std::endl;
    }

    return 0;
}