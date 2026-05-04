/**
 * @file convert_to_csv.cc
 * @author Kevin Scheuer
 * @brief
 * @date 2026-01-30
 *
 */

#include <algorithm>
#include <charconv>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <regex>
#include <string>
#include <vector>

#include "UTILITIES/FitConverter.h"
#include "IUAmpTools/report.h"

const char *kModule = "convert_to_csv";

// forward declarations
std::vector<std::string> sort_files(const std::vector<std::string> &files, int sort_index);
bool are_valid_fit_files(const std::vector<std::string> &files);

int main(int argc, char *argv[])
{
    std::vector<std::string_view> args(argv, argv + argc);

    // ==== PARSE ARGUMENTS ====
    // required arguments
    std::vector<std::string> input_files;
    std::string output_file = "";

    // optional arguments
    bool is_sorted = false;
    int sort_index = -1;
    bool is_acceptance_corrected = false;
    bool create_data_file = false;
    bool create_covariance = false;
    bool create_correlation = false;
    bool verbose = false;
    bool preview = false;
    std::string mass_branch = "M4Pi"; // TODO: change by building up from AmpTools 4-vectors. Instead option for lower vertex indices, and possibly isobar indices

    auto print_help = []()
    {
        std::cout << "Usage: convert_to_csv [-h] [-i INPUT_FILES] [-o OUTPUT_PATH]"
                  << " [-s] [--sort-index INDEX] [-a] [-d] [-m MASS_BRANCH] [-p] [-v]"
                  << " [--correlation] [--covariance] [-d]\n"
                  << "  -i INPUT_FILES:\t\tFull path to the .fit file(s)\n"
                  << "  -o OUTPUT_PATH:\t\tFull path to the output .csv file\n"
                  << "  -s, --sort:\t\t\tSort files by last number in the file name or path (default:true)\n"
                  << "  --sort-index INDEX:\t\tWhat number in file path to sort by (default:-1, so last number in path is used)\n"
                  << "  -a --acceptance-correct:\tWhether to acceptance correct intensities (default:true)\n"
                  << "  -d --data-file:\t\tCreate separate data file (default:false)\n"
                  << "  -m --mass-branch:\t\tBranch name for final mass spectrum of interest (default:M4Pi)\n"
                  << "  -p --preview:\t\t\tPrint files to be processed, but don't run (default:false)\n"
                  << "  -v --verbose:\t\t\tPrint information from converter scripts as they run\n"
                  << "  --correlation:\t\tSave correlation matrix to separate csv file\n"
                  << "  --covariance:\t\t\tSave covariance matrix to separate csv file\n"
                  << "  -h, --help:\t\t\tShow this help message\n";
    };

    // first check if help message is requested
    for (const auto &arg : args)
    {
        if (arg == "-h" || arg == "--help")
        {
            print_help();
            return 0;
        }
    }

    // helper to distinguish flags from negative numbers, since both use "-" character
    auto is_flag = [](std::string_view s)
    {
        return s.starts_with("-") && s.size() > 1 && !std::isdigit(static_cast<unsigned char>(s[1]));
    };

    // parse arguments
    for (size_t i = 1; i < args.size(); ++i)
    {
        auto arg = args[i];

        if (arg == "-i" || arg == "--input")
        {
            // collect all following arguments until next flag
            while (i + 1 < args.size() && !is_flag(args[i + 1]))
            {
                input_files.emplace_back(args[++i]);
            }
        }
        else if (arg == "-o" || arg == "--output")
        {
            if (i + 1 < args.size() && !is_flag(args[i + 1]))
            {
                output_file = std::string(args[++i]);
            }
            else
            {
                std::cerr << "Error: -o/--output requires a value. See help message.\n";
                print_help();
                return 1;
            }
        }
        else if (arg == "-s" || arg == "--sort")
        {
            is_sorted = true;
        }
        else if (arg == "--sort-index")
        {
            if (i + 1 < args.size() && !is_flag(args[i + 1]))
            {
                auto val = args[++i];
                auto [ptr, ec] = std::from_chars(val.data(), val.data() + val.size(), sort_index);
                if (ec != std::errc())
                {
                    std::cerr << "Error: --sort-index requires a valid integer. See help message.\n";
                    print_help();
                    return 1;
                }
            }
            else
            {
                std::cerr << "Error: --sort-index requires a value. See help message.\n";
                print_help();
                return 1;
            }
        }
        else if (arg == "-a" || arg == "--acceptance-correct")
        {
            is_acceptance_corrected = true;
        }
        else if (arg == "-d" || arg == "--data-file")
        {
            create_data_file = true;
        }
        else if (arg == "--covariance")
        {
            create_covariance = true;
        }
        else if (arg == "--correlation")
        {
            create_correlation = true;
        }
        else if (arg == "-v" || arg == "--verbose")
        {
            verbose = true;
        }
        else if (arg == "-p" || arg == "--preview")
        {
            preview = true;
        }
        else if (arg == "-m" || arg == "--mass-branch")
        {
            if (i + 1 < args.size() && !is_flag(args[i + 1]))
            {
                mass_branch = std::string(args[++i]);
            }
            else
            {
                std::cerr << "Error: -m/--mass-branch requires a value. See help message.\n";
                print_help();
                return 1;
            }
        }
        else if (is_flag(arg))
        {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_help();
            return 1;
        }
    }

    // Error checks for required arguments
    if (input_files.empty())
    {
        std::cerr << "Input files must be provided. See help message." << "\n";
        print_help();
        return 1;
    }
    if (output_file.empty())
    {
        std::cerr << "Output file must be provided. See help message." << "\n";
        print_help();
        return 1;
    }

    // ==== VALIDATE INPUTS ====

    // ensure only .fit files are provided
    if (!are_valid_fit_files(input_files))
    {
        std::cerr << "One or more input files are not valid .fit files.\n";
        return 1;
    }

    // guarantee that output file ends with .csv
    if (std::filesystem::path(output_file).extension() != ".csv")
    {
        output_file += ".csv";
    }

    // sort files if requested
    if (is_sorted)
        input_files = sort_files(input_files, sort_index);

    if (preview)
    {
        std::cout << "Preview mode enabled. Printing files and exiting\n";
        for (const auto &file : input_files)
        {
            std::cout << "\t" << file << "\n";
        }
        return 0;
    }

    // ==== PROCESS FILES AND WRITE TO CSV ====
    std::stringstream csv_result_data, csv_cov_data, csv_corr_data;
    bool header_written = false;
    std::string first_file_result_header, first_file_cov_header, first_file_corr_header;

    // extract fit results, and optionally the cov/corr matrices
    for (const auto &file : input_files)
    {
        if (verbose)
            report(INFO, kModule) << "Processing file: " << file << "\n";
        FitConverter converter(file, is_acceptance_corrected, !verbose);

        if (!header_written)
        {
            std::string header = converter.getCSVHeader();
            csv_result_data << header << "\n";
            header_written = true;
            first_file_result_header = header;

            if (create_covariance)
            {
                std::string cov_header = converter.getCSVCovarianceMatrixHeader();
                csv_cov_data << cov_header << "\n";
                first_file_cov_header = cov_header;
            }

            if (create_correlation)
            {
                std::string corr_header = converter.getCSVCorrelationMatrixHeader();
                csv_corr_data << corr_header << "\n";
                first_file_corr_header = corr_header;
            }
        }
        else
        { 
            // this block ensures all files have the same csv header format. If this is
            // not the case, then the csv files will be malformed and difficult to work 
            // with, so we enforce that here before processing any files
            std::string header = converter.getCSVHeader();
            if (header != first_file_result_header)
            {
                report(ERROR, kModule) << "File "
                                       << file
                                       << " has a different CSV header than the first file. Aborting\n";
                return 1;
            }

            if (create_covariance)
            {
                std::string cov_header = converter.getCSVCovarianceMatrixHeader();
                if (cov_header != first_file_cov_header)
                {
                    report(ERROR, kModule) << "File "
                                           << file
                                           << " has a different covariance CSV header than the first file. Aborting\n";
                    return 1;
                }
            }
            if (create_correlation)
            {
                std::string corr_header = converter.getCSVCorrelationMatrixHeader();
                if (corr_header != first_file_corr_header)
                {
                    report(ERROR, kModule) << "File "
                                           << file
                                           << " has a different correlation CSV header than the first file. Aborting\n";
                    return 1;
                }
            }
        }
        csv_result_data << converter.getCSVRow() << "\n";

        if (create_covariance)
            csv_cov_data << converter.getCSVCovarianceMatrix() << "\n";
        if (create_correlation)
            csv_corr_data << converter.getCSVCorrelationMatrix() << "\n";
    }

    // write the csv data to the output file
    std::ofstream result_file(output_file);
    result_file << csv_result_data.str();
    result_file.close();

    // covariance and correlation matrix files are written to separate files in the same
    // directory as the main output file, with suffixes _covariance.csv and 
    // _correlation.csv
    if (create_covariance)
    {        
        std::filesystem::path output_path(output_file);
        std::string covariance_file = (output_path.parent_path() / (output_path.stem().string() + "_covariance.csv")).string();
        std::ofstream cov_file(covariance_file);
        std::cout << covariance_file << "\n";
        cov_file << csv_cov_data.str();
        cov_file.close();
    }
    if (create_correlation)
    {
        std::filesystem::path output_path(output_file);
        std::string correlation_file = (output_path.parent_path() / (output_path.stem().string() + "_correlation.csv")).string();
        std::ofstream corr_file(correlation_file);
        corr_file << csv_corr_data.str();
        corr_file.close();
    }

    if (create_data_file)
        ; // TODO: execute extract_data to create data file. Use output_file_data.csv as output name

    return 0;
}

/**
 * @brief Sort files based on a number extracted from the file path
 *
 * Extracts all numbers (including decimals) from each file path and sorts based on
 * the number at the specified index position.
 *
 * @param[in] files Vector of file paths to sort
 * @param[in] sort_index Index of the number to sort by (supports negative indexing).
 *  Default -1 means use the last number found in the path.
 * @return std::vector<std::string> Sorted vector of file paths
 */
std::vector<std::string> sort_files(const std::vector<std::string> &files, int sort_index)
{
    auto extract_number = [sort_index](const std::string &path) -> double
    {
        // Regex pattern to match integers and decimal numbers
        std::regex number_pattern(R"(\d*\.?\d+)");
        std::vector<double> numbers;

        // Find all numbers in the path
        auto begin = std::sregex_iterator(path.begin(), path.end(), number_pattern);
        auto end = std::sregex_iterator();

        for (auto it = begin; it != end; ++it)
        {
            numbers.push_back(std::stod(it->str()));
        }

        if (numbers.empty())
        {
            return std::numeric_limits<double>::infinity();
        }

        int index = sort_index;

        // handle negative indexing
        if (index < 0)
        {
            index = static_cast<int>(numbers.size()) + index;
        }

        // Return the number at the specified index, or infinity if out of bounds
        if (index >= 0 && index < static_cast<int>(numbers.size()))
        {
            return numbers[index];
        }

        return std::numeric_limits<double>::infinity();
    };

    // Create a copy and sort it
    std::vector<std::string> sorted_files = files;
    std::sort(sorted_files.begin(), sorted_files.end(),
              [&extract_number](const std::string &a, const std::string &b)
              {
                  return extract_number(a) < extract_number(b);
              });

    return sorted_files;
}

/**
 * @brief Check if all provided files are .fit files
 *
 * @param[in] files Vector of file paths to check
 * @return true if all files have a .fit extension and exist, false otherwise
 * @return false if any file is not a .fit file or does not exist
 */
bool are_valid_fit_files(const std::vector<std::string> &files)
{
    for (const auto &file : files)
    {
        if (std::filesystem::path(file).extension() != ".fit" || !std::filesystem::exists(file))
        {
            return false;
        }
    }
    return true;
}