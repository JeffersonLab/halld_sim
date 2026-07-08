/**
 * @file convert_to_csv.cc
 * @author Kevin Scheuer
 * @brief A utility for converting AmpTools .fit files to CSV format
 * @date 2026-05-08
 *
 * This program converts AmpTools fit results in .fit files to CSV format. It can also
 * optionally create separate CSV files for the covariance, correlation, and
 * normalization integral matrices, as well as a CSV file containing the data used in
 * the fit. The program is designed to be flexible and can handle various amplitude
 * naming schemes, which it uses to group amplitudes into coherent sums based on shared
 * quantum numbers. The output CSV files can then be easily imported into just about any
 * data analysis or plotting software.
 *
 * Print the help message with
 * @verbatim
 * user@host:~$ convert_to_csv -h
 * @endverbatim
 * for more details.
 */

#include <algorithm>
#include <charconv>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <memory>
#include <sstream>
#include <limits>
#include <regex>
#include <string>
#include <vector>

#include "UTILITIES/FitConverter.h"
#include "UTILITIES/RootDataConverter.h"
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
    std::vector<int> lower_vertex_indices;
    std::string naming_scheme = "auto";
    bool create_data_file = false;
    bool create_covariance = false;
    bool create_correlation = false;
    bool create_norm_int = false;
    bool verbose = false;
    bool preview = false;

    auto print_help = []()
    {
        report(INFO, kModule) << "Usage: convert_to_csv [-h] [-i INPUT_FILES] [-o OUTPUT_PATH] [-s] [--sort-index INDEX] [-a] [-d] [-l LOWER_VERTEX_INDICES] [-n NAMING_SCHEME] [-p] [-v] [--correlation] [--covariance] [--norm-int] \n";
        report(INFO, kModule) << "  -i INPUT_FILES:\t\tFull path to the .fit file(s)\n";
        report(INFO, kModule) << "  -o OUTPUT_PATH:\t\tFull path to the output .csv file\n";
        report(INFO, kModule) << "  -s, --sort:\t\t\tSort files by last number in the file name or path (default:true)\n";
        report(INFO, kModule) << "  --sort-index INDEX:\t\tWhat number in file path to sort by (default:-1, so last number in path is used)\n";
        report(INFO, kModule) << "  -a --acceptance-correct:\tAcceptance correct the intensities\n";
        report(INFO, kModule) << "  -d --data-file:\t\tCreate separate data file\n";
        report(INFO, kModule) << "  -l --lower-vertex-indices:\tSpecify indices of lower vertex particles in the data 4-vectors\n";
        report(INFO, kModule) << "  -n --naming-scheme:\t\tAmp naming scheme for coherent sum groupings. Default is auto (infers scheme from the amp names)\n";
        report(INFO, kModule) << "  -p --preview:\t\t\tPrint files to be processed, but don't run\n";
        report(INFO, kModule) << "  -v --verbose:\t\t\tPrint information from converter scripts as they run\n";
        report(INFO, kModule) << "  --correlation:\t\tSave correlation matrix to separate csv file\n";
        report(INFO, kModule) << "  --covariance:\t\t\tSave covariance matrix to separate csv file\n";
        report(INFO, kModule) << "  --norm-int:\t\t\tSave normalization integral matrix to separate csv file\n";
        report(INFO, kModule) << "  -h, --help:\t\t\tShow this help message\n";
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
                report(ERROR, kModule) << "Error: -o/--output requires a value. See help message.\n";
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
                    report(ERROR, kModule) << "Error: --sort-index requires a valid integer. See help message.\n";
                    print_help();
                    return 1;
                }
            }
            else
            {
                report(ERROR, kModule) << "Error: --sort-index requires a value. See help message.\n";
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
        else if (arg == "--norm-int")
        {
            create_norm_int = true;
        }
        else if (arg == "-v" || arg == "--verbose")
        {
            verbose = true;
        }
        else if (arg == "-p" || arg == "--preview")
        {
            preview = true;
        }
        else if (arg == "-l" || arg == "--lower-vertex-indices")
        {
            // collect all following arguments until next flag
            while (i + 1 < args.size() && !is_flag(args[i + 1]))
            {
                auto val = args[++i];
                int index;
                auto [ptr, ec] = std::from_chars(val.data(), val.data() + val.size(), index);
                if (ec != std::errc())
                {
                    report(ERROR, kModule) << "Error: -l/--lower-vertex-indices requires valid integers. See help message.\n";
                    print_help();
                    return 1;
                }
                lower_vertex_indices.push_back(index);
            }
        }
        else if (arg == "-n" || arg == "--naming-scheme")
        {
            if (i + 1 < args.size() && !is_flag(args[i + 1]))
            {
                naming_scheme = std::string(args[++i]);
            }
            else
            {
                report(ERROR, kModule) << "Error: -n/--naming-scheme requires a value. See help message.\n";
                print_help();
                return 1;
            }
        }
        else if (is_flag(arg))
        {
            report(ERROR, kModule) << "Unknown argument: " << arg << "\n";
            print_help();
            return 1;
        }
    }

    // ==== VALIDATE INPUTS ====

    // Error checks for required arguments
    if (input_files.empty())
    {
        report(ERROR, kModule) << "Input files must be provided. See help message." << "\n";
        print_help();
        return 1;
    }
    if (output_file.empty())
    {
        report(ERROR, kModule) << "Output file must be provided. See help message." << "\n";
        print_help();
        return 1;
    }

    // ensure only .fit files are provided
    if (!are_valid_fit_files(input_files))
    {
        report(ERROR, kModule) << "One or more input files are not valid .fit files.\n";
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
        report(INFO, kModule) << "Preview mode enabled. Printing files and exiting\n";
        for (const auto &file : input_files)
        {
            report(INFO, kModule) << "\t" << file << "\n";
        }
        return 0;
    }

    // ==== PROCESS FILES AND WRITE TO CSV ====
    std::stringstream csv_result_data, csv_cov_data, csv_corr_data, csv_norm_int_data, csv_root_data;
    bool header_written = false;
    std::string first_file_result_header, first_file_cov_header, first_file_corr_header, first_file_norm_int_header, first_file_data_header;

    for (const auto &file : input_files)
    {
        if (verbose)
            report(INFO, kModule) << "Processing file: " << file << "\n";
        FitConverter converter(file, is_acceptance_corrected, !verbose, naming_scheme);
        std::unique_ptr<RootDataConverter> data_converter;
        if (create_data_file)
        {
            if (!lower_vertex_indices.empty())
            {
                data_converter = std::make_unique<RootDataConverter>(file, lower_vertex_indices, !verbose);
            }
            else
            {
                data_converter = std::make_unique<RootDataConverter>(file, !verbose);
            }
        }

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

            if (create_norm_int)
            {
                std::string norm_int_header = converter.getCSVNormIntMatrixHeader();
                csv_norm_int_data << norm_int_header << "\n";
                first_file_norm_int_header = norm_int_header;
            }

            if (create_data_file)
            {
                std::string data_header = data_converter->getCSVHeader();
                csv_root_data << data_header << "\n";
                first_file_data_header = data_header;
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
            if (create_norm_int)
            {
                std::string norm_int_header = converter.getCSVNormIntMatrixHeader();
                if (norm_int_header != first_file_norm_int_header)
                {
                    report(ERROR, kModule) << "File "
                                           << file
                                           << " has a different normalization integral CSV header than the first file. Aborting\n";
                    return 1;
                }
            }
            if (create_data_file)
            {
                std::string data_header = data_converter->getCSVHeader();
                if (data_header != first_file_data_header)
                {
                    report(ERROR, kModule) << "File "
                                           << file
                                           << " has a different data CSV header than the first file. Aborting\n";
                    return 1;
                }
            }
        }
        csv_result_data << converter.getCSVRow() << "\n";

        if (create_covariance)
            csv_cov_data << converter.getCSVCovarianceMatrix() << "\n";
        if (create_correlation)
            csv_corr_data << converter.getCSVCorrelationMatrix() << "\n";
        if (create_norm_int)
            csv_norm_int_data << converter.getCSVNormIntMatrix() << "\n";
        if (create_data_file)
            csv_root_data << data_converter->getCSVRow() << "\n";

        // loop to next input file to make a new row
    }

    // write the csv data to the output file
    std::ofstream result_file(output_file);
    result_file << csv_result_data.str();
    result_file.close();

    // supplemental files are written to separate files in the same directory as the
    // main output file, with appropriate suffixes e.g. "_data.csv", "correlation.csv"
    if (create_covariance)
    {
        std::filesystem::path output_path(output_file);
        std::string covariance_file = (output_path.parent_path() / (output_path.stem().string() + "_covariance.csv")).string();
        std::ofstream cov_file(covariance_file);
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
    if (create_norm_int)
    {
        std::filesystem::path output_path(output_file);
        std::string norm_int_file = (output_path.parent_path() / (output_path.stem().string() + "_norm_int.csv")).string();
        std::ofstream norm_int_out_file(norm_int_file);
        norm_int_out_file << csv_norm_int_data.str();
        norm_int_out_file.close();
    }
    if (create_data_file)
    {
        std::filesystem::path output_path(output_file);
        std::string data_file = (output_path.parent_path() / (output_path.stem().string() + "_data.csv")).string();
        std::ofstream data_out_file(data_file);
        data_out_file << csv_root_data.str();
        data_out_file.close();
    }

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