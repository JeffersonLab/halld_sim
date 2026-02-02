/**
 * @file FitConverter.h
 * @author Kevin Scheuer
 * @brief Class for converting AmpTools fit results to any format
 * @date 2026-01-31
 *
 * Note that this currently is built to convert to CSV, but the class can be extended
 * to convert to any other format by accessing the various containers that store
 * the fit results.
 */

#ifndef FIT_CONVERTER_H
#define FIT_CONVERTER_H

#include <complex>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "IUAmpTools/FitResults.h"

/**
 * @class FitConverter
 * @brief Converts AmpTools FitResults to any format
 *
 * This class extracts fit results from a single AmpTools .fit file. It can then output
 * the data in various formats, such as CSV. The extracted data includes:
 * - Standard fit outputs (likelihood, event counts, status codes)
 * - AmpTools parameters
 * - Real/Imaginary parts of production coefficients
 * - Coherent sums of amplitudes by quantum numbers
 * - Phase differences
 *
 * Example usage in convert_to_csv.cc:
 * @code
 * std::ofstream csv_file(output_file);
 * std::stringstream csv_data;
 * bool header_written = false;
 *
 * for (const auto& file : input_files) {
 *     FitConverter converter(file, acceptance_corrected);
 *
 *     if (!header_written) {
 *         csv_data << converter.getCSVHeader() << "\n";
 *         header_written = true;
 *     }
 *     csv_data << converter.getCSVRow() << "\n";
 * }
 *
 * csv_file << csv_data.str();
 * csv_file.close();
 * @endcode
 */
class FitConverter
{
public:
    /**
     * @brief Construct a FitConverter for a single fit file
     * @param[in] filename Path to the .fit file
     * @param[in] acceptance_correct Whether to use acceptance-corrected intensities
     * @throw ERROR If the fit file cannot be read or is invalid
     */
    FitConverter(const std::string &filename, const bool &acceptance_correct = false);

    /**
     * @brief Extract fit data from the .fit file
     */
    void extract();

    // TODO: implement these methods
    std::string getCSVHeader() const;
    std::string getCSVRow() const;

    /**
     * @brief return set of all unique amplitude names in a fit result
     *
     * A "full" amplitude is a "reaction::sum::ampName" string, so this extracts all
     * unique "ampName" strings
     */
    std::set<std::string> uniqueAmpNames() const;

    bool isValid() const { return m_fit_results.valid(); }

    /**
     * @brief Map of AmpTool standard fit outputs and their values
     *
     * Standard outputs include values such as:
     * - eMatrixStatus
     * - lastMinuitCommandStatus
     * - likelihood
     * - intensity
     */
    const std::map<std::string, double> &standardResults() const { return m_standard_results; }

    /**
     * @brief Map of AmpTools parameters and their values
     *
     * Ignores production coefficients, which typically have "::" in their names
     */
    const std::map<std::string, double> &parameters() const { return m_parameters; }

    /**
     * @brief Map of unique amplitude names and their related constrained amplitudes
     *
     * A "full" amplitude is a "reaction::sum::ampName" string, and often amplitudes
     * are constrained across reactions and sums. This map groups all such constrained
     * amplitudes by their common "ampName".
     * "myAmpName" -> {
     *      "reaction1::sum1::myAmpName",
     *      "reaction1::sum2::myAmpName",
     *      "reaction2::sum1::myAmpName",
     * ...}
     * Note that the self-constraint is included in the list, so even amplitudes with
     * no constraints will have a vector of size one (itself).
     */
    const std::map<std::string, std::vector<std::string>> &
    constrainedAmps() const { return m_constrained_amps; }

    /**
     * @brief Map of single amplitude intensities to their values and errors
     *
     * Contains the intensity (and error) for each unique amplitude name,
     * calculated using all the constrained amplitudes for that name.
     */
    const std::map<std::string, std::pair<double, double>> &
    singleAmpIntensities() const { return m_single_amp_intensities; }

    /**
     * @brief Map of unique amplitude names to their production coefficients
     *
     * Contains complex production coefficients for each unique amplitude name
     * extracted from the fit results.
     */
    const std::map<std::string, std::complex<double>> &
    productionCoefficients() const { return m_production_coefficients; }

    /**
     * @brief Map of full amplitude pairs to their phase differences and errors
     *
     * Only includes pairs of unique amplitude names that share a reaction and sum. For
     * 2 amplitudes "amp1" and "amp2", we search through every combination in
     * m_constrained_amps[amp1] and m_constrained_amps[amp2] to find pairs that share
     * the same "reaction::sum". The key is a vector of the two full amplitude strings,
     * and the value is a pair of the phase difference and its error.
     */
    const std::map<std::pair<std::string, std::string>, std::pair<double, double>> &
    phaseDifferences() const { return m_phase_differences; }

    ~FitConverter() = default;

private:
    const std::string m_fit_file;         ///< Path to the .fit file
    const FitResults m_fit_results;       ///< AmpTools FitResults object
    const ConfigurationInfo *m_cfg_info;  ///< ConfigurationInfo from the fit
    const bool m_is_acceptance_corrected; ///< Whether to use acceptance-corrected intensities

    // see accessor methods for descriptions
    std::map<std::string, double> m_standard_results;
    std::map<std::string, double> m_parameters;
    std::map<std::string, std::vector<std::string>> m_constrained_amps;
    std::map<std::string, std::pair<double, double>> m_single_amp_intensities;
    std::map<std::string, std::complex<double>> m_production_coefficients;
    std::map<std::pair<std::string, std::string>, std::pair<double, double>> m_phase_differences;

    /**
     * @brief Finds and returns a map of unique amplitude names to their constrained full amplitudes
     * 
     * This assumes that all terms sharing the same amplitude name are constrained to each other.
     * 
     * @throw ERROR If the above assumption is violated
     */
    std::map<std::string, std::vector<std::string>> findConstrainedAmps() const;

    static const char *kModule;
};

#endif // FIT_CONVERTER_H