/**
 * @file FitConverter.h
 * @author Kevin Scheuer
 * @brief Class for converting AmpTools fit results to any format
 * @date 2026-05-08
 *
 * @note that this currently is built to convert to CSV, but the class can be extended
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

#include "AmplitudeParser.h"
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
 * - Intensities of unique amplitudes. See constrainedAmps() for how unique amplitudes
 *   groups are defined.
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
     * @param[in] mute_warning Suppress warning about free parameters if true
     * @param[in] naming_scheme Amplitude naming scheme: auto, JLme, eJPmL, or Lme
     * @throw ERROR If the fit file cannot be read or is invalid
     */
    FitConverter(
        const std::string &filename,
        const bool &acceptance_correct = false,
        bool mute_warning = false,
        const std::string &naming_scheme = "auto");

    /**
     * @brief Extract fit data from the .fit file
     */
    void extract();

    /**
     * @brief Get CSV header string for the fit results
     *
     * Constructs a comma-separated header string based on the keys of the various
     * containers storing the fit results.
     *
     * @note It is crucial that the order of entries in the header matches the order of
     * values in getCSVRow().
     *
     * @return Comma-separated header string
     */
    std::string getCSVHeader() const;

    /**
     * @brief Get CSV row string for the fit results
     *
     * Constructs a comma-separated row string based on the values of the various
     * containers storing the fit results.
     *
     * @note It is crucial that the order of values in this row matches the order of
     * entries in getCSVHeader().
     *
     * @return Comma-separated row string
     */
    std::string getCSVRow() const;

    /**
     * @brief Get the CSV header string for the covariance matrix
     *
     * @return std::string "file, parameter, par1, par2, ..., parN"
     */
    std::string getCSVCovarianceMatrixHeader() const;

    /**
     * @brief Get the CSV representation of the covariance matrix
     *
     * Every row corresponds to a parameter, with the first two columns being the file
     * and parameter name, followed by the covariance values with all other parameters.
     */
    std::string getCSVCovarianceMatrix() const;

    /**
     * @brief Get the CSV header string for the correlation matrix
     *
     * @return std::string "file, parameter, par1, par2, ..., parN"
     */
    std::string getCSVCorrelationMatrixHeader() const;

    /**
     * @brief Get the CSV representation of the correlation matrix
     *
     * Every row corresponds to a parameter, with the first two columns being the file
     * and parameter name, followed by the correlation values with all other parameters.
     */
    std::string getCSVCorrelationMatrix() const;

    /**
     * @brief Get the CSV header string for the normalization integral matrix
     *
     * @return std::string "file, amplitude, amp1, amp2, ..., ampN"
     */
    std::string getCSVNormIntMatrixHeader() const;

    /**
     * @brief Get the CSV representation of the normalization integral matrix
     *
     * Every row corresponds to an amplitude, with the first two columns being the file
     * and amplitude names, followed by complex normalization integral values with all
     * other amplitudes. Complex values are formatted as pandas-readable strings (e.g.,
     * "1+2j").
     */
    std::string getCSVNormIntMatrix() const;

    /**
     * @brief return set of all unique amplitude names in a fit result
     *
     * A "full" amplitude is a "reaction::sum::ampName" string, so this extracts all
     * unique "ampName" strings
     */
    std::set<std::string> uniqueAmpNames() const;

    /**
     * @brief return set of all unique coherent sum + amp names in a fit result
     *
     * A "full" amplitude is a "reaction::sum::ampName" string, so this extracts all
     * unique "sum::ampName" strings
     *
     * @return std::set<std::string>
     */
    std::set<std::string> uniqueSumAmpNames() const;

    bool isValid() const { return m_fit_results.valid(); }

    /**
     * @brief Map of AmpTools standard fit outputs and their values
     *
     * Standard outputs include values such as:
     * - eMatrixStatus
     * - lastMinuitCommandStatus
     * - likelihood
     * - intensity
     */
    const std::map<std::string, double> &standardResults() const { return m_standard_results; }

    /**
     * @brief Map of AmpTools parameters and their values and errors
     *
     * Ignores production coefficients, which typically have "::" in their names
     */
    const std::map<std::string, std::pair<double, double>> &parameters() const { return m_parameters; }

    /**
     * @brief Map of unique amplitude names and their related constrained amplitudes
     *
     * A "full" amplitude is a "reaction::sum::ampName" string, and often amplitudes
     * are constrained across reactions and sums. The goal is to identify the shortest
     * identifying string that can be used to group amplitudes that are constrained to
     * each other. It will first attempt to group by "ampName", then by "sum::ampName",
     * and if neither of those assumptions hold, it will just group by the full
     * amplitude string.
     *
     * For example, if we have the following full amplitude strings:
     * - "reaction1::sum1::myAmpName"
     * - "reaction1::sum2::myAmpName"
     * - "reaction2::sum1::myAmpName"
     *
     * We first check if all amplitudes with the same "ampName" (myAmpName) are
     * constrained to each other, and if so, we group them under the unique amplitude
     * name "myAmpName" and save the intensity + error for that group. If that
     * assumption does not hold, we check if all amplitudes with the same "sum::ampName"
     * are constrained to each other (i.e. across reactions) and if so, we group them
     * under the unique amplitude name "sum::ampName". If neither assumption holds, we
     * just group by the full amplitude name.
     *
     * @note that the self-constraint is included in the list, so even amplitudes with
     * no constraints will have a vector of size one (itself).
     */
    const std::map<std::string, std::vector<std::string>> &
    constrainedAmps() const { return m_constrained_amps; }

    /**
     * @brief Map of unique amplitude intensities to their values and errors
     *
     * Contains the intensity (and error) for each unique amplitude grouping,
     * calculated using all the constrained amplitudes for that group.
     */
    const std::map<std::string, std::pair<double, double>> &
    uniqueAmpIntensities() const { return m_unique_amp_intensities; }

    /**
     * @brief Map of coherent-sum labels to their values and errors
     *
     * The labels are built from the selected amplitude naming scheme and the
     * configured coherent-sum groupings
     */
    const std::map<std::string, std::pair<double, double>> &
    coherentSumIntensities() const { return m_coherent_sum_intensities; }

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
     * the same "reaction::sum". The key is a size 2 vector of the full amplitude
     * strings, and the value is a pair of the phase difference and its error.
     */
    const std::map<std::pair<std::string, std::string>, std::pair<double, double>> &
    phaseDifferences() const { return m_phase_differences; }

    /**
     * @brief Get the covariance matrix of the fit parameters.
     *
     * The error matrix is a 2D vector where each entry [i][j] corresponds to
     * the covariance between parameter i and parameter j. The order of the rows and
     * columns matches the order of parameters in FitResults::parNameList(). Note
     * that we do include production parameters in this matrix, as opposed to the
     * m_parameters map which excludes them.
     */
    const std::vector<std::vector<double>> &
    errorMatrix() const { return m_error_matrix; }

    ~FitConverter() = default;

private:
    const std::string m_fit_file;             ///< Path to the .fit file
    const FitResults m_fit_results;           ///< AmpTools FitResults object
    const ConfigurationInfo *m_cfg_info;      ///< ConfigurationInfo from the fit
    const bool m_is_acceptance_corrected;     ///< Whether to use acceptance-corrected intensities
    const AmplitudeParser m_amplitude_parser; ///< Parser for amplitude labels and sums

    // see accessor methods for descriptions
    std::map<std::string, double> m_standard_results;
    std::map<std::string, std::pair<double, double>> m_parameters;
    std::map<std::string, std::vector<std::string>> m_constrained_amps;
    std::map<std::string, std::pair<double, double>> m_unique_amp_intensities;
    std::map<std::string, std::pair<double, double>> m_coherent_sum_intensities;
    std::map<std::string, std::complex<double>> m_production_coefficients;
    std::map<std::pair<std::string, std::string>, std::pair<double, double>> m_phase_differences;
    std::vector<std::vector<double>> m_error_matrix;

    /**
     * @brief Get map of unique amplitude groups to their constrained full amplitudes
     *
     * This attempts to identify the unique strings that can be a coupled to a set of
     * constrained amplitudes. For a "full" amplitude string of the form
     * "reaction::sum::ampName", it first assumes that amplitudes with the same ampName
     * are constrained to each other, across reactions and sums. If not, it peels back a
     * layer and assumes that amplitudes with the same "sum::ampName" are constrained to
     * each other, across reactions. If that fails, it simply saves a mapping of the
     * full amplitude to itself. The map keys are whats extracted as the unique
     * amplitude names (ampName, sum::ampName, or full amplitude), and the values are
     * vectors of the full amplitude strings that are constrained to each other.
     *
     * @note the keys are whats written to the CSV, and is why it first attempts to
     * extract just the ampName for simplicity.
     */
    std::map<std::string, std::vector<std::string>> findConstrainedAmps() const;

    /**
     * @brief Checks if amplitude names are constrained to each other
     *
     * A "full" amplitude is a "reaction::sum::ampName" string, so this checks if all
     * amplitudes with the same "ampName" are constrained to each other across reactions
     * and sums.
     */
    bool ampNamesAreConstrained() const;

    /**
     * @brief Checks if coherent sum + amplitude names are constrained to each other
     *
     * A "full" amplitude is a "reaction::sum::ampName" string, so this checks if all
     * amplitudes with the same "sum::ampName" are constrained to each other across
     * reactions.
     */
    bool sumAmpNamesAreConstrained() const;

    /**
     * @brief Get the reaction string from a full amplitude name
     *
     * This extracts the reaction part from a "reaction::sum::ampName" string.
     *
     * @param full_amp_name The full amplitude name
     * @return The reaction string
     */
    static std::string getReactionString(const std::string &full_amp_name);

    static const char *kModule;
};

#endif // FIT_CONVERTER_H