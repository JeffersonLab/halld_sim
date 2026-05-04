/**
 * @file FitConverter.cc
 * @author Kevin Scheuer
 * @brief Implementation of FitConverter class for converting AmpTools fit results to CSV
 * @date 2026-02-01
 *
 * TODO:
 * - Add coherent sum extraction. This is amplitude name dependent, so may need custom
 *   AmplitudeParsers for different naming schemes.
 */

#include <cassert>

#include "FitConverter.h"
#include "AmplitudeParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/report.h"

const char *FitConverter::kModule = "FitConverter";

FitConverter::FitConverter(const std::string &filename, const bool &acceptance_correct, bool mute_warning) : m_fit_file(filename),
                                                                                                             m_fit_results(filename, mute_warning),
                                                                                                             m_cfg_info(m_fit_results.configInfo()),
                                                                                                             m_is_acceptance_corrected(acceptance_correct),
                                                                                                             m_error_matrix(m_fit_results.errorMatrix())
{
    if (!m_fit_results.valid())
    {
        report(ERROR, kModule) << "FitConverter ERROR: Unable to read fit results from file " + m_fit_file << "\n";
        assert(false);
    }
    extract();
}

void FitConverter::extract()
{
    // Store the standard AmpTools fit outputs
    m_standard_results["eMatrixStatus"] = m_fit_results.eMatrixStatus();
    m_standard_results["lastMinuitCommandStatus"] = m_fit_results.lastMinuitCommandStatus();
    m_standard_results["likelihood"] = m_fit_results.likelihood();
    m_standard_results["intensity"] = m_fit_results.intensity(false).first;
    m_standard_results["intensity_err"] = m_fit_results.intensity(false).second;
    m_standard_results["ac_intensity"] = m_fit_results.intensity(true).first;
    m_standard_results["ac_intensity_err"] = m_fit_results.intensity(true).second;
    for (const auto &pair : m_standard_results)
    {
        report(DEBUG, kModule) << "Standard result: " << pair.first << " = " << pair.second << "\n";
    }

    // Store parameters, but ignore production coefficients (typically have "::" in name)
    for (const auto &par_name : m_fit_results.parNameList())
    {
        if (par_name.find("::") == std::string::npos)
            m_parameters[par_name] = {m_fit_results.parValue(par_name), m_fit_results.parError(par_name)};
    }
    for (const auto &pair : m_parameters)
    {
        report(DEBUG, kModule) << "Parameter: " << pair.first << " = " << pair.second.first << " +/- " << pair.second.second << "\n";
    }

    // Build map of unique amplitude names and their constrained amplitudes
    m_constrained_amps = findConstrainedAmps();
    for (const auto &pair : m_constrained_amps)
    {
        report(DEBUG, kModule) << "Unique amplitude: "
                               << pair.first
                               << " has constrained amplitudes:\n";
        for (const auto &amp : pair.second)
            report(DEBUG, kModule) << "\t" << amp << "\n";
    }

    // Calculate unique amplitude intensities
    for (const auto &pair : m_constrained_amps)
    {
        const std::string &unique_amp = pair.first;
        const std::vector<std::string> &full_amps = pair.second;

        std::pair<double, double> intensity = m_fit_results.intensity(full_amps, m_is_acceptance_corrected);
        m_unique_amp_intensities[unique_amp] = intensity;

        report(DEBUG, kModule) << "Unique amplitude intensity for "
                               << unique_amp
                               << " = "
                               << intensity.first
                               << " +/- "
                               << intensity.second
                               << "\n";
    }

    // Extract production coefficients for each unique amplitude name
    for (const auto &unique_amp : uniqueAmpNames())
    {
        // any full amplitude name will do, they are all constrained so their
        // production coefficients are identical
        std::string full_amp_name = m_constrained_amps[unique_amp][0];
        m_production_coefficients[unique_amp] = m_fit_results.scaledProductionParameter(full_amp_name);
    }

    // Calculate phase differences between all pairs of unique amplitude names that
    // share the same "reaction::sum"
    for (const auto &pair : m_constrained_amps)
    {
        const std::string &amp1_name = pair.first;
        const std::vector<std::string> &amp1_full_names = pair.second;

        for (const auto &amp2_pair : m_constrained_amps)
        {
            const std::string &amp2_name = amp2_pair.first;
            if (amp1_name >= amp2_name)
                continue; // avoid self-comparison and symmetric duplicates
            const std::vector<std::string> &amp2_full_names = amp2_pair.second;

            // this will ensure we only calculate one phase difference per unique amplitude pair
            bool found_phase_diff = false;

            // search for shared reaction::sum between amp1 and amp2
            for (const auto &full_amp1 : amp1_full_names)
            {
                for (const auto &full_amp2 : amp2_full_names)
                {
                    // check if they share reaction::sum
                    size_t pos1 = full_amp1.find("::", 0);
                    size_t pos2 = full_amp1.find("::", pos1 + 2);
                    std::string reaction_sum1 = full_amp1.substr(0, pos2);

                    pos1 = full_amp2.find("::", 0);
                    pos2 = full_amp2.find("::", pos1 + 2);
                    std::string reaction_sum2 = full_amp2.substr(0, pos2);

                    if (reaction_sum1 == reaction_sum2)
                    {
                        // they share a coherent sum, so calculate phase difference
                        std::pair<double, double> phase_diff = m_fit_results.phaseDiff(full_amp1, full_amp2);
                        m_phase_differences[{full_amp1, full_amp2}] = phase_diff;

                        report(DEBUG, kModule) << "Phase difference between "
                                               << full_amp1
                                               << " and "
                                               << full_amp2
                                               << " = "
                                               << phase_diff.first
                                               << " +/- "
                                               << phase_diff.second
                                               << "\n";

                        found_phase_diff = true;
                        break; // break inner amp2 loop
                    }
                }
                if (found_phase_diff)
                    break; // break outer amp1 loop
            }
        }
    }
}

std::set<std::string> FitConverter::uniqueAmpNames() const
{
    std::set<std::string> unique_names;
    for (const auto &amplitude : m_fit_results.ampList())
    {
        size_t last_colon = amplitude.rfind("::");
        if (last_colon != std::string::npos)
        {
            std::string amp_name = amplitude.substr(last_colon + 2);
            unique_names.insert(amp_name);
        }
    }
    return unique_names;
}

std::set<std::string> FitConverter::uniqueSumAmpNames() const
{
    std::set<std::string> unique_names;
    for (const auto &amplitude : m_fit_results.ampList())
    {
        size_t first_colon = amplitude.find("::");
        size_t last_colon = amplitude.rfind("::");
        if (first_colon != std::string::npos && last_colon != std::string::npos && first_colon != last_colon)
        {
            std::string sum_amp_name = amplitude.substr(first_colon + 2, last_colon - first_colon - 2);
            unique_names.insert(sum_amp_name);
        }
    }
    return unique_names;
}

std::string FitConverter::getCSVHeader() const
{
    std::string header = "file,";
    // standard results
    for (const auto &pair : m_standard_results)
    {
        header += pair.first + ",";
    }
    // parameters
    for (const auto &pair : m_parameters)
    {
        header += pair.first + "," + pair.first + "_err,";
    }
    // production coefficients
    for (const auto &pair : m_production_coefficients)
    {
        header += pair.first + "_re," + pair.first + "_im,";
    }
    // unique amplitude intensities
    for (const auto &pair : m_unique_amp_intensities)
    {
        header += pair.first + "," + pair.first + "_err,";
    }
    // phase differences
    for (const auto &pair : m_phase_differences)
    {
        std::string full_amp1 = pair.first.first;
        std::string full_amp2 = pair.first.second;

        // extract unique amplitude names from full amplitude strings
        size_t last_colon1 = full_amp1.rfind("::");
        size_t last_colon2 = full_amp2.rfind("::");
        std::string amp1 = (last_colon1 != std::string::npos) ? full_amp1.substr(last_colon1 + 2) : full_amp1;
        std::string amp2 = (last_colon2 != std::string::npos) ? full_amp2.substr(last_colon2 + 2) : full_amp2;

        header += amp1 + "_" + amp2 + "," + amp1 + "_" + amp2 + "_err,";
    }
    // remove trailing comma
    if (!header.empty() && header.back() == ',')
    {
        header.pop_back();
    }
    return header;
}

std::string FitConverter::getCSVRow() const
{
    std::string row = m_fit_file + ",";
    // standard results
    for (const auto &pair : m_standard_results)
    {
        row += std::to_string(pair.second) + ",";
    }
    // parameters
    for (const auto &pair : m_parameters)
    {
        row += std::to_string(pair.second.first) + "," + std::to_string(pair.second.second) + ",";
    }
    // production coefficients
    for (const auto &pair : m_production_coefficients)
    {
        row += std::to_string(pair.second.real()) + "," + std::to_string(pair.second.imag()) + ",";
    }
    // unique amplitude intensities
    for (const auto &pair : m_unique_amp_intensities)
    {
        row += std::to_string(pair.second.first) + "," + std::to_string(pair.second.second) + ",";
    }
    // phase differences
    for (const auto &pair : m_phase_differences)
    {
        row += std::to_string(pair.second.first) + "," + std::to_string(pair.second.second) + ",";
    }
    // remove trailing comma
    if (!row.empty() && row.back() == ',')
    {
        row.pop_back();
    }
    return row;
}

std::string FitConverter::getCSVCovarianceMatrixHeader() const
{
    std::string header = "file,parameter";
    for (const auto &par : m_fit_results.parNameList())
    {
        header += "," + par;
    }
    return header;
}

std::string FitConverter::getCSVCorrelationMatrixHeader() const
{
    return getCSVCovarianceMatrixHeader(); // same header format
}

std::string FitConverter::getCSVCovarianceMatrix() const
{
    std::string csv;
    for (long unsigned int row = 0; row < m_error_matrix.size(); row++)
    {
        const std::string row_par = m_fit_results.parNameList()[row];
        csv += m_fit_file + "," + row_par;

        for (long unsigned int col = 0; col < m_error_matrix[row].size(); col++)
        {
            const std::string col_par = m_fit_results.parNameList()[col];

            // double check that our indexing is correct
            if (m_error_matrix[row][col] != m_fit_results.covariance(row_par, col_par))
            {
                report(ERROR, kModule) << "Mismatch in covariance values between parameters "
                                       << row_par << " and " << col_par << ".\n";
                assert(false);
            }

            csv += "," + std::to_string(m_error_matrix[row][col]);
        }

        // don't add newline to last row, to avoid extra blank line at end of file
        // and to match the style of the other csv functions
        if (row != m_error_matrix.size() - 1)
            csv += "\n";
    }
    return csv;
}

std::string FitConverter::getCSVCorrelationMatrix() const
{
    std::string csv;
    for (long unsigned int row = 0; row < m_error_matrix.size(); row++)
    {
        const std::string row_par = m_fit_results.parNameList()[row];
        const double row_par_error = m_fit_results.parError(row_par);
        csv += m_fit_file + "," + row_par;

        for (long unsigned int col = 0; col < m_error_matrix[row].size(); col++)
        {
            const std::string col_par = m_fit_results.parNameList()[col];
            const double col_par_error = m_fit_results.parError(col_par);

            // double check that our indexing is correct
            if (m_error_matrix[row][col] != m_fit_results.covariance(row_par, col_par))
            {
                report(ERROR, kModule) << "Mismatch in covariance values between parameters "
                                       << row_par << " and " << col_par << ".\n";
                assert(false);
            }

            // avoid division by zero
            if (row_par_error == 0.0 || col_par_error == 0.0)
            {
                csv += ",0.0"; // set correlation to 0 if error is zero
                continue;
            }

            double correlation = m_error_matrix[row][col] / (row_par_error * col_par_error);
            csv += "," + std::to_string(correlation);
        }

        // don't add newline to last row, to avoid extra blank line at end of file
        // and to match the style of the other csv functions
        if (row != m_error_matrix.size() - 1)
            csv += "\n";
    }
    return csv;
}

std::map<std::string, std::vector<std::string>> FitConverter::findConstrainedAmps() const
{
    std::map<std::string, std::vector<std::string>> constrained_amp_map;

    if (ampNamesAreConstrained())
    {
        for (const auto &unique_amp : uniqueAmpNames())
        {
            // the bool assures us that all terms with this amp name are constrained to
            // each other
            std::vector<TermInfo *> shared_terms = m_cfg_info->termList("", "", unique_amp);

            // store the full names of the constrained amplitudes
            for (const auto &term : shared_terms)
                constrained_amp_map[unique_amp].push_back(term->fullName());
        }
        return constrained_amp_map;
    }

    if (sumAmpNamesAreConstrained())
    {
        for (const auto &unique_sum_amp : uniqueSumAmpNames())
        {
            std::string sum = unique_sum_amp.substr(0, unique_sum_amp.find("::"));
            std::string amp = unique_sum_amp.substr(unique_sum_amp.find("::") + 2);
            // the bool assures us that all terms with this sum::amp name are
            // constrained to each other
            std::vector<TermInfo *> shared_terms = m_cfg_info->termList("", sum, amp);

            // store the full names of the constrained amplitudes
            for (const auto &term : shared_terms)
                constrained_amp_map[unique_sum_amp].push_back(term->fullName());
        }
        return constrained_amp_map;
    }

    // if we reach here, neither assumption held, so just map full amplitude names to 
    // themselves
    for (const auto &amplitude : m_fit_results.ampList())
    {
        constrained_amp_map[amplitude] = {amplitude};
    }
    return constrained_amp_map;
}


bool FitConverter::ampNamesAreConstrained() const
{
    for (const auto &unique_amp : uniqueAmpNames())
    {
        std::vector<TermInfo *> shared_terms = m_cfg_info->termList("", "", unique_amp);

        // we assume terms with common amplitude names are all constrained to each
        // other, so use first term to check this
        std::vector<TermInfo *> constrained_terms = shared_terms[0]->constraints();

        // filter constrained terms to only those with the same amplitude name
        for (auto it = constrained_terms.begin(); it != constrained_terms.end();)
        {
            if ((*it)->fullName().find(unique_amp) == std::string::npos)
                it = constrained_terms.erase(it);
            else
                ++it;
        }

        // The constrained terms does not include the self term, so add it back in for 
        // comparison
        constrained_terms.push_back(shared_terms[0]);

        // Compare the two sets of terms to ensure they are the same
        std::sort(shared_terms.begin(), shared_terms.end());
        std::sort(constrained_terms.begin(), constrained_terms.end());
        if (shared_terms != constrained_terms)
        {
            return false;
        }
    }
    return true;
}

bool FitConverter::sumAmpNamesAreConstrained() const
{
    for (const auto &unique_sum_amp : uniqueSumAmpNames())
    {
        std::string sum = unique_sum_amp.substr(0, unique_sum_amp.find("::"));
        std::string amp = unique_sum_amp.substr(unique_sum_amp.find("::") + 2);

        std::vector<TermInfo *> shared_terms = m_cfg_info->termList("", sum, amp);

        // we assume terms with common sum::amp names are all constrained to each
        // other, so use first term to check this
        std::vector<TermInfo *> constrained_terms = shared_terms[0]->constraints();

        // filter constrained terms to only those with the same sum::amp name
        for (auto it = constrained_terms.begin(); it != constrained_terms.end();)
        {
            if ((*it)->fullName().find(unique_sum_amp) == std::string::npos)
                it = constrained_terms.erase(it);
            else
                ++it;
        }

        // The constrained terms does not include the self term, so add it back in for 
        // comparison
        constrained_terms.push_back(shared_terms[0]);

        // Compare the two sets of terms to ensure they are the same
        std::sort(shared_terms.begin(), shared_terms.end());
        std::sort(constrained_terms.begin(), constrained_terms.end());
        if (shared_terms != constrained_terms)
        {
            return false;
        }
    }
    return true;
}