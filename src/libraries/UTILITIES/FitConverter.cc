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

FitConverter::FitConverter(const std::string &filename, const bool &acceptance_correct) : m_fit_file(filename),
                                                                                          m_fit_results(filename),
                                                                                          m_cfg_info(m_fit_results.configInfo()),
                                                                                          m_is_acceptance_corrected(acceptance_correct)
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
            m_parameters[par_name] = m_fit_results.parValue(par_name);
    }
    for (const auto &pair : m_parameters)
    {
        report(DEBUG, kModule) << "Parameter: " << pair.first << " = " << pair.second << "\n";
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

    // Calculate single amplitude intensities
    for (const auto &pair : m_constrained_amps)
    {
        const std::string &unique_amp = pair.first;
        const std::vector<std::string> &full_amps = pair.second;

        std::pair<double, double> intensity = m_fit_results.intensity(full_amps, m_is_acceptance_corrected);
        m_single_amp_intensities[unique_amp] = intensity;

        report(DEBUG, kModule) << "Single amplitude intensity for "
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
                        break; // break inner loop
                    }
                }
                if (found_phase_diff)
                    break; // break outer loop
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

std::map<std::string, std::vector<std::string>> FitConverter::findConstrainedAmps() const
{
    std::map<std::string, std::vector<std::string>> constrained_amp_map;

    for (const auto &unique_amp : uniqueAmpNames())
    {
        std::vector<TermInfo *> shared_terms = m_cfg_info->termList("", "", unique_amp);

        // we assume terms with common amplitude names are all constrained to each other
        // so use first term to check this
        std::vector<TermInfo *> constrained_terms = shared_terms[0]->constraints();

        // filter constrained terms to only those with the same amplitude name
        for (auto it = constrained_terms.begin(); it != constrained_terms.end();)
        {
            if ((*it)->fullName().find(unique_amp) == std::string::npos)
                it = constrained_terms.erase(it);
            else
                ++it;
        }

        // The constrained terms does not include the self term, so add it back in for comparison
        constrained_terms.push_back(shared_terms[0]);

        // Compare the two sets of terms to ensure they are the same
        std::sort(shared_terms.begin(), shared_terms.end());
        std::sort(constrained_terms.begin(), constrained_terms.end());
        if (shared_terms != constrained_terms)
        {
            report(ERROR, kModule) << "amplitude name "
                                   << unique_amp
                                   << " has inconsistent constraints among terms sharing this name."
                                   << "\n";
            report(ERROR, kModule) << "shared terms:\n";
            for (const auto &term : shared_terms)
                report(ERROR, kModule) << "\t" << term->fullName() << "\n";
            report(ERROR, kModule) << "constrained terms:\n";
            for (const auto &term : constrained_terms)
                report(ERROR, kModule) << "\t" << term->fullName() << "\n";

            assert(false);
        }

        // store the full names of the constrained amplitudes
        for (const auto &term : constrained_terms)
            constrained_amp_map[unique_amp].push_back(term->fullName());
    }

    return constrained_amp_map;
}