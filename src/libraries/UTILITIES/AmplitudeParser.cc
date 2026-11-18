/**
 * @file AmplitudeParser.cc
 * @author Kevin Scheuer
 * @brief Implementation of the amplitude naming parser and sum grouping logic.
 * @date 2026-05-08
 */

#include <algorithm>
#include <cctype>
#include <cassert>

#include "AmplitudeParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/report.h"

const char *AmplitudeParser::kModule = "AmplitudeParser";

AmplitudeParser::AmplitudeParser(const std::string &naming_scheme) : m_requested_scheme(namingSchemeFromString(naming_scheme)) {}

AmplitudeParser::NamingScheme
AmplitudeParser::namingSchemeFromString(const std::string &scheme_name)
{
    if (scheme_name == "JLme")
        return NamingScheme::JLme;
    if (scheme_name == "eJPmL")
        return NamingScheme::eJPmL;
    if (scheme_name == "Lme")
        return NamingScheme::Lme;
    return NamingScheme::Auto;
}

AmplitudeParser::NamingScheme
AmplitudeParser::inferNamingScheme(const std::string &amp_name)
{
    if (amp_name.empty())
    {
        report(ERROR, kModule) << "Amplitude name is empty, cannot infer naming scheme.\n";
        assert(false);
    }

    if (std::isdigit(static_cast<unsigned char>(amp_name.front())))
        return NamingScheme::JLme;

    if (std::isupper(static_cast<unsigned char>(amp_name.front())) &&
        (amp_name.back() == '+' || amp_name.back() == '-'))
        return NamingScheme::Lme;

    if ((amp_name.front() == 'p' || amp_name.front() == 'm') &&
        std::isupper(static_cast<unsigned char>(amp_name.back())))
        return NamingScheme::eJPmL;

    return NamingScheme::Auto; // default if we can't infer
}

std::vector<std::string>
AmplitudeParser::filterAmplitudesByScheme(
    const std::vector<std::string> &full_amplitudes,
    NamingScheme scheme) const
{
    std::vector<std::string> filtered;
    for (const auto &full_amp : full_amplitudes)
    {
        const std::string amp_name = ampName(full_amp);
        if (inferNamingScheme(amp_name) == scheme)
            filtered.push_back(full_amp);
        else if (!m_suppress_warnings)
        {
            report(WARNING, kModule) << "Amplitude '"
                                     << full_amp
                                     << "' does not match the inferred naming scheme and will be excluded from sum groups.\n";
        }
    }
    return filtered;
}

std::string AmplitudeParser::ampName(const std::string &full_name)
{
    const size_t pos = full_name.rfind("::");
    return (pos == std::string::npos) ? full_name : full_name.substr(pos + 2);
}

AmplitudeParser::ParsedAmplitude
AmplitudeParser::parseAmplitudeName(const std::string &amp_name, NamingScheme scheme) const
{
    ParsedAmplitude parsed;
    parsed.amp_name = amp_name;

    if (scheme == NamingScheme::Auto)
    {
        scheme = inferNamingScheme(amp_name);
    }

    // Parse quantum numbers according to the structure of the naming scheme
    switch (scheme)
    {
    case NamingScheme::JLme:
        parsed.J = amp_name.substr(0, 1);
        parsed.L = amp_name.substr(1, 1);
        parsed.m = amp_name.substr(2, amp_name.size() - 3);
        parsed.e = amp_name.substr(amp_name.size() - 1, 1);
        break;
    case NamingScheme::Lme:
        parsed.L = amp_name.substr(0, 1);
        parsed.m = amp_name.substr(1, amp_name.size() - 2);
        parsed.e = amp_name.substr(amp_name.size() - 1, 1);
        break;
    case NamingScheme::eJPmL:
        parsed.e = amp_name.substr(0, 1);
        parsed.J = amp_name.substr(1, 1);
        parsed.P = amp_name.substr(2, 1);
        parsed.m = amp_name.substr(3, amp_name.size() - 4);
        parsed.L = amp_name.substr(amp_name.size() - 1, 1);
        break;
    default:
        if (!m_suppress_warnings)
            report(WARNING, kModule) << "Unknown naming scheme, cannot parse amplitude name '" << amp_name << "'.\n";
        break;
    }

    return parsed;
}

std::vector<AmplitudeParser::SumGroupDef>
AmplitudeParser::sumsForScheme(NamingScheme scheme) const
{
    if (scheme == NamingScheme::Lme)
    {
        return {
            {{'e'}},
            {{'L'}},
            {{'L', 'e'}},
            {{'L', 'm'}},
        };
    }

    else if (scheme == NamingScheme::JLme)
    {
        return {
            {{'e'}},
            {{'J'}},
            {{'J', 'e'}},
            {{'J', 'L'}},
            {{'J', 'L', 'e'}},
            {{'J', 'L', 'm'}},
        };
    }
    else if (scheme == NamingScheme::eJPmL)
    {
        return {
            {{'e'}},
            {{'J', 'P'}},
            {{'e', 'J', 'P'}},
            {{'J', 'P', 'L'}},
            {{'e', 'J', 'P', 'L'}},
            {{'J', 'P', 'm', 'L'}},
        };
    }
    else
    {
        if(!m_suppress_warnings)
            report(WARNING, kModule) << "Unknown naming scheme, cannot determine sum groups.\n";
        return {};
    }
}

std::string AmplitudeParser::buildSumLabel(
    const SumGroupDef &sum_group,
    const ParsedAmplitude &parsed) const
{
    std::string label;
    for (const auto &qn : sum_group.quantum_numbers)
    {
        std::string value = "";
        switch (qn)
        {
        case 'e':
            value = parsed.e;
            break;
        case 'J':
            value = parsed.J;
            break;
        case 'P':
            value = parsed.P;
            break;
        case 'm':
            value = parsed.m;
            break;
        case 'L':
            value = parsed.L;
            break;
        default:
            if(!m_suppress_warnings)
                report(WARNING, kModule) << "Unknown quantum number '" << qn << "' in sum group definition.\n";
            break;
        }
        label += value;
    }
    return label;
}

std::map<std::string, std::vector<std::string>>
AmplitudeParser::buildSumGroups(const std::vector<std::string> &full_amplitudes) const
{
    std::map<std::string, std::vector<std::string>> groups;

    NamingScheme effective_scheme = m_requested_scheme;

    if (full_amplitudes.empty())
    {
        if(!m_suppress_warnings)
            report(WARNING, kModule) << "No amplitudes found in fit results, so no sums will be built.\n";
        return groups;
    }

    // "auto" means we need to infer the naming scheme from the amplitude labels.
    // We'll use the first amplitude label as a sample to infer the scheme.
    if (effective_scheme == NamingScheme::Auto)
        effective_scheme = inferNamingScheme(ampName(full_amplitudes.front()));

    // If we still can't infer a scheme, warn the user and don't attempt to build sums
    // NOTE: if you have gotten this message and would like to build your own scheme,
    // see the header file for instructions to do so.
    if (effective_scheme == NamingScheme::Auto)
    {
        if(!m_suppress_warnings)
        {
            report(WARNING, kModule) << "Amplitude naming scheme could not be inferred, and so no sums will be built.\n";
            report(WARNING, kModule) << "This likely means your amplitude names don't follow a recognized pattern. Supported patterns include:\n";
            report(WARNING, kModule) << "\t1) JLme: e.g. 1P-1+ \n";
            report(WARNING, kModule) << "\t2) Lme: e.g. P1+ \n";
            report(WARNING, kModule) << "\t3) eJPmL: e.g. p1ppS \n";
            report(WARNING, kModule) << "You can also specify the naming scheme pattern (from the choices above) directly in the FitConverter constructor to avoid inference.\n";
        }
        return groups;
    }

    // Ensure that all amplitudes we attempt to group follow the same naming scheme,
    // otherwise we will get undefined behavior when parsing the amplitude names. If we
    // find any amplitudes that don't match the effective scheme, we will skip them and
    // warn the user.
    std::vector<std::string> amplitudes_in_scheme = filterAmplitudesByScheme(full_amplitudes, effective_scheme);

    // determine what sum groups will be made
    const std::vector<SumGroupDef> sum_groups = sumsForScheme(effective_scheme);

    // loop over all amplitudes across all reactions and sums, parse the amplitude name,
    // and group according to the requested scheme
    for (const auto &full_amp : amplitudes_in_scheme)
    {
        const ParsedAmplitude parsed = parseAmplitudeName(ampName(full_amp), effective_scheme);

        for (const auto &sum_group : sum_groups)
        {
            // create the label that will be used as the CSV header
            const std::string label = buildSumLabel(sum_group, parsed);
            if (label.empty())
                continue;

            // add the full amplitude name that belongs to this sum group in the map
            groups[label].push_back(full_amp);
        }
    }

    // remove duplicates from the sum groups
    for (auto &pair : groups)
    {
        std::sort(pair.second.begin(), pair.second.end());
        pair.second.erase(std::unique(pair.second.begin(), pair.second.end()), pair.second.end());
    }

    return groups;
}


std::set<std::string> 
AmplitudeParser::uniqueAmpNames(const FitResults &result)
{
    std::set<std::string> unique_names;
    for (const auto &amplitude : result.ampList())
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


std::set<std::string> 
AmplitudeParser::uniqueSumAmpNames(const FitResults &result)
{
    std::set<std::string> unique_names;
    for (const auto &amplitude : result.ampList())
    {
        size_t first_colon = amplitude.find("::");
        size_t last_colon = amplitude.rfind("::");
        if (first_colon != std::string::npos && last_colon != std::string::npos && first_colon != last_colon)
        {
            std::string sum_amp_name = amplitude.substr(first_colon + 2);
            unique_names.insert(sum_amp_name);
        }
    }
    return unique_names;
}

std::map<std::string, std::vector<std::string>> 
AmplitudeParser::findConstrainedAmps(const FitResults &result)
{
    std::map<std::string, std::vector<std::string>> constrained_amp_map;

    if (ampNamesAreConstrained(result))
    {
        for (const auto &unique_amp : uniqueAmpNames(result))
        {
            // the bool assures us that all terms with this amp name are constrained to
            // each other
            std::vector<TermInfo *> shared_terms = result.configInfo()->termList("", "", unique_amp);

            // store the full names of the constrained amplitudes
            for (const auto &term : shared_terms)
                constrained_amp_map[unique_amp].push_back(term->fullName());
        }
        return constrained_amp_map;
    }

    if (sumAmpNamesAreConstrained(result))
    {
        for (const auto &unique_sum_amp : uniqueSumAmpNames(result))
        {
            std::string sum = unique_sum_amp.substr(0, unique_sum_amp.find("::"));
            std::string amp = unique_sum_amp.substr(unique_sum_amp.find("::") + 2);
            // the bool assures us that all terms with this sum::amp name are
            // constrained to each other
            std::vector<TermInfo *> shared_terms = result.configInfo()->termList("", sum, amp);

            // store the full names of the constrained amplitudes
            for (const auto &term : shared_terms)
                constrained_amp_map[unique_sum_amp].push_back(term->fullName());
        }
        return constrained_amp_map;
    }

    // if we reach here, neither assumption held, so just map full amplitude names to
    // themselves
    for (const auto &amplitude : result.ampList())
    {
        constrained_amp_map[amplitude] = {amplitude};
    }
    return constrained_amp_map;
}

bool AmplitudeParser::ampNamesAreConstrained(const FitResults &result)
{
    for (const auto &unique_amp : uniqueAmpNames(result))
    {
        std::vector<TermInfo *> shared_terms = result.configInfo()->termList("", "", unique_amp);

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

bool AmplitudeParser::sumAmpNamesAreConstrained(const FitResults &result)
{
    for (const auto &unique_sum_amp : uniqueSumAmpNames(result))
    {
        std::string sum = unique_sum_amp.substr(0, unique_sum_amp.find("::"));
        std::string amp = unique_sum_amp.substr(unique_sum_amp.find("::") + 2);

        std::vector<TermInfo *> shared_terms = result.configInfo()->termList("", sum, amp);

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

std::string AmplitudeParser::getReactionString(const std::string &full_amplitude)
{
    size_t first_colon = full_amplitude.find("::");
    size_t second_colon = full_amplitude.find("::", first_colon + 2);
    if (first_colon == std::string::npos || second_colon == std::string::npos)
    {
        report(ERROR, kModule) << "Amplitude name " << full_amplitude << " does not match expected format 'reaction::sum::amp'.\n";
        assert(false);
    }
    std::string reaction = full_amplitude.substr(0, first_colon);
    return reaction;
}