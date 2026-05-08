/**
 * @file AmplitudeParser.cc
 * @author Kevin Scheuer
 * @brief Implementation of the amplitude naming parser and sum grouping logic.
 * @date 2026-05-06
 */

#include <algorithm>
#include <cctype>
#include <cassert>

#include "AmplitudeParser.h"
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
        else
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
            {{'e', 'L'}},
            {{'L', 'm'}},
        };
    }

    if (scheme == NamingScheme::JLme)
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
    if (scheme == NamingScheme::eJPmL)
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
        report(WARNING, kModule) << "Amplitude naming scheme could not be inferred, and so no sums will be built.\n";
        report(WARNING, kModule) << "This likely means your amplitude names don't follow a recognized pattern. Supported patterns include:\n";
        report(WARNING, kModule) << "\t1) JLme: e.g. 1P-1+ \n";
        report(WARNING, kModule) << "\t2) Lme: e.g. P1+ \n";
        report(WARNING, kModule) << "\t3) eJPmL: e.g. p1ppS \n";
        report(WARNING, kModule) << "You can also specify the naming scheme pattern (from the choices above) directly in the FitConverter constructor to avoid inference.\n";
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
    for (const auto &full_amp : full_amplitudes)
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