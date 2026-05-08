/**
 * @file AmplitudeParser.h
 * @author Kevin Scheuer
 * @brief Helpers for parsing amplitude labels and building sum groups.
 * @date 2026-05-06
 *
 * @note All amplitudes that match the criteria for a sum group will be included,
 * meaning amplitudes will be gathered across reactions and sums as long as they share
 * the relevant quantum numbers based on the selected naming scheme. For example, if we
 * have the following amplitudes:
 * - "reaction1::sum1::p1ppS"
 * - "reaction1::sum2::p1ppS"
 * - "reaction2::sum1::p1ppS"
 * Then all three of these amplitudes would be included in the same sum group "JPL"
 * because they all share the same J, P, and L values.
 *
 * @page How to Build Your A Custom Naming Scheme
 * The AmplitudeParser is designed to be flexible and support different amplitude naming
 * schemes. If your amplitudes follow a different naming convention than the ones
 * currently supported, you can easily extend the AmplitudeParser to accommodate your
 * scheme.
 *
 * Add the label for your naming scheme in the NamingScheme enum, and then implement
 * the parsing logic in the inferNamingScheme() and parseAmplitudeName() methods. Build
 * your desired sum groups in sumsForScheme(), and make sure your label will be
 * properly made in buildSumLabel(). The parser will then be able to extract the
 * relevant quantum numbers from your amplitude names and group them into sums according
 * to your definitions.
 *
 * It's very important that your scheme has unique identifiers for the sums you want to
 * build. For example, if you want to build sums over all amplitudes with the same J,
 * and another sum just over amplitudes with the same m, you need to ensure the J and
 * m values can be uniquely extracted from the amplitude name. If both use digits,
 * then the parser will not know which digits correspond to J and which correspond to m.
 * In that case, you would need to modify the amplitude parser to look for specific
 * delimiters or patterns in the amplitude name to extract the quantum numbers
 * correctly.
 *
 */

#ifndef FIT_AMPLITUDE_PARSER_H
#define FIT_AMPLITUDE_PARSER_H

#include <map>
#include <string>
#include <vector>

class AmplitudeParser
{
public:
    enum class NamingScheme
    {
        Auto,
        JLme,
        eJPmL,
        Lme
    };

    /**
     * @brief Construct a new Amplitude Parser object
     *
     * @param naming_scheme what amplitude naming scheme to use for parsing amplitude names and building sum groups.
     * Supported schemes include:
     * - "auto": the parser will attempt to infer the naming scheme from the amplitude labels
     * - "JLme": amp names like "1P-1+" (the one recommended for most amplitude analyses)
     * - "eJPmL": amp names like "p1ppS" (a common vector-pseudoscalar naming scheme)
     * - "Lme": amp names like "P1+" (a common two-pseudoscalar naming scheme)
     */
    explicit AmplitudeParser(const std::string &naming_scheme = "auto");

    /**
     * @brief
     *
     * @param full_amplitudes
     * @return std::map<std::string, std::vector<std::string>>
     */
    std::map<std::string, std::vector<std::string>>
    buildSumGroups(const std::vector<std::string> &full_amplitudes) const;

private:
    /**
     * @brief struct to contain the parsed quantum numbers from an amplitude name
     *
     * For example:
     * "reaction::sum::p1ppS" would be parsed to:
     * amp_name = "p1ppS"
     * e = "p"
     * J = "1"
     * P = "p"
     * m = "p"
     * L = "S"
     *
     * "reaction::sum::2P1+" would be parsed to:
     * amp_name = "2P1+"
     * e = "+"
     * J = "2"
     * L = "P"
     * m = "1"
     * P is not included in this naming scheme, so it would be empty.
     */
    struct ParsedAmplitude
    {
        std::string amp_name; ///< The "ampName" in the full "reaction::sum::ampName"
        std::string e;        ///< Reflectivity quantum number
        std::string J;        ///< Total angular momentum
        std::string P;        ///< Parity (sometimes excluded as its unnecessary)
        std::string m;        ///< Spin projection
        std::string L;        ///< Orbital angular momentum
    };

    NamingScheme m_requested_scheme; ///< The naming scheme to use, or "Auto" to infer from amplitude labels

    /**
     * @brief Set the naming scheme based on the input string from the user.
     *
     * If the input is "auto" or is an unrecognized string, we will attempt to infer the
     * naming scheme from the amplitude labels later.
     *
     * @param scheme_name user string indicating the naming scheme
     * @return NamingScheme objective naming scheme enum value
     */
    static NamingScheme namingSchemeFromString(const std::string &scheme_name);

    /**
     * @brief Infer the naming scheme from the amplitude name format if "auto" is selected.
     *
     * @param amp_name the part of the full amplitude string after the last "::"
     *
     * @throw ERROR if the name is empty and we cannot infer a scheme
     *
     * @return NamingScheme::JLme if the name starts with a digit
     * NamingScheme::Lme if the name starts with an uppercase letter and ends with + or -
     * NamingScheme::eJPmL if the name starts with "p" or "m" and ends with an uppercase letter
     * NamingScheme::Auto if we cannot infer a scheme
     */
    static NamingScheme inferNamingScheme(const std::string &amp_name);

    /**
     * @brief Filter amplitudes to belong to the chosen naming scheme
     *
     * To avoid undefined behavior when parsing amplitude names, we need to ensure that
     * all amplitudes we attempt to group follow the same naming scheme. Any amps not
     * belonging to the scheme will be reported in a printed warning, and excluded from
     * the returned list.
     *
     * @param full_amplitudes all amplitudes from a fit result, in their full "reaction::sum::ampName" form
     * @param scheme naming scheme to filter by
     * @return std::vector<std::string> full amplitudes that match the naming scheme
     */
    std::vector<std::string> filterAmplitudesByScheme(const std::vector<std::string> &full_amplitudes, NamingScheme scheme) const;

    /**
     * @brief Extract the amplitude name from the full amplitude string.
     *
     * For example, for "reaction::sum::amp", this would return "amp". This is the part
     * of the amplitude label that typically contains the quantum numbers and is used
     * for grouping into sums.
     *
     * @param full_name "reaction::sum::ampName" style amplitude string
     * @return std::string ampName, the part of the full amplitude string after the last "::"
     */
    static std::string ampName(const std::string &full_name);

    /**
     * @brief Extract quantum numbers from the amplitude name according to the scheme
     *
     * @param amp_name ampName in "reaction::sum::ampName"
     * @param scheme naming scheme to use for parsing
     * @return ParsedAmplitude object that contains the original amp_name and the parsed quantum numbers
     */
    ParsedAmplitude parseAmplitudeName(const std::string &amp_name, NamingScheme scheme) const;

    /**
     * @brief Define the grouping of quantum numbers that will become a sum label
     */
    struct SumGroupDef
    {
        std::vector<char> quantum_numbers; ///< Quantum number fields to concatenate for the group label
    };

    /**
     * @brief Get the sum definitions for a given naming scheme
     *
     * A sum label is denoted by an absence of a quantum number in the amplitude
     * naming scheme. For example, the JLme scheme will have a "JL" sum group that
     * includes all amplitudes with the same J and L, but any m and e. Since the
     * ordering of the quantum numbers is dependent on the scheme, the group labels are
     * thus also scheme-dependent.
     *
     * @param scheme naming scheme to get group definitions for
     * @return std::vector<SumGroupDef>
     */
    std::vector<SumGroupDef> sumsForScheme(NamingScheme scheme) const;

    /**
     * @brief concatenates the quantum numbers for a SumGroupDef into a string label
     *
     * @param sum_group a grouping of quantum numbers to concatenate into a sum label, e.g. J and L
     * @param parsed a parsed amplitude name
     * @return std::string the label to be used as the key in sum map e.g. "1P" for a JL group with J=1 and L=P
     */
    std::string buildSumLabel(const SumGroupDef &sum_group, const ParsedAmplitude &parsed) const;

    static const char *kModule;
};

#endif // FIT_AMPLITUDE_PARSER_H