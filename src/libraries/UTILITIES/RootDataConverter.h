/**
 * @file RootDataConverter.h
 * @author Kevin Scheuer
 * @brief Class for converting ROOT tree data associated with an AmpTools fit to any format
 * @date 2026-02-02
 *
 */

#ifndef DATA_CONVERTER_H
#define DATA_CONVERTER_H

#include <map>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/ConfigurationInfo.h"

/**
 * @class RootDataConverter
 * @brief Converts AmpTools fit result data to any format
 *
 * This class extracts data information from a single AmpTools .fit file. It can then
 * output the data in various formats, such as CSV. The extracted data includes:
 * - Bin edges, centers, averages, and RMS values for t, beam energy, and mass
 *   histograms
 * - Total number of events and errors for data and acceptance-corrected data
 * - Efficiency calculation based on genMC and accMC event counts
 *
 * Example usage in convert_to_csv.cc:
 * @code
 * TODO: complete example
 * @endcode
 */
class RootDataConverter
{
public:
    /**
     * @brief Construct a RootDataConverter for a single fit file
     *
     * This assumes the 4-momenta indices are as follows:
     * - index 0 (B): beam photon
     * - index 1: single lower vertex particle, typically the recoil proton
     * - indices 2 to N: final state particles from upper vertex
     *
     * @param[in] filename path to the .fit file
     * @param[in] mute_warnings whether to mute warnings from FitResults and when
     *  background files are not found
     */
    RootDataConverter(const std::string &filename, bool mute_warnings = false);

    /**
     * @brief Construct a new Data Converter object with specified lower vertex indices
     *
     * The 4-momenta indices are as follows:
     * - index 0 (B): beam photon
     * - indices 1 to M: lower vertex particles, e.g. proton, pi+ that decay from Delta++
     * - indices M+1 to N: final state particles from upper vertex
     * 
     * TODO: implement
     *
     * @param[in] filename path to the .fit file
     * @param[in] lower_vertex_indices 4-momenta indices for lower vertex particles
     * @param[in] mute_warnings whether to mute warnings from FitResults and when
     *  background files are not found
     */
    RootDataConverter(const std::string &filename,
                      const std::vector<int> &lower_vertex_indices,
                      bool mute_warnings = false);

    /**
     * @brief Construct a new Data Converter object with specified lower vertex and isobar indices
     *
     * This will add new headers and values for isobar mass and efficiency calculations.
     *
     * The 4-momenta indices are as follows:
     * - index 0 (B): beam photon
     * - indices 1 to M: lower vertex particles, e.g. proton, pi+ that decay from Delta++
     * - indices M+1 to P: isobar particles, e.g. pi+ pi- from rho
     * - indices P+1 to N: remaining final state particles from upper vertex
     *
     * Note that indices M+1 to N are still used to calculate the total mass values.
     *
     * TODO: implement. make sure no conflicts with lower vertex indices. Also this
     * should add isobar mass and event counts to the output.
     *
     * @param[in] filename path to the .fit file
     * @param[in] lower_vertex_indices 4-momenta indices for lower vertex particles
     * @param[in] isobar_indices 4-momenta indices for isobar particles
     * @param[in] mute_warnings whether to mute warnings from FitResults and when
     *  background files are not found
     */
    RootDataConverter(const std::string &filename,
                      const std::vector<int> &lower_vertex_indices,
                      const std::vector<int> &isobar_indices,
                      bool mute_warnings = false);

    /**
     * @brief Extract fit data from the .fit file
     */
    void extract();

    /**
     * @brief Data files associated with the fit results
     * @note It assumes one data file per reaction, and that it is the first argument
     * passed to the data reader.
     * @return std::vector<std::string> where each element is an absolute path to a data file
     */
    std::vector<std::string> dataFiles() const { return findFiles("data"); };

    /**
     * @brief Background files associated with the fit results, if they exist
     * @note It assumes one background file per reaction, and that it is the first
     * argument passed to the data reader.
     * @return std::vector<std::string> where each element is an absolute path to a background file
     */
    std::vector<std::string> backgroundFiles() const { return findFiles("background"); };

    /**
     * @brief Generated MC files associated with the fit results
     * @note It assumes one genMC file per reaction, and that it is the first argument
     * passed to the data reader.
     * @return std::vector<std::string> where each element is an absolute path to a genMC file
     */
    std::vector<std::string> genMCFiles() const { return findFiles("genMC"); };

    /**
     * @brief Accepted MC files associated with the fit results
     * @note It assumes one accMC file per reaction, and that it is the first argument
     * passed to the data reader.
     * @return std::vector<std::string> where each element is an absolute path to an accMC file
     */
    std::vector<std::string> accMCFiles() const { return findFiles("accMC"); };

    // TODO: write these docstrings
    void validateDataFiles() const { validateFiles(m_data_files, "data"); }
    void validateBackgroundFiles() const { validateFiles(m_background_files, "background"); }
    void validateGenMCFiles() const { validateFiles(m_genMC_files, "genMC"); }
    void validateAccMCFiles() const { validateFiles(m_accMC_files, "accMC"); }


    /**
     * @brief data tree name associated with the fit results
     * If one tree exists in the data files, that tree name is returned. Otherwise,
     * the tree name is extracted based on the arguments given to known data readers.
     * If the data reader is unknown, it will warn and default to "kin".
     * @return std::string tree name
     */
    std::string dataTreeName() const { return findTreeName("data"); };

    /**
     * @brief background tree name associated with the fit results
     * If one tree exists in the data files, that tree name is returned. Otherwise,
     * the tree name is extracted based on the arguments given to known data readers.
     * If the data reader is unknown, it will warn and default to "kin".
     * @return std::string tree name
     */
    std::string backgroundTreeName() const { return findTreeName("background"); };

    /**
     * @brief generated MC tree name associated with the fit results
     * If one tree exists in the data files, that tree name is returned. Otherwise,
     * the tree name is extracted based on the arguments given to known data readers.
     * If the data reader is unknown, it will warn and default to "kin".
     * @return std::string tree name
     */
    std::string genMCTreeName() const { return findTreeName("genMC"); };

    /**
     * @brief accepted MC tree name associated with the fit results
     * If one tree exists in the data files, that tree name is returned. Otherwise,
     * the tree name is extracted based on the arguments given to known data readers.
     * If the data reader is unknown, it will warn and default to "kin".
     * @return std::string tree name
     */
    std::string accMCTreeName() const { return findTreeName("accMC"); };

    /**
     * @brief Get the name of the weight branch in the specified tree and file type
     *
     * Some data readers, such as FSRootDataReader, can store the event weights in a
     * friend tree. This function will check for that case first, and return the
     * appropriate weight branch name. If no friend tree is used, it will check for
     * "weight" or "Weight" branches in the main tree. If those branches do not exist,
     * it will return an empty string to indicate no weight branch found.
     *
     * @param[in] file_type The type of file ("data", "background")
     * @param[in] tree_name The name of the tree within the file
     * @return std::string The name of the weight branch, or empty string if none found
     *
     * @todo FSRoot friend tree case untested
     */
    std::string weightBranchName(
        const std::string &file_type,
        const std::string &tree_name) const;

    /**
     * @brief Set indices that belong to the lower vertex final state particles
     *
     * An AmpTools ROOT tree contains 4-vectors for all particles in the event, labelled
     * "EnX", "PxX", "PyX", "PzX" where X is the particle index. The "B" index is always
     * the beam photon, but the other indices are not guaranteed to be in any particular
     * order. The simplest case is a recoil proton at the lower vertex, which would just
     *  be index 1. However, there can be multiple particles associated with the lower
     * vertex like a Delta++ -> p + pi+ decay. 
     * 
     * This function allows the user to specify which indices correspond to the lower 
     * vertex system. The remaining indices are assumed to be final state particles from
     * the upper vertex. This is important for correctly calculating the 4-momentum
     * transfer t and the upper vertex mass.
     */
    void setLowerVertexIndices(const std::vector<int> &indices) { m_lower_vertex_indices = indices; }

    /**
     * @brief Get indices for the lower vertex particles
     *
     *  See setLowerVertexIndices for more details.
     */
    std::vector<int> getLowerVertexIndices() const { return m_lower_vertex_indices; }

    /**
     * @brief Set indices that belong to the upper vertex final state particles
     * 
     * Any particle indices not in m_lower_vertex_indices that are also not 0 are thus
     * assumed to belong to the upper vertex. This saves all leftover indices to
     * m_upper_vertex_indices. This is important for correctly calculating the upper
     * vertex mass.
     *
     * @note This function must be called AFTER setLowerVertexIndices
     */
    void setUpperVertexIndices();

    /**
     * @brief Get indices for the upper vertex particles
     * 
     * See setUpperVertexIndices for more details.
     */
    std::vector<int> getUpperVertexIndices() const { return m_upper_vertex_indices; }


    /**
     * @brief Map of header names to their values
     */
    std::map<std::string, double> getValues() const { return m_values; }

    bool isValid() const { return m_fit_results.valid(); }

    ~RootDataConverter() = default;

private:
    const std::string m_fit_file;            ///< Path to the .fit file
    const FitResults m_fit_results;          ///< AmpTools FitResults object
    const ConfigurationInfo *m_cfg_info;     ///< ConfigurationInfo from the fit
    std::vector<int> m_lower_vertex_indices; ///< Indices of lower vertex particles in the data 4-vectors
    std::vector<int> m_upper_vertex_indices; ///< Indices of upper vertex particles in the data 4-vectors
    std::vector<int> m_isobar_indices;       ///< Indices of isobar particles in the data 4-vectors
    bool m_mute_warnings;                    ///< Whether to mute warnings about missing background files

    const std::vector<std::string> m_data_files;
    const std::vector<std::string> m_background_files;
    const std::vector<std::string> m_genMC_files;
    const std::vector<std::string> m_accMC_files;

    bool m_background_files_exist;

    const std::string m_data_tree_name;
    const std::string m_background_tree_name;
    const std::string m_genMC_tree_name;
    const std::string m_accMC_tree_name;

    // TODO: variables we'll want to implement 
    // "file",
    // "t_low",
    // "t_high",
    // "t_center",
    // "t_avg",
    // "t_rms",
    // "e_low", DONE
    // "e_high", DONE
    // "e_center", DONE
    // "e_avg", DONE
    // "e_rms", DONE
    // "m_low",
    // "m_high",
    // "m_center",
    // "m_avg",
    // "m_rms",
    // "events",
    // "events_err",
    // "ac_events",
    // "ac_events_err",
    // "efficiency",
    

    std::map<std::string, double> m_values; //< Map of headers to their values for the CSV output (except for the "file" header)

    /**
     * @brief Find files of a given type associated with the fit results
     * @param[in] file_type data, background, genMC, or accMC
     * @note It assumes one file per reaction, and that it is the first argument passed
     * to the data reader.
     * @return std::vector<std::string> where each element is an absolute path to a file of the given type
     */
    std::vector<std::string> findFiles(const std::string &file_type) const;


    /**
     * @brief Confirm that files vector is non-empty, unless its background
     * 
     * Generally will throw an error if no files are present for a given type. If 
     * file_type is "background", this can happen, but it warns the user that it assumes
     * any event weights are thus stored directly in the data trees.
     * 
     * @param files absolute paths of files associated with the fit results for a given type
     * @param file_type data, background, genMC, or accMC
     */
    void validateFiles(const std::vector<std::string> &files,
                       const std::string &file_type) const;

    /**
     * @brief Find the tree name associated with the fit results
     *
     * This function will first attempt to see if only one tree exists in the ROOT file.
     * If multiple trees exist, it will extract the tree name based on the arguments
     * given to known data readers in halld_sim. If the data reader is unknown, it will
     * warn and default to "kin".
     *
     * @param[in] file_type data, background, genMC, or accMC
     *
     * @return std::string tree name
     */
    std::string findTreeName(const std::string &file_type) const;

    /**
     * @brief Get the names of all trees in a ROOT file
     *
     * @param filename Path to the ROOT file
     * @return std::vector<std::string> List of tree names in the file
     */
    std::vector<std::string> treesInFile(const std::string &filename) const;

    /**
     * @brief Extract beam energy statistics from the trees
     *
     * Calculates min, max, mean, and RMS values for the "EnPB" branch from the data
     * file. Stores results in m_values with keys: "e_low", "e_high", "e_center",
     * "e_avg", "e_rms"
     *
     * @param[in] weight_branch_name Name of the weight branch (if empty, weights assumed to be 1.0)
     */
    void extractBeamEnergyStats(const std::string &weight_branch_name);

    /**
     * @brief Extract -t 4-momentum transfer statistics from the trees
     * 
     * Calculates min, max, mean, and RMS values from the momentum transfer between 
     * the at-rest proton and the lower-vertex recoil system. This is calculated as 
     * t = (P_proton - P_recoil)^2, where P_recoil is the sum of the 4-vectors of the 
     * lower vertex particles (see setLowerVertexIndices)
     * Stores results in m_values with keys: "t_low", "t_high", "t_center", "t_avg", and
     * "t_rms"
     * 
     * @param weight_branch_name Name of the weight branch (if empty, weights assumed to be 1.0)
     * 
     * @todo The function currently has not been tested for the FSRootFriendTree scenario,
     * and so warns the user and assigns a weight of 1.0 for this case.
     */
    void extractFourMomentumTransferStats(const std::string &weight_branch_name);

    /**
     * @brief Extract mass statistics for the upper vertex system from the trees
     * 
     * Calculates min, max, mean, and RMS values from the total upper-vertex mass
     * system, which is calculated as the invariant mass of the sum of the 4-vectors of
     * the upper vertex particles (see setUpperVertexIndices). Stores the results in 
     * m_values with keys: "m_low", "m_high", "m_center", "m_avg", and "m_rms"
     * 
     * @param weight_branch_name Name of the weight branch (if empty, weights assumed to be 1.0)
     * 
     * @todo The function currently has not been tested for the FSRootFriendTree scenario,
     * and so warns the user and assigns a weight of 1.0 for this case.
     */     
    void extractUpperVertexMassStats(const std::string &weight_branch_name);

    /**
     * @brief Get the max/min values of a branch for a set of files
     * 
     * @param files set of ROOT files containing a common tree and branch name
     * @param tree_name ROOT tree containing the branch of interest
     * @param branch_name branch name we want to find the max/min values of
     * @return std::pair<double, double> minimum and maximum values of the branch across all files
     */
    std::pair<double, double> findMinMaxOfBranch(const std::vector<std::string> &files,
                                               const std::string &tree_name,
                                               const std::string &branch_name) const;

    static const char *kModule;
};

#endif // DATA_CONVERTER_H