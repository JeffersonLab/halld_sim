/**
 * @file RootDataConverter.cc
 * @author Kevin Scheuer
 * @brief Implementation of RootDataConverter class
 * @date 2026-02-02
 *
 */

#include <cassert>
#include <filesystem>
#include <limits>
#include <tuple>

#include "RootDataConverter.h"
#include "IUAmpTools/report.h"

#include "TFile.h"
#include "TKey.h"
#include "TClass.h"
#include "TCollection.h"
#include "TROOT.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "ROOT/RDataFrame.hxx"

const char *RootDataConverter::kModule = "RootDataConverter";

RootDataConverter::RootDataConverter(const std::string &filename, bool mute_warnings) : m_fit_file(filename),
                                                                                        m_fit_results(filename, mute_warnings),
                                                                                        m_cfg_info(m_fit_results.configInfo()),
                                                                                        m_mute_warnings(mute_warnings),
                                                                                        m_data_files(dataFiles()),
                                                                                        m_background_files(backgroundFiles()),
                                                                                        m_genMC_files(genMCFiles()),
                                                                                        m_accMC_files(accMCFiles()),
                                                                                        m_data_tree_name(findTreeName("data")),
                                                                                        m_background_tree_name(findTreeName("background")),
                                                                                        m_genMC_tree_name(findTreeName("genMC")),
                                                                                        m_accMC_tree_name(findTreeName("accMC"))
{

    report(DEBUG, kModule) << "Constructing RootDataConverter for file: " << m_fit_file << "\n";

    if (!m_fit_results.valid())
    {
        report(ERROR, kModule) << "RootDataConverter ERROR: Unable to read fit results from file " + m_fit_file << "\n";
        assert(false);
    }

    setLowerVertexIndices({1}); // default to index 1 being the lower vertex particle
    setUpperVertexIndices();

    // print debug info detailing which particles are labeled as upper or lower vertex 
    // note that index 0 is always the beam photon
    const std::vector<std::string> particle_list = m_cfg_info->reaction(m_fit_results.reactionList()[0])->particleList();
    report(DEBUG, kModule) << "The following particles have been labeled as lower vertex particles by default (index 1):\n";
    report(DEBUG, kModule) << "i = 1 : " << particle_list[1] << "\n";
    report(DEBUG, kModule) << "The following particles are thus labeled as upper vertex particles by default (indices 2 to N):\n";
    for (const auto &idx : m_upper_vertex_indices)
    {
        report(DEBUG, kModule) << "i = " << idx << " : " << particle_list[idx] << "\n";
    }

    m_background_files_exist = !m_background_files.empty();

    validateDataFiles();
    validateBackgroundFiles(); // if no bkgnd files exist, warns user that weights assumed to be in data tree
    validateGenMCFiles();
    validateAccMCFiles();

    // fill the m_values map with all the relevant information from the fit results and
    // associated files that we want to save to csv later. This includes:
    // - Bin edges, centers, averages, and RMS values
    // - total number of detected and acceptance-corrected events
    // - efficiency based on genMC and accMC event counts
    extract();
}

// TODO: this block should be used for the constructor that specifies custom lv indices
// for (size_t i = 0; i < particle_list.size(); i++)
//     {
//         if (std::find(m_lower_vertex_indices.begin(), m_lower_vertex_indices.end(), i) != m_lower_vertex_indices.end())
//             report(DEBUG, kModule) << "i = " << i << " : " << particle_list[i] << "\n";
//     }

void RootDataConverter::extract()
{
    // get weight branch name. If weight branch is not found, will return empty string
    // and we'll assume weights of 1.0 throughout our stats functions
    std::string weight_branch_name;
    if (m_background_files_exist)
        weight_branch_name = weightBranchName("background", m_background_tree_name);        
    else
        weight_branch_name = weightBranchName("data", m_data_tree_name);

    if (weight_branch_name.empty())
    {
        report(DEBUG, kModule) << "No weight branch found. Assuming weights of 1.0 for all events.\n";
    }
    else
    {
        report(DEBUG, kModule) << "Weight branch found: " << weight_branch_name << "\n";
    }

    // Extract min, max, mean, and RMS of the beam energy, -t, and mass of the upper 
    // vertex system. Background subtraction and proper weighting are all included.
    extractBeamEnergyStats(weight_branch_name);
    extractFourMomentumTransferStats(weight_branch_name);
    // extractUpperVertexMassStats(weight_branch_name);

    // TODO: next step is to calculate the upper vertex mass from the remaining
    // particles. Finally, get number of events, and calculate efficiency from MC. Then
    // use efficiency to get ac_events.
}

std::vector<std::string> RootDataConverter::findFiles(const std::string &file_type) const
{
    std::vector<std::string> files;

    // get directory of the fit file
    std::filesystem::path dir = std::filesystem::path(m_fit_file).parent_path();

    const std::vector<std::string> reaction_list = m_fit_results.reactionList();
    for (std::string reaction : reaction_list)
    {
        std::pair<std::string, std::vector<std::string>> file_pair;

        if (file_type == "data")
            file_pair = m_cfg_info->reaction(reaction)->data();
        else if (file_type == "background")
            file_pair = m_cfg_info->reaction(reaction)->bkgnd();
        else if (file_type == "genMC")
            file_pair = m_cfg_info->reaction(reaction)->genMC();
        else if (file_type == "accMC")
            file_pair = m_cfg_info->reaction(reaction)->accMC();
        else
        {
            report(ERROR, kModule) << "Unknown file type requested: "
                                   << file_type
                                   << "\n";
            assert(false);
        }

        std::string data_reader_name = file_pair.first;
        std::vector<std::string> data_reader_args = file_pair.second;

        std::filesystem::path file_path = data_reader_args[0];

        // assumes that relative paths are relative to the fit file directory
        if (!file_path.is_absolute())
            file_path = dir / file_path;

        // canonical will remove any ../ or ./ from the path, while also checking
        // for existence
        try
        {
            files.push_back(std::filesystem::canonical(file_path).string());
        }
        catch (const std::filesystem::filesystem_error &e)
        {
            report(ERROR, kModule) << "Error accessing "
                                   << file_type
                                   << " file: "
                                   << file_path << "\n"
                                   << e.what() << "\n";
            assert(false);
        }
    }
    return files;
}

std::string RootDataConverter::findTreeName(const std::string &file_type) const
{
    const std::vector<std::string> reaction_list = m_fit_results.reactionList();

    // assume all reactions use the same tree name, so just get the first one
    std::string reaction = reaction_list[0];
    std::pair<std::string, std::vector<std::string>> file_pair;
    std::string file_path;
    if (file_type == "data")
    {
        file_pair = m_cfg_info->reaction(reaction)->data();
        file_path = m_data_files[0];
    }
    else if (file_type == "background")
    {
        file_pair = m_cfg_info->reaction(reaction)->bkgnd();
        file_path = m_background_files[0];
    }
    else if (file_type == "genMC")
    {
        file_pair = m_cfg_info->reaction(reaction)->genMC();
        file_path = m_genMC_files[0];
    }
    else if (file_type == "accMC")
    {
        file_pair = m_cfg_info->reaction(reaction)->accMC();
        file_path = m_accMC_files[0];
    }
    else
    {
        report(ERROR, kModule) << "Unknown file type requested: "
                               << file_type
                               << "\n";
        assert(false);
    }

    // if only one tree exists in the file, return that tree name
    std::vector<std::string> tree_names = treesInFile(file_path);
    if (tree_names.size() == 1)
        return tree_names[0];

    // otherwise, we have to extract the tree name from the data reader arguments
    // which unfortunately differ for each data reader
    std::string data_reader_name = file_pair.first;
    std::vector<std::string> data_reader_args = file_pair.second;

    if (data_reader_name == "ROOTDataReader" && data_reader_args.size() == 2)
    {
        return data_reader_args[1];
    }
    else if (data_reader_name == "ROOTDataReaderBootstrap" && data_reader_args.size() == 3)
    {
        return data_reader_args[2];
    }
    else if (data_reader_name == "FSRootDataReader")
    {
        return data_reader_args[1];
    }
    else if (data_reader_name == "FSRootDataReaderBootstrap")
    {
        return data_reader_args[1];
    }
    else
    {
        if (!m_mute_warnings)
            report(WARNING, kModule) << "Unknown data reader or insufficient arguments to"
                                     << " determine tree name. Defaulting to 'kin'\n";
        return "kin"; // default tree name
    }
}

std::vector<std::string> RootDataConverter::treesInFile(const std::string &filename) const
{
    std::vector<std::string> tree_names;

    TFile *f = TFile::Open(filename.c_str());
    if (!f || f->IsZombie())
    {
        report(ERROR, kModule) << "Error opening ROOT file: " + filename + "\n";
        assert(false);
    }

    TIter nextkey(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)nextkey()))
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (cl->InheritsFrom("TTree"))
        {
            tree_names.push_back(key->GetName());
        }
    }

    f->Close();
    return tree_names;
}

std::string RootDataConverter::weightBranchName(
    const std::string &file_type,
    const std::string &tree_name) const
{
    // TODO: verify that this friend tree works. May have to add friend tree to each
    // file to access the weight branch later

    // FSRoot-based data readers can optionally pass a friend weight branch name as an
    // argument, so check explicitly for that first
    const std::vector<std::string> reaction_list = m_fit_results.reactionList();
    std::string reaction = reaction_list[0]; // assume all reactions use the same weight branch
    std::pair<std::string, std::vector<std::string>> file_pair; // data reader name and args
    std::string root_file;
    if (file_type == "data")
    {
        file_pair = m_cfg_info->reaction(reaction)->data();
        root_file = m_data_files[0];
    }
    else if (file_type == "background")
    {
        file_pair = m_cfg_info->reaction(reaction)->bkgnd();
        root_file = m_background_files[0];
    }
    else
    {
        report(ERROR, kModule) << "Unknown file type requested: "
                               << file_type
                               << "\n";
        assert(false);
    }

    std::string data_reader_name = file_pair.first;
    std::vector<std::string> data_reader_args = file_pair.second;

    // the FSRootDataReader can have special friend trees that hold the event weights
    if (data_reader_name.find("FSRootDataReader") != std::string::npos && data_reader_args.size() >= 6)
    {
        std::string friend_file_name = data_reader_args[3];
        std::string friend_tree_name = data_reader_args[4];
        std::string weight_branch_name = data_reader_args[5];

        // to access the weight branch, we need to add the friend tree
        TFile *f = TFile::Open(root_file.c_str());
        if (!f || f->IsZombie())
        {
            report(ERROR, kModule) << "Error opening ROOT file: " + root_file + "\n";
            assert(false);
        }
        TTree *tree = f->Get<TTree>(tree_name.c_str());
        if (!tree)
        {
            report(ERROR, kModule) << "Error: Tree '" + tree_name + "' not found in file: " + root_file + "\n";
            assert(false);
        }
        tree->AddFriend(friend_tree_name.c_str(), friend_file_name.c_str());
        return friend_tree_name + "." + weight_branch_name;
    }

    // not FSROOTDataReader, so no friend tree weights. Proceed to check for weight
    // branches in the main tree
    TFile *f = TFile::Open(root_file.c_str());
    if (!f || f->IsZombie())
    {
        report(ERROR, kModule) << "Error opening ROOT file: " + root_file + "\n";
        assert(false);
    }
    TTree *tree = f->Get<TTree>(tree_name.c_str());
    if (!tree)
    {
        report(ERROR, kModule) << "Error: Tree '" + tree_name + "' not found in file: " + root_file + "\n";
        assert(false);
    }
    // check if weight branch exists before returning it. Most analyses use "weight" or
    // "Weight" branches
    if (tree->GetBranch("weight"))
    {
        f->Close();
        return "weight";
    }
    else if (tree->GetBranch("Weight"))
    {
        f->Close();
        return "Weight";
    }
    else
    {
        f->Close();
        if (!m_mute_warnings)
            report(WARNING, kModule) << "Error: Weight branch not found in tree '"
                                     << tree_name
                                     << "' in file: "
                                     << root_file
                                     << ". Setting weights to 1.0\n";
        return ""; // return empty string to indicate no weight branch found
    }
}


void RootDataConverter::extractBeamEnergyStats(const std::string &weight_branch_name)
{
    // we will add all the data (and possible) weighted background energy histograms
    // together to get overall beam energy stats
    TH1D *h_energy_total = nullptr;

    double min, max;
    std::tie(min, max) = findMinMaxOfBranch(m_data_files, m_data_tree_name, "EnPB");

    if(m_background_files_exist)
    {
        double bg_min, bg_max;
        std::tie(bg_min, bg_max) = findMinMaxOfBranch(m_background_files, m_background_tree_name, "EnPB");

        // set the min and max for the total histogram to encompass both data and background
        min = std::min(min, bg_min);
        max = std::max(max, bg_max);
    }

    // ==== DATA ENERGY HISTOGRAMS ====
    for (const std::string &data_file : m_data_files)
    {
        TFile *f = TFile::Open(data_file.c_str());
        if (!f || f->IsZombie())
        {
            report(ERROR, kModule) << "Error opening data ROOT file: " << data_file << "\n";
            assert(false);
        }
        TTree *tree = f->Get<TTree>(m_data_tree_name.c_str());
        if (!tree)
        {
            report(ERROR, kModule) << "Error retrieving tree " << m_data_tree_name << " from file " << data_file << "\n";
            f->Close();
            assert(false);
        }

        // Set up branch addresses
        double energy;
        tree->SetBranchAddress("EnPB", &energy);

        // draw into temporary histogram
        const int n_bins = 200;        
        TH1D *h_energy = new TH1D("h_energy", "", n_bins, min, max);        

        // If no background files exist, we assume weights are stored in the data tree
        if (!m_background_files_exist && !weight_branch_name.empty() && tree->GetBranch(weight_branch_name.c_str()))
        {
            double weight;
            tree->SetBranchAddress(weight_branch_name.c_str(), &weight);
            tree->Draw("EnPB>>h_energy", weight_branch_name.c_str(), "goff");
        }
        else if (m_background_files_exist)
        {
            tree->Draw("EnPB>>h_energy", "1.0", "goff"); // if bkgnd files exist, we will subtract them later, so just use weight of 1.0 here
        }
        else
        {
            tree->Draw("EnPB>>h_energy", "1.0", "goff");
            if (!m_mute_warnings)
                report(WARNING, kModule) << "No weight branch found in tree '"
                                         << m_data_tree_name
                                         << "' in file: "
                                         << data_file
                                         << "\n"
                                         << "Assuming weights of 1.0 for all events.\n";
        }        

        // if this is first hist, clone it to the total hist. Otherwise, add to the total hist
        if (!h_energy_total)
        {
            h_energy_total = (TH1D *)h_energy->Clone("h_energy_total");
            h_energy_total->SetDirectory(0); // detach from file
        }
        else
        {
            h_energy_total->Add(h_energy);
        }

        delete h_energy; // clean up temporary histogram
        f->Close();
    }

    // ==== BACKGROUND ENERGY HISTOGRAMS ====
    // if background files exist, we'll subtract the weighted background energy
    // histograms from the total energy histogram, which we will use for our stats. 
    for (const std::string &bg_file : m_background_files)
    {
        TFile *f = TFile::Open(bg_file.c_str());
        if (!f || f->IsZombie())
        {
            report(ERROR, kModule) << "Error opening background ROOT file: " << bg_file << "\n";
            assert(false);
        }
        TTree *tree = f->Get<TTree>(m_background_tree_name.c_str());
        if (!tree)
        {
            report(ERROR, kModule) << "Error retrieving tree " << m_background_tree_name << " from file " << bg_file << "\n";
            f->Close();
            assert(false);
        }

        // Set up branch addresses
        double energy;
        double weight = 1.0; // default weight is 1.0 if no weight branch exists

        tree->SetBranchAddress("EnPB", &energy);

        // If weight branch exists, set its address
        if (!weight_branch_name.empty() && tree->GetBranch(weight_branch_name.c_str()))
        {
            tree->SetBranchAddress(weight_branch_name.c_str(), &weight);
        }

        // draw into temporary histogram
        const int n_bins = 200;
        TH1D *h_energy = new TH1D("h_energy", "", n_bins, min, max);
        tree->Draw("EnPB>>h_energy", weight_branch_name.c_str(), "goff");        

        // subtract background histogram from total histogram
        if (h_energy_total)
            h_energy_total->Add(h_energy, -1.0);

        delete h_energy;
        f->Close();        
    }

    const double center_energy = 0.5 * (min + max);
    const double mean_energy = h_energy_total->GetMean();
    const double rms_energy = h_energy_total->GetRMS();

    // Store values in the map
    m_values["e_low"] = min;
    m_values["e_high"] = max;
    m_values["e_mid"] = center_energy;
    m_values["e_avg"] = mean_energy;
    m_values["e_rms"] = rms_energy;

    delete h_energy_total;

    report(DEBUG, kModule) << "Beam energy statistics:\n";
    report(DEBUG, kModule) << "  Min: " << min << "\n";
    report(DEBUG, kModule) << "  Max: " << max << "\n";
    report(DEBUG, kModule) << "  Center: " << center_energy << "\n";
    report(DEBUG, kModule) << "  Mean: " << mean_energy << "\n";
    report(DEBUG, kModule) << "  RMS: " << rms_energy << "\n";
}

void RootDataConverter::extractFourMomentumTransferStats(const std::string &weight_branch_name)
{
    // target proton mass (GeV)
    const double m_proton = 0.938;

    if (m_lower_vertex_indices.empty()) {
        report(ERROR, kModule) << "No lower vertex indices set; cannot compute -t\n";
        assert(false);
    }

    // Prepare branch expression strings like "PxP1 + PxP2 + ..."
    std::string recoil_px_expr;
    std::string recoil_py_expr;
    std::string recoil_pz_expr;
    std::string recoil_e_expr;
    for (size_t i = 0; i < m_lower_vertex_indices.size(); ++i) {
        int idx = m_lower_vertex_indices[i];
        std::string sep = (i == 0) ? "" : " + ";
        recoil_px_expr += sep + "PxP" + std::to_string(idx);
        recoil_py_expr += sep + "PyP" + std::to_string(idx);
        recoil_pz_expr += sep + "PzP" + std::to_string(idx);
        recoil_e_expr  += sep + "EnP" + std::to_string(idx);
    }

    // ---- DATA histogram via RDataFrame ----
    try {
        ROOT::RDataFrame df_data(m_data_tree_name, m_data_files);

        // Define recoil components and t (Mandelstam t = (target - recoil)^2)
        auto df = df_data.Define("recoil_px", recoil_px_expr)
                         .Define("recoil_py", recoil_py_expr)
                         .Define("recoil_pz", recoil_pz_expr)
                         .Define("recoil_E", recoil_e_expr)
                         .Define("t_value",
                                 [m_proton](double rE, double rpx, double rpy, double rpz) {
                                     double dE = m_proton - rE;
                                     double p2 = rpx * rpx + rpy * rpy + rpz * rpz;
                                     double t = dE * dE - p2;
                                     return std::fabs(t);
                                 },
                                 {"recoil_E", "recoil_px", "recoil_py", "recoil_pz"});

        // Determine min and max from data (unweighted)
        auto data_min_r = df.Min("t_value");
        auto data_max_r = df.Max("t_value");
        double data_min = *data_min_r;
        double data_max = *data_max_r;

        double global_min = data_min;
        double global_max = data_max;

        // Build data histogram (if background exists, don't apply weights to data here)
        const int n_bins = 200;
        ROOT::RDF::RResultPtr<TH1D> h_data;
        if (m_background_files_exist) {
            h_data = df.Histo1D({"h_t_data", "t (data)", n_bins, global_min, global_max}, "t_value");
        }
        else {
            // if no background files, use weight branch from data (if provided)
            if (!weight_branch_name.empty()) {
                // if weight branch references a friend (contains a '.'), unclear if RDataFrame can handle that, so skip weights in that case with a warning
                // TODO: test the friend tree case
                if (weight_branch_name.find('.') != std::string::npos) {
                    report(WARNING, kModule) << "Weight branch appears to reference a friend tree (" << weight_branch_name << ") - skipping weights for RDataFrame histogram\n";
                    h_data = df.Histo1D({"h_t_data", "t (data)", n_bins, global_min, global_max}, "t_value");
                }
                else {
                    h_data = df.Histo1D({"h_t_data", "t (data)", n_bins, global_min, global_max}, "t_value", weight_branch_name);
                }
            }
            else {
                h_data = df.Histo1D({"h_t_data", "t (data)", n_bins, global_min, global_max}, "t_value");
            }
        }

        // ---- BACKGROUND histogram via RDataFrame (if present) ----        
        if (m_background_files_exist && !m_background_files.empty()) {
            ROOT::RDF::RResultPtr<TH1D> h_bkg;
            ROOT::RDataFrame df_bg(m_background_tree_name, m_background_files);
            auto dfb = df_bg.Define("recoil_px", recoil_px_expr)
                          .Define("recoil_py", recoil_py_expr)
                          .Define("recoil_pz", recoil_pz_expr)
                          .Define("recoil_E", recoil_e_expr)
                          .Define("t_value",
                                 [m_proton](double rE, double rpx, double rpy, double rpz) {
                                     double dE = m_proton - rE;
                                     double p2 = rpx * rpx + rpy * rpy + rpz * rpz;
                                     double t = dE * dE - p2;
                                     return std::fabs(t);
                                 },
                                 {"recoil_E", "recoil_px", "recoil_py", "recoil_pz"});

            // Expand global range to include background extremes
            auto bg_min_r = dfb.Min("t_value");
            auto bg_max_r = dfb.Max("t_value");
            double bg_min = *bg_min_r;
            double bg_max = *bg_max_r;
            global_min = std::min(global_min, bg_min);
            global_max = std::max(global_max, bg_max);

            // Recreate histograms with unified ranges
            h_data = df.Histo1D({"h_t_data", "t (data)", n_bins, global_min, global_max}, "t_value");

            // background histogram should use weight branch (weight_branch_name was chosen from background earlier)
            if (!weight_branch_name.empty() && weight_branch_name.find('.') == std::string::npos) {
                h_bkg = dfb.Histo1D({"h_t_total", "t (bg)", n_bins, global_min, global_max}, "t_value", weight_branch_name);
            }
            else {
                if (weight_branch_name.find('.') != std::string::npos)
                    report(WARNING, kModule) << "Background weight branch references friend tree; skipping weights for bg RDataFrame histogram\n";
                h_bkg = dfb.Histo1D({"h_t_total", "t (bg)", n_bins, global_min, global_max}, "t_value");
            }

            // Subtract background from a standalone clone so RResultPtr ownership stays intact.
            TH1D *h_result = (TH1D *)h_data->Clone("h_t_result");
            h_result->SetDirectory(nullptr);
            h_result->Add(h_bkg.GetPtr(), -1.0);

            double mean_t = h_result->GetMean();
            double rms_t = h_result->GetRMS();
            double center_t = 0.5 * (global_min + global_max);

            m_values["t_low"] = global_min;
            m_values["t_high"] = global_max;
            m_values["t_mid"] = center_t;
            m_values["t_avg"] = mean_t;
            m_values["t_rms"] = rms_t;

            report(DEBUG, kModule) << "-t statistics (data - background):" << "\n";
            report(DEBUG, kModule) << "  Min: " << global_min << "\n";
            report(DEBUG, kModule) << "  Max: " << global_max << "\n";
            report(DEBUG, kModule) << "  Center: " << center_t << "\n";
            report(DEBUG, kModule) << "  Mean: " << mean_t << "\n";
            report(DEBUG, kModule) << "  RMS: " << rms_t << "\n";

            delete h_result;
        }
        else {
            // No background files: final histogram is h_data
            double mean_t = h_data->GetMean();
            double rms_t = h_data->GetRMS();
            double center_t = 0.5 * (global_min + global_max);

            m_values["t_low"] = global_min;
            m_values["t_high"] = global_max;
            m_values["t_mid"] = center_t;
            m_values["t_avg"] = mean_t;
            m_values["t_rms"] = rms_t;

            report(DEBUG, kModule) << "-t statistics (data):\n";
            report(DEBUG, kModule) << "  Min: " << global_min << "\n";
            report(DEBUG, kModule) << "  Max: " << global_max << "\n";
            report(DEBUG, kModule) << "  Center: " << center_t << "\n";
            report(DEBUG, kModule) << "  Mean: " << mean_t << "\n";
            report(DEBUG, kModule) << "  RMS: " << rms_t << "\n";
        }
    }
    catch (const std::exception &e) {
        report(ERROR, kModule) << "RDataFrame error computing -t: " << e.what() << "\n";
        assert(false);
    }
}


void RootDataConverter::validateFiles(const std::vector<std::string> &files,
                                       const std::string &file_type) const
{
    if (files.empty() && file_type != "background")
    {
        report(ERROR, kModule) << "No " 
                               << file_type << " files found in fit results for file: "
                               << m_fit_file << "\n";
        assert(false);
    }
    else if (files.empty() && file_type == "background" && !m_mute_warnings)
    {
        report(WARNING, kModule) << "No background files found in fit results for file: "
                                 << m_fit_file
                                 << ". Assuming weights are stored directly in the data"
                                 << " trees.\n";
    }
    else
    {
        return; // files are valid
    }
}

std::pair<double, double> RootDataConverter::findMinMaxOfBranch(const std::vector<std::string> &files,
                                               const std::string &tree_name,
                                               const std::string &branch_name) const
{
    double global_min = std::numeric_limits<double>::max();
    double global_max = std::numeric_limits<double>::lowest();

    for (const std::string &file : files)
    {
        TFile *f = TFile::Open(file.c_str());
        if (!f || f->IsZombie())
        {
            report(ERROR, kModule) << "Error opening ROOT file: " << file << "\n";
            assert(false);
        }
        TTree *tree = f->Get<TTree>(tree_name.c_str());
        if (!tree)
        {
            report(ERROR, kModule) << "Error retrieving tree " << tree_name << " from file " << file << "\n";
            f->Close();
            assert(false);
        }

        double local_min = tree->GetMinimum(branch_name.c_str());
        double local_max = tree->GetMaximum(branch_name.c_str());

        if (local_min < global_min)
            global_min = local_min;
        if (local_max > global_max)
            global_max = local_max;

        f->Close();
    }

    return {global_min, global_max};
}

void RootDataConverter::setUpperVertexIndices()
{
    m_upper_vertex_indices.clear();
    const std::vector<std::string> particle_list = m_cfg_info->reaction(m_fit_results.reactionList()[0])->particleList();

    // Note that we always skip particle i=0, as this will be the beam photon
    for (size_t i = 1; i < particle_list.size(); i++)
    {
        if (std::find(m_lower_vertex_indices.begin(), m_lower_vertex_indices.end(), i) == m_lower_vertex_indices.end())
        {
            m_upper_vertex_indices.push_back(i);
        }
    }
}