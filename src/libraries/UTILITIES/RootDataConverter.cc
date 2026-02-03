/**
 * @file RootDataConverter.cc
 * @author Kevin Scheuer
 * @brief Implementation of RootDataConverter class
 * @date 2026-02-02
 *
 */

#include <cassert>
#include <filesystem>

#include "RootDataConverter.h"
#include "IUAmpTools/report.h"

#include "TFile.h"
#include "TKey.h"
#include "TClass.h"
#include "TTree.h"

const char *RootDataConverter::kModule = "RootDataConverter";

RootDataConverter::RootDataConverter(const std::string &filename, bool mute_warnings) : m_fit_file(filename),
                                                                                        m_fit_results(filename, mute_warnings),
                                                                                        m_cfg_info(m_fit_results.configInfo()),
                                                                                        m_mute_warnings(mute_warnings)
{
    if (!m_fit_results.valid())
    {
        report(ERROR, kModule) << "RootDataConverter ERROR: Unable to read fit results from file " + m_fit_file << "\n";
        assert(false);
    }
    setLowerVertexIndices({1}); // default to index 1 being the lower vertex particle
    extract();
}

void RootDataConverter::extract()
{
    // extract data files, background files, genMC, and accMC files
    std::vector<std::string> data_files = findFiles("data");
    std::vector<std::string> background_files = findFiles("background");
    std::vector<std::string> genMC_files = findFiles("genMC");
    std::vector<std::string> accMC_files = findFiles("accMC");

    if (data_files.empty())
    {
        report(ERROR, kModule) << "No data files found in fit results for file: "
                               << m_fit_file << "\n";
    }
    if (background_files.empty() && !m_mute_warnings)
    {
        report(WARNING, kModule) << "No background files found in fit results for file: "
                                 << m_fit_file
                                 << ". Assuming weights are stored directly in the data"
                                 << " trees.\n";
    }
    if (genMC_files.empty())
    {
        report(ERROR, kModule) << "No genMC files found in fit results for file: "
                               << m_fit_file << "\n";
    }
    if (accMC_files.empty())
    {
        report(ERROR, kModule) << "No accMC files found in fit results for file: "
                               << m_fit_file << "\n";
    }

    // get tree names for each file
    std::string data_tree_name = dataTreeName();
    std::string background_tree_name = (background_files.empty()) ? "" : backgroundTreeName();
    std::string genMC_tree_name = genMCTreeName();
    std::string accMC_tree_name = accMCTreeName();

    // get weight branch name. If weight branch is not found, will return empty string
    // and we'll assume weights of 1.0 later
    std::string weight_branch_name;
    if (background_files.empty())
        weight_branch_name = weightBranchName("data", data_tree_name);
    else
        weight_branch_name = weightBranchName("background", background_tree_name);

    // get beamE branch name
    std::string beamE_branch_name = beamEBranchName("data");

    // TODO: next step is to calculate t from the beam and lower vertex particle
    // 4-vectors. Look at other plotters for reference. Then calculate the upper
    // vertex mass from the remaining particles. Finally, fill histograms and
    // calculate bin info with a better mapping between the header and value storage.
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
            file_pair = m_cfg_info->reaction(reaction)->background();
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
        file_path = findFiles("data")[0];
    }
    else if (file_type == "background")
    {
        file_pair = m_cfg_info->reaction(reaction)->background();
        file_path = findFiles("background")[0];
    }
    else if (file_type == "genMC")
    {
        file_pair = m_cfg_info->reaction(reaction)->genMC();
        file_path = findFiles("genMC")[0];
    }
    else if (file_type == "accMC")
    {
        file_pair = m_cfg_info->reaction(reaction)->accMC();
        file_path = findFiles("accMC")[0];
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
    const std::pair<std::string, std::vector<std::string>> file_pair;
    std::string root_file;
    if (file_type == "data")
    {
        file_pair = m_cfg_info->reaction(reaction)->data();
        root_file = findFiles("data")[0];
    }
    else if (file_type == "background")
    {
        file_pair = m_cfg_info->reaction(reaction)->background();
        root_file = findFiles("background")[0];
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
    else // not FSDataReader, so no friend tree weights
    {
        // proceed to check for weight branches in the main tree
    }

    // otherwise, simply check if "weight" or "Weight" branch exists in the tree
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

std::string beamEBranchName(const std::string &file_type) const
{
    std::string file_path;
    std::string tree_name;
    if (file_type == "data")
    {
        file_path = findFiles("data")[0];
        tree_name = dataTreeName();
    }
    else if (file_type == "background")
    {
        file_path = findFiles("background")[0];
        tree_name = backgroundTreeName();
    }
    else if (file_type == "genMC")
    {
        file_path = findFiles("genMC")[0];
        tree_name = genMCTreeName();
    }
    else if (file_type == "accMC")
    {
        file_path = findFiles("accMC")[0];
        tree_name = accMCTreeName();
    }
    else
    {
        report(ERROR, kModule) << "Unknown file type requested: "
                               << file_type
                               << "\n";
        assert(false);
    }

    TFile *f = TFile::Open(file_path.c_str());
    if (!f || f->IsZombie())
    {
        report(ERROR, kModule) << "Error opening ROOT file: " + file_path + "\n";
        assert(false);
    }
    TTree *tree = f->Get<TTree>(tree_name.c_str());
    if (!tree)
    {
        report(ERROR, kModule) << "Error: Tree '" + tree_name + "' not found in file: " + file_path + "\n";
        assert(false);
    }

    if (tree->GetBranch("E_Beam"))
    {
        f->Close();
        return "E_Beam";
    }
    else if (tree->GetBranch("EnPB"))
    {
        f->Close();
        return "EnPB";
    }
    else
    {
        f->Close();
        report(ERROR, kModule) << "Error: beamE branch not found in tree '"
                               << tree_name
                               << "' in file: "
                               << file_path
                               << "\n";
        assert(false);
    }
}