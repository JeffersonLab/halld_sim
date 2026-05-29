#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "AMPTOOLS_AMPS/ComplexCoeff.h"
#include "AMPTOOLS_AMPS/OmegaDalitz.h"
#include "AMPTOOLS_AMPS/PhaseOffset.h"
#include "AMPTOOLS_AMPS/Piecewise.h"
#include "AMPTOOLS_AMPS/Uniform.h"
#include "AMPTOOLS_AMPS/Vec_ps_moment.h"
#include "AMPTOOLS_AMPS/Vec_ps_refl.h"
#include "AMPTOOLS_DATAIO/FSRootDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderBootstrap.h"
#include "AMPTOOLS_DATAIO/ROOTDataReaderTEM.h"
#include "AMPTOOLS_DATAIO/VecPsPlotGenerator.h"
#include "AmpPlotter/PlotFactory.h"
#include "AmpPlotter/PlotterMainWindow.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/FitResults.h"
#include "MinuitInterface/MinuitMinimizationManager.h"
#include "TApplication.h"
#include "TClass.h"
#include "TFile.h"
#include "TGClient.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"

typedef VecPsPlotGenerator vecps_PlotGen;

void atiSetup() {
    AmpToolsInterface::registerAmplitude(BreitWigner());
    AmpToolsInterface::registerAmplitude(Uniform());
    AmpToolsInterface::registerAmplitude(Vec_ps_refl());
    AmpToolsInterface::registerAmplitude(PhaseOffset());
    AmpToolsInterface::registerAmplitude(ComplexCoeff());
    AmpToolsInterface::registerAmplitude(Piecewise());
    AmpToolsInterface::registerAmplitude(OmegaDalitz());
    AmpToolsInterface::registerAmplitude(Vec_ps_moment());

    AmpToolsInterface::registerDataReader(ROOTDataReader());
    AmpToolsInterface::registerDataReader(FSRootDataReader());
    AmpToolsInterface::registerDataReader(ROOTDataReaderTEM());
}

// ...oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//              Enum Classes to Simplify Boolean Logic
// ...oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

// Classify the sum of the reaction according to reflectivity
enum class ReflSum {
    Neither,     // something like an incoherent background sum
    Positive,
    Negative
};

// Specify the group of amplitudes to be plotted
enum class AmpPlotGroup{
    AllAmps,      // for the total fit result
    AllReflAmps,  // for the total refl 
    JlFamily,     // sum over the m projections 
    SingleAmp     // one amplitude at a time    
};

// Specify the group of sums to be plotted
enum class SumPlotGroup{
    AllSums,         // for the total fit results
    BothReflSums,    // Incoherent sum of both reflectivities
    PosReflSum,      // one coherent sum
    NegReflSum,      // one coherent sum
    NoReflSum        // This could be broken down for multiple non-refl sums
};

// ...oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//              Structures to Gather Relevant Info Together
// ...oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

// Store the information from the command line arguments and configuration file
// and allows to easily change the default values of the options
// TODO: add an option for only GUI mode (aka no root file output)
struct ProgramOptions{
    std::string resultsName;
    std::string outName = "vecps_plot.root";
    bool showGui = false;
    bool mcGen = false;
    bool bkgFile = false;
    bool sFile = false;
    bool skipPlots = false;
    bool txtOutput = false;
    bool csvOutput = false;
};

// Store the information needed to identify a wave 
struct WaveKey {
    std::string ampName;     // full amp name (e.g. 1S+1)
    std::string jFamily;     // Jl family (e.g. 1S)
    std::string m;           // m projection (e.g. +1)
    bool isNumeric;          // non-Jl amps do not start with a number                      
};

// Store all the information that identifies an amplitude
struct AmplitudeKey {
    std::string fullName;        // unparsed AmpTools full name
    std::string reactionName; 
    std::string sumName;         
    ReflSum reflSum;             // Pos/neg refl or neither
    WaveKey key; 

    // index in uniqueAmplitudes()
    size_t ampIndex  = std::numeric_limits<size_t>::max();
    // index in uniqueSums()
    size_t sumIndex  = std::numeric_limits<size_t>::max();
};

// Store all the relation between the unique sums and unique amps
// could be extended to global (between reactions) relations
struct Sisterhood{
    // maps to look up the indices of the unique amplitudes and sums
    // the key is the string/enum name of the amplitude/sum  
    std::map<std::string, std::vector<size_t>> jlIndexMap;
    std::map<ReflSum, std::vector<size_t>> sumIndexMap;

    // Relationship between sums and each amplitude (index-based)
    std::vector<std::vector<size_t>> ampToSums; 
    std::vector<std::vector<size_t>> sumToAmps; 
    // Note: amplitudes that do not have refl are not in the following map
    std::map<ReflSum, std::vector<size_t>> reflToAmps; 

    // Relations for amplitudes based on the full name
    std::vector<std::vector<size_t>> ampToFull;   // uniqueAmp to full indices
    std::vector<std::vector<size_t>> sumToFull;   // uniqueSum to full indices
    std::map<ReflSum, std::vector<size_t>> reflToFull;
    std::map<std::pair<size_t, ReflSum>, std::vector<size_t>> ampReflToFull;
    std::map<std::pair<std::string, ReflSum>, std::vector<size_t>> jlReflToFull;
};

// Store the information that is static when a reaction is read
// this structure avoids having to parse the string information again
// and it becomes a look up table for the indices of the amplitudes and 
// sums when plotting
struct AmplitudeRegistry {
    // Collection of all of the amplitudes in the fit file
    std::vector<AmplitudeKey> ampsList;

    // Parsed unique amplitudes
    std::vector<WaveKey> uniqueAmps;  
    size_t nUnqAmps;
    // Parsed unique sums
    std::vector<ReflSum> uniqueSums;
    size_t nUnqSums;
    // containers to look up the indices of the unique amplitudes and sums
    // and their relation to other sums/amps
    Sisterhood relations;
};

// This stores the specific configuration for a plotting
struct PlotConfig{
    AmpPlotGroup ampGroup;
    SumPlotGroup sumGroup;

    // Store the indices to be enable
    std::vector<size_t> ampIndices;
    std::vector<size_t> sumIndices;
};

struct PlotLabels{
    std::string human;        // human friendly output
    std::string histName;     // better for saving in a ROOT file
    // TODO: add csv format
    // std::string csvFormat; // better for saving in a CSV file
};

// ...oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//            Helper Functions to Fill the Structures, etc
// ...oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

// Parse the name of the sum to classify it according to reflectivity
// Note: this follow standard naming conventions for the vec_ps amplitudes
// this means that:
// [ImagNegSign & RealPosSign] = Negative Reflectivity
// [RealNegSign & ImagPosSign] = Positive Reflectivity
// The 3rd option is "Neither", but it can be easily expanded if needed
// dont forget to update the ReflSum enum if you do!
// General naming convention for amps:
// https://halldweb.jlab.org/doc-private/DocDB/ShowDocument?docid=7036
ReflSum classifyReflSum(const std::string& sumName){
    if (sumName.find("RealNegSign") != std::string::npos ||
        sumName.find("ImagPosSign") != std::string::npos){
        return ReflSum::Positive;  
    }
    else if (sumName.find("RealPosSign") != std::string::npos ||
             sumName.find("ImagNegSign") != std::string::npos){
        return ReflSum::Negative; 
    }
    else {
        return ReflSum::Neither;
    }
}

std::string reflToString(ReflSum r){
    switch(r){
        case ReflSum::Positive: return "PosRefl";
        case ReflSum::Negative: return "NegRefl";
        default: return "Neither";
    }
}

// Store in a vector the different parts of the full amplitude name, given 
// the common delimiter used in AmpTools (::)
std::vector<std::string> breakDownName(const std::string& fullAmpName,
                                       const std::string& delim){
    std::vector<std::string> brokenDownName;
    size_t start = 0; size_t end;

    while((end = fullAmpName.find(delim, start)) != std::string::npos){
        brokenDownName.push_back(fullAmpName.substr(start, end-start));
        start = end + delim.size();
    }

    brokenDownName.push_back(fullAmpName.substr(start));
    return brokenDownName;
}

WaveKey getJlm(const std::string& ampName){
    WaveKey key;
    key.ampName = ampName;
    // For the standard format "Jl.m" (e.g. "1S+1")
    // it splits the m-projection and the Jl family
    key.isNumeric =  std::isdigit(ampName[0]);
    if (!ampName.empty() && key.isNumeric){
            key.jFamily = ampName.substr(0, 2); 
            key.m = ampName.substr(2);        
    }
    else {
        key.jFamily = "";
        key.m = ""; 
    }
    return key;
}

AmplitudeKey parseAmp(const std::string& fullAmpName){
    AmplitudeKey ampKey;
    ampKey.fullName = fullAmpName;

    std::vector<std::string> tokens = breakDownName(fullAmpName, "::");

    // Expect exactly 3 tokens:
    // [reaction] :: [sum] :: [ampName]
    if (tokens.size() != 3){
        std::cerr << "Warning: Unexpected amplitude name format: " 
                  << fullAmpName << "\n";
        ampKey.reactionName = "";
        ampKey.sumName      = "";
        ampKey.key.jFamily  = "";
        ampKey.key.m        = "";
    } 
    else {
        ampKey.reactionName = tokens[0];
        ampKey.sumName      = tokens[1];
        ampKey.key          = getJlm(tokens[2]);
    }
    ampKey.reflSum = classifyReflSum(ampKey.sumName);

    return ampKey;
}

// Just a reminder that the registry is meant to be once per file, and 
// it serves as a look up table for everything related to the amps
AmplitudeRegistry buildRegistry(std::unique_ptr<vecps_PlotGen>& plotResults){    
    AmplitudeRegistry registry;
    
    std::vector<std::string> uniqueAmps = plotResults->uniqueAmplitudes();
    registry.nUnqAmps = uniqueAmps.size();
    std::vector<std::string> uniqueSums = plotResults->uniqueSums();
    registry.nUnqSums = uniqueSums.size();

    // Look up maps for parsed info
    std::map<std::string, size_t> ampNameToIndex;
    std::map<std::string, size_t> sumNameToIndex;

    std::vector<WaveKey> uniqueAmpsParsed;
    for (size_t i = 0; i < uniqueAmps.size(); ++i){
        ampNameToIndex[uniqueAmps[i]] = i;
        WaveKey key = getJlm(uniqueAmps[i]);
        uniqueAmpsParsed.push_back(key);
        if(!key.jFamily.empty()) // non-numeric amps have and empty jFamily
            registry.relations.jlIndexMap[key.jFamily].push_back(i);
    }
    registry.uniqueAmps = uniqueAmpsParsed;
    
    std::vector<ReflSum> uniqueSumsParsed;
    for (size_t i = 0; i < uniqueSums.size(); ++i){
        sumNameToIndex[uniqueSums[i]] = i;
        ReflSum refl = classifyReflSum(uniqueSums[i]);
        uniqueSumsParsed.push_back(refl);
        registry.relations.sumIndexMap[refl].push_back(i);
    }
    registry.uniqueSums = uniqueSumsParsed;
    
    std::vector<AmplitudeKey> parsedAmps;
    for (const auto& amp :  plotResults->fullAmplitudes()){
        AmplitudeKey ampKey = parseAmp(amp);
        ampKey.ampIndex = ampNameToIndex[ampKey.key.ampName];
        ampKey.sumIndex = sumNameToIndex[ampKey.sumName];
        ampKey.reflSum = classifyReflSum(ampKey.sumName);
        parsedAmps.push_back(ampKey);
    }
    registry.ampsList = parsedAmps;

    // Sisterhood vectors
    std::vector<std::vector<size_t>> ampToSums(registry.nUnqAmps);
    std::vector<std::vector<size_t>> sumToAmps(registry.nUnqSums);
    std::map<ReflSum, std::vector<size_t>> reflToAmps;

    // Sisterhood vectors for full reference
    std::vector<std::vector<size_t>> ampToFull(registry.nUnqAmps);
    std::vector<std::vector<size_t>> sumToFull(registry.nUnqSums);
    std::map<ReflSum, std::vector<size_t>> reflToFull;
    std::map<std::pair<size_t, ReflSum>, std::vector<size_t>> ampReflToFull;
    std::map<std::pair<std::string, ReflSum>, std::vector<size_t>> jlReflToFull;

    for(size_t iFull = 0; iFull < registry.ampsList.size(); ++iFull){
        const auto& amp = registry.ampsList[iFull];
        size_t ampIndex = amp.ampIndex;
        size_t sumIndex = amp.sumIndex;

        auto& amp2Sum = ampToSums[ampIndex];
        if(std::find(amp2Sum.begin(), amp2Sum.end(), sumIndex) == amp2Sum.end())
            amp2Sum.push_back(sumIndex);

        auto& sum2Amp = sumToAmps[sumIndex];
        if(std::find(sum2Amp.begin(), sum2Amp.end(), ampIndex) == sum2Amp.end())
            sum2Amp.push_back(ampIndex);

        auto& refl2Amp = reflToAmps[amp.reflSum];
        if(std::find(refl2Amp.begin(), refl2Amp.end(), ampIndex) == refl2Amp.end())
            refl2Amp.push_back(ampIndex);

        ampToFull[amp.ampIndex].push_back(iFull);
        sumToFull[amp.sumIndex].push_back(iFull);
        reflToFull[amp.reflSum].push_back(iFull);    

        ampReflToFull[{ampIndex, amp.reflSum}].push_back(iFull);
        if(!amp.key.jFamily.empty())
            jlReflToFull[{amp.key.jFamily, amp.reflSum}].push_back(iFull);
    }

    registry.relations.ampToSums = ampToSums;
    registry.relations.sumToAmps = sumToAmps;
    registry.relations.reflToAmps = reflToAmps;

    registry.relations.ampToFull = ampToFull;
    registry.relations.sumToFull = sumToFull;
    registry.relations.reflToFull = reflToFull;
    registry.relations.ampReflToFull = ampReflToFull;
    registry.relations.jlReflToFull = jlReflToFull;
    
    return registry;
}

// We want to avoid producing empty histograms. This normally happens because
// an active sum is paired to a set of amplitudes that does not belong to it.
// For example, an isotropic background added incoherently will not have 
// any reflectivity
bool isNonEmptyConfig(const std::vector<size_t>& amps,
                      const std::vector<size_t>& sums,
                      const Sisterhood& rel){
    if(amps.empty() || sums.empty())
        return false;

    std::unordered_set<size_t> sumSet(sums.begin(), sums.end());
    for(size_t a : amps){
        for(size_t s : rel.ampToSums[a]){
            if(sumSet.count(s))
                return true;
        }
    }
    return false; 
}

// The vector that this function outputs contains all the configurations that
// are allowed by the enum classes defined above with the caveat that we 
// screen for "physically meaningful" combinations in isNonEmptyConfig()
std::vector<PlotConfig> generatePlotConfigs(const AmplitudeRegistry& reg){
    std::vector<PlotConfig> configs;
    std::cout << reg.uniqueAmps.size() << " unique amps and "
              << reg.uniqueSums.size() << " unique sums." << std::endl;

    const std::vector<AmpPlotGroup> ampGroups = {
        AmpPlotGroup::AllAmps,
        AmpPlotGroup::AllReflAmps,
        AmpPlotGroup::SingleAmp,
        AmpPlotGroup::JlFamily
    };

    const std::vector<SumPlotGroup> sumGroups = {
        SumPlotGroup::AllSums,
        SumPlotGroup::PosReflSum,
        SumPlotGroup::NegReflSum,
        SumPlotGroup::NoReflSum,
        SumPlotGroup::BothReflSums
    };

    std::vector<size_t> allAmpIndices(reg.nUnqAmps);
    std::iota(allAmpIndices.begin(), allAmpIndices.end(), 0);

    std::vector<size_t> allSumIndices(reg.nUnqSums);
    std::iota(allSumIndices.begin(), allSumIndices.end(), 0);

    // Extract the relation from the registry
    const auto& rel = reg.relations;

    // These vector is meant to be used with AllReflAmps. Note that generally
    // both Positive and Negative share the same wave. But there are situations
    // where only one refl is present. So, it's safer to loop over both and
    // then remove duplicates
    std::vector<size_t> reflAmps;
    if(rel.reflToAmps.count(ReflSum::Positive))
        reflAmps.insert(reflAmps.end(),
            rel.reflToAmps.at(ReflSum::Positive).begin(),
            rel.reflToAmps.at(ReflSum::Positive).end());
    
    if(rel.reflToAmps.count(ReflSum::Negative))
        reflAmps.insert(reflAmps.end(),
            rel.reflToAmps.at(ReflSum::Negative).begin(),
            rel.reflToAmps.at(ReflSum::Negative).end());
    
    // Remove duplicates
    std::sort(reflAmps.begin(), reflAmps.end());
    reflAmps.erase(std::unique(reflAmps.begin(), reflAmps.end()), reflAmps.end());
    
    // To build all possible combinations we need to loop according to the 
    // following logic: (AmpPlotGroup × SumPlotGroup) × (ampSets × sumSets)
    // the first group, (AmpPlotGroup × SumPlotGroup), tell us what type of 
    // plot we want; while the second group, (ampSets × sumSets), tell us which 
    // amplitudes and sums to include in the plot
    for(auto aGroup : ampGroups){

        std::vector<std::vector<size_t>> ampSets;

        // build ampSets
        switch(aGroup){
            case AmpPlotGroup::AllAmps:
                ampSets.push_back(allAmpIndices);
                break;

            case AmpPlotGroup::AllReflAmps:
                 ampSets.push_back(reflAmps);
                break;

            case AmpPlotGroup::SingleAmp:
                for(size_t i : allAmpIndices)
                    ampSets.push_back({i});
                break;

            case AmpPlotGroup::JlFamily:
                for(const auto& [jl, indices] : rel.jlIndexMap)
                    ampSets.push_back(indices);
                break;
        }

        for(auto sGroup : sumGroups){
            // AmpPlotGroup::AllAmps can only be with SumPlotGroup::AllSums
            if ((aGroup == AmpPlotGroup::AllAmps) !=
                (sGroup == SumPlotGroup::AllSums))
                continue;

            // AmpPlotGroup::AllReflAmps can only be with 
            // SumPlotGroup::PosReflSum, SumPlotGroup::NegReflSum, or
            // SumPlotGroup::BothReflSums
            if (aGroup == AmpPlotGroup::AllReflAmps) {
                if (!(sGroup == SumPlotGroup::PosReflSum ||
                      sGroup == SumPlotGroup::NegReflSum ||
                      sGroup == SumPlotGroup::BothReflSums))
                    continue;
            }

            std::vector<std::vector<size_t>> sumSets;

            // build sumSets
            switch(sGroup){
                case SumPlotGroup::AllSums:
                    sumSets.push_back(allSumIndices);
                    break;

                case SumPlotGroup::PosReflSum:
                    if(rel.sumIndexMap.count(ReflSum::Positive))
                        sumSets.push_back(rel.sumIndexMap.at(ReflSum::Positive));
                    break;

                case SumPlotGroup::NegReflSum:
                    if(rel.sumIndexMap.count(ReflSum::Negative))
                        sumSets.push_back(rel.sumIndexMap.at(ReflSum::Negative));
                    break;

                case SumPlotGroup::NoReflSum:
                    if(rel.sumIndexMap.count(ReflSum::Neither))
                        sumSets.push_back(rel.sumIndexMap.at(ReflSum::Neither));
                    break;

                case SumPlotGroup::BothReflSums:{
                    std::vector<size_t> posNeg;

                    if(rel.sumIndexMap.count(ReflSum::Positive)){
                        posNeg.insert(posNeg.end(),
                            rel.sumIndexMap.at(ReflSum::Positive).begin(),
                            rel.sumIndexMap.at(ReflSum::Positive).end());
                    }
                    if(rel.sumIndexMap.count(ReflSum::Negative)){
                        posNeg.insert(posNeg.end(),
                            rel.sumIndexMap.at(ReflSum::Negative).begin(),
                            rel.sumIndexMap.at(ReflSum::Negative).end());
                    }
                    sumSets.push_back(posNeg);
                    break;
                }
            }

            // Use the sets to build amp x sum combinations
            for(const auto& amps : ampSets){
                for(const auto& sums : sumSets){
                    if(!isNonEmptyConfig(amps, sums, rel)) continue;
                    configs.push_back({aGroup, sGroup, amps, sums});
                }
            }
        }
    }
    return configs;
}

// Rather than reseting the active amps/sums,
// this function checks the active status for them and changes what's needed.
void applyConfig(const PlotConfig& newConfig, const AmplitudeRegistry& reg, 
                  std::unique_ptr<vecps_PlotGen>& plotResults){
    // Initialize a vector with all entries as false and flip to true the
    // indices for the amps that need to be active
    std::vector<bool> newAmps(reg.nUnqAmps,false);
    for(auto a : newConfig.ampIndices) newAmps[a] = true;

    for(size_t i=0;i<reg.nUnqAmps;++i){
        if(plotResults->isAmpEnabled(i) && !newAmps[i]) 
            plotResults->disableAmp(i);
        if(!plotResults->isAmpEnabled(i) && newAmps[i]) 
            plotResults->enableAmp(i);
    }

    // Initialize a vector with all entries as false and flip to true the
    // indices for the sums that need to be active
    std::vector<bool> newSums(reg.nUnqSums,false);
    for(auto s : newConfig.sumIndices) newSums[s] = true;

    for(size_t i=0;i<reg.nUnqSums;++i){
        if(plotResults->isSumEnabled(i) && !newSums[i]) 
            plotResults->disableSum(i);
        if(!plotResults->isSumEnabled(i) && newSums[i]) 
            plotResults->enableSum(i);
    }
}

PlotLabels makeLabel(const PlotConfig& cfg, const AmplitudeRegistry& reg){
    PlotLabels labels;
    std::string amp;
    std::string sum;

    switch(cfg.ampGroup){
        case AmpPlotGroup::AllAmps:
            amp = "All";
            break;

        case AmpPlotGroup::AllReflAmps:
            amp = "AllRefl";
            break;

        case AmpPlotGroup::SingleAmp:{
            const auto& ampKey = reg.uniqueAmps[cfg.ampIndices[0]];
            amp = ampKey.ampName;
            break;
        }

        case AmpPlotGroup::JlFamily:{
            const auto& ampKey = reg.uniqueAmps[cfg.ampIndices[0]];
            amp = ampKey.jFamily;
            break;
        }
    }

    switch(cfg.sumGroup){
        case SumPlotGroup::AllSums:
            sum = "Full";
            break;

        case SumPlotGroup::PosReflSum:
            sum = "posRefl";
            break;

        case SumPlotGroup::NegReflSum:
            sum = "negRefl";
            break;

        case SumPlotGroup::NoReflSum:
            sum = "noRefl";
            break;

        case SumPlotGroup::BothReflSums:
            sum = "BothRefl";
            break;
    }

    // Construct the labels for human ond histogram names.
    // TODO: add csv format
    labels.human = amp + " | " + sum;
    labels.histName  = "_" + amp + "_" + sum;

    return labels;
}

// This function creates the plots that are assigned in VecPsPlotGenerator
void generatePlots(ProgramOptions& options, std::unique_ptr<vecps_PlotGen>& plotResults,
                   const  std::string& reactionName,
                   const TFile& plotFile, bool& singleData,
                   const PlotConfig& cfg, const PlotLabels& cfgLabel){
    // std::cout << "Creating plots for: " << cfgLabel.human << std::endl;
    // Loop over data types
    for (unsigned int iType = 0; iType < PlotGenerator::kNumTypes; ++iType){
        if (!options.bkgFile && iType == PlotGenerator::kBkgnd) continue;

        if (iType == PlotGenerator::kData){
            if(!singleData) singleData = true;
            else if(singleData) continue;  // only plot data once
        }

        std::string dataType;
        switch (iType) {
            case PlotGenerator::kData:    dataType = "dat"; break;
            case PlotGenerator::kAccMC:   dataType = "acc"; break;
            case PlotGenerator::kGenMC:   dataType = "gen"; break;
            case PlotGenerator::kBkgnd:   dataType = "bkgnd"; break;
            default: continue;
        }

        // Loop over all histograms defined in VecPsPlotGenerator
        // The loop variable should be defined using the first entry in the 
        // enum, in this case kVecPsMass (change if needed). The last one is 
        // kNumHists and this one should not change
        for (VecPsPlotGenerator::Hist_index ivar = VecPsPlotGenerator::kVecPsMass; 
            ivar < VecPsPlotGenerator::kNumHists; 
            ivar = static_cast<VecPsPlotGenerator::Hist_index>(ivar + 1)){
            Histogram* hist = plotResults->projection(ivar, reactionName, iType);
            if (!hist) continue; // safety skip
            TH1* histRoot = hist->toRoot();

            // Assign proper names to each histogram. Note that the names 
            // and the translation table are defined in VecPsPlotGenerator 
            std::string histName = VecPsPlotGenerator::numToString(ivar);

            if(!(iType==PlotGenerator::kData))
                histName += "_" + dataType + cfgLabel.histName;
            histRoot->SetName(histName.c_str());
            histRoot->Write();
        }
    }      
                       
}

// This function output a semi-human readable text file with the fit fractions
// and phase differences between amplitudes. It has some patterns to help 
// parsing, but csv file might be preferred for that purpose.
// Note: the default behavior is the acceptance corrected intensities
void writeAmpInfo(const AmplitudeRegistry& reg,
                  std::unique_ptr<FitResults>& results){

    std::ofstream parsOutputFile("params.txt");
    // This is the total number of events for the acceptance correction
    auto total = results->intensity();
    double totalVal = total.first;
    double totalErr = total.second;

    parsOutputFile << "################################################\n";
    parsOutputFile << "#       All Values are acceptance corrected    #\n";
    parsOutputFile << "################################################\n";
    parsOutputFile << "Total Events = " << totalVal << " +- " << totalErr << "\n";
    
    // Coherent sums
    for(size_t iAmp = 0; iAmp < reg.nUnqAmps; ++iAmp){
        const auto& wave = reg.uniqueAmps[iAmp];

        // Coherent sums for unique amp based in refl
        for(auto refl : {ReflSum::Positive, ReflSum::Negative, ReflSum::Neither}){
            auto key = std::make_pair(iAmp, refl);
            if(reg.relations.ampReflToFull.count(key) == 0)
                continue;
            const auto& fullIndices = reg.relations.ampReflToFull.at(key);
            if(fullIndices.empty())
                continue;
            // Build full name list
            std::vector<std::string> fullNames;
            for(auto idx : fullIndices)
                fullNames.push_back(reg.ampsList[idx].fullName);

            auto val = results->intensity(fullNames);

            parsOutputFile << "FIT FRACTION (coherent sum) "
                           << reflToString(refl) << " "
                           << wave.ampName << " = "
                           << val.first / totalVal << " +- "
                           << val.second / totalVal << "\n";
        }
    }

    // Coherent sums for Jl families 
    for(const auto& [key, fullIndices] : reg.relations.jlReflToFull){
        auto [jlFamily, refl] = key;
        if(fullIndices.empty()) continue;

        std::vector<std::string> fullNames;
        for(auto idx : fullIndices)
            fullNames.push_back(reg.ampsList[idx].fullName);

        auto val = results->intensity(fullNames);

        parsOutputFile << "FIT FRACTION (coherent sum) "
                       << reflToString(refl) << " "
                       << jlFamily << " = "
                       << val.first / totalVal << " +- "
                       << val.second / totalVal << "\n";
    }

    // Individual amplitudes
    for(const auto& ampKey : reg.ampsList){

        std::vector<std::string> singleAmp = { ampKey.fullName };

        auto val = results->intensity(singleAmp);

        parsOutputFile << "FIT FRACTION " << ampKey.fullName << " = "
                       << val.first / totalVal << " +- "
                       << val.second / totalVal << "\n";
    }
    
    parsOutputFile << "################################################\n";
    parsOutputFile << "################################################\n";
    parsOutputFile << "#              Phase Differences               #\n";
    parsOutputFile << "################################################\n";
    parsOutputFile << "################################################\n";

    // In general, the phase differences are the same for all 
    // reactions. This is because the current use case is that 
    // multiple "reactions" are actually the same physical reaction
    // but for a different data set like different run periods or 
    // beam pol. angle
    const auto& firstReaction = results->reactionList()[0];

    for(size_t i = 0; i < reg.ampsList.size(); ++i){
        const auto& amp1 = reg.ampsList[i];

        // Only amplitudes in the first reaction to avoid duplicates
        if(amp1.reactionName != firstReaction) continue;
        if(amp1.reflSum == ReflSum::Neither) continue;

        for(size_t j = i+1; j < reg.ampsList.size(); ++j){
            const auto& amp2 = reg.ampsList[j];

            if(amp2.reactionName != firstReaction) continue;
            if(amp2.reflSum == ReflSum::Neither) continue;

            // Only compare amplitudes in the same coherent sum
            if(amp1.sumIndex != amp2.sumIndex) continue;

            auto phaseDiff = results->phaseDiff(amp1.fullName, amp2.fullName);

            parsOutputFile << "PHASE DIFF " 
                           << amp1.fullName << " vs " 
                           << amp2.fullName << " = " 
                           << phaseDiff.first << " +- " 
                           << phaseDiff.second << "\n";
        }
    }
    
    parsOutputFile << "################################################\n";
    parsOutputFile << "################################################\n";
    parsOutputFile << "#                Amp Parameters                #\n";
    parsOutputFile << "################################################\n";
    parsOutputFile << "################################################\n";

    // currently plots all parameters, including production amplitudes
    // a helper function could be added to filter the parameters of interest
    std::vector<std::string> parNameList = results->parNameList();
    for(const auto& parName : parNameList){     
        double parVal = results->parValue(parName);
        double parErr = results->parError(parName);
        parsOutputFile << parName << ": " 
                       << parVal <<  " +- "
                       << parErr << "\n";
    }
    parsOutputFile.close();
} 

void writeAmpInfoCSV(const AmplitudeRegistry& reg,
                     std::unique_ptr<FitResults>& results){
    std::ofstream csvOutputFile("params.csv");
    // This is the total number of events for the acceptance correction
    auto total = results->intensity();
    double totalVal = total.first;
    double totalErr = total.second;

    csvOutputFile << "Type,Name,Value,Error\n";
    csvOutputFile << "Total,Acceptance Corrected Events," << totalVal << "," 
        << totalErr << "\n";

    // Coherent sums
    for(size_t iAmp = 0; iAmp < reg.nUnqAmps; ++iAmp){
        const auto& wave = reg.uniqueAmps[iAmp];

        // Coherent sums for unique amp based in refl
        for(auto refl : {ReflSum::Positive, ReflSum::Negative, ReflSum::Neither}){
            auto key = std::make_pair(iAmp, refl);
            if(reg.relations.ampReflToFull.count(key) == 0) continue;
            const auto& fullIndices = reg.relations.ampReflToFull.at(key);
            if(fullIndices.empty()) continue;

            std::vector<std::string> fullNames;
            for(auto idx : fullIndices)
                fullNames.push_back(reg.ampsList[idx].fullName);

            auto val = results->intensity(fullNames);

            csvOutputFile << "CoherentSum," 
                << reflToString(refl) << " " << wave.ampName << ","
                << val.first / totalVal << ","
                << val.second / totalVal << "\n";
        }
    }

    // Coherent sums for Jl families 
    for(const auto& [key, fullIndices] : reg.relations.jlReflToFull){
        auto [jlFamily, refl] = key;
        if(fullIndices.empty()) continue;

        std::vector<std::string> fullNames;
        for(auto idx : fullIndices)
            fullNames.push_back(reg.ampsList[idx].fullName);

        auto val = results->intensity(fullNames);

        csvOutputFile << "CoherentSum," 
            << reflToString(refl) << " " << jlFamily << ","
            << val.first / totalVal << ","
            << val.second / totalVal << "\n";
    }

    // Individual amplitudes
    for(const auto& ampKey : reg.ampsList){
        std::vector<std::string> singleAmp = { ampKey.fullName };
        auto val = results->intensity(singleAmp);

        csvOutputFile << "IndividualAmp," 
            << ampKey.fullName << ","
            << val.first / totalVal << ","
            << val.second / totalVal << "\n";
    }

    // Phase differences
    const auto& firstReaction = results->reactionList()[0];

    for(size_t i = 0; i < reg.ampsList.size(); ++i){
        const auto& amp1 = reg.ampsList[i];
        if(amp1.reactionName != firstReaction) continue;
        if(amp1.reflSum == ReflSum::Neither) continue;

        for(size_t j = i+1; j < reg.ampsList.size(); ++j){
            const auto& amp2 = reg.ampsList[j];
            if(amp2.reactionName != firstReaction) continue;
            if(amp2.reflSum == ReflSum::Neither) continue;
            if(amp1.sumIndex != amp2.sumIndex) continue;

            auto phaseDiff = results->phaseDiff(amp1.fullName, amp2.fullName);

            csvOutputFile << "PhaseDiff,"
                << amp1.fullName << " vs " << amp2.fullName << ","
                << phaseDiff.first << ","
                << phaseDiff.second << "\n";
        }
    }

    // currently plots all parameters, including production amplitudes
    // a helper function could be added to filter the parameters of interest
    std::vector<std::string> parNameList = results->parNameList();
    for(const auto& parName : parNameList){     
        double parVal = results->parValue(parName);
        double parErr = results->parError(parName);
        csvOutputFile << "Parameter,"
                       << parName << "," 
                       << parVal <<  ","
                       << parErr << "\n";
    }
    csvOutputFile.close();
}

// ...oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//                          Main Functions 
// ...oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

bool parseInputs(int argc, char* argv[], ProgramOptions& options){
    // CLI parsing
    if (argc < 2) {
        std::cerr << "Usage: \n";
        std::cerr << "\tvecps_plotter <results file name>" 
            "-o <output file name>\n"<< std::endl;
        return false;
    }

    for (int i = 1; i < argc; i++){
        std::string arg(argv[i]);
        if (arg == "-f"){
            options.resultsName = std::string(argv[++i]);
        }
        else if (arg == "-g"){
            options.showGui = true;
        }
        else if (arg == "-o"){
            options.outName = std::string(argv[++i]);
        }
        else if (arg == "-mcGen"){
            options.mcGen = true;
        }
        else if (arg == "-bkgFile"){
            options.bkgFile = true;
        }
        else if (arg == "-sFile"){
            options.sFile = true;
        }
        else if (arg == "-skipPlots"){
            options.skipPlots = true;
        }
        else if (arg == "-txt"){
            options.txtOutput = true;
        }
        else if (arg == "-csv"){
            options.csvOutput = true;
        }
        else if (arg == "-h"){
            std::cout << "\n Usage for: " << argv[0] << "\n\n";
            std::cout << "\t -f <file>\t input fit results file path\n";
            std::cout << "\t -o <file>\t output file path\n";
            std::cout << "\t -g \t show GUI. Default off\n";
            std::cout << "\t -mcGen\t output genMC. Default off\n";
            std::cout << "\t -sFile\t use a single reaction file set. Default off\n";
            std::cout << "\t -skipPlots\t skip generating plots. Default off\n";
            std::cout << "\t -txt\t output text file. Default off\n";
            std::cout << "\t -csv\t output CSV file. Default off\n";
            std::cout << "\t -h\t show this help message\n"
                      << std::endl;
            return false;
        }
        else {
            std::cerr << "Error: Unknown option '" << arg << "'\n";
            std::cerr << "Run with -h for help." << std::endl;
            return false;
        }
    }
    
    // report used options   
    std::cout << "Fit results file name    = " << options.resultsName << "\n";
    if(!options.skipPlots)
    std::cout << "Output file name         = " << options.outName << "\n";
    std::cout << "Show GUI?                = " 
              << (options.showGui ? "Yes" : "No") << "\n";
    std::cout << "MC Gen?                  = " 
              << (options.mcGen ? "Yes" : "No") << "\n";
    std::cout << "Single File?             = " 
              << (options.sFile ? "Yes" : "No") << "\n";
    std::cout << "Text Output?             = " 
              << (options.txtOutput ? "Yes" : "No") << "\n";
    std::cout << "CSV Output?              = " 
              << (options.csvOutput ? "Yes" : "No") << "\n"
              << std::endl;
    
    return true;
}

bool initializeAmpTools(ProgramOptions& options, 
                       std::unique_ptr<FitResults>& results,
                       std::unique_ptr<vecps_PlotGen>& plotResults,
                       AmplitudeRegistry& registry){
    // load the results and display the configuration info
    std::cout << "Initializing AmpTools interface...\n";
    atiSetup();
    std::cout << "Interface initialized\n";
    std::cout << "Loading Fit results\n";
    results = std::make_unique<FitResults>(options.resultsName);
    if (!results->valid()){
        std::cerr << "Invalid fit results in file: " 
                  << options.resultsName << std::endl;
        return false;
    }
    std::cout << "Fit results loaded\n";
    // set up the plot generator
    if (options.mcGen){
      std::cout << "Plot results (MC gen on) \n";
      plotResults = 
            std::make_unique<vecps_PlotGen >(*results, PlotGenerator::kDefault);
    } 
    else {
      std::cout << "Plot results (MC gen off) \n";
      plotResults = 
            std::make_unique<vecps_PlotGen >(*results, PlotGenerator::kNoGenMC);
    }
    std::cout << "Initialized ati and PlotGen" << std::endl;

    // Build amplitude registry (once per file)
    registry = buildRegistry(plotResults);

    return true;
}

void analyzeFitFile(ProgramOptions& options,
                    std::unique_ptr<FitResults>& results,
                    std::unique_ptr<vecps_PlotGen>& plotResults,
                     AmplitudeRegistry& registry){    
    if(options.skipPlots){
        std::cout << "Skipping plot generation.\n";
        return;
    }
    cout << "\n*** Viewing Results Using AmpPlotter***\n\n";
    // Loop over reactions
    // User will have to add the different root file outputs manually
    auto reactionList = results->reactionList();
    for (size_t fileSet = 0; fileSet < reactionList.size(); ++fileSet){
        if (options.sFile && fileSet > 0) break;

        std::string reactionName = reactionList[fileSet];
        std::cout << "Processing reaction " << reactionName << " ("
                  << fileSet + 1 << "/" << reactionList.size() << ")\n";

        if (!options.sFile){
            options.outName = reactionName + ".root";
            std::cout << "Output file: " << options.outName << std::endl;
        }

        TFile* plotFile = new TFile(options.outName.c_str(), "RECREATE");
        TH1::AddDirectory(kFALSE); // avoid automatic ownership for GUI use

        plotResults->enableReaction(reactionName);
        // Generate plot configurations
        auto plotConfigs = generatePlotConfigs(registry);
        // Loop over plot configurations
        bool first = true;
        bool singleData = false;
        for (const auto& cfg : plotConfigs){

            if (!first){
                applyConfig(cfg, registry, plotResults);
            } 
            else{
                // On first iteration, just enable the relevant amps/sums
                for (auto i : cfg.ampIndices) plotResults->enableAmp(i);
                for (auto i : cfg.sumIndices) plotResults->enableSum(i);
                first = false;
            }

           PlotLabels cfgLabel = makeLabel(cfg, registry);
           generatePlots(options, plotResults, reactionName,
                   *plotFile, singleData,cfg, cfgLabel);
        }
        plotFile->Close();
    }
}

void startGUI(ProgramOptions& options,
              std::unique_ptr<vecps_PlotGen>& plotResults){
    if(options.skipPlots){
        std::cout << "No plots = No GUI\n";
        return;
    }
    auto plotResults_local = std::move(plotResults);
    // start the GUI
    if (options.showGui){
        std::cout << ">> Plot generator ready, starting GUI..." << std::endl;

        int dummy_argc = 0;
        char* dummy_argv[] = {};
        TApplication app("app", &dummy_argc, dummy_argv);
        app.SetReturnFromRun(true);

        gStyle->SetFillColor(10);
        gStyle->SetCanvasColor(10);
        gStyle->SetPadColor(10);
        gStyle->SetFillStyle(1001);
        gStyle->SetPalette(1);
        gStyle->SetFrameFillColor(10);
        gStyle->SetFrameFillStyle(1001);

        std::cout << " Initialized App " << std::endl;
        auto factory = new PlotFactory(*plotResults_local);
        std::cout << " Created Plot Factory " << std::endl;
        auto mainFrame = new PlotterMainWindow(gClient->GetRoot(), *factory);
        std::cout << " Main frame created " << std::endl;
        mainFrame->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
        std::cout << " App running" << std::endl;
        app.Run();
	    std::cout << "GUI closed cleanly\n";
    }

}

int main(int argc, char* argv[]){
    
    ProgramOptions options;
    if (!parseInputs(argc, argv, options)){
        std::cout << "Error parsing inputs" << std::endl;
        return 1;   // Parsing failed
    }

    std::unique_ptr<FitResults> results;
    std::unique_ptr<vecps_PlotGen> plotResults;
    AmplitudeRegistry registry;
    if (!initializeAmpTools(options, results, plotResults, registry)){
        std::cout << "Error initializing AmpTools" << std::endl;
        return 1;   // Parsing failed
    }
    
    analyzeFitFile(options, results, plotResults, registry);
    
    writeAmpInfo(registry,results);
    
    writeAmpInfoCSV(registry,results);
    
    // GUI needs to be started last
    startGUI(options, plotResults);

    return 0;
}