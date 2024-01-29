////////////////////////////////////////////////////////////////////
/// \file BASED ON: PlotHitTimeResiduals.cc
///
/// \brief Functions to get residual hit time PDFs.
///
/// \author James Page <j.page@sussex.ac.uk>
///
/// REVISION HISTORY:\n
///
/// \details EV Calibrated hit times are plotted minus transit times
/// based on the MC position or the fitted position.
/// Multiple PDFs can be made at the same time, and cuts can be performed
/// for all events, or specific to each PDF.
///
/// To compile: g++ -g -std=c++1y getPDF.cpp -o getPDF.exe `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux
///
////////////////////////////////////////////////////////////////////

#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/DU/LightPathCalculator.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include "RAT/DU/Point3D.hh"
#include <RAT/DS/FitResult.hh>
#include <RAT/DataCleaningUtility.hh>

#include <RAT/TrackNav.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/TrackNode.hh>
#include <RAT/DB.hh>

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TVectorD.h>

#include <string>
#include <fstream>


// #define PRINT_T_RES

bool find_prompt_event(RAT::DU::DSReader& dsReader, const int entry, const int evt, std::vector<double>& E, std::vector<double>& times, std::vector<TVector3>& pos, std::vector<double>& Nhits,
                       double& Delta_T, std::ofstream& t_res_file, RAT::DU::TimeResidualCalculator& fTRCalc, const RAT::DU::ReconCalibrator& e_cal, RAT::DU::DetectorStateCorrection& stateCorr,
                       bool is_data, const bool verbose, const bool use_pos_dep_corr, const bool use_E_corr, const bool use_Nhit, const bool cut_R_min);
bool get_recon_info(std::vector<double>& E, std::vector<double>& times, std::vector<TVector3>& pos, std::vector<double>& Nhits, const unsigned int idx, const RAT::DS::EV& evt, const RAT::DU::ReconCalibrator& e_cal,
                    RAT::DU::DetectorStateCorrection& stateCorr, bool is_data, const bool use_pos_dep_corr, const bool use_E_corr);
void get_t_res(std::ofstream& t_res_file, const RAT::DS::EV& evt, RAT::DU::TimeResidualCalculator& fTRCalc, const TVector3& position, const double vertex_time);
void make_hists(const std::vector<std::string>& fileNames, const std::string output_root_address, const std::string output_txt_address, bool is_data, const bool verbose, const std::string fitName = "");
bool pass_prompt_cuts(const double energy, const double Nhit, TVector3 position, const bool use_Nhit, const bool cut_R_min);
bool pass_delayed_cuts(const double energy, const double Nhit, TVector3 position, const bool use_Nhit, const bool cut_R_min);
bool pass_coincidence_cuts(const double delay, TVector3 prompt_pos, TVector3 delayed_pos);


/* ~~~~~~~~~~~~~~~~~~~~~~ MAIN FUNCTION ~~~~~~~~~~~~~~~~~~~~~ */

int main(int argc, char** argv) {
    std::string output_root_address = argv[1];
    std::string output_txt_address = argv[2];
    bool flat_E_prompt = std::stoi(argv[3]);  // Whether to flatten the prompt energy spectrum or not
    bool is_data = std::stoi(argv[4]);
    bool verbose = std::stoi(argv[5]);
    // Addresses of simulation output files to be analysed
    std::vector<std::string> input_files;
    for (unsigned int i = 6; i < argc; ++i) {
        input_files.push_back(argv[i]);
    }

    // Loop through files to get info from every event (including t_res), and print all to text file
    if (verbose) {std::cout << "Getting info..." << std::endl;}
    make_hists(input_files, output_root_address, output_txt_address, is_data, verbose);

    return 0;
}


/* ~~~~~~~~~~~~~~~~~~~~~~ CUT FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~ */

const double R_MIN = 4000.0, R_MAX = 5700.0;
const double MAX_DIST = 1500.0;
const double MIN_DELAY = 500.0, MAX_DELAY = 0.8E6;

const double MIN_PROMPT_E = 0.7, MAX_PROMPT_E = 3.5;
const double MIN_PROMPT_Nhit = 150.0, MAX_PROMPT_Nhit = 750.0;

const double MIN_DELAYED_E = 1.8, MAX_DELAYED_E = 2.5;
const double MIN_DELAYED_Nhit = 400.0, MAX_DELAYED_Nhit = 650.0;

ULong64_t dcAnalysisWord = 36283883733698;  // Converted hex to decimal from 0x2100000042C2


bool pass_prompt_cuts(const double energy, const double Nhit, TVector3 position, const bool use_Nhit, const bool cut_R_min) {
    if (use_Nhit) {
        if (Nhit < MIN_PROMPT_Nhit) return false;  // min Nhit cut
        if (Nhit > MAX_PROMPT_Nhit) return false;  // max Nhit cut
    } else {
        if (energy < MIN_PROMPT_E) return false;  // min energy cut (MeV)
        if (energy > MAX_PROMPT_E) return false;  // max energy cut (MeV)
    }
    if (position.Mag() > R_MAX) return false;  // FV cut (mm)
    if (cut_R_min && (position.Mag() < R_MIN)) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts(const double energy, const double Nhit, TVector3 position, const bool use_Nhit, const bool cut_R_min) {
    if (use_Nhit) {
        if (Nhit < MIN_DELAYED_Nhit) return false;  // min Nhit cut
        if (Nhit > MAX_DELAYED_Nhit) return false;  // max Nhit cut
    } else {
        if (energy < MIN_DELAYED_E) return false;  // min energy cut (MeV)
        if (energy > MAX_DELAYED_E) return false;  // max energy cut (MeV)
    }
    if (position.Mag() > R_MAX) return false;  // FV cut (mm)
    if (cut_R_min && (position.Mag() < R_MIN)) return false;  // FV cut (mm)

    return true;
}

bool pass_coincidence_cuts(const double delay, TVector3 prompt_pos, TVector3 delayed_pos) {
    // double delay = (delayed_time - prompt_time) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
    double distance = (delayed_pos - prompt_pos).Mag();

    if (delay < MIN_DELAY) return false;  // min delay cut (ns)
    if (delay > MAX_DELAY) return false;  // max delay cut (ns)
    if (distance > MAX_DIST) return false;  // max distance cut (mm)

    return true;
}



/* ~~~~~~~~~~~~~~~~~~~~~~ PRIMARY FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~ */

/**
 * @brief Makes various histograms of reconstrcuted quantities of both delayed and prompt events,
 * as well as the time residuals of prompt events. Only makes these for event-pairs that pass all
 * cuts and tagging.
 * This is done by looping through all events, finding an event that passes delayed cuts, recording its
 * information and looping back to find:
 * - Either a second event is found that pass both the prompt event and tagging cuts -> record events.
 * - Or an event is found whose recon event time is futher from the delayed event than the max delta-t
 * cut -> through events away.
 * - Or the file runs out of events -> throw events away.
 * The loop then picks up again on the event just following the previous prompt event.
 * POTENTIAL ISSUE: double counting an event as delayed in one pair AND prompt in another.
 * 
 * @param fileNames  list of root file addresses to loop through with event info (ratds files)
 * @param output_file_address  output root file to save histograms to
 * @param verbose  flag
 * @param fitName  name of fitter to get recon info from (default is "")
 */
void make_hists(const std::vector<std::string>& fileNames, const std::string output_root_address, const std::string output_txt_address, bool is_data, const bool verbose, const std::string fitName) {
    if (verbose) {std::cout << "Running make_hists()" << std::endl;}

    /*********** Set up print file ***********/
    std::ofstream t_res_file;
    t_res_file.open(output_txt_address);

    /*********** Set histograms ***********/
    const double min_delayed_E = 1.2, max_delayed_E = 3.2, min_prompt_E = 0.3, max_prompt_E = 4.5;
    const double min_delayed_Nhit = 200.0, max_delayed_Nhit = 800.0, min_prompt_Nhit = 100.0, max_prompt_Nhit = 900.0;
    const double min_R = 0.0, max_R = R_MAX, min_deltaT = MIN_DELAY, max_deltaT = MAX_DELAY, min_deltaR = 0.0, max_deltaR = MAX_DIST;

    // nhit cuts, no R cut, no pos dep corr, no E corr
    TH1D hist_prompt_E_NhitCut("prompt_E_NhitCut", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_NhitCut("delayed_E_NhitCut", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_NhitCut("prompt_R_NhitCut", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_NhitCut("delayed_R_NhitCut", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_NhitCut("prompt_Nhits_NhitCut", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_NhitCut("delayed_Nhits_NhitCut", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_NhitCut("deltaR_NhitCut", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_NhitCut("deltaT_NhitCut", "deltaT", 1000, min_deltaT, max_deltaT);

    // E cuts, no R cut, no pos dep corr, no E corr
    TH1D hist_prompt_E_ECut("prompt_E_ECut", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_ECut("delayed_E_ECut", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_ECut("prompt_R_ECut", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_ECut("delayed_R_ECut", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_ECut("prompt_Nhits_ECut", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_ECut("delayed_Nhits_ECut", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_ECut("deltaR_ECut", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_ECut("deltaT_ECut", "deltaT", 1000, min_deltaT, max_deltaT);

    // nhit cuts, R cut, no pos dep corr, no E corr
    TH1D hist_prompt_E_NhitCut_Rcut("prompt_E_NhitCut-Rcut", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_NhitCut_Rcut("delayed_E_NhitCut-Rcut", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_NhitCut_Rcut("prompt_R_NhitCut-Rcut", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_NhitCut_Rcut("delayed_R_NhitCut-Rcut", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_NhitCut_Rcut("prompt_Nhits_NhitCut-Rcut", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_NhitCut_Rcut("delayed_Nhits_NhitCut-Rcut", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_NhitCut_Rcut("deltaR_NhitCut-Rcut", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_NhitCut_Rcut("deltaT_NhitCut-Rcut", "deltaT", 1000, min_deltaT, max_deltaT);

    // E cuts, R cut, no pos dep corr, no E corr
    TH1D hist_prompt_E_ECut_Rcut("prompt_E_ECut-Rcut", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_ECut_Rcut("delayed_E_ECut-Rcut", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_ECut_Rcut("prompt_R_ECut-Rcut", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_ECut_Rcut("delayed_R_ECut-Rcut", "delayed_R", 1000, min_R, max_R);
    TH1D hist_prompt_Nhits_ECut_Rcut("prompt_Nhits_ECut-Rcut", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_ECut_Rcut("delayed_Nhits_ECut-Rcut", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_ECut_Rcut("deltaR_ECut-Rcut", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_ECut_Rcut("deltaT_ECut-Rcut", "deltaT", 1000, min_deltaT, max_deltaT);

    /* ---- */

    // nhit cuts, no R cut, pos dep corr, no E corr
    TH1D hist_prompt_E_NhitCut_PosDepCorr("prompt_E_NhitCut-PosDepCorr", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_NhitCut_PosDepCorr("delayed_E_NhitCut-PosDepCorr", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_NhitCut_PosDepCorr("prompt_R_NhitCut-PosDepCorr", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_NhitCut_PosDepCorr("delayed_R_NhitCut-PosDepCorr", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_NhitCut_PosDepCorr("prompt_Nhits_NhitCut-PosDepCorr", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_NhitCut_PosDepCorr("delayed_Nhits_NhitCut-PosDepCorr", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_NhitCut_PosDepCorr("deltaR_NhitCut-PosDepCorr", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_NhitCut_PosDepCorr("deltaT_NhitCut-PosDepCorr", "deltaT", 1000, min_deltaT, max_deltaT);

    // E cuts, no R cut, pos dep corr, no E corr
    TH1D hist_prompt_E_ECut_PosDepCorr("prompt_E_ECut-PosDepCorr", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_ECut_PosDepCorr("delayed_E_ECut-PosDepCorr", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_ECut_PosDepCorr("prompt_R_ECut-PosDepCorr", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_ECut_PosDepCorr("delayed_R_ECut-PosDepCorr", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_ECut_PosDepCorr("prompt_Nhits_ECut-PosDepCorr", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_ECut_PosDepCorr("delayed_Nhits_ECut-PosDepCorr", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_ECut_PosDepCorr("deltaR_ECut-PosDepCorr", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_ECut_PosDepCorr("deltaT_ECut-PosDepCorr", "deltaT", 1000, min_deltaT, max_deltaT);

    // nhit cuts, R cut, pos dep corr, no E corr
    TH1D hist_prompt_E_NhitCut_Rcut_PosDepCorr("prompt_E_NhitCut-Rcut-PosDepCorr", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_NhitCut_Rcut_PosDepCorr("delayed_E_NhitCut-Rcut-PosDepCorr", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_NhitCut_Rcut_PosDepCorr("prompt_R_NhitCut-Rcut-PosDepCorr", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_NhitCut_Rcut_PosDepCorr("delayed_R_NhitCut-Rcut-PosDepCorr", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_NhitCut_Rcut_PosDepCorr("prompt_Nhits_NhitCut-Rcut-PosDepCorr", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_NhitCut_Rcut_PosDepCorr("delayed_Nhits_NhitCut-Rcut-PosDepCorr", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_NhitCut_Rcut_PosDepCorr("deltaR_NhitCut-Rcut-PosDepCorr", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_NhitCut_Rcut_PosDepCorr("deltaT_NhitCut-Rcut-PosDepCorr", "deltaT", 1000, min_deltaT, max_deltaT);

    // E cuts, R cut, pos dep corr, no E corr
    TH1D hist_prompt_E_ECut_Rcut_PosDepCorr("prompt_E_ECut-Rcut-PosDepCorr", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_ECut_Rcut_PosDepCorr("delayed_E_ECut-Rcut-PosDepCorr", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_ECut_Rcut_PosDepCorr("prompt_R_ECut-Rcut-PosDepCorr", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_ECut_Rcut_PosDepCorr("delayed_R_ECut-Rcut-PosDepCorr", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_ECut_Rcut_PosDepCorr("prompt_Nhits_ECut-Rcut-PosDepCorr", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_ECut_Rcut_PosDepCorr("delayed_Nhits_ECut-Rcut-PosDepCorr", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_ECut_Rcut_PosDepCorr("deltaR_ECut-Rcut-PosDepCorr", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_ECut_Rcut_PosDepCorr("deltaT_ECut-Rcut-PosDepCorr", "deltaT", 1000, min_deltaT, max_deltaT);

    /* ±±±±±±±±±±±±±±±±±±±±±±±± */

    // nhit cuts, no R cut, no pos dep corr, E corr
    TH1D hist_prompt_E_NhitCut_Ecorr("prompt_E_NhitCut-Ecorr", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_NhitCut_Ecorr("delayed_E_NhitCut-Ecorr", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_NhitCut_Ecorr("prompt_R_NhitCut-Ecorr", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_NhitCut_Ecorr("delayed_R_NhitCut-Ecorr", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_NhitCut_Ecorr("prompt_Nhits_NhitCut-Ecorr", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_NhitCut_Ecorr("delayed_Nhits_NhitCut-Ecorr", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_NhitCut_Ecorr("deltaR_NhitCut-Ecorr", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_NhitCut_Ecorr("deltaT_NhitCut-Ecorr", "deltaT", 1000, min_deltaT, max_deltaT);

    // E cuts, no R cut, no pos dep corr, E corr
    TH1D hist_prompt_E_ECut_Ecorr("prompt_E_ECut-Ecorr", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_ECut_Ecorr("delayed_E_ECut-Ecorr", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_ECut_Ecorr("prompt_R_ECut-Ecorr", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_ECut_Ecorr("delayed_R_ECut-Ecorr", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_ECut_Ecorr("prompt_Nhits_ECut-Ecorr", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_ECut_Ecorr("delayed_Nhits_ECut-Ecorr", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_ECut_Ecorr("deltaR_ECut-Ecorr", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_ECut_Ecorr("deltaT_ECut-Ecorr", "deltaT", 1000, min_deltaT, max_deltaT);

    // nhit cuts, R cut, no pos dep corr, E corr
    TH1D hist_prompt_E_NhitCut_Rcut_Ecorr("prompt_E_NhitCut-Rcut-Ecorr", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_NhitCut_Rcut_Ecorr("delayed_E_NhitCut-Rcut-Ecorr", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_NhitCut_Rcut_Ecorr("prompt_R_NhitCut-Rcut-Ecorr", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_NhitCut_Rcut_Ecorr("delayed_R_NhitCut-Rcut-Ecorr", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_NhitCut_Rcut_Ecorr("prompt_Nhits_NhitCut-Rcut-Ecorr", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_NhitCut_Rcut_Ecorr("delayed_Nhits_NhitCut-Rcut-Ecorr", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_NhitCut_Rcut_Ecorr("deltaR_NhitCut-Rcut-Ecorr", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_NhitCut_Rcut_Ecorr("deltaT_NhitCut-Rcut-Ecorr", "deltaT", 1000, min_deltaT, max_deltaT);

    // E cuts, R cut, no pos dep corr, E corr
    TH1D hist_prompt_E_ECut_Rcut_Ecorr("prompt_E_ECut-Rcut-Ecorr", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_ECut_Rcut_Ecorr("delayed_E_ECut-Rcut-Ecorr", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_ECut_Rcut_Ecorr("prompt_R_ECut-Rcut-Ecorr", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_ECut_Rcut_Ecorr("delayed_R_ECut-Rcut-Ecorr", "delayed_R", 1000, min_R, max_R);
    TH1D hist_prompt_Nhits_ECut_Rcut_Ecorr("prompt_Nhits_ECut-Rcut-Ecorr", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_ECut_Rcut_Ecorr("delayed_Nhits_ECut-Rcut-Ecorr", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_ECut_Rcut_Ecorr("deltaR_ECut-Rcut-Ecorr", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_ECut_Rcut_Ecorr("deltaT_ECut-Rcut-Ecorr", "deltaT", 1000, min_deltaT, max_deltaT);

    /* ---- */

    // nhit cuts, no R cut, pos dep corr, E corr
    TH1D hist_prompt_E_NhitCut_PosDepCorr_Ecorr("prompt_E_NhitCut-PosDepCorr-Ecorr", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_NhitCut_PosDepCorr_Ecorr("delayed_E_NhitCut-PosDepCorr-Ecorr", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_NhitCut_PosDepCorr_Ecorr("prompt_R_NhitCut-PosDepCorr-Ecorr", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_NhitCut_PosDepCorr_Ecorr("delayed_R_NhitCut-PosDepCorr-Ecorr", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_NhitCut_PosDepCorr_Ecorr("prompt_Nhits_NhitCut-PosDepCorr-Ecorr", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_NhitCut_PosDepCorr_Ecorr("delayed_Nhits_NhitCut-PosDepCorr-Ecorr", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_NhitCut_PosDepCorr_Ecorr("deltaR_NhitCut-PosDepCorr-Ecorr", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_NhitCut_PosDepCorr_Ecorr("deltaT_NhitCut-PosDepCorr-Ecorr", "deltaT", 1000, min_deltaT, max_deltaT);

    // E cuts, no R cut, pos dep corr, E corr
    TH1D hist_prompt_E_ECut_PosDepCorr_Ecorr("prompt_E_ECut-PosDepCorr-Ecorr", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_ECut_PosDepCorr_Ecorr("delayed_E_ECut-PosDepCorr-Ecorr", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_ECut_PosDepCorr_Ecorr("prompt_R_ECut-PosDepCorr-Ecorr", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_ECut_PosDepCorr_Ecorr("delayed_R_ECut-PosDepCorr-Ecorr", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_ECut_PosDepCorr_Ecorr("prompt_Nhits_ECut-PosDepCorr-Ecorr", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_ECut_PosDepCorr_Ecorr("delayed_Nhits_ECut-PosDepCorr-Ecorr", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_ECut_PosDepCorr_Ecorr("deltaR_ECut-PosDepCorr-Ecorr", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_ECut_PosDepCorr_Ecorr("deltaT_ECut-PosDepCorr-Ecorr", "deltaT", 1000, min_deltaT, max_deltaT);

    // nhit cuts, R cut, pos dep corr, E corr
    TH1D hist_prompt_E_NhitCut_Rcut_PosDepCorr_Ecorr("prompt_E_NhitCut-Rcut-PosDepCorr-Ecorr", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_NhitCut_Rcut_PosDepCorr_Ecorr("delayed_E_NhitCut-Rcut-PosDepCorr-Ecorr", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_NhitCut_Rcut_PosDepCorr_Ecorr("prompt_R_NhitCut-Rcut-PosDepCorr-Ecorr", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_NhitCut_Rcut_PosDepCorr_Ecorr("delayed_R_NhitCut-Rcut-PosDepCorr-Ecorr", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_NhitCut_Rcut_PosDepCorr_Ecorr("prompt_Nhits_NhitCut-Rcut-PosDepCorr-Ecorr", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_NhitCut_Rcut_PosDepCorr_Ecorr("delayed_Nhits_NhitCut-Rcut-PosDepCorr-Ecorr", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_NhitCut_Rcut_PosDepCorr_Ecorr("deltaR_NhitCut-Rcut-PosDepCorr-Ecorr", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_NhitCut_Rcut_PosDepCorr_Ecorr("deltaT_NhitCut-Rcut-PosDepCorr-Ecorr", "deltaT", 1000, min_deltaT, max_deltaT);

    // E cuts, R cut, pos dep corr, E corr
    TH1D hist_prompt_E_ECut_Rcut_PosDepCorr_Ecorr("prompt_E_ECut-Rcut-PosDepCorr-Ecorr", "prompt_E", 100, min_prompt_E, max_prompt_E);
    TH1D hist_delayed_E_ECut_Rcut_PosDepCorr_Ecorr("delayed_E_ECut-Rcut-PosDepCorr-Ecorr", "delayed_E", 100, min_delayed_E, max_delayed_E);
    TH1D hist_prompt_R_ECut_Rcut_PosDepCorr_Ecorr("prompt_R_ECut-Rcut-PosDepCorr-Ecorr", "prompt_R", 100, min_R, max_R);
    TH1D hist_delayed_R_ECut_Rcut_PosDepCorr_Ecorr("delayed_R_ECut-Rcut-PosDepCorr-Ecorr", "delayed_R", 100, min_R, max_R);
    TH1D hist_prompt_Nhits_ECut_Rcut_PosDepCorr_Ecorr("prompt_Nhits_ECut-Rcut-PosDepCorr-Ecorr", "prompt_Nhits", 100, min_prompt_Nhit, max_prompt_Nhit);
    TH1D hist_delayed_Nhits_ECut_Rcut_PosDepCorr_Ecorr("delayed_Nhits_ECut-Rcut-PosDepCorr-Ecorr", "delayed_Nhits", 100, min_delayed_Nhit, max_delayed_Nhit);
    TH1D hist_deltaR_ECut_Rcut_PosDepCorr_Ecorr("deltaR_ECut-Rcut-PosDepCorr-Ecorr", "deltaR", 100, min_deltaR, max_deltaR);
    TH1D hist_deltaT_ECut_Rcut_PosDepCorr_Ecorr("deltaT_ECut-Rcut-PosDepCorr-Ecorr", "deltaT", 1000, min_deltaT, max_deltaT);


    // RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();
    // ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );

    /*********** Loop through all files, entries, events, and PMTs ***********/
    
    // Set up empty vectors for {prompt, delayed} recon quatities
    std::vector<double> E = {0.0, 0.0};
    std::vector<double> times = {0.0, 0.0};
    std::vector<double> Nhits = {0.0, 0.0};
    std::vector<TVector3> pos = {TVector3(0.0, 0.0, 0.0), TVector3(0.0, 0.0, 0.0)};
    double Delta_T = 0.0;

    // loops through files
    for (unsigned int i = 0; i < fileNames.size(); ++i) {
        std::cout << "Reading in file: " << fileNames.at(i) << std::endl;
        RAT::DU::DSReader dsReader(fileNames.at(i));
        fTRCalc.BeginOfRun();  // Re-initialize time residual calculator (light-path calculator) after it gets geo info from DSReader

        // Initialise DetectorStateCorrection (assume only one run in each file)
        RAT::DU::Utility::Get()->BeginOfRun();
        RAT::DU::DetectorStateCorrection stateCorr = RAT::DU::Utility::Get()->GetDetectorStateCorrection();
        RAT::DU::ReconCalibrator e_cal = RAT::DU::Utility::Get()->GetReconCalibrator();
        
        if (verbose) {std::cout << "Looping through entries..." << std::endl;}
        // loops through entries
        for (unsigned int iEntry = 0; iEntry < dsReader.GetEntryCount(); ++iEntry) {
            // if (verbose) std::cout << "iEntry = " << iEntry << std::endl;
            RAT::DS::Entry rDS = dsReader.GetEntry(iEntry);

            // loops through events
            for (unsigned int iEV = 0; iEV < rDS.GetEVCount(); ++iEV) {
                if (verbose) std::cout << "iEV = " << iEV << std::endl;
                RAT::DS::EV rEV = rDS.GetEV(iEV);

                if (is_data && !RAT::EventIsClean(rEV, dcAnalysisWord)) continue;

                // Get recon info, if it exists

                if (get_recon_info(E, times, pos, Nhits, 1, rEV, e_cal, stateCorr, is_data, false, false)) {
                    // nhit cuts, no R cut, no pos dep corr, no E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], true, false)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, false, false, true, false)) {
                            hist_prompt_E_NhitCut.Fill(E[0]);
                            hist_delayed_E_NhitCut.Fill(E[1]);
                            hist_prompt_R_NhitCut.Fill(pos[0].Mag());
                            hist_delayed_R_NhitCut.Fill(pos[1].Mag());
                            hist_prompt_Nhits_NhitCut.Fill(Nhits[0]);
                            hist_delayed_Nhits_NhitCut.Fill(Nhits[1]);
                            hist_deltaR_NhitCut.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_NhitCut.Fill(Delta_T);

                            if ((pos[0].Mag() < 4000.0) && (pos[1].Mag() < 4000.0)) {
                                std::cout << "Inner AV event from file " << fileNames.at(i) << ", delayed entry # " << iEntry << " and delayed event # " << iEV << "." << std::endl;
                            }
                        }
                    }

                    // E cuts, no R cut, no pos dep corr, no E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], false, false)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, false, false, false, false)) {
                            hist_prompt_E_ECut.Fill(E[0]);
                            hist_delayed_E_ECut.Fill(E[1]);
                            hist_prompt_R_ECut.Fill(pos[0].Mag());
                            hist_delayed_R_ECut.Fill(pos[1].Mag());
                            hist_prompt_Nhits_ECut.Fill(Nhits[0]);
                            hist_delayed_Nhits_ECut.Fill(Nhits[1]);
                            hist_deltaR_ECut.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_ECut.Fill(Delta_T);
                        }
                    }

                    // nhit cuts, R cut, no pos dep corr, no E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], true, true)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, false, false, true, true)) {
                            hist_prompt_E_NhitCut_Rcut.Fill(E[0]);
                            hist_delayed_E_NhitCut_Rcut.Fill(E[1]);
                            hist_prompt_R_NhitCut_Rcut.Fill(pos[0].Mag());
                            hist_delayed_R_NhitCut_Rcut.Fill(pos[1].Mag());
                            hist_prompt_Nhits_NhitCut_Rcut.Fill(Nhits[0]);
                            hist_delayed_Nhits_NhitCut_Rcut.Fill(Nhits[1]);
                            hist_deltaR_NhitCut_Rcut.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_NhitCut_Rcut.Fill(Delta_T);
                        }
                    }

                    // E cuts, R cut, no pos dep corr, no E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], false, true)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, false, false, false, true)) {
                            hist_prompt_E_ECut_Rcut.Fill(E[0]);
                            hist_delayed_E_ECut_Rcut.Fill(E[1]);
                            hist_prompt_R_ECut_Rcut.Fill(pos[0].Mag());
                            hist_delayed_R_ECut_Rcut.Fill(pos[1].Mag());
                            hist_prompt_Nhits_ECut_Rcut.Fill(Nhits[0]);
                            hist_delayed_Nhits_ECut_Rcut.Fill(Nhits[1]);
                            hist_deltaR_ECut_Rcut.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_ECut_Rcut.Fill(Delta_T);
                        }
                    }
                }

                /* ---- */

                if (get_recon_info(E, times, pos, Nhits, 1, rEV, e_cal, stateCorr, is_data, true, false)) {
                    // nhit cuts, no R cut, pos dep corr, no E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], true, false)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, true, false, true, false)) {
                            hist_prompt_E_NhitCut_PosDepCorr.Fill(E[0]);
                            hist_delayed_E_NhitCut_PosDepCorr.Fill(E[1]);
                            hist_prompt_R_NhitCut_PosDepCorr.Fill(pos[0].Mag());
                            hist_delayed_R_NhitCut_PosDepCorr.Fill(pos[1].Mag());
                            hist_prompt_Nhits_NhitCut_PosDepCorr.Fill(Nhits[0]);
                            hist_delayed_Nhits_NhitCut_PosDepCorr.Fill(Nhits[1]);
                            hist_deltaR_NhitCut_PosDepCorr.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_NhitCut_PosDepCorr.Fill(Delta_T);
                        }
                    }

                    // E cuts, no R cut, pos dep corr, no E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], false, false)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, true, false, false, false)) {
                            hist_prompt_E_ECut_PosDepCorr.Fill(E[0]);
                            hist_delayed_E_ECut_PosDepCorr.Fill(E[1]);
                            hist_prompt_R_ECut_PosDepCorr.Fill(pos[0].Mag());
                            hist_delayed_R_ECut_PosDepCorr.Fill(pos[1].Mag());
                            hist_prompt_Nhits_ECut_PosDepCorr.Fill(Nhits[0]);
                            hist_delayed_Nhits_ECut_PosDepCorr.Fill(Nhits[1]);
                            hist_deltaR_ECut_PosDepCorr.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_ECut_PosDepCorr.Fill(Delta_T);
                        }
                    }

                    // nhit cuts, R cut, pos dep corr, no E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], true, true)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, true, false, true, true)) {
                            hist_prompt_E_NhitCut_Rcut_PosDepCorr.Fill(E[0]);
                            hist_delayed_E_NhitCut_Rcut_PosDepCorr.Fill(E[1]);
                            hist_prompt_R_NhitCut_Rcut_PosDepCorr.Fill(pos[0].Mag());
                            hist_delayed_R_NhitCut_Rcut_PosDepCorr.Fill(pos[1].Mag());
                            hist_prompt_Nhits_NhitCut_Rcut_PosDepCorr.Fill(Nhits[0]);
                            hist_delayed_Nhits_NhitCut_Rcut_PosDepCorr.Fill(Nhits[1]);
                            hist_deltaR_NhitCut_Rcut_PosDepCorr.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_NhitCut_Rcut_PosDepCorr.Fill(Delta_T);
                        }
                    }

                    // E cuts, R cut, pos dep corr, no E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], false, true)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, true, false, false, true)) {
                            hist_prompt_E_ECut_Rcut_PosDepCorr.Fill(E[0]);
                            hist_delayed_E_ECut_Rcut_PosDepCorr.Fill(E[1]);
                            hist_prompt_R_ECut_Rcut_PosDepCorr.Fill(pos[0].Mag());
                            hist_delayed_R_ECut_Rcut_PosDepCorr.Fill(pos[1].Mag());
                            hist_prompt_Nhits_ECut_Rcut_PosDepCorr.Fill(Nhits[0]);
                            hist_delayed_Nhits_ECut_Rcut_PosDepCorr.Fill(Nhits[1]);
                            hist_deltaR_ECut_Rcut_PosDepCorr.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_ECut_Rcut_PosDepCorr.Fill(Delta_T);
                        }
                    }
                    
                }

                /* ±±±±±±±±±±±±±±±±±±±±±±±± */

                if (get_recon_info(E, times, pos, Nhits, 1, rEV, e_cal, stateCorr, is_data, false, true)) {
                    // nhit cuts, no R cut, no pos dep corr, E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], true, false)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, false, true, true, false)) {
                            hist_prompt_E_NhitCut_Ecorr.Fill(E[0]);
                            hist_delayed_E_NhitCut_Ecorr.Fill(E[1]);
                            hist_prompt_R_NhitCut_Ecorr.Fill(pos[0].Mag());
                            hist_delayed_R_NhitCut_Ecorr.Fill(pos[1].Mag());
                            hist_prompt_Nhits_NhitCut_Ecorr.Fill(Nhits[0]);
                            hist_delayed_Nhits_NhitCut_Ecorr.Fill(Nhits[1]);
                            hist_deltaR_NhitCut_Ecorr.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_NhitCut_Ecorr.Fill(Delta_T);

                            if ((pos[0].Mag() < 4000.0) && (pos[1].Mag() < 4000.0)) {
                                std::cout << "Inner AV event from file " << fileNames.at(i) << ", delayed entry # " << iEntry << " and delayed event # " << iEV << "." << std::endl;
                            }
                        }
                    }

                    // E cuts, no R cut, no pos dep corr, E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], false, false)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, false, true, false, false)) {
                            hist_prompt_E_ECut_Ecorr.Fill(E[0]);
                            hist_delayed_E_ECut_Ecorr.Fill(E[1]);
                            hist_prompt_R_ECut_Ecorr.Fill(pos[0].Mag());
                            hist_delayed_R_ECut_Ecorr.Fill(pos[1].Mag());
                            hist_prompt_Nhits_ECut_Ecorr.Fill(Nhits[0]);
                            hist_delayed_Nhits_ECut_Ecorr.Fill(Nhits[1]);
                            hist_deltaR_ECut_Ecorr.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_ECut_Ecorr.Fill(Delta_T);
                        }
                    }

                    // nhit cuts, R cut, no pos dep corr, E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], true, true)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, false, true, true, true)) {
                            hist_prompt_E_NhitCut_Rcut_Ecorr.Fill(E[0]);
                            hist_delayed_E_NhitCut_Rcut_Ecorr.Fill(E[1]);
                            hist_prompt_R_NhitCut_Rcut_Ecorr.Fill(pos[0].Mag());
                            hist_delayed_R_NhitCut_Rcut_Ecorr.Fill(pos[1].Mag());
                            hist_prompt_Nhits_NhitCut_Rcut_Ecorr.Fill(Nhits[0]);
                            hist_delayed_Nhits_NhitCut_Rcut_Ecorr.Fill(Nhits[1]);
                            hist_deltaR_NhitCut_Rcut_Ecorr.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_NhitCut_Rcut_Ecorr.Fill(Delta_T);
                        }
                    }

                    // E cuts, R cut, no pos dep corr, E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], false, true)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, false, true, false, true)) {
                            hist_prompt_E_ECut_Rcut_Ecorr.Fill(E[0]);
                            hist_delayed_E_ECut_Rcut_Ecorr.Fill(E[1]);
                            hist_prompt_R_ECut_Rcut_Ecorr.Fill(pos[0].Mag());
                            hist_delayed_R_ECut_Rcut_Ecorr.Fill(pos[1].Mag());
                            hist_prompt_Nhits_ECut_Rcut_Ecorr.Fill(Nhits[0]);
                            hist_delayed_Nhits_ECut_Rcut_Ecorr.Fill(Nhits[1]);
                            hist_deltaR_ECut_Rcut_Ecorr.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_ECut_Rcut_Ecorr.Fill(Delta_T);
                        }
                    }
                }

                /* ---- */

                if (get_recon_info(E, times, pos, Nhits, 1, rEV, e_cal, stateCorr, is_data, true, true)) {
                    // nhit cuts, no R cut, pos dep corr, E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], true, false)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, true, true, true, false)) {
                            hist_prompt_E_NhitCut_PosDepCorr_Ecorr.Fill(E[0]);
                            hist_delayed_E_NhitCut_PosDepCorr_Ecorr.Fill(E[1]);
                            hist_prompt_R_NhitCut_PosDepCorr_Ecorr.Fill(pos[0].Mag());
                            hist_delayed_R_NhitCut_PosDepCorr_Ecorr.Fill(pos[1].Mag());
                            hist_prompt_Nhits_NhitCut_PosDepCorr_Ecorr.Fill(Nhits[0]);
                            hist_delayed_Nhits_NhitCut_PosDepCorr_Ecorr.Fill(Nhits[1]);
                            hist_deltaR_NhitCut_PosDepCorr_Ecorr.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_NhitCut_PosDepCorr_Ecorr.Fill(Delta_T);
                        }
                    }

                    // E cuts, no R cut, pos dep corr, E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], false, false)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, true, true, false, false)) {
                            hist_prompt_E_ECut_PosDepCorr_Ecorr.Fill(E[0]);
                            hist_delayed_E_ECut_PosDepCorr_Ecorr.Fill(E[1]);
                            hist_prompt_R_ECut_PosDepCorr_Ecorr.Fill(pos[0].Mag());
                            hist_delayed_R_ECut_PosDepCorr_Ecorr.Fill(pos[1].Mag());
                            hist_prompt_Nhits_ECut_PosDepCorr_Ecorr.Fill(Nhits[0]);
                            hist_delayed_Nhits_ECut_PosDepCorr_Ecorr.Fill(Nhits[1]);
                            hist_deltaR_ECut_PosDepCorr_Ecorr.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_ECut_PosDepCorr_Ecorr.Fill(Delta_T);
                        }
                    }

                    // nhit cuts, R cut, pos dep corr, E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], true, true)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, true, true, true, true)) {
                            hist_prompt_E_NhitCut_Rcut_PosDepCorr_Ecorr.Fill(E[0]);
                            hist_delayed_E_NhitCut_Rcut_PosDepCorr_Ecorr.Fill(E[1]);
                            hist_prompt_R_NhitCut_Rcut_PosDepCorr_Ecorr.Fill(pos[0].Mag());
                            hist_delayed_R_NhitCut_Rcut_PosDepCorr_Ecorr.Fill(pos[1].Mag());
                            hist_prompt_Nhits_NhitCut_Rcut_PosDepCorr_Ecorr.Fill(Nhits[0]);
                            hist_delayed_Nhits_NhitCut_Rcut_PosDepCorr_Ecorr.Fill(Nhits[1]);
                            hist_deltaR_NhitCut_Rcut_PosDepCorr_Ecorr.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_NhitCut_Rcut_PosDepCorr_Ecorr.Fill(Delta_T);
                        }
                    }

                    // E cuts, R cut, pos dep corr, E corr
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1], false, true)) {
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, true, true, false, true)) {
                            hist_prompt_E_ECut_Rcut_PosDepCorr_Ecorr.Fill(E[0]);
                            hist_delayed_E_ECut_Rcut_PosDepCorr_Ecorr.Fill(E[1]);
                            hist_prompt_R_ECut_Rcut_PosDepCorr_Ecorr.Fill(pos[0].Mag());
                            hist_delayed_R_ECut_Rcut_PosDepCorr_Ecorr.Fill(pos[1].Mag());
                            hist_prompt_Nhits_ECut_Rcut_PosDepCorr_Ecorr.Fill(Nhits[0]);
                            hist_delayed_Nhits_ECut_Rcut_PosDepCorr_Ecorr.Fill(Nhits[1]);
                            hist_deltaR_ECut_Rcut_PosDepCorr_Ecorr.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT_ECut_Rcut_PosDepCorr_Ecorr.Fill(Delta_T);
                        }
                    }
                    
                }
            }
        }
        dsReader.Delete();
    }
    t_res_file.close();

    /*********** Open output root file, and write histograms ***********/
    TFile rootfile(output_root_address.c_str(), "RECREATE");

    //now write everything
    rootfile.cd();

    // nhit cuts, no R cut, no pos dep corr
    hist_prompt_E_NhitCut.Write();
    hist_delayed_E_NhitCut.Write();
    hist_prompt_R_NhitCut.Write();
    hist_delayed_R_NhitCut.Write();
    hist_prompt_Nhits_NhitCut.Write();
    hist_delayed_Nhits_NhitCut.Write();
    hist_deltaR_NhitCut.Write();
    hist_deltaT_NhitCut.Write();

    // E cuts, no R cut, no pos dep corr
    hist_prompt_E_ECut.Write();
    hist_delayed_E_ECut.Write();
    hist_prompt_R_ECut.Write();
    hist_delayed_R_ECut.Write();
    hist_prompt_Nhits_ECut.Write();
    hist_delayed_Nhits_ECut.Write();
    hist_deltaR_ECut.Write();
    hist_deltaT_ECut.Write();

    // nhit cuts, R cut, no pos dep corr
    hist_prompt_E_NhitCut_Rcut.Write();
    hist_delayed_E_NhitCut_Rcut.Write();
    hist_prompt_R_NhitCut_Rcut.Write();
    hist_delayed_R_NhitCut_Rcut.Write();
    hist_prompt_Nhits_NhitCut_Rcut.Write();
    hist_delayed_Nhits_NhitCut_Rcut.Write();
    hist_deltaR_NhitCut_Rcut.Write();
    hist_deltaT_NhitCut_Rcut.Write();

    // E cuts, R cut, no pos dep corr
    hist_prompt_E_ECut_Rcut.Write();
    hist_delayed_E_ECut_Rcut.Write();
    hist_prompt_R_ECut_Rcut.Write();
    hist_delayed_R_ECut_Rcut.Write();
    hist_prompt_Nhits_ECut_Rcut.Write();
    hist_delayed_Nhits_ECut_Rcut.Write();
    hist_deltaR_ECut_Rcut.Write();
    hist_deltaT_ECut_Rcut.Write();

    /* ---- */

    // nhit cuts, no R cut, pos dep corr
    hist_prompt_E_NhitCut_PosDepCorr.Write();
    hist_delayed_E_NhitCut_PosDepCorr.Write();
    hist_prompt_R_NhitCut_PosDepCorr.Write();
    hist_delayed_R_NhitCut_PosDepCorr.Write();
    hist_prompt_Nhits_NhitCut_PosDepCorr.Write();
    hist_delayed_Nhits_NhitCut_PosDepCorr.Write();
    hist_deltaR_NhitCut_PosDepCorr.Write();
    hist_deltaT_NhitCut_PosDepCorr.Write();

    // E cuts, no R cut, pos dep corr
    hist_prompt_E_ECut_PosDepCorr.Write();
    hist_delayed_E_ECut_PosDepCorr.Write();
    hist_prompt_R_ECut_PosDepCorr.Write();
    hist_delayed_R_ECut_PosDepCorr.Write();
    hist_prompt_Nhits_ECut_PosDepCorr.Write();
    hist_delayed_Nhits_ECut_PosDepCorr.Write();
    hist_deltaR_ECut_PosDepCorr.Write();
    hist_deltaT_ECut_PosDepCorr.Write();

    // nhit cuts, R cut, pos dep corr
    hist_prompt_E_NhitCut_Rcut_PosDepCorr.Write();
    hist_delayed_E_NhitCut_Rcut_PosDepCorr.Write();
    hist_prompt_R_NhitCut_Rcut_PosDepCorr.Write();
    hist_delayed_R_NhitCut_Rcut_PosDepCorr.Write();
    hist_prompt_Nhits_NhitCut_Rcut_PosDepCorr.Write();
    hist_delayed_Nhits_NhitCut_Rcut_PosDepCorr.Write();
    hist_deltaR_NhitCut_Rcut_PosDepCorr.Write();
    hist_deltaT_NhitCut_Rcut_PosDepCorr.Write();

    // E cuts, R cut, pos dep corr
    hist_prompt_E_ECut_Rcut_PosDepCorr.Write();
    hist_delayed_E_ECut_Rcut_PosDepCorr.Write();
    hist_prompt_R_ECut_Rcut_PosDepCorr.Write();
    hist_delayed_R_ECut_Rcut_PosDepCorr.Write();
    hist_prompt_Nhits_ECut_Rcut_PosDepCorr.Write();
    hist_delayed_Nhits_ECut_Rcut_PosDepCorr.Write();
    hist_deltaR_ECut_Rcut_PosDepCorr.Write();
    hist_deltaT_ECut_Rcut_PosDepCorr.Write();

    /* ±±±±±±±±±±±±±±±±±±±±±±±± */

    // nhit cuts, no R cut, no pos dep corr
    hist_prompt_E_NhitCut_Ecorr.Write();
    hist_delayed_E_NhitCut_Ecorr.Write();
    hist_prompt_R_NhitCut_Ecorr.Write();
    hist_delayed_R_NhitCut_Ecorr.Write();
    hist_prompt_Nhits_NhitCut_Ecorr.Write();
    hist_delayed_Nhits_NhitCut_Ecorr.Write();
    hist_deltaR_NhitCut_Ecorr.Write();
    hist_deltaT_NhitCut_Ecorr.Write();

    // E cuts, no R cut, no pos dep corr
    hist_prompt_E_ECut_Ecorr.Write();
    hist_delayed_E_ECut_Ecorr.Write();
    hist_prompt_R_ECut_Ecorr.Write();
    hist_delayed_R_ECut_Ecorr.Write();
    hist_prompt_Nhits_ECut_Ecorr.Write();
    hist_delayed_Nhits_ECut_Ecorr.Write();
    hist_deltaR_ECut_Ecorr.Write();
    hist_deltaT_ECut_Ecorr.Write();

    // nhit cuts, R cut, no pos dep corr
    hist_prompt_E_NhitCut_Rcut_Ecorr.Write();
    hist_delayed_E_NhitCut_Rcut_Ecorr.Write();
    hist_prompt_R_NhitCut_Rcut_Ecorr.Write();
    hist_delayed_R_NhitCut_Rcut_Ecorr.Write();
    hist_prompt_Nhits_NhitCut_Rcut_Ecorr.Write();
    hist_delayed_Nhits_NhitCut_Rcut_Ecorr.Write();
    hist_deltaR_NhitCut_Rcut_Ecorr.Write();
    hist_deltaT_NhitCut_Rcut_Ecorr.Write();

    // E cuts, R cut, no pos dep corr
    hist_prompt_E_ECut_Rcut_Ecorr.Write();
    hist_delayed_E_ECut_Rcut_Ecorr.Write();
    hist_prompt_R_ECut_Rcut_Ecorr.Write();
    hist_delayed_R_ECut_Rcut_Ecorr.Write();
    hist_prompt_Nhits_ECut_Rcut_Ecorr.Write();
    hist_delayed_Nhits_ECut_Rcut_Ecorr.Write();
    hist_deltaR_ECut_Rcut_Ecorr.Write();
    hist_deltaT_ECut_Rcut_Ecorr.Write();

    /* ---- */

    // nhit cuts, no R cut, pos dep corr
    hist_prompt_E_NhitCut_PosDepCorr_Ecorr.Write();
    hist_delayed_E_NhitCut_PosDepCorr_Ecorr.Write();
    hist_prompt_R_NhitCut_PosDepCorr_Ecorr.Write();
    hist_delayed_R_NhitCut_PosDepCorr_Ecorr.Write();
    hist_prompt_Nhits_NhitCut_PosDepCorr_Ecorr.Write();
    hist_delayed_Nhits_NhitCut_PosDepCorr_Ecorr.Write();
    hist_deltaR_NhitCut_PosDepCorr_Ecorr.Write();
    hist_deltaT_NhitCut_PosDepCorr_Ecorr.Write();

    // E cuts, no R cut, pos dep corr
    hist_prompt_E_ECut_PosDepCorr_Ecorr.Write();
    hist_delayed_E_ECut_PosDepCorr_Ecorr.Write();
    hist_prompt_R_ECut_PosDepCorr_Ecorr.Write();
    hist_delayed_R_ECut_PosDepCorr_Ecorr.Write();
    hist_prompt_Nhits_ECut_PosDepCorr_Ecorr.Write();
    hist_delayed_Nhits_ECut_PosDepCorr_Ecorr.Write();
    hist_deltaR_ECut_PosDepCorr_Ecorr.Write();
    hist_deltaT_ECut_PosDepCorr_Ecorr.Write();

    // nhit cuts, R cut, pos dep corr
    hist_prompt_E_NhitCut_Rcut_PosDepCorr_Ecorr.Write();
    hist_delayed_E_NhitCut_Rcut_PosDepCorr_Ecorr.Write();
    hist_prompt_R_NhitCut_Rcut_PosDepCorr_Ecorr.Write();
    hist_delayed_R_NhitCut_Rcut_PosDepCorr_Ecorr.Write();
    hist_prompt_Nhits_NhitCut_Rcut_PosDepCorr_Ecorr.Write();
    hist_delayed_Nhits_NhitCut_Rcut_PosDepCorr_Ecorr.Write();
    hist_deltaR_NhitCut_Rcut_PosDepCorr_Ecorr.Write();
    hist_deltaT_NhitCut_Rcut_PosDepCorr_Ecorr.Write();

    // E cuts, R cut, pos dep corr
    hist_prompt_E_ECut_Rcut_PosDepCorr_Ecorr.Write();
    hist_delayed_E_ECut_Rcut_PosDepCorr_Ecorr.Write();
    hist_prompt_R_ECut_Rcut_PosDepCorr_Ecorr.Write();
    hist_delayed_R_ECut_Rcut_PosDepCorr_Ecorr.Write();
    hist_prompt_Nhits_ECut_Rcut_PosDepCorr_Ecorr.Write();
    hist_delayed_Nhits_ECut_Rcut_PosDepCorr_Ecorr.Write();
    hist_deltaR_ECut_Rcut_PosDepCorr_Ecorr.Write();
    hist_deltaT_ECut_Rcut_PosDepCorr_Ecorr.Write();


    rootfile.Write();
    rootfile.Close();
}


/* ~~~~~~~~~~~~~~~~~~~~~~ GET INFO FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~ */

/**
 * @brief Loop through events (before evt) to find a suitable prompt event, if it exists:
 * - Either an event is found that pass both the prompt event and tagging cuts -> record event info and return true.
 * - Or an event is found whose recon event time is futher from the delayed event than the max delta-t cut -> return false.
 * - Or the file runs out of events, or the max number of loops in reached -> return false.
 * 
 * @param dsReader  file info currently being analysed
 * @param entry  entry number delayed event is from
 * @param evt  event number delayed event is from
 * @param E  list of prompt and delayed pair recon energies (delayed energy is alread there, prompt will be saved there if found)
 * @param times  list of prompt and delayed pair recon times (delayed energy is alread there, prompt will be saved there if found)
 * @param pos  list of prompt and delayed pair recon positions (delayed energy is alread there, prompt will be saved there if found)
 * @param verbose  flag
 * @param fitName  name of fitter to get recon info from (default is "")
 * @return true 
 * @return false 
 */
bool find_prompt_event(RAT::DU::DSReader& dsReader, const int entry, const int evt, std::vector<double>& E, std::vector<double>& times, std::vector<TVector3>& pos, std::vector<double>& Nhits,
                       double& Delta_T, std::ofstream& t_res_file, RAT::DU::TimeResidualCalculator& fTRCalc, const RAT::DU::ReconCalibrator& e_cal, RAT::DU::DetectorStateCorrection& stateCorr,
                       bool is_data, const bool verbose, const bool use_pos_dep_corr, const bool use_E_corr, const bool use_Nhit, const bool cut_R_min) {

    RAT::DS::Entry rDS = dsReader.GetEntry(entry);
    RAT::DS::EV rEV = rDS.GetEV(evt);
    ULong64_t delayed_50MHz_time = rEV.GetClockCount50();
    if (verbose) std::cout << "delayed_50MHz_time = " << delayed_50MHz_time << std::endl;

    unsigned int k = 0;
    for (int iEntry = entry; iEntry >= 0; --iEntry) {
        if (verbose) std::cout << "iEntry = " << iEntry << std::endl;
        rDS = dsReader.GetEntry(iEntry);

        // loops through events
        int start_evt = rDS.GetEVCount() - 1;
        if (iEntry == entry) start_evt = evt - 1;
        for (int iEV = start_evt; iEV >= 0; --iEV) {
            if (verbose) std::cout << "iEV = " << iEV << std::endl;
            rEV = rDS.GetEV(iEV);

            if (!RAT::EventIsClean(rEV, dcAnalysisWord)) continue;

            // Get recon info, if it exists
            if (get_recon_info(E, times, pos, Nhits, 0, rEV, e_cal, stateCorr, is_data, use_pos_dep_corr, use_E_corr)) {
                if (verbose) std::cout << "rEV.GetClockCount50() = " << rEV.GetClockCount50() << std::endl;
                Delta_T = ((int64_t(delayed_50MHz_time) - int64_t(rEV.GetClockCount50())) & 0x7FFFFFFFFFF) * 20.0;
                if (verbose) std::cout << "Delta_T = " << Delta_T << std::endl;

                if (Delta_T > MAX_DELAY) return false;

                if (pass_prompt_cuts(E[0], Nhits[0], pos[0], use_Nhit, cut_R_min)) {
                    if (verbose) std::cout << "Passed prompt cuts! Checking coincidence cuts..." << std::endl;
                    if (pass_coincidence_cuts(Delta_T, pos[0], pos[1])) {
                        #ifdef PRINT_T_RES
                            // Also get time residuals while at it
                            get_t_res(t_res_file, rEV, fTRCalc, pos[0], times[0]);
                        #endif
                        return true;
                    }
                }
            }
            if (k > 100) return false;  // cut off looping after a certain point anyway
            k++;
        }
    }

    return false;
}

/**
 * @brief Save recon info of event, if it exists and is valid
 * 
 * @param E  list of prompt and delayed pair recon energies
 * @param times  list of prompt and delayed pair recon times
 * @param pos  list of prompt and delayed pair recon positions
 * @param idx  indicates which event (0 = prompt, 1 = delayed) to save info to
 * @param evt  event object
 * @param fitName  name of fitter to get recon info from (default is "")
 * @return true 
 * @return false 
 */
bool get_recon_info(std::vector<double>& E, std::vector<double>& times, std::vector<TVector3>& pos, std::vector<double>& Nhits,
                    const unsigned int idx, const RAT::DS::EV& evt, const RAT::DU::ReconCalibrator& e_cal, RAT::DU::DetectorStateCorrection& stateCorr,
                    bool is_data, const bool use_pos_dep_corr, const bool use_E_corr) {

    try {
        // Get recon info
        const RAT::DS::FitResult fitResult = evt.GetFitResult(evt.GetDefaultFitName());
        if (!fitResult.GetValid()) return false; // fit invalid
        const RAT::DS::FitVertex& rVertex = fitResult.GetVertex(0);
        if (!(rVertex.ValidPosition() && rVertex.ValidTime() && rVertex.ValidEnergy())) return false; // fit invalid
        times[idx] = rVertex.GetTime();
        pos[idx] = rVertex.GetPosition();
        E[idx] = rVertex.GetEnergy();
        Nhits[idx] = (double)evt.GetNhitsCleaned();

        if (use_E_corr) {
            // Data vs MC energy correction (Tony's)
            const double corrected_energy = e_cal.CalibrateEnergyRTF(is_data, E[idx], std::sqrt(pos[idx].X()*pos[idx].X() + pos[idx].Y()*pos[idx].Y()), pos[idx].Z()); // gives the new E
            Nhits[idx] *= corrected_energy / E[idx];
            E[idx] = corrected_energy;
        }

        if (use_pos_dep_corr) {
            // Correct for position coverage dependence (Logan's)
            RAT::DU::Point3D position(0, pos[idx]);  // position of event [mm] (as Point3D in PSUP coordinates, see system_id in POINT3D_SHIFTS tables)
            // std::cout << "pos[idx] = (" << pos[idx].X() << ", " << pos[idx].Y() << ", " << pos[idx].Z() << "), "
            //           << "position = (" << position.X() << ", " << position.Y() << ", " << position.Z() << ")." << std::endl;
            double E_state_corr = stateCorr.GetCorrectionPos(position, 0, 0) / stateCorr.GetCorrection(9394, 0.75058); // a correction factor (divide E by it)
            Nhits[idx] /= E_state_corr;
            E[idx] /= E_state_corr;
        }
    }
    catch (const RAT::DS::DataNotFound&) {return false;}  // no fit data
    catch (const RAT::DS::FitCollection::NoResultError&) {return false;} // no fit result by the name of fitName
    catch (const RAT::DS::FitResult::NoVertexError&) {return false;} // no fit vertex
    catch (const RAT::DS::FitVertex::NoValueError&) {return false;} // position or time missing
    // catch (const RAT::DS::ClassifierResult::NoClassificationError&) {return false;} // classifier result error

    return true;
}

/**
 * @brief Fill histogram with time residual information from event.
 * 
 * @param evt  event object
 * @param fTRCalc  time residual calculator object
 * @param position  prompt event recon position
 * @param vertex_time  prompt event recon time
 * @param hist  time residual histogram to be filled
 */
void get_t_res(std::ofstream& t_res_file, const RAT::DS::EV& evt, RAT::DU::TimeResidualCalculator& fTRCalc, const TVector3& position, const double vertex_time) {

    RAT::DU::Point3D pos(0, position);  // position of event [mm] (as Point3D in PSUP coordinates, see system_id in POINT3D_SHIFTS tables)

    // Compute time residuals
    const RAT::DS::CalPMTs& calibratedPMTs = evt.GetCalPMTs();
    for (size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); ++iPMT) {
        const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT(iPMT);
        // Use new time residual calculator
        t_res_file << fTRCalc.CalcTimeResidual(pmtCal, pos, vertex_time) << ", ";
    }
}

// bin = 0;       underflow bin
// bin = 1;       first bin with low-edge xlow INCLUDED
// bin = nbins;   last bin with upper-edge xup EXCLUDED
// bin = nbins+1; overflow bin
