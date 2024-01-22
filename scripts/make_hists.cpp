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


bool find_prompt_event(RAT::DU::DSReader& dsReader, const int entry, const int evt, std::vector<double>& E, std::vector<double>& times, std::vector<TVector3>& pos, std::vector<double>& Nhits,
                       double& Delta_T, std::ofstream& t_res_file, RAT::DU::TimeResidualCalculator& fTRCalc, const RAT::DU::ReconCalibrator& e_cal, RAT::DU::DetectorStateCorrection& stateCorr, bool is_data, const bool verbose, const std::string fitName);
bool get_recon_info(std::vector<double>& E, std::vector<double>& times, std::vector<TVector3>& pos, std::vector<double>& Nhits, const unsigned int idx, const RAT::DS::EV& evt, const RAT::DU::ReconCalibrator& e_cal, RAT::DU::DetectorStateCorrection& stateCorr, bool is_data, std::string fitName);
void get_t_res(std::ofstream& t_res_file, const RAT::DS::EV& evt, RAT::DU::TimeResidualCalculator& fTRCalc, const TVector3& position, const double vertex_time);
void make_hists(const std::vector<std::string>& fileNames, const std::string output_root_address, const std::string output_txt_address, bool is_data, const bool verbose, const std::string fitName = "");
bool pass_prompt_cuts(const double energy, const double Nhit, TVector3 position);
bool pass_delayed_cuts(const double energy, const double Nhit, TVector3 position);
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

const double R_MIN = 0.0, R_MAX = 5700.0;
const double MAX_DIST = 1500.0;
// const double MIN_DELAY = 5500.0, MAX_DELAY = 0.8E6;
const double MIN_DELAY = 400.0, MAX_DELAY = 0.8E6;

// const double MIN_PROMPT_E = 0.7, MAX_PROMPT_E = 3.5;
double MIN_PROMPT_E = 0.3, MAX_PROMPT_E = 4.5;
const double MIN_PROMPT_Nhit = 120.0, MAX_PROMPT_Nhit = 680.0;

// const double MIN_DELAYED_E = 1.4, MAX_DELAYED_E = 2.8;
const double MIN_DELAYED_E = 1.0, MAX_DELAYED_E = 3.5;
const double MIN_DELAYED_Nhit = 350.0, MAX_DELAYED_Nhit = 600.0;

const bool USE_NHIT = true;

bool pass_prompt_cuts(const double energy, const double Nhit, TVector3 position) {
    if (USE_NHIT) {
        if (Nhit < MIN_PROMPT_Nhit) return false;  // min Nhit cut
        if (Nhit > MAX_PROMPT_Nhit) return false;  // max Nhit cut
    } else {
        if (energy < MIN_PROMPT_E) return false;  // min energy cut (MeV)
        if (energy > MAX_PROMPT_E) return false;  // max energy cut (MeV)
    }
    if (position.Mag() > R_MAX) return false;  // FV cut (mm)
    if (position.Mag() < R_MIN) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts(const double energy, const double Nhit, TVector3 position) {
    if (USE_NHIT) {
        if (Nhit < MIN_DELAYED_Nhit) return false;  // min Nhit cut
        if (Nhit > MAX_DELAYED_Nhit) return false;  // max Nhit cut
    } else {
        if (energy < MIN_DELAYED_E) return false;  // min energy cut (MeV)
        if (energy > MAX_DELAYED_E) return false;  // max energy cut (MeV)
    }
    if (position.Mag() > R_MAX) return false;  // FV cut (mm)
    if (position.Mag() < R_MIN) return false;  // FV cut (mm)

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
    TH1D hist_prompt_E("prompt_E", "prompt_E", 100, MIN_PROMPT_E, MAX_PROMPT_E);
    TH1D hist_delayed_E("delayed_E", "delayed_E", 100, MIN_DELAYED_E, MAX_DELAYED_E);
    TH1D hist_prompt_R("prompt_R", "prompt_R", 100, R_MIN, R_MAX);
    TH1D hist_delayed_R("delayed_R", "delayed_R", 100, R_MIN, R_MAX);
    TH1D hist_prompt_Nhits("prompt_Nhits", "prompt_Nhits", 100, MIN_PROMPT_Nhit, MAX_PROMPT_Nhit);
    TH1D hist_delayed_Nhits("delayed_Nhits", "delayed_Nhits", 100, MIN_DELAYED_Nhit, MAX_DELAYED_Nhit);
    TH1D hist_deltaR("deltaR", "deltaR", 100, 0.0, MAX_DIST);
    TH1D hist_deltaT("deltaT", "deltaT", 100, MIN_DELAY, MAX_DELAY);

    TH1D hist_prompt_E_54R("prompt_E_5.4<R", "prompt_E_5.4<R", 100, MIN_PROMPT_E, MAX_PROMPT_E);
    TH1D hist_delayed_E_54R("delayed_E_5.4<R", "delayed_E_5.4<R", 100, MIN_DELAYED_E, MAX_DELAYED_E);
    TH1D hist_prompt_Nhits_54R("prompt_Nhits_5.4<R", "prompt_Nhits_5.4<R", 100, MIN_PROMPT_Nhit, MAX_PROMPT_Nhit);
    TH1D hist_delayed_Nhits_54R("delayed_Nhits_5.4<R", "delayed_Nhits_5.4<R", 100, MIN_PROMPT_Nhit, MAX_PROMPT_Nhit);
    TH1D hist_deltaR_54R("deltaR_5.4<R", "deltaR_5.4<R", 100, 0.0, MAX_DIST);
    TH1D hist_deltaT_54R("deltaT_5.4<R", "deltaT_5.4<R", 100, MIN_DELAY, MAX_DELAY);

    TH1D hist_prompt_E_4R54("prompt_E_4<R<5.4", "prompt_E_4<R<5.4", 100, MIN_PROMPT_E, MAX_PROMPT_E);
    TH1D hist_delayed_E_4R54("delayed_E_4<R<5.4", "delayed_E_4<R<5.4", 100, MIN_DELAYED_E, MAX_DELAYED_E);
    TH1D hist_prompt_Nhits_4R54("prompt_Nhits_4<R<5.4", "prompt_Nhits_4<R<5.4", 100, MIN_PROMPT_Nhit, MAX_PROMPT_Nhit);
    TH1D hist_delayed_Nhits_4R54("delayed_Nhits_4<R<5.4", "delayed_Nhits_4<R<5.4", 100, MIN_PROMPT_Nhit, MAX_PROMPT_Nhit);
    TH1D hist_deltaR_4R54("deltaR_4<R<5.4", "deltaR_4<R<5.4", 100, 0.0, MAX_DIST);
    TH1D hist_deltaT_4R54("deltaT_4<R<5.4", "deltaT_4<R<5.4", 100, MIN_DELAY, MAX_DELAY);

    TH1D hist_prompt_E_R4("prompt_E_R<4", "prompt_E_R<4", 100, MIN_PROMPT_E, MAX_PROMPT_E);
    TH1D hist_delayed_E_R4("delayed_E_R<4", "delayed_E_R<4", 100, MIN_DELAYED_E, MAX_DELAYED_E);
    TH1D hist_prompt_Nhits_R4("prompt_Nhits_R<4", "prompt_Nhits_R<4", 100, MIN_PROMPT_Nhit, MAX_PROMPT_Nhit);
    TH1D hist_delayed_Nhits_R4("delayed_Nhits_R<4", "delayed_Nhits_R<4", 100, MIN_PROMPT_Nhit, MAX_PROMPT_Nhit);
    TH1D hist_deltaR_R4("deltaR_R<4", "deltaR_R<4", 100, 0.0, MAX_DIST);
    TH1D hist_deltaT_R4("deltaT_R<4", "deltaT_R<4", 100, MIN_DELAY, MAX_DELAY);

    TH1D hist_prompt_E_rest("prompt_E_rest", "prompt_E_rest", 100, MIN_PROMPT_E, MAX_PROMPT_E);
    TH1D hist_delayed_E_rest("delayed_E_rest", "delayed_E_rest", 100, MIN_DELAYED_E, MAX_DELAYED_E);
    TH1D hist_prompt_Nhits_rest("prompt_Nhits_rest", "prompt_Nhits_rest", 100, MIN_PROMPT_Nhit, MAX_PROMPT_Nhit);
    TH1D hist_delayed_Nhits_rest("delayed_Nhits_rest", "delayed_Nhits_rest", 100, MIN_PROMPT_Nhit, MAX_PROMPT_Nhit);
    TH1D hist_deltaR_rest("deltaR_rest", "deltaR_rest", 100, 0.0, MAX_DIST);
    TH1D hist_deltaT_rest("deltaT_rest", "deltaT_rest", 100, MIN_DELAY, MAX_DELAY);

    TH2D hist_prompt_E_vs_R("prompt_E_vs_R", "prompt_E_vs_R", 100, MIN_PROMPT_E, MAX_PROMPT_E, 100, R_MIN, R_MAX);
    TH2D hist_delayed_E_vs_R("delayed_E_vs_R", "delayed_E_vs_R", 100, MIN_DELAYED_E, MAX_DELAYED_E, 100, R_MIN, R_MAX);
    TH2D hist_prompt_Nhits_vs_R("prompt_Nhits_vs_R", "prompt_Nhits_vs_R", 100, MIN_PROMPT_Nhit, MAX_PROMPT_Nhit, 100, R_MIN, R_MAX);
    TH2D hist_delayed_Nhits_vs_R("delayed_Nhits_vs_R", "delayed_Nhits_vs_R", 100, MIN_PROMPT_Nhit, MAX_PROMPT_Nhit, 100, R_MIN, R_MAX);

    // RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();
    // RAT::DU::ReconCalibrator e_cal;
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
        // e_cal.Get();
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

                // if (RAT::EventIsClean(rEV, dcAnalysisWord)) {}

                // Get recon info, if it exists
                if (get_recon_info(E, times, pos, Nhits, 1, rEV, e_cal, stateCorr, is_data, fitName)) {
                    // Check if event passes delayed cuts
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1])) {
                        if (verbose) std::cout << "Passed delayed cuts! Checking for prompt event..." << std::endl;
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, fitName)) {
                            if (verbose) std::cout << "Passed coincidence cuts! Filling histograms..." << std::endl;
                            // original hists
                            hist_prompt_E.Fill(E[0]);
                            hist_delayed_E.Fill(E[1]);
                            hist_prompt_R.Fill(pos[0].Mag());
                            hist_delayed_R.Fill(pos[1].Mag());
                            hist_prompt_Nhits.Fill(Nhits[0]);
                            hist_delayed_Nhits.Fill(Nhits[1]);
                            hist_deltaR.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT.Fill(Delta_T);

                            // 2D hists
                            hist_prompt_E_vs_R.Fill(E[0], pos[0].Mag());
                            hist_delayed_E_vs_R.Fill(E[1], pos[1].Mag());
                            hist_prompt_Nhits_vs_R.Fill(Nhits[0], pos[0].Mag());
                            hist_delayed_Nhits_vs_R.Fill(Nhits[1], pos[1].Mag());

                            // FV cut hists
                            if ((pos[0].Mag() > 5400.0) && (pos[1].Mag() > 5400.0)) {
                                hist_prompt_E_54R.Fill(E[0]);
                                hist_delayed_E_54R.Fill(E[1]);
                                hist_prompt_Nhits_54R.Fill(Nhits[0]);
                                hist_delayed_Nhits_54R.Fill(Nhits[1]);
                                hist_deltaR_54R.Fill((pos[1] - pos[0]).Mag());
                                hist_deltaT_54R.Fill(Delta_T);
                            } else if ((pos[0].Mag() <= 5400.0) && (pos[1].Mag() <= 5400.0) && (pos[0].Mag() > 4000.0) && (pos[1].Mag() > 4000.0)) {
                                hist_prompt_E_4R54.Fill(E[0]);
                                hist_delayed_E_4R54.Fill(E[1]);
                                hist_prompt_Nhits_4R54.Fill(Nhits[0]);
                                hist_delayed_Nhits_4R54.Fill(Nhits[1]);
                                hist_deltaR_4R54.Fill((pos[1] - pos[0]).Mag());
                                hist_deltaT_4R54.Fill(Delta_T);
                            } else if ((pos[0].Mag() <= 4000.0) && (pos[1].Mag() <= 4000.0)) {
                                hist_prompt_E_R4.Fill(E[0]);
                                hist_delayed_E_R4.Fill(E[1]);
                                hist_prompt_Nhits_R4.Fill(Nhits[0]);
                                hist_delayed_Nhits_R4.Fill(Nhits[1]);
                                hist_deltaR_R4.Fill((pos[1] - pos[0]).Mag());
                                hist_deltaT_R4.Fill(Delta_T);
                            } else {
                                hist_prompt_E_rest.Fill(E[0]);
                                hist_delayed_E_rest.Fill(E[1]);
                                hist_prompt_Nhits_rest.Fill(Nhits[0]);
                                hist_delayed_Nhits_rest.Fill(Nhits[1]);
                                hist_deltaR_rest.Fill((pos[1] - pos[0]).Mag());
                                hist_deltaT_rest.Fill(Delta_T);
                            }
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

    hist_prompt_E.Write();
    hist_delayed_E.Write();
    hist_prompt_Nhits.Write();
    hist_delayed_Nhits.Write();
    hist_prompt_R.Write();
    hist_delayed_R.Write();
    hist_deltaR.Write();
    hist_deltaT.Write();

    hist_prompt_E_54R.Write();
    hist_delayed_E_54R.Write();
    hist_prompt_Nhits_54R.Write();
    hist_delayed_Nhits_54R.Write();
    hist_deltaR_54R.Write();
    hist_deltaT_54R.Write();
    
    hist_prompt_E_4R54.Write();
    hist_delayed_E_4R54.Write();
    hist_prompt_Nhits_4R54.Write();
    hist_delayed_Nhits_4R54.Write();
    hist_deltaR_4R54.Write();
    hist_deltaT_4R54.Write();

    hist_prompt_E_R4.Write();
    hist_delayed_E_R4.Write();
    hist_prompt_Nhits_R4.Write();
    hist_delayed_Nhits_R4.Write();
    hist_deltaR_R4.Write();
    hist_deltaT_R4.Write();

    hist_prompt_E_rest.Write();
    hist_delayed_E_rest.Write();
    hist_prompt_Nhits_rest.Write();
    hist_delayed_Nhits_rest.Write();
    hist_deltaR_rest.Write();
    hist_deltaT_rest.Write();

    hist_prompt_E_vs_R.Write();
    hist_delayed_E_vs_R.Write();
    hist_prompt_Nhits_vs_R.Write();
    hist_delayed_Nhits_vs_R.Write();

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
                       bool is_data, const bool verbose, const std::string fitName) {

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

            // Get recon info, if it exists
            if (get_recon_info(E, times, pos, Nhits, 0, rEV, e_cal, stateCorr, is_data, fitName)) {
                if (verbose) std::cout << "rEV.GetClockCount50() = " << rEV.GetClockCount50() << std::endl;
                Delta_T = ((int64_t(delayed_50MHz_time) - int64_t(rEV.GetClockCount50())) & 0x7FFFFFFFFFF) * 20.0;
                if (verbose) std::cout << "Delta_T = " << Delta_T << std::endl;

                if (Delta_T > MAX_DELAY) return false;

                if (pass_prompt_cuts(E[0], Nhits[0], pos[0])) {
                    if (verbose) std::cout << "Passed prompt cuts! Checking coincidence cuts..." << std::endl;
                    if (pass_coincidence_cuts(Delta_T, pos[0], pos[1])) {
                        // Also get time residuals while at it
                        get_t_res(t_res_file, rEV, fTRCalc, pos[0], times[0]);
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
                    bool is_data, std::string fitName) {

    // Grab the fit information
    if (fitName == "") fitName = evt.GetDefaultFitName();
    try {
        // Get recon info
        const RAT::DS::FitResult fitResult = evt.GetFitResult(fitName);
        if (!fitResult.GetValid()) return false; // fit invalid
        const RAT::DS::FitVertex& rVertex = fitResult.GetVertex(0);
        if (!(rVertex.ValidPosition() && rVertex.ValidTime() && rVertex.ValidEnergy())) return false; // fit invalid
        times[idx] = rVertex.GetTime();
        pos[idx] = rVertex.GetPosition();
        E[idx] = rVertex.GetEnergy();
        Nhits[idx] = (double)evt.GetNhitsCleaned();

        // Data vs MC energy correction (Tony's)
        const double corrected_energy = e_cal.CalibrateEnergyRTF(is_data, E[idx], std::sqrt(pos[idx].X()*pos[idx].X() + pos[idx].Y()*pos[idx].Y()), pos[idx].Z()); // gives the new E

        // Correct for position coverage dependence (Logan's)
        RAT::DU::Point3D position(0, pos[idx]);  // position of event [mm] (as Point3D in PSUP coordinates, see system_id in POINT3D_SHIFTS tables)
        double E_corr = stateCorr.GetCorrectionPos(position, 0, 0) / stateCorr.GetCorrection(7733, 0.768972); // a correction factor (divide E by it)

        // Apply corrections: also scale the Nhits by the same amount, for consistency (and in case Nhits are used for cuts)
        Nhits[idx] *= (corrected_energy / E[idx]) / E_corr;
        E[idx] = corrected_energy / E_corr;
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
