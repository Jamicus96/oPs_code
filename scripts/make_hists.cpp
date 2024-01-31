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
    bool flat_E_prompt = std::stoi(argv[3]);  // Old legacy entry, left for backwards compatibility with wrapper code (unused here)
    bool is_data = std::stoi(argv[4]);
    bool verbose = std::stoi(argv[5]);
    // Addresses of simulation output files to be analysed
    std::vector<std::string> input_files;
    for (unsigned int i = 6; i < argc; ++i) {
        input_files.push_back(argv[i]);
    }

    // Loop through files to get info from every event (including t_res) and write to files (t_res to txt, others to root).
    if (verbose) {std::cout << "Getting info..." << std::endl;}
    make_hists(input_files, output_root_address, output_txt_address, is_data, verbose);

    return 0;
}


/* ~~~~~~~~~~~~~~~~~~~~~~ CUT FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~ */

const double R_MIN = 0.0, R_MAX = 5700.0;
const double MAX_DIST = 1500.0;
const double MIN_DELAY = 500.0, MAX_DELAY = 0.8E6;
const double MIN_PROMPT_E = 0.7, MAX_PROMPT_E = 3.5;
const double MIN_DELAYED_E = 1.8, MAX_DELAYED_E = 2.5;

ULong64_t dcAnalysisWord = 36283883733698;  // Converted hex to decimal from 0x2100000042C2


bool pass_prompt_cuts(const double energy, const double Nhit, TVector3 position) {
    if (energy < MIN_PROMPT_E) return false;  // min energy cut (MeV)
    if (energy > MAX_PROMPT_E) return false;  // max energy cut (MeV)
    if (position.Mag() > R_MAX) return false;  // FV cut (mm)
    if (position.Mag() < R_MIN) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts(const double energy, const double Nhit, TVector3 position) {
    if (energy < MIN_DELAYED_E) return false;  // min energy cut (MeV)
    if (energy > MAX_DELAYED_E) return false;  // max energy cut (MeV)
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
 * @param fileNames  list of root file addresses to loop through with event info (ratds files).
 * @param output_root_address  output root file to save histograms to.
 * @param output_txt_address   output text file to save t_res of every PMT hit (of passed prompt events) to.
 * @param is_data  flag for energy correction and data cleaning.
 * @param verbose  flag.
 * @param fitName  name of fitter to get recon info from (default is "").
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
    TH1D hist_prompt_Nhits("prompt_Nhits", "prompt_Nhits", 100, 100, 900);
    TH1D hist_delayed_Nhits("delayed_Nhits", "delayed_Nhits", 100, 200, 800);
    TH1D hist_deltaR("deltaR", "deltaR", 100, 0.0, MAX_DIST);
    TH1D hist_deltaT("deltaT", "deltaT", 100, MIN_DELAY, MAX_DELAY);

    // RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();

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
                if (get_recon_info(E, times, pos, Nhits, 1, rEV, e_cal, stateCorr, is_data, fitName)) {
                    // Check if event passes delayed cuts
                    if (pass_delayed_cuts(E[1], Nhits[1], pos[1])) {
                        if (verbose) std::cout << "Passed delayed cuts! Checking for prompt event..." << std::endl;
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Nhits, Delta_T, t_res_file, fTRCalc, e_cal, stateCorr, is_data, verbose, fitName)) {
                            if (verbose) std::cout << "Passed coincidence cuts! Filling histograms..." << std::endl;
                            hist_prompt_E.Fill(E[0]);
                            hist_delayed_E.Fill(E[1]);
                            hist_prompt_R.Fill(pos[0].Mag());
                            hist_delayed_R.Fill(pos[1].Mag());
                            hist_prompt_Nhits.Fill(Nhits[0]);
                            hist_delayed_Nhits.Fill(Nhits[1]);
                            hist_deltaR.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT.Fill(Delta_T);
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
 * @param dsReader  file info currently being analysed.
 * @param entry  entry number delayed event is from.
 * @param evt  event number delayed event is from.
 * @param E  list of prompt and delayed pair recon energies (delayed energy is alread there, prompt will be saved there if found).
 * @param times  list of prompt and delayed pair recon times (delayed time is alread there, prompt will be saved there if found).
 * @param pos  list of prompt and delayed pair recon positions (delayed pos is alread there, prompt will be saved there if found).
 * @param Nhits  list of prompt and delayed pair Nhits (delayed Nhit is alread there, prompt will be saved there if found).
 * @param Delta_T  reference to Delta_T that will be computed here.
 * @param t_res_file  text file that t_res will be written to.
 * @param fTRCalc  time residual calculator.
 * @param e_cal  energy correction calculator.
 * @param stateCorr  position dependent energy correction calculator.
 * @param is_data  flag for energy correction.
 * @param verbose  flag.
 * @param fitName  name of fitter to get recon info from (default is "").
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
                        get_t_res(t_res_file, rEV, fTRCalc, pos[0], times[0]);
                        return true;
                    }
                }
            }
            if (k > 200) return false;  // cut off looping after a certain point anyway
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
 * @param e_cal  energy correction calculator.
 * @param stateCorr  position dependent energy correction calculator.
 * @param is_data  flag for energy correction.
 * @param fitName  name of fitter to get recon info from (default is "")
 * @return true 
 * @return false 
 */
bool get_recon_info(std::vector<double>& E, std::vector<double>& times, std::vector<TVector3>& pos, std::vector<double>& Nhits,
                    const unsigned int idx, const RAT::DS::EV& evt, const RAT::DU::ReconCalibrator& e_cal, RAT::DU::DetectorStateCorrection& stateCorr,
                    bool is_data, std::string fitName) {

    // Grab the fit information
    if (fitName == "") fitName = evt.GetDefaultFitName();
    Nhits[idx] = evt.GetNhitsCleaned();
    try {
        // Get recon info
        const RAT::DS::FitResult fitResult = evt.GetFitResult(fitName);
        if (!fitResult.GetValid()) return false; // fit invalid
        const RAT::DS::FitVertex& rVertex = fitResult.GetVertex(0);
        if (!(rVertex.ValidPosition() && rVertex.ValidTime() && rVertex.ValidEnergy())) return false; // fit invalid
        times[idx] = rVertex.GetTime();
        pos[idx] = rVertex.GetPosition();
        E[idx] = rVertex.GetEnergy();

        // Data vs MC energy correction (Tony's)
        E[idx] = e_cal.CalibrateEnergyRTF(is_data, E[idx], std::sqrt(pos[idx].X()*pos[idx].X() + pos[idx].Y()*pos[idx].Y()), pos[idx].Z()); // gives the new E

        // Correct for position coverage dependence (Logan's)
        RAT::DU::Point3D position(0, pos[idx]);  // position of event [mm] (as Point3D in PSUP coordinates, see system_id in POINT3D_SHIFTS tables)
        E[idx] /= stateCorr.GetCorrectionPos(position, 0, 0) / stateCorr.GetCorrection(9394, 0.75058); // a correction factor (divide E by it)
    }
    catch (const RAT::DS::DataNotFound&) {return false;}  // no fit data
    catch (const RAT::DS::FitCollection::NoResultError&) {return false;} // no fit result by the name of fitName
    catch (const RAT::DS::FitResult::NoVertexError&) {return false;} // no fit vertex
    catch (const RAT::DS::FitVertex::NoValueError&) {return false;} // position or time missing
    // catch (const RAT::DS::ClassifierResult::NoClassificationError&) {return false;} // classifier result error

    return true;
}

/**
 * @brief Print ime residual information from event to text file.
 * 
 * @param t_res_file Text file.
 * @param evt  event object.
 * @param fTRCalc  time residual calculator object.
 * @param position  prompt event recon position.
 * @param vertex_time  prompt event recon time.
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
