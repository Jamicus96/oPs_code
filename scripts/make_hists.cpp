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

#define LPC_uses_Point3D  // Uncomment if rat version is >= 7.0.11 (the light path calculator uses Point3D instead of Vector3 after this point)

bool find_prompt_event(RAT::DU::DSReader& dsReader, const int entry, const int evt, std::vector<double>& E, std::vector<double>& times, std::vector<TVector3>& pos, double& Delta_T, const bool verbose, const std::string fitName);
bool get_recon_info(std::vector<double>& E, std::vector<double>& times, std::vector<TVector3>& pos, const unsigned int idx, const RAT::DS::EV& evt, std::string fitName);
void get_t_res(const RAT::DS::EV& evt, RAT::DU::TimeResidualCalculator& fTRCalc, const TVector3& position, const double vertex_time, TH1D& hist);
void make_hists(const std::vector<std::string>& fileNames, const std::string output_file_address, const bool verbose, const std::string fitName = "");
bool pass_prompt_cuts(const double energy, TVector3 position);
bool pass_delayed_cuts(const double energy, TVector3 position);
bool pass_coincidence_cuts(const double delay, TVector3 prompt_pos, TVector3 delayed_pos);


/* ~~~~~~~~~~~~~~~~~~~~~~ MAIN FUNCTION ~~~~~~~~~~~~~~~~~~~~~ */

int main(int argc, char** argv) {
    std::string output_file_address = argv[1];
    bool verbose = std::stoi(argv[2]);
    // Addresses of simulation output files to be analysed
    std::vector<std::string> input_files;
    for (unsigned int i = 3; i < argc; ++i) {
        input_files.push_back(argv[i]);
    }

    // Loop through files to get info from every event (including t_res), and print all to text file
    if (verbose) {std::cout << "Getting info..." << std::endl;}
    make_hists(input_files, output_file_address, verbose);

    return 0;
}


/* ~~~~~~~~~~~~~~~~~~~~~~ CUT FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~ */

double MAX_DELAY = 1.1E6;
double R_CUT = 5700.0;

bool pass_prompt_cuts(const double energy, TVector3 position) {
    if (energy < 0.9) return false;  // min energy cut (MeV)
    if (energy > 3.5) return false;  // max energy cut (MeV)
    if (position.Mag() > R_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts(const double energy, TVector3 position) {
    if (energy < 1.7) return false;  // min energy cut (MeV)
    if (energy > 2.5) return false;  // max energy cut (MeV)
    if (position.Mag() > R_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_coincidence_cuts(const double delay, TVector3 prompt_pos, TVector3 delayed_pos) {
    // double delay = (delayed_time - prompt_time) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
    double distance = (delayed_pos - prompt_pos).Mag();

    if (delay < 200.0) return false;  // min delay cut (ns)
    if (delay > MAX_DELAY) return false;  // max delay cut (ns)
    if (distance > 1500.0) return false;  // max distance cut (mm)

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
void make_hists(const std::vector<std::string>& fileNames, const std::string output_file_address, const bool verbose, const std::string fitName) {
    if (verbose) {std::cout << "Running make_hists()" << std::endl;}

    /*********** Open output file, to print info to ***********/
    TFile rootfile(output_file_address.c_str(), "RECREATE");

    /*********** Setup hists ***********/
    TH1D hist_prompt_t_res("prompt_t_res", "prompt_t_res", 1300, -300.5, 999.5);
    TH1D hist_prompt_E("prompt_E", "prompt_E", 100, 0.9, 3.5);
    TH1D hist_delayed_E("delayed_E", "delayed_E", 100, 1.7, 2.5);
    TH1D hist_prompt_R("prompt_R", "prompt_R", 100, 0.0, R_CUT);
    TH1D hist_delayed_R("delayed_R", "delayed_R", 100, 0.0, R_CUT);
    TH1D hist_deltaR("deltaR", "deltaR", 100, 0.0, 1500.0);
    TH1D hist_deltaT("deltaT", "deltaT", 100, 200.0, MAX_DELAY);

    // RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();


    /*********** Loop through all files, entries, events, and PMTs ***********/
    
    // Set up empty vectors for {prompt, delayed} recon quatities
    std::vector<double> E = {0.0, 0.0};
    std::vector<double> times = {0.0, 0.0};
    std::vector<TVector3> pos = {TVector3(0.0, 0.0, 0.0), TVector3(0.0, 0.0, 0.0)};
    double Delta_T = 0.0;

    // loops through files
    for (unsigned int i = 0; i < fileNames.size(); ++i) {
        RAT::DU::DSReader dsReader(fileNames.at(i));
        fTRCalc.BeginOfRun();  // Re-initialize time residual calculator (light-path calculator) after it gets geo info from DSReader
        
        if (verbose) {std::cout << "Looping through entries..." << std::endl;}
        // loops through entries
        for (unsigned int iEntry = 0; iEntry < dsReader.GetEntryCount(); ++iEntry) {
            if (verbose) std::cout << "iEntry = " << iEntry << std::endl;
            RAT::DS::Entry rDS = dsReader.GetEntry(iEntry);

            // loops through events
            for (unsigned int iEV = 0; iEV < rDS.GetEVCount(); ++iEV) {
                if (verbose) std::cout << "iEV = " << iEV << std::endl;
                RAT::DS::EV rEV = rDS.GetEV(iEV);

                // Get recon info, if it exists
                if (get_recon_info(E, times, pos, 1, rEV, fitName)) {
                    // Check if event passes delayed cuts
                    if (pass_delayed_cuts(E[1], pos[1])) {
                        if (verbose) std::cout << "Passed delayed cuts! Checking for delayed event..." << std::endl;
                        // Loop over previous events for one that passes prompt and tag cuts (stops after time diff exceeds cut, or max loops)
                        if (find_prompt_event(dsReader, iEntry, iEV, E, times, pos, Delta_T, verbose, fitName)) {
                            if (verbose) std::cout << "Passed coincidence cuts! Filling histograms..." << std::endl;

                            // Fill histograms!
                            hist_prompt_E.Fill(E[0]);
                            hist_delayed_E.Fill(E[1]);
                            hist_prompt_R.Fill(pos[0].Mag());
                            hist_delayed_R.Fill(pos[1].Mag());
                            hist_deltaR.Fill((pos[1] - pos[0]).Mag());
                            hist_deltaT.Fill(Delta_T);
                            get_t_res(rEV, fTRCalc, pos[0], times[0], hist_prompt_t_res);
                        }
                    }
                }
            }
        }
        dsReader.Delete();
    }

    //now write everything
    rootfile.cd();
    hist_prompt_E.Write();
    hist_delayed_E.Write();
    hist_prompt_R.Write();
    hist_delayed_R.Write();
    hist_deltaR.Write();
    hist_deltaT.Write();
    hist_prompt_t_res.Write();

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
bool find_prompt_event(RAT::DU::DSReader& dsReader, const int entry, const int evt, std::vector<double>& E,
                        std::vector<double>& times, std::vector<TVector3>& pos, double& Delta_T, const bool verbose, const std::string fitName) {

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
            if (get_recon_info(E, times, pos, 0, rEV, fitName)) {
                if (verbose) std::cout << "rEV.GetClockCount50() = " << rEV.GetClockCount50() << std::endl;
                Delta_T = ((int64_t(delayed_50MHz_time) - int64_t(rEV.GetClockCount50())) & 0x7FFFFFFFFFF) * 20.0;
                if (verbose) std::cout << "Delta_T = " << Delta_T << std::endl;

                if (Delta_T > MAX_DELAY) return false;

                if (pass_prompt_cuts(E[0], pos[0])) {
                    if (verbose) std::cout << "Passed prompt cuts! Checking coincidence cuts..." << std::endl;
                    if (pass_coincidence_cuts(Delta_T, pos[0], pos[1])) return true;
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
bool get_recon_info(std::vector<double>& E, std::vector<double>& times, std::vector<TVector3>& pos, const unsigned int idx, const RAT::DS::EV& evt, std::string fitName) {

    // Grab the fit information
    if (fitName == "") fitName = evt.GetDefaultFitName();
    try {
        // Get recon info
        const RAT::DS::FitVertex& rVertex = evt.GetFitResult(fitName).GetVertex(0);
        if (!(rVertex.ValidPosition() && rVertex.ValidTime() && rVertex.ValidEnergy())) {return false;} // fit invalid
        E[idx] = rVertex.GetEnergy();
        times[idx] = rVertex.GetTime();
        pos[idx] = rVertex.GetPosition();
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
void get_t_res(const RAT::DS::EV& evt, RAT::DU::TimeResidualCalculator& fTRCalc, const TVector3& position, const double vertex_time, TH1D& hist) {

    // Deal with RAT backwards compatifility
    #ifdef LPC_uses_Point3D
        RAT::DU::Point3D pos(0, position);  // position of event [mm] (as Point3D in PSUP coordinates, see system_id in POINT3D_SHIFTS tables)
    #else
        TVector3 pos = position;
    #endif

    // Compute time residuals
    const RAT::DS::CalPMTs& calibratedPMTs = evt.GetCalPMTs();
    for (size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); ++iPMT) {
        const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT(iPMT);
        // Use new time residual calculator
        hist.Fill(fTRCalc.CalcTimeResidual(pmtCal, pos, vertex_time));
    }
}

// bin = 0;       underflow bin
// bin = 1;       first bin with low-edge xlow INCLUDED
// bin = nbins;   last bin with upper-edge xup EXCLUDED
// bin = nbins+1; overflow bin