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

void print_info_to_file(const std::vector<std::string>& fileNames, const std::string output_file_address, unsigned int starting_fileNum, bool verbose, std::string fitName="");
std::vector<double> get_MC_info(const RAT::DS::Entry& entry, const RAT::DS::EV& evt, unsigned int evt_idx, RAT::DU::TimeResidualCalculator fTRCalc, unsigned int Nbins, double lower_lim, double upper_lim);
std::vector<double> get_recon_info(const RAT::DS::EV& evt, RAT::DU::TimeResidualCalculator fTRCalc, unsigned int Nbins, double lower_lim, double upper_lim, std::string fitName);



/* ~~~~~~~~~~~~~~~~~~~~~~ MAIN FUNCTION ~~~~~~~~~~~~~~~~~~~~~ */

int main(int argc, char** argv) {
    std::string output_file_address = argv[1];
    unsigned int starting_fileNum = std::stoi(argv[2]);
    bool verbose = std::stoi(argv[3]);
    // Addresses of simulation output files to be analysed
    std::vector<std::string> input_files;
    for (unsigned int i = 4; i < argc; ++i) {
        input_files.push_back(argv[i]);
    }

    // Loop through files to get info from every event (including t_res), and print all to text file
    if (verbose) {std::cout << "Getting info..." << std::endl;}
    print_info_to_file(input_files, output_file_address, starting_fileNum, verbose);


    return 0;
}


/* ~~~~~~~~~~~~~~~~~~~~~~ CUT FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~ */

double MAX_DELAY = 1.1E6;

bool pass_prompt_cuts(double energy, TVector3 position, double R_cut = 5700.0) {
    if (energy < 0.9) return false;  // min energy cut (MeV)
    if (energy > 3.5) return false;  // max energy cut (MeV)
    if (position.Mag() > R_cut) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts(double energy, TVector3 position, double R_cut = 5700.0) {
    if (energy < 1.7) return false;  // min energy cut (MeV)
    if (energy > 2.5) return false;  // max energy cut (MeV)
    if (position.Mag() > R_cut) return false;  // FV cut (mm)

    return true;
}

bool pass_coincidence_cuts(double delay, TVector3 prompt_pos, TVector3 delayed_pos) {
    // double delay = (delayed_time - prompt_time) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
    double distance = (delayed_pos - prompt_pos).Mag();

    if (delay < 200) return false;  // min delay cut (ns)
    if (delay > MAX_DELAY) return false;  // max delay cut (ns)
    if (distance > 1500) return false;  // max distance cut (mm)

    return true;
}



/* ~~~~~~~~~~~~~~~~~~~~~~ PRIMARY FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~ */


void print_info_to_file(const std::vector<std::string>& fileNames, const std::string output_file_address, unsigned int starting_fileNum, bool verbose, std::string fitName) {
    if (verbose) {std::cout << "Running print_info_to_file()" << std::endl;}

    /*********** Open output file, to print info to ***********/
    std::ofstream output_file;
    output_file.open(output_file_address);

    /*********** Setup hists ***********/
    double R_cut = 5700.0;
    TH1D* hist_t_res = new TH1D("t_res", "t_res", 1300, -300.5, 999.5);
    TH1D* hist_prompt_E = new TH1D("prompt_E", "prompt_E", 100, 0.9, 3.5);
    TH1D* hist_delayed_E = new TH1D("delayed_E", "delayed_E", 100, 1.7, 2.5);
    TH1D* hist_prompt_R = new TH1D("prompt_R", "prompt_R", 100, 0.0, R_cut);
    TH1D* hist_delayed_R = new TH1D("delayed_R", "delayed_R", 100, 0.0, R_cut);
    TH1D* hist_deltaR = new TH1D("deltaR", "deltaR", 100, 0.0, 1500.0);
    TH1D* hist_deltaT = new TH1D("deltaT", "deltaT", 100, 0.0, MAX_DELAY);

    // RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();


    /*********** Loop through all files, entries, events, and PMTs ***********/
    // loops through files
    unsigned int entry_num = starting_fileNum;
    ULong64_t previous_50MHz_time = 0;
    for (unsigned int i = 0; i < fileNames.size(); ++i) {
        RAT::DU::DSReader dsReader(fileNames.at(i));
        fTRCalc.BeginOfRun();  // Re-initialize time residual calculator (light-path calculator) after it gets geo info from DSReader
        
        if (verbose) {std::cout << "Looping through entries..." << std::endl;}
        // loops through entries
        for (size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); ++iEntry) {
            if (verbose) {std::cout << "iEntry = " << iEntry << std::endl;}
            const RAT::DS::Entry& rDS = dsReader.GetEntry(iEntry);

            // loops through events
            previous_50MHz_time = 0;
            for (size_t iEV = 0; iEV < rDS.GetEVCount(); ++iEV) {
                if (verbose) {std::cout << "iEV = " << iEV << std::endl;}
                const RAT::DS::EV& rEV = rDS.GetEV(iEV);

                double Delta_T = 0.0;
                if (previous_50MHz_time == 0) {
                    previous_50MHz_time = rEV.GetClockCount50();
                } else {
                    // Compute Delta_T in ns (from last event in entry). The "& 0x7FFFFFFFFFF" deals with the 50MHz clock roll-over,
                    // by doing a bitwise comparison with the largest 8-bit integer
                    Delta_T = ((int64_t(rEV.GetClockCount50()) - int64_t(previous_50MHz_time)) & 0x7FFFFFFFFFF) * 20.0;
                }

                // Get recon info, and (if it exists) print to file
                std::vector<double> info_recon = get_recon_info(rEV, fTRCalc, Nbins, lower_lim, upper_lim, fitName);
                if (info_recon.size() > 0) {
                    output_file  << "recon " << entry_num << " " << iEV << " " << info_recon.at(0) << " "
                    << info_recon.at(1) << " " << info_recon.at(2) << "," << info_recon.at(3) << "," << info_recon.at(4)
                    << " " << Delta_T << " " << rEV.GetNhitsCleaned() << " "; // MC entry evt alphaN_classifier_result E x,y,z Delta_T Nhits t_res_1,t_res_2,...,t_res_n
                    for (unsigned int j = 5; j < (info_recon.size() - 1); ++j) {
                        output_file << info_recon.at(j) << ",";
                    }
                    output_file << info_recon.at(info_recon.size() - 1) << std::endl;
                }
            }
            ++entry_num;
        }
        dsReader.Delete();
    }
}


/* ~~~~~~~~~~~~~~~~~~~~~~ GET INFO FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~ */


/**
 * @brief Returns vector of recon info for event (including entries of bins from time residual histogram):
 * {alphaN_classifier_result, E, x, y, z, t_res_0, t_res_1, ..., t_res_n}
 * 
 * @param evt 
 * @param fTRCalc 
 * @param Nbins 
 * @param lower_lim 
 * @param upper_lim 
 * @param fitname 
 * @return std::vector<double> 
 */
std::vector<double> get_recon_info(const RAT::DS::EV& evt, RAT::DU::TimeResidualCalculator fTRCalc, unsigned int Nbins,
                                    double lower_lim, double upper_lim, std::string fitName) {

    // Get event info (grab the fit information)
    std::vector<double> output;
    if (fitName == "")
        fitName = evt.GetDefaultFitName();
    try {
        // Get recon info
        const RAT::DS::FitVertex& rVertex = evt.GetFitResult(fitName).GetVertex(0);
        if (!(rVertex.ValidPosition() && rVertex.ValidTime() && rVertex.ValidEnergy())) {return output;} // fit invalid
        double energy = rVertex.GetEnergy();
        RAT::DU::Point3D pos(0, rVertex.GetPosition());  // position of fibre [mm] (as Point3D in PSUP coordinates, see system_id in POINT3D_SHIFTS tables)
        double vertex_time = rVertex.GetTime();

        // Get classifier info
        double alphaNreactor_classier_result;
        if (!evt.ClassifierResultExists("AlphaNReactorIBDClassifier")) {
            std::cout << "No AlphaNReactorIBDClassifier results." << std::endl;
            alphaNreactor_classier_result = -999.;
        } else if (!evt.GetClassifierResult("AlphaNReactorIBDClassifier").GetValid()) {
            std::cout << "No valid AlphaNReactorIBDClassifier result." << std::endl;
            alphaNreactor_classier_result = -999.;
        } else {
            RAT::DS::ClassifierResult alphaNreactor_result = evt.GetClassifierResult("AlphaNReactorIBDClassifier");
            alphaNreactor_classier_result = alphaNreactor_result.GetClassification("AlphaNReactorIBDClassifier");
        }

        // Package info
        output = {alphaNreactor_classier_result, energy, pos.X(), pos.Y(), pos.Z()};

        // calculate time residuals (loop through PMTs) and dump them in a histogram
        const RAT::DS::CalPMTs& calibratedPMTs = evt.GetCalPMTs();
        TH1D* hist = new TH1D("name", "title", Nbins, lower_lim, upper_lim);
        for (size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++) {
            const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT(iPMT);
            // Use new time residual calculator
            hist->Fill(fTRCalc.CalcTimeResidual(pmtCal, pos, vertex_time));
        }

        // Now loop through histogram bins, and save entries in output vector
        unsigned int N_bins = hist->GetNbinsX();
        for (unsigned int i = 0; i < N_bins+2; ++i) { // (first (0) and last (nbins+1) bins are overflow bins)
            output.push_back(hist->GetBinContent(i));
        }
        delete hist;
    }
    catch (const RAT::DS::DataNotFound&) {return output;}  // no fit data
    catch (const RAT::DS::FitCollection::NoResultError&) {return output;} // no fit result by the name of fitName
    catch (const RAT::DS::FitResult::NoVertexError&) {return output;} // no fit vertex
    catch (const RAT::DS::FitVertex::NoValueError&) {return output;} // position or time missing
    catch (const RAT::DS::ClassifierResult::NoClassificationError&) {return output;} // classifier result error

    return output;
}


// bin = 0;       underflow bin
// bin = 1;       first bin with low-edge xlow INCLUDED
// bin = nbins;   last bin with upper-edge xup EXCLUDED
// bin = nbins+1; overflow bin