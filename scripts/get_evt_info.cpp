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

void print_info_to_file(const std::vector<std::string>& fileNames, const std::string output_file_address, bool verbose, std::string fitName="");
std::vector<double> get_MC_info(const RAT::DS::Entry& entry, const RAT::DS::EV& evt, unsigned int evt_idx, RAT::DU::TimeResidualCalculator fTRCalc, unsigned int Nbins, double lower_lim, double upper_lim);
std::vector<double> get_recon_info(const RAT::DS::EV& evt, RAT::DU::TimeResidualCalculator fTRCalc, unsigned int Nbins, double lower_lim, double upper_lim, std::string fitname);



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
    print_info_to_file(input_files, output_file_address, verbose);


    return 0;
}


/* ~~~~~~~~~~~~~~~~~~~~~~ PRIMARY FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~ */

/**
 * @brief Loop through simulation output root files, and print MC and recon info for all events
 * to a text file. Format as follows:
 * 
 * entry:number                             (0 to total number of entries in all files)
 * evt:number                               (0 to total number of events in entry)
 * MC:
 * KE:number                                (kinetic energy of first particle in event track)
 * PDG:number                               (PDG code of first particle in event track)
 * pos:number,number,number                 (3-D coordinates of event: x,y,z)
 * time:number                              (event time)
 * t_res:number,number,number,...,number    (entries from time residual histogram)
 * recon:
 * E:number                                 (event energy)
 * pos:number,number,number                 (3-D coordinates of event: x,y,z)
 * time:number                              (event time)
 * t_res:number,number,number,...,number    (entries from time residual histogram)
 * 
 * This repeats for every event in an entry, and every entry. If recon failed, no recon-related lines
 * will be present for that event.
 * 
 * @param fileNames Input simulation file addresses
 * @param output_file_address Output text file address
 * @param verbose 
 * @param fitName 
 */
void print_info_to_file(const std::vector<std::string>& fileNames, const std::string output_file_address, bool verbose, std::string fitName) {
    if (verbose) {std::cout << "Running print_info_to_file()" << std::endl;}

    /*********** Open output file, to print info to ***********/
    std::ofstream output_file;
    output_file.open(output_file_address);

    /*********** Setup ***********/
    unsigned int Nbins = 1300;
    double lower_lim = -300.5;
    double upper_lim = 999.5;
    // Print these to first line of file
    output_file << Nbins << ", " << lower_lim << ", " <<  upper_lim << std::endl;
    // Create output histograms we will get PDFs from
    std::vector<TH1D*> histograms;
    for (unsigned int i = 0; i < 3; ++i) {
        std::string hist_name = "hHitTimeResiduals_" + type + "_" + to_string(i);
        std::string hist_title = "Hit time residuals using the " + type + " position" + ", " + to_string(i);
        histograms.push_back(new TH1D(hist_name.c_str(), hist_title.c_str(), 1300, -300.5, 999.5));
    }

    RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();


    /*********** Loop through all files, entries, events, and PMTs ***********/
    // loops through files
    unsigned int entry_num = 0;
    for (unsigned int i = 0; i < fileNames.size(); ++i) {
        RAT::DU::DSReader dsReader(fileNames.at(i));
        fTRCalc.BeginOfRun();  // Re-initialize time residual calculator (light-path calculator) after it gets geo info from DSReader
        
        if (verbose) {std::cout << "Looping through entries..." << std::endl;}
        // loops through entries
        for (size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); ++iEntry) {
            if (verbose) {std::cout << "iEntry = " << iEntry << std::endl;}
            const RAT::DS::Entry& rDS = dsReader.GetEntry(iEntry);
            output_file << "entry:" << entry_num << std::endl;
            ++entry_num;

            // loops through events
            for (size_t iEV = 0; iEV < rDS.GetEVCount(); ++iEV) {
                if (verbose) {std::cout << "iEV = " << iEV << std::endl;}
                const RAT::DS::EV& rEV = rDS.GetEV(iEV);
                output_file << "evt:" << iEV << std::endl;

                // Get MC info, print it to file
                std::vector<double> info_MC = get_MC_info(rDS, rEV, iEV, fTRCalc, Nbins, lower_lim, upper_lim);
                output_file << "MC:" << std::endl;
                output_file << "KE:" << info_MC.at(0) << std::endl;
                output_file << "PDG:" << info_MC.at(1) << std::endl;
                output_file << "pos:" << info_MC.at(2) << "," << info_MC.at(3) << "," << info_MC.at(4) << std::endl;
                output_file << "time:" << info_MC.at(5) << std::endl;
                output_file << "t_res:";
                for (unsigned int j = 6; j < (info_MC.size() - 1); ++j) {
                    output_file << info_MC.at(j) << ",";
                }
                output_file << info_MC.at(info_MC.size() - 1) << std::endl;

                // Get recon info, and (if it exists) print to file
                std::vector<double> info_recon = get_recon_info(rEV, fTRCalc, Nbins, lower_lim, upper_lim, fitname);
                if (info_recon.size() > 0) {
                    output_file << "recon:" << std::endl;
                    output_file << "E:" << info_recon.at(0) << std::endl;
                    output_file << "pos:" << info_recon.at(1) << "," << info_recon.at(2) << "," << info_recon.at(3) << std::endl;
                    output_file << "time:" << info_recon.at(4) << std::endl;
                    output_file << "t_res:";
                    for (unsigned int j = 5; j < (info_recon.size() - 1); ++j) {
                        output_file << info_recon.at(j) << ",";
                    }
                    output_file << info_recon.at(info_recon.size() - 1) << std::endl;
                }
            }
        }
        dsReader.Delete();
    }
}


/* ~~~~~~~~~~~~~~~~~~~~~~ GET INFO FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~ */

/**
 * @brief Returns vector of MC info for event (including entries of bins from time residual histogram):
 * {KE, PDG, x, y, z, time, t_res_0, t_res_1, ..., t_res_n}
 * 
 * @param entry 
 * @param evt 
 * @param evt_idx 
 * @param fTRCalc 
 * @param Nbins 
 * @param lower_lim 
 * @param upper_lim 
 * @return std::vector<double> 
 */
std::vector<double> get_MC_info(const RAT::DS::Entry& entry, const RAT::DS::EV& evt, unsigned int evt_idx,
                                RAT::DU::TimeResidualCalculator fTRCalc, unsigned int Nbins, double lower_lim, double upper_lim) {

    // Get event info
    double KE = entry.GetMC().GetMCParticle(evt_idx).GetKineticEnergy();
    double PDG_code = entry.GetMC().GetMCParticle(evt_idx).GetPDGCode();
    TVector3 pos = entry.GetMC().GetMCParticle(evt_idx).GetPosition();
    double time = 390 - entry.GetMCEV(evt_idx).GetGTTime();  // event time is 390ns - GT time.

    // Create output vector, with general event info first
    std::vector<double> output = {KE, PDG_code, pos.X(), pos.Y(), pos.Z(), time};

    // calculate time residuals (loop through PMTs) and dump them in a histogram
    const RAT::DS::CalPMTs& calibratedPMTs = evt.GetCalPMTs();
    TH1D* hist = new TH1D("name", "title", Nbins, lower_lim, upper_lim);
    for (size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++) {
        const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT(iPMT);
        // Use new time residual calculator
        hist->Fill(fTRCalc.CalcTimeResidual(pmtCal, evt_pos, evt_time));
    }

    // Now loop through histogram bins, and save entries in output vector
    unsigned int N_bins = hist->GetNbinsX();
    for (unsigned int i = 0; i < N_bins+1; ++i) { // (first and last bins are overflow bins)
        output.push_back(hist->GetBinContent(i));
    }

    return output;
}

/**
 * @brief Returns vector of recon info for event (including entries of bins from time residual histogram):
 * {E, x, y, z, time, t_res_0, t_res_1, ..., t_res_n}
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
                                    double lower_lim, double upper_lim, std::string fitname) {

    // Get event info (grab the fit information)
    std::vector<double> output = {}
    if(fitName == "")
        fitName = evt.GetDefaultFitName();
    try {
        // Get recon info
        const RAT::DS::FitVertex& rVertex = evt.GetFitResult(fitName).GetVertex(0);
        if (!(rVertex.ValidPosition() && rVertex.ValidTime() && rVertex.ValidEnergy())) {return output;} // fit invalid
        double energy = rVertex.GetEnergy();
        TVector3 pos = rVertex.GetPosition();
        double time = rVertex.GetTime();

        // Package info 
        output = {energy, pos.X(), pos.Y(), pos.Z(), time};
    }
    catch (const RAT::DS::FitCollection::NoResultError&) {return output;} // no fit result by the name of fitName
    catch (const RAT::DS::FitResult::NoVertexError&) {return output;} // no fit vertex
    catch (const RAT::DS::FitVertex::NoValueError&) {return output;} // position or time missing

    // calculate time residuals (loop through PMTs) and dump them in a histogram
    const RAT::DS::CalPMTs& calibratedPMTs = evt.GetCalPMTs();
    TH1D* hist = new TH1D("name", "title", Nbins, lower_lim, upper_lim);
    for (size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++) {
        const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT(iPMT);
        // Use new time residual calculator
        hist->Fill(fTRCalc.CalcTimeResidual(pmtCal, evt_pos, evt_time));
    }

    // Now loop through histogram bins, and save entries in output vector
    unsigned int N_bins = hist->GetNbinsX();
    for (unsigned int i = 0; i < N_bins+1; ++i) { // (first and last bins are overflow bins)
        output.push_back(hist->GetBinContent(i));
    }

    return output;
}


// bin = 0;       underflow bin
// bin = 1;       first bin with low-edge xlow INCLUDED
// bin = nbins;   last bin with upper-edge xup EXCLUDED
// bin = nbins+1; overflow bin