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
void print_info_to_file(const std::vector<std::string>& fileNames, const std::string output_file_address, unsigned int starting_fileNum, bool verbose, std::string fitName) {
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

    RAT::DB::Get()->SetAirplaneModeStatus(true);
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
                    // Compute Delta_T in ns (from last event in entry). The "& 0x7FFFFFFFFFF" just removes the sign
                    // (so we get the absolute value), by doing a bitwise comparison with the largest 8-bit integer
                    Delta_T = ((int64_t(rEV.GetClockCount50()) - int64_t(previous_50MHz_time)) & 0x7FFFFFFFFFF) * 20.0;
                }

                // // Get MC info, and (if it exists) print to file
                // std::vector<double> info_MC = get_MC_info(rDS, rEV, iEV, fTRCalc, Nbins, lower_lim, upper_lim);
                // if (info_MC.size() > 0) {
                //     output_file  << "MC " << entry_num << " " << iEV << " " << info_MC.at(0) << " " << info_MC.at(1) << " "
                //     << info_MC.at(2) << " " << info_MC.at(3) << "," << info_MC.at(4) << "," << info_MC.at(5)
                //     << " " << Delta_T << " " << rEV.GetNhitsCleaned() << " "; // MC entry evt alphaN_classifier_result KE PDG x,y,z Delta_t Nhits t_res_1,t_res_2,...,t_res_n
                //     for (unsigned int j = 6; j < (info_MC.size() - 1); ++j) {
                //         output_file << info_MC.at(j) << ",";
                //     }
                //     output_file << info_MC.at(info_MC.size() - 1) << std::endl;
                // }

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
 * @brief Returns vector of MC info for event (including entries of bins from time residual histogram):
 * {alphaN_classifier_result, KE, PDG, x, y, z, t_res_0, t_res_1, ..., t_res_n}
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

    std::vector<double> output;
    if ((evt_idx < entry.GetMC().GetMCParticleCount()) && evt_idx < entry.GetMCEVCount()) {
        // Get event info
        double KE = entry.GetMC().GetMCParticle(evt_idx).GetKineticEnergy();
        double PDG_code = entry.GetMC().GetMCParticle(evt_idx).GetPDGCode();
        RAT::DU::Point3D pos(0, entry.GetMC().GetMCParticle(evt_idx).GetPosition());  // position of fibre [mm] (as Point3D in PSUP coordinates, see system_id in POINT3D_SHIFTS tables)
        double GT_time = 390 - entry.GetMCEV(evt_idx).GetGTTime();  // event time is 390ns - GT time.

        try {
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

            // Package output info
            output = {alphaNreactor_classier_result, KE, PDG_code, pos.X(), pos.Y(), pos.Z()};

            // calculate time residuals (loop through PMTs) and dump them in a histogram
            const RAT::DS::CalPMTs& calibratedPMTs = evt.GetCalPMTs();
            TH1D* hist = new TH1D("name", "title", Nbins, lower_lim, upper_lim);
            for (size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++) {
                const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT(iPMT);
                // Use new time residual calculator
                hist->Fill(fTRCalc.CalcTimeResidual(pmtCal, pos, GT_time));
            }

            // Now loop through histogram bins, and save entries in output vector
            unsigned int N_bins = hist->GetNbinsX();
            for (unsigned int i = 0; i < N_bins+2; ++i) { // (first and last bins are overflow bins)
                output.push_back(hist->GetBinContent(i));
            }
            delete hist;
        }
        catch (const RAT::DS::ClassifierResult::NoClassificationError&) {return output;} // classifier result error
    }


    return output;
}

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