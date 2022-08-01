////////////////////////////////////////////////////////////////////
/// \file BASED ON: PlotHitTimeResiduals.cc
///
/// \brief Functions to get o-Ps residual hit time pdf.
///
/// \author James Page <j.page@sussex.ac.uk>
///
/// REVISION HISTORY:\n
///
/// \details EV Calibrated hit times are plotted minus transit times
/// based on the MC position or the fitted position.
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

#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFile.h>
#include <TVectorD.h>

#include <string>
#include <fstream>

std::vector<std::vector<double> > findPositronDelays(const std::vector<std::string>& filenames, bool verbose);
void printPDF(const std::string& output_filename, TH1D* MC_hist, TH1D* Fitted_hist);
TH1D* HitTimeResidualsMCPosition(const std::vector<std::string>& fileNames, std::vector<double> delays, double vol_cut, bool is_oPs, bool verbose);
TH1D* HitTimeResidualsFitPosition( const std::vector<std::string>& fileNames, std::vector<double> delays, double vol_cut, bool is_oPs, bool verbose, std::string fitName = "");

int main(int argc, char** argv) {
    std::string output_file = argv[1];
    bool is_oPs = std::stoi(argv[2]);
    bool verbose = std::stoi(argv[3]);
    // Addresses of sim files to be analysed
    std::vector<std::string> input_files;
    for (unsigned int i = 4; i < argc; ++i) {
        input_files.push_back(argv[i]);
    }

    // Get e+ delays
    std::vector<std::vector<double> > delays_vecs = findPositronDelays(input_files, verbose);
    std::vector<double> delays = delays_vecs.at(0);
    std::vector<double> delays_sansZero = delays_vecs.at(1);

    // Create time residual histograms (copied from rat/example/root/PlotHitTimeResiduals.cc)
    if (verbose) {std::cout << "Getting hists..." << std::endl;}
    TH1D* MC_summed_hist = HitTimeResidualsMCPosition(input_files, delays, 5700, is_oPs, verbose);
    TH1D* Fitted_summed_hist = HitTimeResidualsFitPosition(input_files, delays, 5700, is_oPs, verbose);

    // Get o-Ps pdf and print to file
    printPDF(output_file, MC_summed_hist, Fitted_summed_hist);

    // Create output root file if recordeing extra info
    if (verbose) {
        // Creat file name
        std::cout << "Creating output histogram file name" << std::endl;
        std::size_t botDirPos = output_file.find_last_of("/");
        std::string filename = output_file.substr(botDirPos+1, output_file.length() - 5);
        std::string path = output_file.substr(0, botDirPos+1);
        std::string saveroot = path + "Hists_" + filename + ".root";

        // Save root file
        TFile *rootfile = new TFile(saveroot.c_str(),"RECREATE");

        // Now write everything to the file and close
        if (verbose) {std::cout << "Writing everything to file and closing" << std::endl;}
        rootfile->cd();
        MC_summed_hist->Write();
        Fitted_summed_hist->Write();

        // raw_hist->Write();
        rootfile->Write();
        rootfile->Close();
    }

    return 0;
}

/**
 * @brief Returns a list of the delays imparted to positron decays (emulating oPs)
 * 
 * @param filenames
 * @param verbose
 * @return std::vector<double> 
 */
std::vector<std::vector<double> > findPositronDelays(const std::vector<std::string>& filenames, bool verbose) {
    if (verbose) {std::cout << "Finding e+ delays..." << std::endl;}

    std::vector<double> delays;
    std::vector<double> delays_sansZero;
    for (unsigned int i = 0; i < filenames.size(); ++i) {
        RAT::DU::DSReader dsReader(filenames.at(i));

        // Loop through events. Each one should have one primary track (first child) in MC
        if (verbose) {std::cout << "Looping through events..." << std::endl;}
        for (size_t iEv =0; iEv<dsReader.GetEntryCount(); iEv++) {
            const RAT::DS::Entry& rDS = dsReader.GetEntry(iEv);
            RAT::TrackNav nav(&rDS);
            RAT::TrackCursor cursor = nav.Cursor(false);

            // Check there is an event in this entry (won't get associated t_res plot if not)
            if (rDS.GetEVCount() > 0) {
                // Should only go through this loop once in MC.
                for (size_t iCh = 0; iCh<(size_t)cursor.ChildCount(); iCh++) {
                    cursor.GoChild(iCh);

                    // Go to the end of the e+ track
                    cursor.GoTrackEnd();
                    RAT::TrackNode* parent_node = cursor.Here();
                    double start_time = parent_node->GetGlobalTime();
                    if (verbose) {std::cout << "Parent particle: " << parent_node->GetParticleName() << std::endl;}
                    if (verbose) {std::cout << "Last step process: " << parent_node->GetProcess() << std::endl;}

                    // Go to the start of the first child track (gamma)
                    cursor.GoChild(0);
                    RAT::TrackNode* child_node = cursor.Here();
                    double end_time = child_node->GetGlobalTime();
                    delays.push_back(end_time - start_time);
                    if (end_time - start_time) {
                        delays_sansZero.push_back(end_time - start_time);
                    }
                    if (verbose) {std::cout << "Child particle: " << child_node->GetParticleName() << std::endl;}
                    if (verbose) {std::cout << "e+ delay: " << end_time - start_time << std::endl;}

                    // Go back to e+ track
                    cursor.GoParent();

                    // Go back to parent node to redo loop
                    cursor.GoTrackStart();
                    cursor.GoParent();
                } //Primary Particle Tracks
            }
        } //event
        dsReader.Delete();
    }

    if (verbose) {std::cout << "Num delays: " << delays.size() << std::endl;}
    std::vector<std::vector<double> > results = {delays, delays_sansZero};
    return results;
}

/**
 * @brief Print pdf to text file based on inputted histogram.
 * 
 * @param output_filename
 * @param MC_hist 
 * @param Fitted_hist 
 */
void printPDF(const std::string& output_filename, TH1D* MC_hist, TH1D* Fitted_hist) {

    // Open file to print results to
    std::ofstream datafile;
    datafile.open(output_filename, std::ofstream::out | std::ofstream::trunc);
    datafile << "[";

    double bin_width = MC_hist->GetBinCenter(2) - MC_hist->GetBinCenter(1);
    unsigned int N_bins = MC_hist->GetNbinsX();

    // Checks
    if (bin_width != Fitted_hist->GetBinCenter(2) - Fitted_hist->GetBinCenter(1)) {
        std::cout << "Bin widths different" << std::endl;
        exit(1);
    }
    if (N_bins != Fitted_hist->GetNbinsX()) {
        std::cout << "Number of bins different" << std::endl;
        exit(1);
    }

    double tot_MC_hits = 0.0;
    double tot_Fitted_hits = 0.0;
    std::vector<double> times;
    std::vector<double> MC_probs;
    std::vector<double> Fitted_probs;
    // Get histogram bin values
    for (unsigned int i = 1; i < N_bins+1; ++i) { //loop over histogram bins
        times.push_back(MC_hist->GetBinCenter(i));
        MC_probs.push_back(MC_hist->GetBinContent(i));
        tot_MC_hits += MC_probs.at(i-1);
        Fitted_probs.push_back(Fitted_hist->GetBinContent(i));
        tot_Fitted_hits += Fitted_probs.at(i-1);

        // print times to file while at it
        datafile << times.at(i-1);
        if (i < N_bins) {
            datafile << ", ";
        }
    }
    datafile << "]" << std::endl;

    std::cout << "Total number of MC PMT hits = " << tot_MC_hits << std::endl;
    std::cout << "Total number of Fitted PMT hits = " << tot_Fitted_hits << std::endl;

    // Normalise MC to pdf
    datafile << "[";
    for (unsigned int i = 0; i < MC_probs.size(); ++i) {
        MC_probs.at(i) /= (tot_MC_hits * bin_width);
        datafile << MC_probs.at(i);
        if (i < MC_probs.size() - 1) {
            datafile << ", ";
        }
    }
    datafile << "]" << std::endl;

    // Normalise Fitted to pdf
    datafile << "[";
    for (unsigned int i = 0; i < Fitted_probs.size(); ++i) {
        Fitted_probs.at(i) /= (tot_Fitted_hits * bin_width);
        datafile << Fitted_probs.at(i);
        if (i < Fitted_probs.size() - 1) {
            datafile << ", ";
        }
    }
    datafile << "]" << std::endl;

    datafile.close();
}

/**
 * @brief Make hit time residual histogram using MC position and event time (summed over all events).
 * 
 * @param fileNames Simulation root file names.
 * @param delays Vector of o-Ps decay times.
 * @param vol_cut Volume cut, given by radius in mm.
 * @param is_oPs Are events we wish to look at only o-Ps events? If so, ignore delay=0 events (some get simulated anyway for some reason, so cut them out here).
 * @param verbose
 */
TH1D* HitTimeResidualsMCPosition(const std::vector<std::string>& fileNames, std::vector<double> delays, double vol_cut, bool is_oPs, bool verbose) {
    if (verbose) {std::cout << "Running HitTimeResidualsMCPosition()" << std::endl;}

    TH1D* histTimeResiduals = new TH1D("pdfTimeResidualsMC", "PDF for Hit time residuals using the MC position", 1300, -300.5, 999.5);
    // If this is being done on data that does not require remote database connection
    // eg.: a simple simulation with default run number (0)
    // We can disable the remote connections:
    //
    // NOTE: Don't do this if you are using real data!!!
    RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();

    unsigned int num_evts = 0;
    for (unsigned int i = 0; i < fileNames.size(); ++i) {
        RAT::DU::DSReader dsReader(fileNames.at(i));

        if (verbose) {std::cout << "Looping through entries..." << std::endl;}
        unsigned int evt_idx = 0;
        for(size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++) {
            const RAT::DS::Entry& rDS = dsReader.GetEntry(iEntry);
            const TVector3 eventPosition = rDS.GetMC().GetMCParticle(0).GetPosition(); // At least 1 is somewhat guaranteed
            for(size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++) {
                const RAT::DS::EV& rEV = rDS.GetEV(iEV);
                const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
                if (verbose) {std::cout << "evt_idx = " << evt_idx << std::endl;}
                if (evt_idx >= delays.size()) {
                    if (verbose) {std::cout << "evt_idx out of range of delays vector" << std::endl;}
                    ++evt_idx;
                    continue;
                } else if (is_oPs && delays.at(evt_idx) == 0.0) {  // Filter out non o-Ps events
                    if (verbose) {std::cout << "Delay = 0, ignoring event." << std::endl;}
                    ++evt_idx;
                    continue;
                } else if (eventPosition.Mag() > vol_cut) {
                    if (verbose) {std::cout << "Outside volume cut, ignoring event." << std::endl;}
                    ++evt_idx;
                    continue;
                } else {
                    if (verbose) {std::cout << "Adding event to histogram." << std::endl;}
                    // Use new time residual calculator
                    for(size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++) {
                        const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT(iPMT);
                        histTimeResiduals->Fill(fTRCalc.CalcTimeResidual(pmtCal, eventPosition, 390 - rDS.GetMCEV(iEV).GetGTTime()));  // event time is 390ns - GT time.
                    }
                    if (verbose) {std::cout << "...added." << std::endl;}
                    ++num_evts;
                    ++evt_idx;
                }
            }
        }
        dsReader.Delete();
    }
    std::cout << "Number of events recorded = " << num_evts << std::endl;

    histTimeResiduals->GetYaxis()->SetTitle( "Count per 1 ns bin" );
    histTimeResiduals->GetXaxis()->SetTitle( "Hit time residuals [ns]" );
    histTimeResiduals->Draw();
    return histTimeResiduals;
}


/**
 * @brief Make hit time residual histogram using Fit position and event time (summed over all events).
 * 
 * @param fileNames Simulation root file .
 * @param delays Vector of o-Ps decay times.
 * @param vol_cut Volume cut, given by radius in mm.
 * @param is_oPs Are events we wish to look at only o-Ps events? If so, ignore delay=0 events (some get simulated anyway for some reason, so cut them out here).
 * @param verbose
 * @param fitName Which fitter to use if not default.
 */
TH1D* HitTimeResidualsFitPosition( const std::vector<std::string>& fileNames, std::vector<double> delays, double vol_cut, bool is_oPs, bool verbose, std::string fitName) {
    if (verbose) {std::cout << "Running HitTimeResidualsFitPosition()" << std::endl;}

    TH1D* histTimeResiduals = new TH1D("hHitTimeResidualsFit", "Hit time residuals using the Fit position", 1300, -300.5, 999.5);
    // If this is being done on data that does not require remote database connection
    // eg.: a simple simulation with default run number (0)
    // We can disable the remote connections:
    //
    // NOTE: Don't do this if you are using real data!!!
    RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();

    unsigned int num_evts = 0;
    for (unsigned int i = 0; i < fileNames.size(); ++i) {
        RAT::DU::DSReader dsReader(fileNames.at(i));

        if (verbose) {std::cout << "Looping through entries..." << std::endl;}
        unsigned int evt_idx = 0;
        for (size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++) {
            const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
            for (size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++) {
                const RAT::DS::EV& rEV = rDS.GetEV( iEV );

                // grab the fit information
                if(fitName == "")
                    fitName = rEV.GetDefaultFitName();

                TVector3 eventPosition;
                double   eventTime;

                try{
                    const RAT::DS::FitVertex& rVertex = rEV.GetFitResult(fitName).GetVertex(0);
                    if(!(rVertex.ValidPosition() && rVertex.ValidTime()))
                        continue; // fit invalid

                    eventPosition = rVertex.GetPosition();
                    eventTime = rVertex.GetTime();
                }
                catch(const RAT::DS::FitCollection::NoResultError&){
                    // no fit result by the name of fitName
                    continue;
                }
                catch (const RAT::DS::FitResult::NoVertexError&){
                    // no fit vertex
                    continue;
                }
                catch(const RAT::DS::FitVertex::NoValueError&){
                    // position or time missing
                    continue;
                }
                // DataNotFound --> implies no fit results are present, don't catch.

                // calculate time residuals
                const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
                unsigned int nhits = rEV.GetNhitsCleaned();
                if (verbose) {std::cout << "evt_idx = " << evt_idx << std::endl;}
                if (evt_idx >= delays.size()) {
                    if (verbose) {std::cout << "evt_idx out of range of delays vector" << std::endl;}
                    ++evt_idx;
                    continue;
                } else if (is_oPs && delays.at(evt_idx) == 0.0) {  // Filter out non o-Ps events
                    if (verbose) {std::cout << "Delay = 0, ignoring event." << std::endl;}
                    ++evt_idx;
                    continue;
                } else if (eventPosition.Mag() > vol_cut) {
                    if (verbose) {std::cout << "Outside volume cut, ignoring event." << std::endl;}
                    ++evt_idx;
                    continue;
                } else if (nhits < 200 || nhits > 6000) {
                    if (verbose) {std::cout << "Outside nhit cut, ignoring event." << std::endl;}
                    ++evt_idx;
                    continue;
                } else {
                    if (verbose) {std::cout << "Adding event to histogram." << std::endl;}
                    // Use new time residual calculator
                    for (size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++) {
                        const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT(iPMT);
                        histTimeResiduals->Fill(fTRCalc.CalcTimeResidual(pmtCal, eventPosition, eventTime));
                        
                    }
                }
                if (verbose) {std::cout << "...added." << std::endl;}
                ++num_evts;
                ++evt_idx;
            }
        }
        dsReader.Delete();
    }
    std::cout << "Number of events recorded = " << num_evts << std::endl;

    histTimeResiduals->GetYaxis()->SetTitle( "Count per 1 ns bin" );
    histTimeResiduals->GetXaxis()->SetTitle( "Hit time residuals [ns]" );
    histTimeResiduals->Draw();
    return histTimeResiduals;
}