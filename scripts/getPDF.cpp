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
TH1D* HitTimeResiduals(std::string type, const std::vector<std::string>& fileNames, std::vector<double> delays, std::string cuts_name, bool verbose, std::string fitName = "");
std::vector<unsigned int> apply_cuts(std::string cuts_name, unsigned int evt_idx, const RAT::DS::Entry& entry, const RAT::DS::EV& evt, std::string fitName = "");
std::vector<unsigned int> alphaN_cuts(std::string cuts_name, unsigned int evt_idx, const RAT::DS::Entry& entry, const RAT::DS::EV& evt, std::string fitName = "");
std::vector<double> getReconInfo(unsigned int evt_idx, const RAT::DS::Entry& entry, std::string fitName = "");
std::vector<double> getMCInfo(unsigned int evt_idx, const RAT::DS::Entry& entry);



int main(int argc, char** argv) {
    std::string output_file = argv[1];
    std::string cuts_name = argv[2];
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
    TH1D* MC_summed_hist = HitTimeResiduals("true", input_files, delays, cuts_name, verbose);
    TH1D* Fitted_summed_hist = HitTimeResiduals("recon", input_files, delays, cuts_name, verbose);

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
 * @brief Make hit time residual histogram (summed over all events).
 * 
 * @param type = "true" (MC info), or "recon" (Fitted info).
 * @param fileNames list of simulation output root files.
 * @param delays vector with list of e+ decay times.
 * @param cuts_name name of cuts method to be applied. Example "alphaN", if unknown cuts_name, defaults to all events pass cuts
 * @param verbose
 * @param fitName
 * @return TH1D*
 */
TH1D* HitTimeResiduals(std::string type, const std::vector<std::string>& fileNames, std::vector<double> delays, std::string cuts_name, bool verbose, std::string fitName) {
    if (verbose) {std::cout << "Running HitTimeResidualsFitPosition()" << std::endl;}

    if (type != "true" && type != "recon") {
        std::cout << "info_type should be 'true' or 'recon', not '" << type << "'." << std::endl;
        exit(1);
    }

    std::string hist_name = "hHitTimeResiduals_" + type;
    std::string hist_title = "Hit time residuals using the " + type + " position";
    TH1D* histTimeResiduals = new TH1D(hist_name.c_str(), hist_title.c_str(), 1300, -300.5, 999.5);
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
        for (size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++) {
            const RAT::DS::Entry& rDS = dsReader.GetEntry(iEntry);
            for (size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++) {
                const RAT::DS::EV& rEV = rDS.GetEV(iEV);

                // Get indices of events that pass cuts
                std::vector<unsigned int> passed_evt_indices = apply_cuts(cuts_name, iEV, rDS, rEV, fitName);
                // Add passed events to histograms
                for (unsigned int j = 0; j < passed_evt_indices.size(); ++j) {
                    const RAT::DS::EV& passed_evt = rDS.GetEV(passed_evt_indices.at(j));

                    // Get event info
                    TVector3 evt_pos;
                    double evt_time;
                    if (type == 'recon') {
                        std::vector<double> evt_recon_info = getReconInfo(iEV, rDS, fitName);
                        evt_pos = TVector3(evt_recon_info.at(1), evt_recon_info.at(2), evt_recon_info.at(3));
                        evt_time = evt_recon_info.at(4);
                    } else {
                        std::vector<double> evt_true_info = getMCInfo(iEV, rDS);
                        evt_pos = TVector3(evt_true_info.at(1), evt_true_info.at(2), evt_true_info.at(3));
                        evt_time = evt_true_info.at(4);
                    }

                    // calculate time residuals
                    const RAT::DS::CalPMTs& calibratedPMTs = passed_evt.GetCalPMTs();
                    if (verbose) {std::cout << "Adding event to histogram." << std::endl;}
                    // Use new time residual calculator
                    for (size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++) {
                        const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT(iPMT);
                        histTimeResiduals->Fill(fTRCalc.CalcTimeResidual(pmtCal, evt_pos, evt_time));
                    }
                    if (verbose) {std::cout << "...added." << std::endl;}
                    ++num_evts;
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



/* ~~~~~~~~~~~~~~~~~~~~~~ CUT FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~ */


/**
 * @brief Returns indices of events that pass cut, for a specific entry.
 * (In the alpha-n case only prompt events from combined prompt-delayed events that pass cuts will be returned).
 * 
 * @param cuts_name 
 * @param evt_idx 
 * @param entry 
 * @param evt 
 * @param info_type 
 * @param fitName 
 * @return std::vector<unsigned int> 
 */
std::vector<unsigned int> apply_cuts(std::string cuts_name, unsigned int evt_idx, const RAT::DS::Entry& entry, const RAT::DS::EV& evt, std::string info_type, std::string fitName) {
    if (info_type != "true" && info_type != "recon") {
        std::cout << "info_type should be 'true' or 'recon', not '" << info_type << "'." << std::endl;
        exit(1);
    }
    if (cuts_name == "alphaN") {
        return alphaN_cuts(cuts_name, evt_idx, entry, evt, info_type, fitName);
    } else {
        std::vector<unsigned int> automatic_pass = {evt_idx};
        return automatic_pass;
    }
}

/**
 * @brief Applies alphaN cut. If delayed event is inputted, and both prompt and delayed events pass cuts,
 * the prompt event index only will be returned.
 * 
 * @param cuts_name 
 * @param evt_idx 
 * @param entry 
 * @param evt 
 * @return std::vector<unsigned int> 
 */
std::vector<unsigned int> alphaN_cuts(std::string cuts_name, unsigned int evt_idx, const RAT::DS::Entry& entry, const RAT::DS::EV& evt, std::string info_type, std::string fitName) {
    std::vector<unsigned int> passed_evt_indices;

    // Only pass the prompt event this is the delayed event, so we know it has a delayed event
    if (evt_idx != 1) {return passed_evt_indices;}

    // Get relevent info for cuts
        // For prompt event
        const RAT::DS::EV& prompt_evt = entry.GetEV(0);
        // unsigned int prompt_nhits = prompt_evt.GetNhitsCleaned();
        double prompt_time = prompt_evt.GetClockCount50() * 20.0;  // nanoseconds
        std::vector<double> prompt_info;
        if (info_type == "recon") {
            prompt_info =  getReconInfo(0, entry, fitName);
        } else {
            prompt_info =  getMCInfo(0, entry);  // FIXME: Add energy to event energy if e+? 1.02MeV
        }
        if (prompt_info.size() != 4) {return passed_evt_indices;}
        double prompt_energy = prompt_info.at(0);
        TVector3 prompt_pos = TVector3(prompt_info.at(1), prompt_info.at(2), prompt_info.at(3));
        double prompt_R = prompt_pos.Mag();

        // For delayed event
        // unsigned int delayed_nhits = delayed_evt.GetNhitsCleaned();
        double delayed_time = evt.GetClockCount50() * 20.0;  // nanoseconds
        std::vector<double> delayed_info;
        if (info_type == "recon") {
            delayed_info =  getReconInfo(evt_idx, entry, fitName);
        } else {
            delayed_info =  getMCInfo(evt_idx, entry);  // FIXME: Add energy to event energy if e+? 1.02MeV
        }
        if (delayed_info.size() != 4) {return passed_evt_indices;}
        double delayed_energy = delayed_info.at(0);
        TVector3 delayed_pos = TVector3(delayed_info.at(1), delayed_info.at(2), delayed_info.at(3));
        double delayed_R = delayed_pos.Mag();

        // Combined info
        double deltaT = delayed_time - prompt_time;
        double delta_R = (delayed_pos - prompt_pos).Mag();

    // Apply cuts
    if (prompt_energy < 0.9 or prompt_energy > 8.0) {return passed_evt_indices;}
    if (delayed_energy < 1.85 or delayed_energy > 2.4) {return passed_evt_indices;}
    if (prompt_R > 5700 or delayed_R > 5700) {return passed_evt_indices;}  // in mm
    if (delta_R > 1500) {return passed_evt_indices;}  // in mm
    if (delayed_time < 400 or delayed_time > 0.8E6) {return passed_evt_indices;}  // in ns
    if (prompt_energy > 3.5) {return passed_evt_indices;}

    // Passed all cuts, return prompt event index
    passed_evt_indices = {0};
    return passed_evt_indices;
}

std::vector<unsigned int> oPs_cuts(std::string cuts_name, unsigned int evt_idx, const RAT::DS::Entry& entry, const RAT::DS::EV& evt, std::string info_type, std::string fitName) {
    // if (verbose) {std::cout << "evt_idx = " << evt_idx << std::endl;}
    // if (evt_idx >= delays.size()) {
    //     if (verbose) {std::cout << "evt_idx out of range of delays vector" << std::endl;}
    //     ++evt_idx;
    //     continue;
    // } else if (is_oPs && delays.at(evt_idx) == 0.0) {  // Filter out non o-Ps events
    //     if (verbose) {std::cout << "Delay = 0, ignoring event." << std::endl;}
    //     ++evt_idx;
    //     continue;
    // } else if (eventPosition.Mag() > vol_cut) {
    //     if (verbose) {std::cout << "Outside volume cut, ignoring event." << std::endl;}
    //     ++evt_idx;
    //     continue;
    // } else if (nhits < 200 || nhits > 6000) {
    //     if (verbose) {std::cout << "Outside nhit cut, ignoring event." << std::endl;}
    //     ++evt_idx;
    //     continue;
}


/**
 * @brief Get the reconstructed energy and position of an event as {E, x, y, z, t}.
 * If error in reconstruction, returns empty vector.
 * 
 * @param evt 
 * @return std::vector<double> 
 */
std::vector<double> getReconInfo(unsigned int evt_idx, const RAT::DS::Entry& entry, std::string fitName) {
    const RAT::DS::EV& evt = entry.GetEV(evt_idx);
    std::vector<double> output_vector;
    // grab the fit information
    if(fitName == "")
        fitName = evt.GetDefaultFitName();
    try {
        const RAT::DS::FitVertex& rVertex = evt.GetFitResult(fitName).GetVertex(0);
        if (!(rVertex.ValidPosition() && rVertex.ValidTime() && rVertex.ValidEnergy())) {return output_vector;} // fit invalid
        double energy = rVertex.GetEnergy();
        TVector3 pos = rVertex.GetPosition();
        double time = rVertex.GetTime();
        output_vector = {energy, pos.X(), pos.Y(), pos.Z(), time};
    }
    catch (const RAT::DS::FitCollection::NoResultError&) {return output_vector;} // no fit result by the name of fitName
    catch (const RAT::DS::FitResult::NoVertexError&) {return output_vector;} // no fit vertex
    catch (const RAT::DS::FitVertex::NoValueError&) {return output_vector;} // position or time missing

    return output_vector;
}

/**
 * @brief Get the MC energy and position of an event as {E, x, y, z, t}.
 * 
 * @param evt 
 * @return std::vector<double> 
 */
std::vector<double> getMCInfo(unsigned int evt_idx, const RAT::DS::Entry& entry) {
    double energy = entry.GetMC().GetMCParticle(evt_idx).GetKineticEnergy();
    if (entry.GetMC().GetMCParticle(evt_idx).GetPDGCode() == -11) {  // If it's a e+ (See PDG encoding in geant4.10.00.p04/source/particle/leptons/src/G4Positron.cc)
        energy += 1.02;  // Add e+e- annihilation energy to kinetic energy to get event energy.
    }
    const TVector3 pos = entry.GetMC().GetMCParticle(evt_idx).GetPosition();
    double time = 390 - entry.GetMCEV(evt_idx).GetGTTime();  // event time is 390ns - GT time.

    std::vector<double> output_vector = {energy, pos.X(), pos.Y(), pos.Z(), time};
    return output_vector;
}
