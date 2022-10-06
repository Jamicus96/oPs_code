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

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFile.h>
#include <TVectorD.h>

#include <string>
#include <fstream>

std::vector<std::vector<double> > findPositronDelays(const std::vector<std::string>& filenames, bool verbose);
void printPDF(const std::string& output_filename, std::vector<TH1D*> MC_hist, std::vector<TH1D*> Fitted_hist, bool verbose);
std::vector<TH1D*> HitTimeResiduals(std::string type, const std::vector<std::string>& fileNames, std::vector<double> delays, bool verbose, std::string fitName = "");
std::vector<unsigned int> apply_cuts(std::string cuts_name, unsigned int evt_idx, const RAT::DS::Entry& entry, const RAT::DS::EV& evt, std::string info_type, TH2F* alpha_scaling_hist, const std::vector<unsigned int>& Nbins, const std::vector<double>& Limits, std::string fitName = "");
std::vector<unsigned int> alphaN_cuts(std::string cuts_name, unsigned int evt_idx, const RAT::DS::Entry& entry, const RAT::DS::EV& evt, std::string info_type, TH2F* alpha_scaling_hist, const std::vector<unsigned int>& Nbins, const std::vector<double>& Limits, std::string fitName = "");
std::vector<double> getReconInfo(unsigned int evt_idx, const RAT::DS::Entry& entry, TH2F* alpha_scaling_hist, const std::vector<unsigned int>& Nbins, const std::vector<double>& Limits, std::string fitName = "");
std::vector<double> getMCInfo(unsigned int evt_idx, const RAT::DS::Entry& entry);
double get_Energy_Scaling(TH2F* alpha_scaling_hist, const std::vector<unsigned int>& Nbins, const std::vector<double>& Limits, double R, double Z);



int main(int argc, char** argv) {
    std::string output_file = argv[1];
    bool verbose = std::stoi(argv[2]);
    // Addresses of sim files to be analysed
    std::vector<std::string> input_files;
    for (unsigned int i = 3; i < argc; ++i) {
        input_files.push_back(argv[i]);
    }

    // Get e+ delays
    // std::vector<std::vector<double> > delays_vecs = findPositronDelays(input_files, verbose);
    // std::vector<double> delays = delays_vecs.at(0);
    // std::vector<double> delays_sansZero = delays_vecs.at(1);
    std::vector<double> delays;

    // Create time residual histograms (copied from rat/example/root/PlotHitTimeResiduals.cc)
    if (verbose) {std::cout << "Getting hists..." << std::endl;}
    std::vector<TH1D*> MC_summed_hists = HitTimeResiduals("true", input_files, delays, verbose);
    std::vector<TH1D*> Fitted_summed_hists = HitTimeResiduals("recon", input_files, delays, verbose);

    // Get o-Ps pdf and print to file
    if (verbose) {std::cout << "Printing PDFs to file..." << std::endl;}
    printPDF(output_file, MC_summed_hists, Fitted_summed_hists, verbose);

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
        for (unsigned int i = 0; i < MC_summed_hists.size(); ++i) {
            MC_summed_hists.at(i)->Write();
            Fitted_summed_hists.at(i)->Write();
        }

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
                    if ((size_t)cursor.ChildCount() == 0) {
                        delays.push_back(0.0);
                        if (verbose) {std::cout << "No child particles. Skipping." << std::endl;}
                        continue;
                    }
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
void printPDF(const std::string& output_filename, std::vector<TH1D*> MC_hists, std::vector<TH1D*> Fitted_hists, bool verbose) {

    // Open file to print results to
    if (verbose) {std::cout << "Open file" << std::endl;}
    std::ofstream datafile;
    datafile.open(output_filename, std::ofstream::out | std::ofstream::trunc);

    if (verbose) {std::cout << "Get bin width and number of bins" << std::endl;}
    double bin_width = MC_hists.at(0)->GetBinCenter(2) - MC_hists.at(0)->GetBinCenter(1);
    unsigned int N_bins = MC_hists.at(0)->GetNbinsX();

    // Checks
    if (verbose) {std::cout << "Perform checks" << std::endl;}
    if (bin_width != Fitted_hists.at(0)->GetBinCenter(2) - Fitted_hists.at(0)->GetBinCenter(1)) {
        std::cout << "Bin widths different" << std::endl;
        exit(1);
    }
    if (N_bins != Fitted_hists.at(0)->GetNbinsX()) {
        std::cout << "Number of bins different" << std::endl;
        exit(1);
    }

    std::vector<double> tot_MC_hits = {0.0, 0.0, 0.0};
    std::vector<double> tot_Fitted_hits = {0.0, 0.0, 0.0};
    std::vector<double> times;
    std::vector<std::vector<double>> MC_probs;
    std::vector<std::vector<double>> Fitted_probs;
    // Get histogram bin values
    if (verbose) {std::cout << "Get histogram bin values, and print times to file" << std::endl;}
    datafile << "[";
    for (unsigned int i = 1; i < N_bins+1; ++i) { //loop over histogram bins
        times.push_back(MC_hists.at(0)->GetBinCenter(i));
        std::vector<double> MC_prob;
        std::vector<double> Fitted_prob;
        for (unsigned int j = 0; j < MC_hists.size(); ++j) {
            MC_prob.push_back(MC_hists.at(j)->GetBinContent(i));
            tot_MC_hits.at(j) += MC_prob.at(j);
            Fitted_prob.push_back(Fitted_hists.at(j)->GetBinContent(i));
            tot_Fitted_hits.at(j) += Fitted_prob.at(j);
        }
        MC_probs.push_back(MC_prob);
        Fitted_probs.push_back(Fitted_prob);

        // print times to file while at it
        datafile << times.at(i-1);
        if (i < N_bins) {
            datafile << ", ";
        }
    }
    datafile << "]" << std::endl;

    for (unsigned int i = 0; i < MC_hists.size(); ++i) {
        std::cout << "Total number of MC PMT hits in PDF " << i << " = " << tot_MC_hits.at(i) << std::endl;
        std::cout << "Total number of Fitted PMT hits in PDF " << i << " = " << tot_Fitted_hits.at(i) << std::endl;
    }

    // Normalise MC to pdf
    if (verbose) {std::cout << "Normalising and writing true PDF" << std::endl;}
    for (unsigned int j = 0; j < MC_hists.size(); ++j) {
        datafile << "[";
        for (unsigned int i = 0; i < MC_probs.size(); ++i) {
            MC_probs.at(i).at(j) /= (tot_MC_hits.at(j) * bin_width);
            datafile << MC_probs.at(i).at(j);
            if (i < MC_probs.size() - 1) {
                datafile << ", ";
            }
        }
        datafile << "]" << std::endl;
    }

    // Normalise Fitted to pdf
    if (verbose) {std::cout << "Normalising and writing recon PDF" << std::endl;}
    for (unsigned int j = 0; j < Fitted_probs.size(); ++j) {
        datafile << "[";
        for (unsigned int i = 0; i < Fitted_probs.size(); ++i) {
            Fitted_probs.at(i).at(j) /= (tot_Fitted_hits.at(j) * bin_width);
            datafile << Fitted_probs.at(i).at(j);
            if (i < Fitted_probs.size() - 1) {
                datafile << ", ";
            }
        }
        datafile << "]" << std::endl;
    }

    if (verbose) {std::cout << "Closing file" << std::endl;}
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
std::vector<TH1D*> HitTimeResiduals(std::string type, const std::vector<std::string>& fileNames, std::vector<double> delays, bool verbose, std::string fitName) {
    if (verbose) {std::cout << "Running HitTimeResidualsFitPosition()" << std::endl;}

    if (type != "true" && type != "recon") {
        std::cout << "info_type should be 'true' or 'recon', not '" << type << "'." << std::endl;
        exit(1);
    }

    // Read in Iwan's energy correction histogram and compute constant needed to use it loop
    // Read in root file
    std::string E_corr_filename = "/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/escaling_data_mc_mc_out_NTUPLE_SEG_6000_xyznhits_erecon_tag_nhits_PartialScintBipo214_ScintRun_257635_264716.root";
    TFile* fin = new TFile(E_corr_filename.c_str());
    if (!fin->IsOpen()) {
        std::cout << "Cannot open input file " << E_corr_filename << std::endl;
        exit(1);
    }
    // Get histogram (Po used to calibrate alpha events. In this case for Energy scaling)
    TH2F* alpha_scaling_hist = (TH2F*)fin->Get("h2_Po_scale_z_R");
    // Get info
    std::vector<unsigned int> Nbins = {(unsigned int)(alpha_scaling_hist->GetNbinsX()), (unsigned int)(alpha_scaling_hist->GetNbinsY())};
    std::vector<double> Limits = {alpha_scaling_hist->GetXaxis()->GetBinCenter(1), alpha_scaling_hist->GetXaxis()->GetBinCenter(Nbins.at(0)),
                                  alpha_scaling_hist->GetYaxis()->GetBinCenter(1), alpha_scaling_hist->GetYaxis()->GetBinCenter(Nbins.at(1))};

    // Create output histograms we will get PDFs from
    std::vector<TH1D*> histograms;
    for (unsigned int i = 0; i < 3; ++i) {
        std::string hist_name = "hHitTimeResiduals_" + type + "_" + to_string(i);
        std::string hist_title = "Hit time residuals using the " + type + " position" + ", " + to_string(i);
        histograms.push_back(new TH1D(hist_name.c_str(), hist_title.c_str(), 1300, -300.5, 999.5));
    }
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
            if (verbose) {std::cout << "iEntry = " << iEntry << std::endl;}
            const RAT::DS::Entry& rDS = dsReader.GetEntry(iEntry);
            for (size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++) {
                if (verbose) {std::cout << "iEV = " << iEV << std::endl;}
                const RAT::DS::EV& rEV = rDS.GetEV(iEV);

                // Get indices of events that pass cuts
                std::vector<std::vector<unsigned int>> passed_evt_indices_list;
                if (verbose) {std::cout << "passed_evt_indices_list.size() = " << passed_evt_indices_list.size() << std::endl;}
                passed_evt_indices_list.push_back(apply_cuts("alphaN_1", iEV, rDS, rEV, type, alpha_scaling_hist, Nbins, Limits, fitName));
                passed_evt_indices_list.push_back(apply_cuts("alphaN_2", iEV, rDS, rEV, type, alpha_scaling_hist, Nbins, Limits, fitName));
                passed_evt_indices_list.push_back(apply_cuts("alphaN_3", iEV, rDS, rEV, type, alpha_scaling_hist, Nbins, Limits, fitName));
                // Add passed events to histograms
                for (unsigned int j = 0; j < passed_evt_indices_list.size(); ++j) {
                    std::vector<unsigned int> passed_evt_indices = passed_evt_indices_list.at(j);
                    if (verbose) {std::cout << "Going through events that pass alphaN_" << j+1 << " cuts..." << std::endl;}
                    for (unsigned int k = 0; k < passed_evt_indices.size(); ++k) {
                        if (verbose) {std::cout << "Passed event number: " << passed_evt_indices.at(k) << std::endl;}
                        const RAT::DS::EV& passed_evt = rDS.GetEV(passed_evt_indices.at(k));

                        // Get event info
                        TVector3 evt_pos;
                        double evt_time;
                        if (type == "recon") {
                            std::vector<double> evt_recon_info = getReconInfo(iEV, rDS, alpha_scaling_hist, Nbins, Limits, fitName);
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
                            histograms.at(j)->Fill(fTRCalc.CalcTimeResidual(pmtCal, evt_pos, evt_time));
                        }
                        if (verbose) {std::cout << "...added." << std::endl;}
                        ++num_evts;
                    }
                }
            }
        }
        dsReader.Delete();
    }
    std::cout << "Number of events recorded = " << num_evts << std::endl;

    for (unsigned int i = 0; i < histograms.size(); ++i) {
        histograms.at(i)->GetYaxis()->SetTitle( "Count per 1 ns bin" );
        histograms.at(i)->GetXaxis()->SetTitle( "Hit time residuals [ns]" );
        histograms.at(i)->Draw();
    }

    return histograms;
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
 * @param info_type = 'true' or 'recon'
 * @param alpha_scaling_hist Energy correction histogram (find bin matching R and Z, then multiply energy by the value in this bin).
 * @param Nbins {Number of bins in x-direction (R), Number of bins in y-direction (Z)}, for alpha_scaling_hist.
 * @param Limits {R_min, R_max, Z_min, Z_max}, from alpha_scaling_hist.
 * @param fitName 
 * @return std::vector<unsigned int> 
 */
std::vector<unsigned int> apply_cuts(std::string cuts_name, unsigned int evt_idx, const RAT::DS::Entry& entry, const RAT::DS::EV& evt, std::string info_type, TH2F* alpha_scaling_hist, const std::vector<unsigned int>& Nbins, const std::vector<double>& Limits, std::string fitName) {
    if (info_type != "true" && info_type != "recon") {
        std::cout << "info_type should be 'true' or 'recon', not '" << info_type << "'." << std::endl;
        exit(1);
    }
    if (cuts_name == "alphaN_1" || cuts_name == "alphaN_2" || cuts_name == "alphaN_3") {
        return alphaN_cuts(cuts_name, evt_idx, entry, evt, info_type, alpha_scaling_hist, Nbins, Limits, fitName);
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
 * @param info_type = 'true' or 'recon'
 * @param alpha_scaling_hist Energy correction histogram (find bin matching R and Z, then multiply energy by the value in this bin).
 * @param Nbins {Number of bins in x-direction (R), Number of bins in y-direction (Z)}, for alpha_scaling_hist.
 * @param Limits {R_min, R_max, Z_min, Z_max}, from alpha_scaling_hist.
 * @param fitName 
 * @return std::vector<unsigned int> 
 */
std::vector<unsigned int> alphaN_cuts(std::string cuts_name, unsigned int evt_idx, const RAT::DS::Entry& entry, const RAT::DS::EV& evt, std::string info_type, TH2F* alpha_scaling_hist, const std::vector<unsigned int>& Nbins, const std::vector<double>& Limits, std::string fitName) {
    std::vector<unsigned int> passed_evt_indices;

    // Only pass the prompt event if this is the delayed event, so we know it has a delayed event
    if (evt_idx != 1) {return passed_evt_indices;}

    // Get relevent info for cuts
        // For prompt event
        const RAT::DS::EV& prompt_evt = entry.GetEV(0);
        // unsigned int prompt_nhits = prompt_evt.GetNhitsCleaned();
        double prompt_time = prompt_evt.GetClockCount50() * 20.0;  // nanoseconds
        std::vector<double> prompt_info;
        if (info_type == "recon") {
            prompt_info =  getReconInfo(0, entry, alpha_scaling_hist, Nbins, Limits, fitName);
        } else {
            prompt_info =  getMCInfo(0, entry);  // FIXME: Add energy to event energy if e+? 1.02MeV
        }
        if (prompt_info.size() != 5) {return passed_evt_indices;}
        double prompt_energy = prompt_info.at(0);
        TVector3 prompt_pos = TVector3(prompt_info.at(1), prompt_info.at(2), prompt_info.at(3));
        double prompt_R = prompt_pos.Mag();

        // For delayed event
        // unsigned int delayed_nhits = delayed_evt.GetNhitsCleaned();
        double delayed_time = evt.GetClockCount50() * 20.0;  // nanoseconds
        std::vector<double> delayed_info;
        if (info_type == "recon") {
            delayed_info =  getReconInfo(evt_idx, entry, alpha_scaling_hist, Nbins, Limits, fitName);
        } else {
            delayed_info =  getMCInfo(evt_idx, entry);  // FIXME: Add energy to event energy if e+? 1.02MeV
        }
        if (delayed_info.size() != 5) {return passed_evt_indices;}
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
    if (deltaT < 400 or deltaT > 0.8E6) {return passed_evt_indices;}  // in ns

    // Extra cuts
        // index: "te_diol_dda_0p5_labppo_scintillator_bisMSB_Jan2018",
        // nhit_threshold: [1500, 2000],
        // energy_threshold: [4.0, 5.4],
        //
        // index: "labppo_1p1_berkeley_scintillator",
        // nhit_threshold: [1300, 1700],
        // energy_threshold: [3.5, 5.4],
        //
        // index: "labppo_0p5_scintillator",
        // nhit_threshold: [1200, 1550],
        // energy_threshold: [3.5, 5.4],
    if (cuts_name == "alphaN_1" && prompt_energy >= 3.5) {return passed_evt_indices;}
    if (cuts_name == "alphaN_2" && (prompt_energy < 3.5 || prompt_energy > 5.4)) {return passed_evt_indices;}
    if (cuts_name == "alphaN_3" && prompt_energy <= 5.4) {return passed_evt_indices;}

    // Passed all cuts, return prompt event index
    passed_evt_indices = {0};
    return passed_evt_indices;
}

// std::vector<unsigned int> oPs_cuts(std::string cuts_name, unsigned int evt_idx, const RAT::DS::Entry& entry, const RAT::DS::EV& evt, std::string info_type, std::string fitName) {
//     if (verbose) {std::cout << "evt_idx = " << evt_idx << std::endl;}
//     if (evt_idx >= delays.size()) {
//         if (verbose) {std::cout << "evt_idx out of range of delays vector" << std::endl;}
//         ++evt_idx;
//         continue;
//     } else if (is_oPs && delays.at(evt_idx) == 0.0) {  // Filter out non o-Ps events
//         if (verbose) {std::cout << "Delay = 0, ignoring event." << std::endl;}
//         ++evt_idx;
//         continue;
//     } else if (eventPosition.Mag() > vol_cut) {
//         if (verbose) {std::cout << "Outside volume cut, ignoring event." << std::endl;}
//         ++evt_idx;
//         continue;
//     } else if (nhits < 200 || nhits > 6000) {
//         if (verbose) {std::cout << "Outside nhit cut, ignoring event." << std::endl;}
//         ++evt_idx;
//         continue;
// }


/**
 * @brief Get the reconstructed energy and position of an event as {E, x, y, z, t}.
 * If error in reconstruction, returns empty vector.
 * 
 * @param evt_idx 
 * @param entry 
 * @param alpha_scaling_hist Energy correction histogram (find bin matching R and Z, then multiply energy by the value in this bin).
 * @param Nbins {Number of bins in x-direction (R), Number of bins in y-direction (Z)}, for alpha_scaling_hist.
 * @param Limits {R_min, R_max, Z_min, Z_max}, from alpha_scaling_hist.
 * @param fitName 
 * @return std::vector<double> 
 */
std::vector<double> getReconInfo(unsigned int evt_idx, const RAT::DS::Entry& entry, TH2F* alpha_scaling_hist, const std::vector<unsigned int>& Nbins, const std::vector<double>& Limits, std::string fitName) {
    const RAT::DS::EV& evt = entry.GetEV(evt_idx);
    std::vector<double> output_vector;
    // grab the fit information
    if(fitName == "")
        fitName = evt.GetDefaultFitName();
    try {
        // Get recon info
        const RAT::DS::FitVertex& rVertex = evt.GetFitResult(fitName).GetVertex(0);
        if (!(rVertex.ValidPosition() && rVertex.ValidTime() && rVertex.ValidEnergy())) {return output_vector;} // fit invalid
        double energy = rVertex.GetEnergy();
        TVector3 pos = rVertex.GetPosition();
        double time = rVertex.GetTime();

        // Apply energy scaling (Iwan's energy scaling for partial-fill)
        energy *= get_Energy_Scaling(alpha_scaling_hist, Nbins, Limits, pos.Mag(), pos.Z());

        // Package info 
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



/* ~~~~~~~~~~~~~~~~~~~~~~ ENERGY SCALING ~~~~~~~~~~~~~~~~~~~~~ */


/**
 * @brief Get the Energy Scaling factor from Iwan's energy correction, for partial-fill
 * 
 * @param alpha_scaling_hist Energy correction histogram (find bin matching R and Z, then multiply energy by the value in this bin).
 * @param Nbins {Number of bins in x-direction (R), Number of bins in y-direction (Z)}, for alpha_scaling_hist.
 * @param Limits {R_min, R_max, Z_min, Z_max}, from alpha_scaling_hist.
 * @param R Distance of event from the AV centre.
 * @param Z Projection of distance of event from the AV centre onto the Z-axis.
 * @return double 
 */
double get_Energy_Scaling(TH2F* alpha_scaling_hist, const std::vector<unsigned int>& Nbins, const std::vector<double>& Limits, double R, double Z) {
    // Unpack values
    unsigned int NbinsR = Nbins.at(0);
    unsigned int NbinsZ = Nbins.at(1);
    double R_min = Limits.at(0);
    double R_max = Limits.at(1);
    double Z_min = Limits.at(2);
    double Z_max = Limits.at(3);

    // Get x and y bin index from R and Z respectively
    double ratio_R = (R - R_min) / (R_max - R_min);
    unsigned int n_bins_R = (unsigned int)(ratio_R + 0.5 - (ratio_R < 0.0)); // Round ratio to nearest integer, since type cast always truncates
    unsigned int R_idx = 1 + NbinsR * n_bins_R;

    double ratio_Z = (Z - Z_min) / (Z_max - Z_min);
    unsigned int n_bins_Z = (unsigned int)(ratio_Z + 0.5 - (ratio_Z < 0.0)); // Round ratio to nearest integer, since type cast always truncates
    unsigned int Z_idx = 1 + NbinsZ * n_bins_Z;

    return alpha_scaling_hist->GetBinContent(R_idx, Z_idx);
}


// bin = 0;       underflow bin
// bin = 1;       first bin with low-edge xlow INCLUDED
// bin = nbins;   last bin with upper-edge xup EXCLUDED
// bin = nbins+1; overflow bin