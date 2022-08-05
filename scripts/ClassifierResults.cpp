////////////////////////////////////////////////////////////////////
/// \file BASED ON: PlotHitTimeResiduals.cc
///
/// \brief Functions to get o-Ps classification and compare to decay times.
///
/// \author James Page <j.page@sussex.ac.uk>
///
/// REVISION HISTORY:\n
///
/// \details Functions to get o-Ps classification and compare to decay times.
///
/// To compile: g++ -g -std=c++1y ClassifierResults.cpp -o ClassifierResults.exe `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux
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
#include <RAT/DS/FitClassifierCollection.hh>

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

std::vector<std::vector<double> > findPositronDelays_andClassification(const std::vector<std::string>& filenames, const std::string& out_file_address, double vol_cut, bool is_oPs, bool make_hists, bool verbose);
void printResults(const std::string& output_filename, const std::vector<double>& delays, const std::vector<double>& oPs_classier_results, const std::vector<double>& alphaNreactor_classier_results, const std::vector<double>& nhits_vec, const std::vector<double>& energies);

int main(int argc, char** argv) {
    std::string output_file = argv[1];
    bool is_oPs = std::stoi(argv[2]);
    bool make_hists = std::stoi(argv[3]);
    bool verbose = std::stoi(argv[4]);
    // Addresses of sim files to be analysed
    std::vector<std::string> input_files;
    for (unsigned int i = 5; i < argc; ++i) {
        input_files.push_back(argv[i]);
    }

    // Get e+ delays
    std::vector<std::vector<double> > results = findPositronDelays_andClassification(input_files, output_file, 5700, is_oPs, make_hists, verbose);
    const std::vector<double> delays = results.at(0);
    const std::vector<double> oPs_classier_results = results.at(1);
    const std::vector<double> alphaNreactor_classier_results = results.at(2);
    const std::vector<double> nhits_vec = results.at(3);
    const std::vector<double> energies = results.at(4);

    // Get o-Ps pdf and print to file
    if (verbose) {std::cout << "Printing results to file..." << std::endl;}
    printResults(output_file, delays, oPs_classier_results, alphaNreactor_classier_results, nhits_vec, energies);

    return 0;
}


/**
 * @brief Returns a list of the delays imparted to positron decays (emulating oPs), classifications, nhits and recon energies.
 * 
 * @param filenames List of simulation filenames to analyse .
 * @param out_file_address Text filename to write results to. 
 * @param vol_cut  Radius of volume cut (mm).
 * @param is_oPs 
 * @param make_hists 
 * @param verbose 
 * @return std::vector<double> = {delays, oPs_classier_results, alphaNreactor_classier_results, nhits_vec, energies}
 */
std::vector<std::vector<double> > findPositronDelays_andClassification(const std::vector<std::string>& filenames, const std::string& out_file_address, double vol_cut, bool is_oPs, bool make_hists, bool verbose) {
    if (verbose) {std::cout << "Finding e+ delays..." << std::endl;}

    RAT::DB::Get()->SetAirplaneModeStatus(true);
    const RAT::DU::PMTCalStatus& PMTCalStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();  // Needed for GetHitStatus() later

    if (make_hists) {
        // Create output file names
        if (verbose) {std::cout << "Creating output root file" << std::endl;}
        std::size_t botDirPos = out_file_address.find_last_of("/");
        std::string filename = out_file_address.substr(botDirPos+1, out_file_address.length() - 5);
        std::string path = out_file_address.substr(0, botDirPos+1);
        std::string saveroot = path + "Hists_" + filename + ".root";
        // output root file
        TFile *rootfile = new TFile(saveroot.c_str(), "RECREATE");
        rootfile->cd();
    }

    // sumed over all events
    TH1D* summed_hist = new TH1D("SummedTimeResidualsMC", "Summed hit time residuals using the MC position, and everage o-Ps delay = ? ns", 1000, -10.0, 500.0);
    double mean_delay = 0.0;
    unsigned int num_evts = 0;

    double delay = 0.0;
    std::vector<double> delays;
    double energy;
    std::vector<double> energies;
    double oPs_classier_result = 0.0;
    std::vector<double> oPs_classier_results;
    double alphaNreactor_classier_result = 0.0;
    std::vector<double> alphaNreactor_classier_results;
    double nhits = 0.0;
    std::vector<double> nhits_vec;

    for (unsigned int i = 0; i < filenames.size(); ++i) {
        RAT::DU::DSReader dsReader(filenames.at(i));
        // Loop through events. Each one should have one primary track (first child) in MC
        if (verbose) {std::cout << "Looping through events..." << std::endl;}
        for (size_t iEv = 0; iEv < dsReader.GetEntryCount(); iEv++) {
            const RAT::DS::Entry& rDS = dsReader.GetEntry(iEv);
            RAT::TrackNav nav(&rDS);
            RAT::TrackCursor cursor = nav.Cursor(false);

            // Check there is an event in this entry (won't get associated t_res plot if not)
                if (rDS.GetEVCount() > 0) {
                const RAT::DS::EV& rEV = rDS.GetEV(0);

                double posRad = 0.0;
                double energy = -1.0;
                try {
                    const RAT::DS::FitVertex& rVertex = rEV.GetFitResult(rEV.GetDefaultFitName()).GetVertex(0);
                    if (!(rVertex.ValidPosition() && rVertex.ValidTime() && rVertex.ValidEnergy())) {continue;} // fit invalid
                    posRad = rVertex.GetPosition().Mag();
                    energy = rVertex.GetEnergy();
                }
                catch (const RAT::DS::FitCollection::NoResultError&) {
                    // no fit result by the name of fitName
                    continue;
                }
                catch (const RAT::DS::FitResult::NoVertexError&) {
                    // no fit vertex
                    continue;
                }
                catch (const RAT::DS::FitVertex::NoValueError&) {
                    // position or time missing
                    continue;
                }
                // DataNotFound --> implies no fit results are present, don't catch.

                // Apply volume cut
                if (posRad > vol_cut) {
                    if (verbose) {std::cout << "Outside volume cut, ignoring event." << std::endl;}
                    continue;
                }

                // Get classifier results
                try{
                    if (!rEV.ClassifierResultExists("PositroniumClassifier")) {
                        std::cout << "No PositroniumClassifier results for entry " << iEv << std::endl;
                        continue;
                    }
                    if (!rEV.GetClassifierResult("PositroniumClassifier").GetValid()) {
                        std::cout << "No valid PositroniumClassifier result for entry " << iEv << std::endl;
                        continue;
                    }
                    RAT::DS::ClassifierResult oPs_result = rEV.GetClassifierResult("PositroniumClassifier");
                    oPs_classier_result = oPs_result.GetClassification("PositroniumClassifier");
                    if (verbose) {std::cout << "POSITRONIUM classifier result (for next particle) = " << oPs_classier_result << std::endl;}
                } catch (RAT::DS::ClassifierResult::NoClassificationError&) {
                    std::cout << "Error in PositroniumClassifier" << std::endl;
                    continue;
                }

                try {
                    //RAT::DS::FitClassifierCollection<RAT::DS::FitResult> fitResult = rEV.GetFitResults(0);
                    //RAT::DS::ClassifierResult alphaNreactor_result = fitResult.GetResult("AlphaNReactorIBDClassifier");

                    if (!rEV.ClassifierResultExists("AlphaNReactorIBDClassifier")) {
                        std::cout << "No AlphaNReactorIBDClassifier results for entry " << iEv << std::endl;
                        continue;
                    }
                    if (!rEV.GetClassifierResult("AlphaNReactorIBDClassifier").GetValid()) {
                        std::cout << "No valid AlphaNReactorIBDClassifier result for entry " << iEv << std::endl;
                        continue;
                    }
                    RAT::DS::ClassifierResult alphaNreactor_result = rEV.GetClassifierResult("AlphaNReactorIBDClassifier");
                    alphaNreactor_classier_result = alphaNreactor_result.GetClassification("AlphaNReactorIBDClassifier");
                    if (verbose) {std::cout << "alphaNreactor classifier result (for next particle) = " << alphaNreactor_classier_result << std::endl;}
                } catch (RAT::DS::ClassifierResult::NoClassificationError&) {
                    std::cout << "Error in AlphaNReactorIBDClassifier" << std::endl;
                    continue;
                }

                nhits = rEV.GetNhitsCleaned();

                // Should only go through this loop once in MC.
                if (verbose) {std::cout << "Getting track history..." << std::endl;}
                //for (size_t iCh = 0; iCh<(size_t)cursor.ChildCount(); iCh++) {
                    //cursor.GoChild(iCh);
                    cursor.GoChild(0);

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
                    delay = end_time - start_time;
                    if (verbose) {std::cout << "Child particle: " << child_node->GetParticleName() << std::endl;}
                    if (verbose) {std::cout << "e+ delay: " << delay << std::endl;}

                    // Go back to e+ track
                    cursor.GoParent();

                    // Go back to parent node to redo loop
                    cursor.GoTrackStart();
                    cursor.GoParent();
                //} //Primary Particle Tracks

                // Make sure that if it's only o-Ps sims to leave out the few bugged events with no delay
                if (!is_oPs || delay != 0.0) {
                    if (verbose) {std::cout << "Writing info..." << std::endl;}
                    delays.push_back(delay);
                    oPs_classier_results.push_back(oPs_classier_result);
                    alphaNreactor_classier_results.push_back(alphaNreactor_classier_result);
                    nhits_vec.push_back(nhits);
                    energies.push_back(energy);

                    /* ~~~~~~ Make time residual hist to check results make sense ~~~~~ */
                    if (make_hists) {
                        if (verbose) {std::cout << "Making histograms... (hist #" << num_evts << ")" << std::endl;}
                        // update mean delay
                        mean_delay += (delay - mean_delay) / (num_evts + 1);

                        // for each individual event
                        std::string hist_name = "hHitTimeResidualsMC_" + std::to_string(num_evts);
                        std::string title = "Hit time residuals using the MC position, o-Ps delay = " + std::to_string(delay)
                                            + " ns, E = " + std::to_string(energy) + " MeV, Classifier result = " + std::to_string(oPs_classier_result);
                        TH1D* evt_hist = new TH1D(hist_name.c_str(), title.c_str(), 1000, -10.0, 500.0);

                        const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
                        const TVector3 eventPosition = rDS.GetMC().GetMCParticle(0).GetPosition(); // At least 1 is somewhat guaranteed
                        double event_time = 390 - rDS.GetMCEV(0).GetGTTime();  // event time is 390ns - GT time.
                        for(size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++) {
                            const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT(iPMT);
                            // Check if PMT passes data cleaning
                            if (PMTCalStatus.GetHitStatus(pmtCal) != 0) {continue;}
                            // Get time residuals
                            RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();
                            evt_hist->Fill(fTRCalc.CalcTimeResidual(pmtCal, eventPosition, event_time));
                            summed_hist->Fill(fTRCalc.CalcTimeResidual(pmtCal, eventPosition, event_time));
                        }
                        // Write event residual hit time to root file
                        evt_hist->Write();
                    }
                    if (verbose) {std::cout << "End of loop." << std::endl;}
                    ++num_evts;
                }
            }
        } //event
        //dsReader.Delete();
        dsReader.~DSReader();
    }

    if (make_hists) {
        // Write summed histogram to root file
        std::string summed_title = "Hit time residuals using the MC position, summed over " + std::to_string(num_evts)
                                    + " events, and mean o-Ps delay = " + std::to_string(mean_delay) + " ns";
        summed_hist->SetTitle(summed_title.c_str());
        summed_hist->Write();
    }

    if (verbose) {std::cout << "Num delays: " << delays.size() << std::endl;}
    if (verbose) {std::cout << "Num oPs_classifier results: " << oPs_classier_results.size() << std::endl;}
    std::vector<std::vector<double> > results = {delays, oPs_classier_results, alphaNreactor_classier_results, nhits_vec, energies};
    return results;
}


/**
 * @brief Print pdf to text file based on inputted histogram.
 * Each line is:
 * oPs_classifier_result alphaNreactor_classifier_result delay nhits energy
 * 
 * @param output_filename 
 * @param hist 
 */
void printResults(const std::string& output_filename, const std::vector<double>& delays, const std::vector<double>& oPs_classier_results, const std::vector<double>& alphaNreactor_classier_results, const std::vector<double>& nhits_vec, const std::vector<double>& energies) {

    // Check results make sense
    if (!(delays.size() == oPs_classier_results.size() && delays.size() == alphaNreactor_classier_results.size() && delays.size() == nhits_vec.size() && delays.size() == energies.size())) {
        std::cout << "ERROR: number of elements not equal:" << std::endl;
        std::cout << "delays: " << delays.size() << ", oPs_classier_results: " << oPs_classier_results.size() << ", alphaNreactor_classier_results: " << alphaNreactor_classier_results.size() << ", nhits_vec: " << nhits_vec.size() << ", energies: " << energies.size() << std::endl;
        exit(1);
    }

    // Open file and print results
    std::ofstream datafile;
    datafile.open(output_filename, std::ofstream::out | std::ofstream::trunc);
    for (unsigned int i = 0; i < delays.size(); ++i) {
        datafile << oPs_classier_results.at(i) << " " << alphaNreactor_classier_results.at(i) << " " << delays.at(i) << " " << nhits_vec.at(i) << " " << energies.at(i) << std::endl;
    }

    datafile.close();
}