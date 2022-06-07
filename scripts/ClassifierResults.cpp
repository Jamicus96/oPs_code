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

std::vector<std::vector<double> > findPositronDelays_andClassification(const std::string& filename, bool verbose);
void printResults(const std::string& output_filename, TH1D* hist);

int main(int argc, char** argv) {
    std::string file = argv[1];
    bool is_oPs = std::stoi(argv[2]);
    bool verbose = std::stoi(argv[3]);

    // Get e+ delays
    std::vector<std::vector<double> > results = findPositronDelays_andClassification(file, verbose);
    std::vector<double> delays = results.at(0);
    std::vector<double> classier_results = results.at(1);

    // Create output file names
    if (verbose) {std::cout << "Creating output file" << std::endl;}
    std::size_t botDirPos = file.find_last_of("/");
    std::string filename = file.substr(botDirPos+1, file.length() - 5);
    std::string classiderRes_filename = "ClassifierRes_" + filename + ".txt";

    // Get o-Ps pdf and print to file
    if (verbose) {std::cout << "Printing results to file..." << std::endl;}
    printResults(classiderRes_filename, delays, classier_results);

    return 0;
}

/**
 * @brief Returns a list of the delays imparted to positron decays (emulating oPs)
 * 
 * @param filename 
 * @return std::vector<double> 
 */
std::vector<std::vector<double> > findPositronDelays_andClassification(const std::string& filename, bool verbose) {
    if (verbose) {std::cout << "Finding e+ delays..." << std::endl;}

    RAT::DU::DSReader dsReader(filename);
    std::vector<double> delays;
    std::vector<double> classier_results;

    // Loop through events. Each one should have one primary track (first child) in MC
    if (verbose) {std::cout << "Looping through events..." << std::endl;}
    for (size_t iEv =0; iEv<dsReader.GetEntryCount(); iEv++) {
        const RAT::DS::Entry& rDS = dsReader.GetEntry(iEv);
        RAT::TrackNav nav(&rDS);
        RAT::TrackCursor cursor = nav.Cursor(false);

        // Check there is an event in this entry (won't get associated t_res plot if not)
        if (rDS.GetEVCount() > 0) {
            const RAT::DS::EV& rEV = rDS.GetEV(iEV);
            RAT::DS::ClassifierResult cResult = rEV.GetClassifierResult("PositroniumClassifier");     // Get classifier result
            classier_results.push_back(cResult.GetClassification("PositroniumClassifier"));
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

    if (verbose) {std::cout << "Num delays: " << delays.size() << std::endl;}
    if (verbose) {std::cout << "Num classifier results: " << classier_results.size() << std::endl;}
    std::vector<std::vector<double> > results = {delays, classier_results};
    return results;
}

/**
 * @brief Print pdf to text file based on inputted histogram.
 * 
 * @param output_filename 
 * @param hist 
 */
void printResults(const std::string& output_filename, std::vector<double> delays, std::vector<double> classier_results) {

    // Open file to print results to
    std::ofstream datafile;
    datafile.open(output_filename, std::ofstream::out | std::ofstream::trunc);

    // print delays
    datafile << "[";
    for (unsigned int i = 0; i < delays.size(); ++i) {
        datafile << delays.at(i);
        if (i < hist->GetNbinsX()) {
            datafile << ", ";
        }
    }
    datafile << "]" << std::endl;

    // print classifier results
    datafile << "[";
    for (unsigned int i = 0; i < classier_results.size(); ++i) {
        datafile << classier_results.at(i);
        if (i < hist->GetNbinsX()) {
            datafile << ", ";
        }
    }
    datafile << "]" << std::endl;

    datafile.close();
}