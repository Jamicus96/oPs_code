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

std::vector<double> findPositronDelays(const std::string& filename, bool verbose);
void printPDF(const std::string& output_filename, TH1D* hist);
TH1D* PlotHitTimeResidualsMCPosition(const std::string& fileName, std::vector<double> delays, bool is_oPs, bool verbose);

int main(int argc, char** argv) {
    std::string file = argv[1];
    bool is_oPs = argv[2];
    bool verbose = argv[3];

    // Get e+ delays
    std::vector<double> delays = findPositronDelays(file, verbose);

    // Create time residual histograms (copied from rat/example/root/PlotHitTimeResiduals.cc)
    if (verbose) {std::cout << "Getting hists..." << std::endl;}
    TH1D* MC_summed_hist = PlotHitTimeResidualsMCPosition(file, delays, is_oPs, verbose);

    // Create output file names
    if (verbose) {std::cout << "Creating output file" << std::endl;}
    std::size_t botDirPos = file.find_last_of("/");
    std::string filename = file.substr(botDirPos+1, file.length() - 5);
    std::string saveroot = "Hists_" + filename + ".root";
    std::string pdf_filename = "pdf_" + filename + ".txt";

    // Get o-Ps pdf and print to file
    printPDF(pdf_filename, MC_summed_hist);

    // Save root file
    TFile *rootfile = new TFile(saveroot.c_str(),"RECREATE");

    // Now write everything to the file and close
    if (verbose) {std::cout << "Writing everything to file and closing" << std::endl;}
    rootfile->cd();
    MC_summed_hist->Write();

    // raw_hist->Write();
    rootfile->Write();
    rootfile->Close();

    return 0;
}

/**
 * @brief Returns a list of the delays imparted to positron decays (emulating oPs)
 * 
 * @param filename 
 * @return std::vector<double> 
 */
std::vector<double> findPositronDelays(const std::string& filename, bool verbose) {
  if (verbose) {std::cout << "Finding e+ delays..." << std::endl;}

  RAT::DU::DSReader dsReader(filename);
  std::vector<double> delays;

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

  return delays;
}

/**
 * @brief Print pdf to text file based on inputted histogram.
 * 
 * @param output_filename 
 * @param hist 
 */
void printPDF(const std::string& output_filename, TH1D* hist) {

    // Open file to print results to
    std::ofstream datafile;
    datafile.open(output_filename, std::ofstream::out | std::ofstream::trunc);
    datafile << "[";

    double bin_width = hist->GetBinCenter(2) - hist->GetBinCenter(1);
    double tot_hits = 0.0;
    std::vector<double> times;
    std::vector<double> probabilities;
    // Get histogram bin values
    for (unsigned int i = 1; i < hist->GetNbinsX()+1; ++i) { //loop over histogram bins
        times.push_back(hist->GetBinCenter(i));
        probabilities.push_back(hist->GetBinContent(i));
        tot_hits += hist->GetBinContent(i);

        // print times to file while as it
        datafile << hist->GetBinCenter(i);
        if (i < hist->GetNbinsX()) {
            datafile << ", ";
        }
    }
    datafile << "]" << std::endl;

    std::cout << "Total number of PMT hits = " << tot_hits << std::endl;

    // Normalise to pdf
    datafile << "[";
    for (unsigned int i = 0; i < probabilities.size(); ++i) {
        probabilities.at(i) /= (tot_hits * bin_width);
        datafile << probabilities.at(i);
        if (i < probabilities.size() - 1) {
            datafile << ", ";
        }
    }
    datafile << "]" << std::endl;
    datafile.close();
}

/// Plot the hit time residuals for the MC position
///
/// @param[in] fileName of the RAT::DS root file to analyse
/// @return the histogram plot
TH1D* PlotHitTimeResidualsMCPosition(const std::string& fileName, std::vector<double> delays, bool is_oPs, bool verbose) {
    if (verbose) {std::cout << "Running pdfMCPosition()" << std::endl;}

    TH1D* histTimeResiduals = new TH1D( "pdfTimeResidualsMC", "PDF for Hit time residuals using the MC position", 1300, -300.5, 999.5 );
    // If this is being done on data that does not require remote database connection
    // eg.: a simple simulation with default run number (0)
    // We can disable the remote connections:
    //
    // NOTE: Don't do this if you are using real data!!!
    RAT::DB::Get()->SetAirplaneModeStatus(true);

    RAT::DU::DSReader dsReader( fileName );

    // RAT::DU::Utility::Get()->GetLightPathCalculator() must be called *after* the RAT::DU::DSReader constructor.
    RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
    const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
    const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo(); // The PMT positions etc...

    if (verbose) {std::cout << "Looping through entries..." << std::endl;}
    unsigned int evt_idx = 0;
    unsigned int num_evts = 0;
    for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ ) {
        const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
        const TVector3 eventPosition = rDS.GetMC().GetMCParticle(0).GetPosition(); // At least 1 is somewhat guaranteed
        for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ ) {
            const RAT::DS::EV& rEV = rDS.GetEV( iEV );
            const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
            for( size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++ ) {
                if (is_oPs && delays.at(evt_idx) == 0.0) {  // Filter out non o-Ps events
                    continue;
                } else {
                    ++num_evts;
                    const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT( iPMT );

                    // lightPath.CalcByPosition( eventPosition, pmtInfo.GetPosition( pmtCal.GetID() ) );
                    // double distInInnerAV = lightPath.GetDistInInnerAV();
                    // double distInAV = lightPath.GetDistInAV();
                    // double distInWater = lightPath.GetDistInWater();

                    // const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
                    // // Time residuals estimate the photon emission time relative to the event start so subtract off the transit time
                    // // hit times are relative to the trigger time, which will depend on event time and detector position so correct for that to line up events
                    // // The 390ns corrects for the electronics delays and places the pulse in the middle of the window

                    // histTimeResiduals->Fill( pmtCal.GetTime() - transitTime - 390 + rDS.GetMCEV(iEV).GetGTTime());

                    RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();
                    histTimeResiduals->Fill(fTRCalc.CalcTimeResidual(pmtCal, eventPosition, rDS.GetMCEV(iEV).GetGTTime()));
                }
            }
        }
    }

    histTimeResiduals->GetYaxis()->SetTitle( "Count per 1 ns bin" );
    histTimeResiduals->GetXaxis()->SetTitle( "Hit time residuals [ns]" );
    histTimeResiduals->Draw();
    return histTimeResiduals;
}