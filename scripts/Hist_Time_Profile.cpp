////////////////////////////////////////////////////////////////////
/// \file COPIED FROM: PlotHitTimeResiduals.cc
///
/// \brief Functions to plot hit time residuals.
///
/// \author P G Jones <p.g.jones@qmul.ac.uk>
///
/// REVISION HISTORY:\n
///     2014-03-27 : P G Jones - First Revision.\n
///
/// \details EV Calibrated hit times are plotted minus transit times
/// based on the MC position or the fitted position.
///
/// To compile: g++ -g -std=c++1y Hist_Time_Profile.cpp -o Hist_Time_Profile.exe `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux
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

std::vector<double> findPositronDelays(const std::string& filename, bool verbose);
TH1D* PlotHitTimeResidualsMCPosition(const std::string& fileName, bool verbose);
std::vector<TH1D*> PlotHitTimeResidualsMCPosition_individual(const std::string& fileName, bool verbose);
TH1D* PlotHitTimeResidualsFitPosition(const std::string& fileName, bool verbose, std::string fitName = "");
// void PlotHitTimeResiduals(const std::string& fileName);

int main(int argc, char** argv) {
  std::string file = argv[1];
  bool verbose = argv[2];

  // Create time residual histograms (copied from rat/example/root/PlotHitTimeResiduals.cc)
  if (verbose) {std::cout << "Getting hists..." << std::endl;}
  TH1D* MC_summed_hist = PlotHitTimeResidualsMCPosition(file, verbose);
  std::vector<TH1D*> evt_hists = PlotHitTimeResidualsMCPosition_individual(file, verbose);
  // TH1D* raw_hist = PlotHitTimeResidualsFitPosition(file, verbose);

  // Get e+ delays
  std::vector<double> delays = findPositronDelays(file, verbose);

  // Check list of hists and delays are the same length
  if (delays.size() == evt_hists.size()) {
    std::string title;
    double mean_delay = 0.0;
    // Write delay to appropriate title
    for (unsigned int i = 0; i < evt_hists.size(); ++i) {
      title = "Hit time residuals using the MC position, and o-Ps delay = " + std::to_string(delays.at(i)) + " ns";
      evt_hists.at(i)->SetTitle(title.c_str());
      // While here, compute mean delay (just for extra info)
      mean_delay += (delays.at(i) - mean_delay) / (i + 1);
    }
    // Set summed hist title too
    title = "Hit time residuals using the MC position, summed over " + std::to_string(evt_hists.size())
                      + " events, and mean o-Ps delay = " + std::to_string(mean_delay) + " ns";
    MC_summed_hist->SetTitle(title.c_str());
  } else {
    std::cout << "delay and hist lists of different lengths: delays_len = " <<  delays.size()
              << ", hists_len = " << evt_hists.size() << std::endl;
  }

  // Create output file
  if (verbose) {std::cout << "Creating output file" << std::endl;}
  std::size_t botDirPos = file.find_last_of("/");
  std::string filename = file.substr(botDirPos+1, file.length());
  std::string saveroot = "Hists_" + filename;
  TFile *rootfile = new TFile(saveroot.c_str(),"RECREATE");

  // Now write everything to the file and close
  if (verbose) {std::cout << "Writing everything to file and closing" << std::endl;}
  rootfile->cd();
  MC_summed_hist->Write();
  for (int i = 0; i < evt_hists.size(); ++i) {
    evt_hists.at(i)->Write();
  }
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

/// Plot the hit time residuals for the MC position
///
/// @param[in] fileName of the RAT::DS root file to analyse
/// @return the histogram plot
TH1D* PlotHitTimeResidualsMCPosition( const std::string& fileName, bool verbose) {
  if (verbose) {std::cout << "Running PlotHitTimeResidualsMCPosition()" << std::endl;}

  TH1D* hHitTimeResiduals = new TH1D( "hHitTimeResidualsMC", "Hit time residuals using the MC position", 1000, -10.0, 500.0 );
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
  for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
    {
      const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
      const TVector3 eventPosition = rDS.GetMC().GetMCParticle(0).GetPosition(); // At least 1 is somewhat guaranteed
      for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ )
        {
          const RAT::DS::EV& rEV = rDS.GetEV( iEV );
          const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
          for( size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++ )
            {
              const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT( iPMT );

              lightPath.CalcByPosition( eventPosition, pmtInfo.GetPosition( pmtCal.GetID() ) );
              double distInInnerAV = lightPath.GetDistInInnerAV();
              double distInAV = lightPath.GetDistInAV();
              double distInWater = lightPath.GetDistInWater();

              const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
              // Time residuals estimate the photon emission time relative to the event start so subtract off the transit time
              // hit times are relative to the trigger time, which will depend on event time and detector position so correct for that to line up events
              // The 390ns corrects for the electronics delays and places the pulse in the middle of the window
              hHitTimeResiduals->Fill( pmtCal.GetTime() - transitTime - 390 + rDS.GetMCEV(iEV).GetGTTime());
            }
        }
    }
  hHitTimeResiduals->GetYaxis()->SetTitle( "Count per 1 ns bin" );
  hHitTimeResiduals->GetXaxis()->SetTitle( "Hit time residuals [ns]" );
  hHitTimeResiduals->Draw();
  return hHitTimeResiduals;
}

/// Plot the hit time residuals for the fit position
///
/// @param[in] fileName of the RAT::DS root file to analyse
/// @return the histogram plot
TH1D* PlotHitTimeResidualsFitPosition( const std::string& fileName, bool verbose, std::string fitName) {
  if (verbose) {std::cout << "Running PlotHitTimeResidualsFitPosition()" << std::endl;}

  TH1D* hHitTimeResiduals = new TH1D( "hHitTimeResidualsFit", "Hit time residuals using the Fit position", 1000, -10.0, 500.0 );
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
  for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
    {
      const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
      for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ )
        {
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
          for( size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++ )
            {
              const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT( iPMT );

              lightPath.CalcByPosition( eventPosition, pmtInfo.GetPosition( pmtCal.GetID() ) );
              double distInInnerAV = lightPath.GetDistInInnerAV();
              double distInAV = lightPath.GetDistInAV();
              double distInWater = lightPath.GetDistInWater();
              const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
              // Time residuals estimate the photon emission time relative to the event start so subtract off the transit time and eventTime
              hHitTimeResiduals->Fill( pmtCal.GetTime() - transitTime - eventTime);
            }
        }
    }
  hHitTimeResiduals->GetYaxis()->SetTitle( "Count per 1 ns bin" );
  hHitTimeResiduals->GetXaxis()->SetTitle( "Hit time residuals [ns]" );
  hHitTimeResiduals->Draw();
  return hHitTimeResiduals;
}

// /// Plot both the MC and Fitted position residuals
// ///
// /// @param[in] fileName of the RAT::DS root file to analyse
// void PlotHitTimeResiduals(const std::string& fileName) {

//   gStyle->SetFillColor( kWhite );
//   TCanvas* c1 = new TCanvas();
//   TH1D* mc = PlotHitTimeResidualsMCPosition(fileName);
//   TH1D* fit = PlotHitTimeResidualsFitPosition(fileName);
//   mc->Draw();
//   fit->SetLineColor( kGreen + 2 );
//   fit->Draw("SAME");
//   TLegend* t1 = new TLegend( 0.7, 0.7, 0.9, 0.9 );
//   t1->AddEntry( mc, "MC Position", "l" );
//   t1->AddEntry( fit, "Fit Position", "l" );
//   t1->Draw();
//   c1->Update();
// }

/// Plot the hit time residuals for the MC position, each event individually
///
/// @param[in] fileName of the RAT::DS root file to analyse
/// @return the a vector of histogram histogram plots
std::vector<TH1D*> PlotHitTimeResidualsMCPosition_individual( const std::string& fileName , bool verbose) {
  if (verbose) {std::cout << "Running PlotHitTimeResidualsMCPosition_individual()" << std::endl;}

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

  // Create vector of histograms: one for each event
  std::vector<TH1D*> evt_hists;
  unsigned int counter = 0;

  if (verbose) {std::cout << "Looping through entries..." << std::endl;}
  for(size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++) {
    /* ~~~~~~~ MCEV info: position ~~~~~~~ */

    const RAT::DS::Entry& rDS = dsReader.GetEntry(iEntry);
    const TVector3 eventPosition = rDS.GetMC().GetMCParticle(0).GetPosition(); // At least 1 is somewhat guaranteed
    std::vector<unsigned int> track_ids = rDS.GetMC().GetMCTrackIDs();

    /* ~~~~~~~ EV info: events and PMT hits ~~~~~~~ */

    for(size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++) {
      const RAT::DS::EV& rEV = rDS.GetEV(iEV);
      const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
      std::string hist_name = "hHitTimeResidualsMC_" + std::to_string(counter);
      std::string title = "Hit time residuals using the MC position";
      
      evt_hists.push_back(new TH1D(hist_name.c_str(), title.c_str(), 1000, -10.0, 500.0));
      for(size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++) {
          const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT(iPMT);

          lightPath.CalcByPosition(eventPosition, pmtInfo.GetPosition(pmtCal.GetID()));
          double distInInnerAV = lightPath.GetDistInInnerAV();
          double distInAV = lightPath.GetDistInAV();
          double distInWater = lightPath.GetDistInWater();

          const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
          // Time residuals estimate the photon emission time relative to the event start so subtract off the transit time
          // hit times are relative to the trigger time, which will depend on event time and detector position so correct for that to line up events
          // The 390ns corrects for the electronics delays and places the pulse in the middle of the window
          evt_hists.at(counter)->Fill( pmtCal.GetTime() - transitTime - 390 + rDS.GetMCEV(iEV).GetGTTime());
      }
      evt_hists.at(counter)->GetYaxis()->SetTitle( "Count per 1 ns bin" );
      evt_hists.at(counter)->GetXaxis()->SetTitle( "Hit time residuals [ns]" );
      evt_hists.at(counter)->Draw();
      ++counter;
    }
  }
  return evt_hists;
}