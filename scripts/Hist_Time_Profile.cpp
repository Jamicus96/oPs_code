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

#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFile.h>
#include <TVectorD.h>

#include <string>

TH1D* PlotHitTimeResidualsMCPosition(const std::string& fileName, bool verbose);
std::vector<TH1D*> PlotHitTimeResidualsMCPosition_individual(const std::string& fileName, bool verbose);
TH1D* PlotHitTimeResidualsFitPosition(const std::string& fileName, bool verbose, std::string fitName = "");
// void PlotHitTimeResiduals(const std::string& fileName);

int main(int argc, char** argv){
  std::string file = argv[1];
  bool verbose = argv[2];

  // Create time residual histograms (copied from rat/example/root/PlotHitTimeResiduals.cc)
  if (verbose) {std::cout << "Getting hists..." << std::endl;}
  TH1D* MC_hist = PlotHitTimeResidualsMCPosition(file, verbose);
  std::vector<TH1D*> evt_hists = PlotHitTimeResidualsMCPosition_individual(file, verbose);
  // TH1D* raw_hist = PlotHitTimeResidualsFitPosition(file, verbose);

  // Create output file
  if (verbose) {std::cout << "Creating output file" << std::endl;}
  std::size_t botDirPos = file.find_last_of("/");
  std::string filename = file.substr(botDirPos+1, file.length());
  std::string saveroot = "Hists_" + filename;
  TFile *rootfile = new TFile(saveroot.c_str(),"RECREATE");

  // Now write everything to the file and close
  if (verbose) {std::cout << "Writing everything to file and closing" << std::endl;}
  rootfile->cd();
  MC_hist->Write();
  for (int i = 0; i < evt_hists.size(); ++i) {
    evt_hists.at(i)->Write();
  }
  // raw_hist->Write();
  rootfile->Write();
  rootfile->Close();

  return 0;
}

/// Plot the hit time residuals for the MC position
///
/// @param[in] fileName of the RAT::DS root file to analyse
/// @return the histogram plot
TH1D* PlotHitTimeResidualsMCPosition( const std::string& fileName , bool verbose) {
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

      /* ~~~~~~~ MCEV info: position and delay ~~~~~~~ */

      const RAT::DS::Entry& rDS = dsReader.GetEntry(iEntry);
      const TVector3 eventPosition = rDS.GetMC().GetMCParticle(0).GetPosition(); // At least 1 is somewhat guaranteed
      std::vector<unsigned int> track_ids = rDS.GetMC().GetMCTrackIDs();

      // Find first track that is a positron
      Double_t annihilTime = -1.0;
      int posi_track = -1;
      for (unsigned int itrack = 0; itrack < track_ids.size(); itrack++) {
        RAT::DS::MCTrack mctrack = rDS.GetMC().GetMCTrack(track_ids.at(itrack));
        std::string name = mctrack.GetParticleName();

        if (name == "e+") {  // is a positron
          if (verbose) {std::cout << "Found e+ in track " << itrack << std::endl;}
          if (mctrack.GetParentID() == 0) {  // has no parent track/particle

            // Loop through steps in track
            size_t numSteps = mctrack.GetMCTrackStepCount();
            for (unsigned int istep = 0; istep < (unsigned int)numSteps; istep++) {
              RAT::DS::MCTrackStep& mcTrackStep = mctrack.GetMCTrackStep(istep);
              if (verbose) {std::cout << mcTrackStep.GetProcess() << std::endl;}

              // Check e+ dies in annihilation
              if (mcTrackStep.GetProcess() == "annihil") {
                // Time of last step of e+ track (approx time of annihilation, I hope)
                if (verbose) {std::cout<< "got annihil process for step " << istep << std::endl;}
                annihilTime = mcTrackStep.GetGlobalTime();
                posi_track = track_ids.at(itrack);
              }
            }
          }
        }
      }

      Double_t delay = 0.0;
      // If a positron that annihilated was found...
      if (posi_track != -1) {
        // Loop through the tracks to find child gamma track with latest creation time
        Double_t creationTime = annihilTime;
        Double_t tempTime;
        for (unsigned int itrack = 0; itrack < track_ids.size(); itrack++) {
          if (rDS.GetMC().GetMCTrack(track_ids.at(itrack)).GetParticleName() == "gamma") {  // check it's a gamma
            if (verbose) {std::cout << "gamma parent track id: " << rDS.GetMC().GetMCTrack(itrack).GetParentID() << std::endl;}
            if (rDS.GetMC().GetMCTrack(track_ids.at(itrack)).GetParentID() == posi_track) {  // check its parent was the initial track/particle
              // Get time of first step in track (approx track creation time, I hope)
              tempTime = rDS.GetMC().GetMCTrack(track_ids.at(itrack)).GetMCTrackStep(0).GetGlobalTime();
              if (tempTime > creationTime) {
                creationTime = tempTime;
              }
            }
          }
        }
        delay = creationTime - annihilTime;
      } 

        /* ~~~~~~~ EV info: events and PMT hits ~~~~~~~ */

      for(size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++) {
        const RAT::DS::EV& rEV = rDS.GetEV(iEV);
        const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
        std::string hist_name = "hHitTimeResidualsMC_" + std::to_string(counter);
        std::string title = "Hit time residuals using the MC position, and o-Ps delay = " + std::to_string(delay);
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