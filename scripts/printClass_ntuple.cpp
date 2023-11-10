// ./printClass_ntuple.exe /mnt/lustre/scratch/epp/jp643/antinu/Fisher/ReactorIBD/ReactorIBD_Fisher.ntuple.root /mnt/lustre/scratch/epp/jp643/antinu/Fisher/alphaN/alphaN_Fisher.ntuple.root /mnt/lustre/scratch/epp/jp643/antinu/Fisher/Class_results.txt
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <string>
#include <fstream>
#include <TLine.h>
#include <TTree.h>
#include <TVector3.h>
#include <RAT/DB.hh>


// define max delay, since it gets used twice (for consistency)
double MIN_PROMPT_E = 0.9, MAX_PROMPT_E = 3.5;
double MIN_DELAYED_E = 1.85, MAX_DELAYED_E = 2.4;
double MAX_DELAY = 0.8E6;
double R_CUT = 5700;


bool pass_prompt_cuts(double energy, TVector3 position) {
    if (energy < MIN_PROMPT_E) return false;  // min energy cut (MeV)
    if (energy > MAX_PROMPT_E) return false;  // max energy cut (MeV)
    if (position.Mag() > R_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts(double energy, TVector3 position) {
    if (energy < MIN_DELAYED_E) return false;  // min energy cut (MeV)
    if (energy > MAX_DELAYED_E) return false;  // max energy cut (MeV)
    if (position.Mag() > R_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_coincidence_cuts(double delay, TVector3 prompt_pos, TVector3 delayed_pos) {
    // double delay = (delayed_time - prompt_time) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
    double distance = (delayed_pos - prompt_pos).Mag();

    if (delay < 400) return false;  // min delay cut (ns)
    if (delay > MAX_DELAY) return false;  // max delay cut (ns)
    if (distance > 1500) return false;  // max distance cut (mm)

    return true;
}


/**
 * @brief Takes in an event tree, and compiles classifier histogram.
 * 
 * @param EventInfo  event tree from ntuple.
 * @return TH1D*
 */
void Apply_tagging_and_cuts(TTree* EventInfo, std::vector<double>& class_results, std::vector<TH1D*>& hists, std::string hist_name_start) {

    /* ~~~~~~~ Set up histograms ~~~~~~~ */

    hists.push_back(new TH1D((hist_name_start + "prompt_E").c_str(), "Prompt E", 100, MIN_PROMPT_E, MAX_PROMPT_E));
    hists.push_back(new TH1D((hist_name_start + "delayed_E").c_str(), "Delayed E", 100, MIN_DELAYED_E, MAX_DELAYED_E));
    hists.push_back(new TH1D((hist_name_start + "prompt_R").c_str(), "Prompt R", 100, 0.0, R_CUT));
    hists.push_back(new TH1D((hist_name_start + "delayed_R").c_str(), "Delayed R", 100, 0.0, R_CUT));
    hists.push_back(new TH1D((hist_name_start + "prompt_NhitsCleaned").c_str(), "Prompt Nhits (cleaned)", 1200, 99.5, 1299.5));
    hists.push_back(new TH1D((hist_name_start + "delayed_NhitsCleaned").c_str(), "Delayed Nhits (cleaned)", 1200, 99.5, 1299.5));
    hists.push_back(new TH1D((hist_name_start + "deltaR").c_str(), "delta_R", 100, 0.0, 1500.0));
    hists.push_back(new TH1D((hist_name_start + "deltaT").c_str(), "delta_T", 100, 200.0, MAX_DELAY));
    hists.push_back(new TH1D((hist_name_start + "prompt_Class").c_str(), "Prompt Classifier Results", 200, -50.0, 30.0));

    // Set branch addresses to unpack TTree
    TString *originReactor = NULL;
    Double_t reconEnergy;
    Double_t reconX;
    Double_t reconY;
    Double_t reconZ;
    ULong64_t eventTime;
    Bool_t valid;
    Int_t mcIndex;
    Int_t NhitsCleaned;
    Double_t classResult;

    // EventInfo->SetBranchAddress("parentKE1", &parentKE1);
    EventInfo->SetBranchAddress("energy", &reconEnergy);
    EventInfo->SetBranchAddress("posx", &reconX);
    EventInfo->SetBranchAddress("posy", &reconY);
    EventInfo->SetBranchAddress("posz", &reconZ);
    EventInfo->SetBranchAddress("clockCount50", &eventTime);
    EventInfo->SetBranchAddress("fitValid", &valid);
    EventInfo->SetBranchAddress("mcIndex", &mcIndex);
    EventInfo->SetBranchAddress("nhitsCleaned", &NhitsCleaned);
    EventInfo->SetBranchAddress("alphaNReactorIBD", &classResult);


    /* ~~~~~~~ loop through events, tag and cut ~~~~~~~ */

    unsigned int nentries = EventInfo->GetEntries();
    unsigned int nvalid = 0;
    unsigned int nvaliddelayed = 0;
    TVector3 promptPos;
    TVector3 delayedPos;
    double delayedEnergy;
    double delayedTime;
    unsigned int delayedMCIndex;
    unsigned int delayedNhitsCleaned;
    double delay;
    unsigned int a = 0;

    // Loop through all events:
    // First, find prompt event that pass all cuts.
    // Then go through next 3 events to try to find a delayed event that also passes all cuts.
    while (a < nentries) {
        if (a % 100 == 0) std::cout << "Done " << ((float)a / (float)nentries) * 100.0 << "%" << std::endl;

        EventInfo->GetEntry(a);
        delayedEnergy = reconEnergy;
        delayedPos = TVector3(reconX, reconY, reconZ);
        delayedTime = eventTime;
        delayedMCIndex = mcIndex;
        delayedNhitsCleaned = NhitsCleaned;

        if (valid and pass_delayed_cuts(delayedEnergy, delayedPos)) {
            nvaliddelayed++;

            // Delayed event is valid, check through the previous 10 events for event that passes prompt + tagging cuts
            for (unsigned int b = 1; b <= 10; ++b) {
                EventInfo->GetEntry(a - b);
                if (mcIndex != delayedMCIndex) continue;  // check just in case

                delay = (delayedTime - eventTime) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
                if (delay > MAX_DELAY) break;  // If delay becomes larger than cut, stop looking for new events

                promptPos = TVector3(reconX, reconY, reconZ);
                if (valid and pass_prompt_cuts(reconEnergy, promptPos) and pass_coincidence_cuts(delay, promptPos, delayedPos)) {
                    // Event pair survived analysis cuts
                    nvalid++;

                    // Fill hists
                    hists.at(0)->Fill(reconEnergy);
                    hists.at(1)->Fill(delayedEnergy);
                    hists.at(2)->Fill(promptPos.Mag());
                    hists.at(3)->Fill(promptPos.Mag());
                    hists.at(4)->Fill(NhitsCleaned);
                    hists.at(5)->Fill(delayedNhitsCleaned);
                    hists.at(6)->Fill((promptPos - delayedPos).Mag());
                    hists.at(7)->Fill(delay);
                    hists.at(8)->Fill(classResult);

                    // Fill vector (to print to txt file, for easier use)
                    class_results.push_back(classResult);
                }
            }
        }
        a++;
    }
    std::cout << "From " << nentries << " entries, number of valid delayed events: "<< nvaliddelayed << " number of valid event pairs surviving all cuts: " << nvalid << std::endl;
}


int main(int argv, char** argc) {
    std::string reactor_events_address = argc[1];
    std::string alphaN_events_address = argc[2];
    std::string outRoot_address = argc[3];
    std::string outText_address = argc[4];

    // Read in files and get their TTrees
    TFile *reactorFile = TFile::Open(reactor_events_address.c_str());
    TFile *alphaNFile = TFile::Open(alphaN_events_address.c_str());

    TTree *reactorEventTree = (TTree *) reactorFile->Get("output");
    TTree *alphaNEventTree = (TTree *) alphaNFile->Get("output");

    // Loop through and apply tagging + cuts
    std::cout << "Looping through reactor IBD events..." << std::endl; 
    std::vector<double> IBD_class_res;
    std::vector<TH1D*> IBD_hists;
    Apply_tagging_and_cuts(reactorEventTree, IBD_class_res, IBD_hists, "IBD_");
    std::cout << "Looping through alpha-n events..." << std::endl; 
    std::vector<double> alphaN_class_res;
    std::vector<TH1D*> alphaN_hists;
    Apply_tagging_and_cuts(alphaNEventTree, alphaN_class_res, alphaN_hists, "alphaN_");

    // Write hists to root file
    TFile rootfile(outRoot_address.c_str(), "RECREATE");
    rootfile.cd();
    for (unsigned int i = 0; i < IBD_hists.size(); ++i) {
        IBD_hists.at(i)->Write();
    }
    for (unsigned int i = 0; i < alphaN_hists.size(); ++i) {
        alphaN_hists.at(i)->Write();
    }
    rootfile.Write();
    rootfile.Close();

    // Print class results to a text file
    std::ofstream output_file;
    output_file.open(outText_address);

    output_file << IBD_class_res.at(0);
    for (unsigned int i = 1; i < IBD_class_res.size(); ++i) {
        output_file << ", " << IBD_class_res.at(i);
    }
    output_file << std::endl;
    output_file << alphaN_class_res.at(0);
    for (unsigned int i = 1; i < alphaN_class_res.size(); ++i) {
        output_file << ", " << alphaN_class_res.at(i);
    }
    output_file << std::endl;

    return 0;
}