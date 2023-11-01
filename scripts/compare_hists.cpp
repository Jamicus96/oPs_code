//g++ -g -std=c++1y -o compare_hists.exe compare_hists.cpp `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux

#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TText.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
#include <RAT/DB.hh>
#include <TMath.h>
#include <TVector3.h>



int main(){

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DEFAULT SNO+ SETTINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

    TStyle *snoStyle= new TStyle("snoplus","SNO+ plots style for publications");

    // use plain black on white colors
    snoStyle->SetFrameBorderMode(0);
    snoStyle->SetCanvasBorderMode(0);
    snoStyle->SetPadBorderMode(0);
    snoStyle->SetPadBorderSize(0);
    snoStyle->SetPadColor(0);
    snoStyle->SetCanvasColor(0);
    snoStyle->SetTitleColor(0);
    snoStyle->SetStatColor(0);
    snoStyle->SetFillColor(0);

    // use bold lines 
    snoStyle->SetHistLineWidth(2);
    snoStyle->SetLineWidth(2);

    // no title, stats box or fit as default
    snoStyle->SetOptTitle(0);
    //snoStyle->SetOptStat(0);
    //snoStyle->SetOptFit(0);
    
    // postscript dashes
    snoStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

    // text style and size
    //snoStyle->SetTextFont(132);
    //snoStyle->SetTextSize(0.24);
    snoStyle->SetLabelOffset(0.01,"x");
    snoStyle->SetTickLength(0.015,"x");
    snoStyle->SetTitleOffset(1.5,"x");
    snoStyle->SetLabelOffset(0.01,"y");
    snoStyle->SetTickLength(0.015,"y");
    snoStyle->SetTitleOffset(1.5,"y");
    snoStyle->SetLabelOffset(0.01,"z");
    snoStyle->SetTickLength(0.015,"z");
    snoStyle->SetTitleOffset(1.5,"z");
    snoStyle->SetLabelFont(132,"x");
    snoStyle->SetLabelFont(132,"y");
    snoStyle->SetLabelFont(132,"z");
    snoStyle->SetTitleFont(132,"x");
    snoStyle->SetTitleFont(132,"y");
    snoStyle->SetTitleFont(132,"z");
    snoStyle->SetLabelSize(0.04,"x");
    snoStyle->SetTitleSize(0.05,"x");
    snoStyle->SetTitleColor(1,"x");
    snoStyle->SetLabelSize(0.04,"y");
    snoStyle->SetTitleSize(0.05,"y");
    snoStyle->SetTitleColor(1,"y");
    snoStyle->SetLabelSize(0.04,"z");
    snoStyle->SetTitleSize(0.05,"z");
    snoStyle->SetTitleColor(1,"z");
    snoStyle->SetPadTickX(1);
    snoStyle->SetPadTickY(1);

    // AXIS OFFSETS
    snoStyle->SetTitleOffset(0.8,"x");
    snoStyle->SetTitleOffset(0.8,"y");
    snoStyle->SetTitleOffset(0.8,"z");

    // Legends
    snoStyle->SetLegendBorderSize(0);
    snoStyle->SetLegendFont(132);
    snoStyle->SetLegendFillColor(0);
        
    // graphs - set default marker to cross, rather than .
    snoStyle->SetMarkerStyle(21);  // filled square not .
        
    // SNO+ Preliminary label
    snoStyle->SetTextFont(132);
    snoStyle->SetTextSize(0.06);

    gROOT->SetStyle("snoplus");

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~ OUTPUT INFO ~~~~~~~~~~~~~~~~~~~~~~~~~ //

    TFile fout("/mnt/lustre/scratch/epp/jp643/antinu/AmBe/new_alpha/compare_AmBe.root", "recreate"); // Create the output file

    // Things needed later
    RAT::DB *db = RAT::DB::Get();
    db->SetAirplaneModeStatus(true);
    db->LoadDefaults();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~ ADD ONE HIST ~~~~~~~~~~~~~~~~~~~~~~~~~ //
    
    // File to read
    std::string file_address1 = "/mnt/lustre/scratch/epp/jp643/antinu/AmBe/new_alpha/tot_hists/tot_hists.root";

    // Open root file and read in hist
    TFile* f1 = TFile::Open(file_address1.c_str());
    TH1D* h1 = (TH1D*)f1->Get("prompt_t_res");
    h1->SetName("prompt_t_res_MC");

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~ ADD SECOND HIST ~~~~~~~~~~~~~~~~~~~~~~~~~ //
    
    // File to read
    std::string file_address2 = "/mnt/lustre/scratch/epp/jp643/antinu/AmBe/data_tothist/tot_hists.root";

    // Open root file and read in hist
    TFile* f2 = TFile::Open(file_address2.c_str());
    TH1D* h2 = (TH1D*)f2->Get("prompt_t_res");
    h2->SetName("prompt_t_res_data");

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~ MAKE DIFFERENCE HIST ~~~~~~~~~~~~~~~~~~~~~~~~~ //

    // Rescale histogram (Assuming they started with the same bins)
    double x_min = -30.0, x_max = 160.0;
    unsigned int bin_min = (unsigned int)(h1->GetXaxis()->GetNbins() * (x_min - h1->GetXaxis()->GetXmin()) / (h1->GetXaxis()->GetXmax() - h1->GetXaxis()->GetXmin()));
    unsigned int bin_max = (unsigned int)(h1->GetXaxis()->GetNbins() * (x_max - h1->GetXaxis()->GetXmin()) / (h1->GetXaxis()->GetXmax() - h1->GetXaxis()->GetXmin()));
    h1->Scale(h2->Integral(bin_min, bin_max) / h1->Integral(bin_min, bin_max));

    TH1D hdiff(*h2); // TH1F or TH1D, same as h1
    hdiff.SetName("prompt_t_res_data");
    if (!(hdiff.GetSumw2N() > 0)) hdiff.Sumw2(kTRUE); // ensure proper error propagation
    hdiff.Add(h1, -1.0);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~ FROMAT HISTOGRAMS ~~~~~~~~~~~~~~~~~~~~~~~~~ //

    // Reset limits
    h1->GetXaxis()->SetRangeUser(x_min, x_max);
    h2->GetXaxis()->SetRangeUser(x_min, x_max);
    hdiff.GetXaxis()->SetRangeUser(x_min, x_max);

    // Labels
    std::string y_label = "Counts";
    std::string x_label = "t_res (ns)";
    double y_axis_label_offset = 1.0, x_axis_label_offset = 1.0;

    h1->SetTitle("Data vs MC t_res Comparison");
    h1->SetLineColor(kBlue);
    h1->GetYaxis()->SetTitle(y_label.c_str());
    h1->GetYaxis()->SetTitleOffset(y_axis_label_offset);
    h1->GetXaxis()->SetTitle(x_label.c_str());
    h1->GetXaxis()->SetTitleOffset(x_axis_label_offset);
    h1->SetStats(0);

    h1->SetTitle("Data vs MC t_res Comparison");
    h2->SetLineColor(kRed);
    h2->GetYaxis()->SetTitle(y_label.c_str());
    h2->GetYaxis()->SetTitleOffset(y_axis_label_offset);
    h2->GetXaxis()->SetTitle(x_label.c_str());
    h2->GetXaxis()->SetTitleOffset(x_axis_label_offset);
    h2->SetStats(0);

    hdiff.SetTitle("");
    hdiff.SetLineColor(kBlack);
    hdiff.GetYaxis()->SetTitle("Difference");
    hdiff.GetYaxis()->SetTitleOffset(y_axis_label_offset);
    hdiff.GetXaxis()->SetTitle(x_label.c_str());
    hdiff.GetXaxis()->SetTitleOffset(x_axis_label_offset);
    hdiff.SetStats(0);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~ CREATE CANVAS AND DRAW ON IT ~~~~~~~~~~~~~~~~~~~~~~~~~ //

    TCanvas c1("c1","c1");  // Create output canvas to be saved in output file
    c1.Divide(1, 2, 0.01, 0.01);  // n_x, n_y, x_margin, y_margin

    // top canvas
    c1.cd(1);
    h1->Draw("HIST E");
    h2->Draw("E1 SAME");

    TLegend leg(0.6, 0.5, 0.8, 0.8);
    leg.AddEntry(h1, "MC", "L");
    leg.AddEntry(h2, "Data", "L");
    leg.SetTextSize(0.045);  // Default is 0.
    leg.SetTextFont(62);  // This makes it look bold
    leg.Draw();

    // bottom canvas
    c1.cd(2);
    hdiff.Draw("E1");

    // Embed text
    // TText tPrelim(3.5, 0.5, "SNO+ Preliminary");
    // tPrelim.Draw();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~ WRITE AND CLOSE ~~~~~~~~~~~~~~~~~~~~~~~~~ //

    fout.cd();
    c1.Write();  // Write canvas to root file
    h1->Write(); h2->Write(); hdiff.Write();  // Write histograms to root file too
    fout.Close();
    
    return 0;
}