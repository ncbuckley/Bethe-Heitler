#include <TLeaf.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TBrowser.h>


int EnergyLookup(double E0){

  //  int E0_bin = 0;
  Int_t bins = 200;
  Int_t E0_bin = 0;
  Int_t E0_count = 0;
  
  //get branch from given root file
  TFile* oFile= TFile::Open("GetBremDistFromHere.root");
  TTree* kin = (TTree* )oFile->Get("kin");
  TBranch* E_Beam = kin->GetBranch("E_Beam");

  Int_t total_entries = E_Beam->GetEntries();
  
  // make new histogram
  TFile* Beam_Hist = TFile::Open("Beam_Hist.root","recreate");
  TH1D* ebeamhist = new TH1D("ebeamhist","ebeamhistogram",bins,0,13);
  
  //fill histogram with data from branch
    for(Int_t j=0; j<total_entries; j++){
      kin->GetEntry(j);
      Double_t entry = kin->GetLeaf("E_Beam")->GetValue(0);
      ebeamhist->Fill(entry);
    }

    for(Int_t k=0; k<bins-1; k++){
      double binlow = ebeamhist->GetBinLowEdge(k);
      double binhigh = ebeamhist->GetBinLowEdge(k+1);

       if(binlow <= E0 and E0 < binhigh){
	 int E0_bin = ebeamhist->GetBin(k);
	 int E0_count = ebeamhist->GetBinContent(E0_bin);
         Beam_Hist->Write();
         Beam_Hist->Close();
         oFile->Close();

         return(E0_count);

 	 break;
       }
    }

  Beam_Hist->Write();
  Beam_Hist->Close();
  oFile->Close();

  return(E0_count);
}
