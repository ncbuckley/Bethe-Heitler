#include "file2.h"
#include "sigmaborn.h"
#include "EnergyLookup.h"
#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <TLeaf.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TBrowser.h>
using namespace std;

int main(){

  clock_t start = clock();
  double N = 1000;
  int bins = 100;
  double me = 5.109989461e-4;
  //get largest value of y
  
  //Create a new TFile with TTree ar containing double values acc
  TFile *bh = new TFile("sigma_born.root","recreate");
  TTree *sb = new TTree("sb", "Tree with generated values");

  Double_t theta1;
  Double_t theta2;
  Double_t phi1;
  Double_t phi2;
  Double_t x;
  Double_t E0;
  Double_t y;

  Double_t atheta1;
  Double_t atheta2;
  Double_t aphi1;
  Double_t aphi2;
  Double_t ax;

  sb ->Branch("atheta1",&atheta1);
  sb ->Branch("atheta2",&atheta2);
  sb ->Branch("aphi1",&aphi1);
  sb ->Branch("aphi2",&aphi2);
  sb ->Branch("ax",&ax);

    sb ->Branch("theta1",&theta1);
    sb ->Branch("theta2",&theta2);
    sb ->Branch("phi1",&phi1);
    sb ->Branch("phi2",&phi2);
    sb ->Branch("x",&x);
    sb ->Branch("E0",&E0);


  //Create a 1d histogram 
  TFile* oFile = TFile::Open("Hist_sb.root","recreate");
  TH1D* Hist_sb = new TH1D("Hist_sb","Accepted Numbers", bins, 0, 1);
  
  //Generate random values x and y N times, if y value is below or on e^(-x/2) fill to acc
  for(Int_t i = 0; i <= N; i++){
    E0 = RandomReal(8.1,8.9);
    x = RandomReal(me/E0,1-(me/E0));
    phi1 = RandomReal(0,2*M_PI);
    phi2 = RandomReal(0,2*M_PI);
    theta1 = 1/(RandomReal( 1/(13.12*(M_PI/180)) , 1/(0.75*(M_PI/180)) ));
    theta2 = 1/(RandomReal( 1/(13.12*(M_PI/180)) , 1/(0.75*(M_PI/180)) ));
    y = RandomReal(0,0.00011111);


    if( y <= dSigmaBorn(x,phi1,phi2,theta1,theta2,E0) ){
      ax = x;
      aphi1 = phi1;
      aphi2 = phi2;
      atheta1 = theta1;
      atheta2 = theta2;
    }
    sb -> Fill();
  }

  //Take accepted x values from TTree ar and fill histogram
  for(Int_t j=0; j<N; j++){
    sb->GetEntry(j);
    Double_t entry = sb->GetLeaf("ax")->GetValue(0);
    Hist_sb->Fill(entry);
  }
  oFile->Write();
  oFile->Close();

  bh->Write();
  bh->Close();

  //Stop clock and calculate time to run
  clock_t end = clock();
  double run_time = ((end - start) / double(CLOCKS_PER_SEC));

  std::cout << std::fixed;
  std::cout << std::setprecision(5);
  std::cout <<"Time: "<< run_time << endl;
  
  return 0;
  
}
