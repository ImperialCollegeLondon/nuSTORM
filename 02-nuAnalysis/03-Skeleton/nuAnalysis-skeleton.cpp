// Read the previously produced N-Tuple and print on screen
// its content

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

#include "TVector.h"
#include "TLorentzVector.h"

#include "RunControl.hpp"

Double_t nueErest(Double_t *xin, Double_t *par) {
  Double_t mmu  = 0.1056583745;
  Double_t s    = par[0];
  Double_t x    = (2.*xin[0]/mmu);
  Double_t scl  = (12.*s);
  Double_t dx   = (2.*par[1]/mmu);
  Double_t y    = scl*x*x*(1.-x)*dx;
  if ( y < 0. or x > 1.) y = 0.; 
  return y;
}

Double_t numuErest(Double_t *xin, Double_t *par) {
  Double_t mmu  = 0.1056583745;
  Double_t s    = par[0];
  Double_t x    = (2.*xin[0]/mmu);
  Double_t scl  = (2.*s);
  Double_t dx   = (2.*par[1]/mmu);
  Double_t y    = scl*x*x*(3. - 2.*x)*dx;
  if ( y < 0. or x > 1.) y = 0.; 
  std::cout << "numuErest: " << y << std::endl;
  return y;
}

void read_ntuple_from_file(){

  // Open ROOT file:
  // TFile in_file("conductivity_experiment.root");
}
 
int main(int nArgs, char **ArgV){

  // Initialse and parse input arguments:
  RunControl* RC = RunControl::getInstance();
  RC->ParseArgs(nArgs, ArgV);
  if ( RC->getDebug() )
    RC->print();
  
  // Open the ROOT file:
  TFile iRf( RC->getROOTfilename().c_str() );
  if ( RC->getDebug() ) {
    std::cout << " Print TFile: "
	      << RC->getROOTfilename().c_str() << std::endl;
    iRf.Print();
    std::cout << " List contents of TFile: "
	      << RC->getROOTfilename().c_str() << std::endl;
    iRf.ls();
  }

  // Get pointers to the Trees one by one:
  TTree *TrunInfo = (TTree *)iRf.Get("runInfo");
  TTree *Tbeam    = (TTree *)iRf.Get("beam");
  TTree *Tflux    = (TTree *)iRf.Get("flux");
  if ( RC->getDebug() ) {
    std::cout << " ----> Run info tree: " << std::endl;
    TrunInfo->Print();
    std::cout << " ----> Beam tree: " << std::endl;
    Tbeam->Print();
    std::cout << " Flux tree: " << std::endl;
    Tflux->Print();
  }

  // Loop over neutrinos in flux ntuple:
  int nEvt = Tbeam->GetEntries();
  if ( RC->getDebug() ) {
    std::cout << " ----> Beam ntuple has "
	      << nEvt << " entries" << std::endl;
  }

  Float_t nueP[4], numuP[4], muP[4];
  Tbeam->SetBranchAddress("Nue",  &nueP);
  Tbeam->SetBranchAddress("NuMu", &numuP);
  Tbeam->SetBranchAddress("muon", &muP);

   TH1F *hnueE      = new TH1F("hnueE",     "nue energy", 100, 0., 6.);
   TH1F *hmumass    = new TH1F("hmumass",   "muon mass", 100, 0.105, 0.107);
   TH1F *hnueErest  = new TH1F("hnueErest", "nue energy in rest frame",
			      100, 0., 0.07);
   TH1F *hnumuErest = new TH1F("hnumuErest", "numu energy in rest frame",
			      100, 0., 0.07);

   for (int i=0 ; i<nEvt ; i++) {
     Tbeam->GetEntry(i);
     if ( RC->getDebug() and i<10) {
       std::cout << " ----> Event: " << i << std::endl;
       std::cout << "    ----> Nue.P: ";
       for(int ii=0;ii<4;ii++){
	 std::cout << nueP[ii] << ", ";
       }
       std::cout << std::endl;
       std::cout << "    ----> NueMu.P: ";
       for(int ii=0;ii<4;ii++){
	 std::cout << numuP[ii] << ", ";
       }
       std::cout << std::endl;
       std::cout << "    ----> muon.P: ";
       for(int ii=0;ii<4;ii++){
	 std::cout << muP[ii] << ", ";
       }
       std::cout << std::endl;
     }

     hnueE->Fill(nueP[0]);

     TLorentzVector Pmu;
     Pmu.SetPxPyPzE(muP[1], muP[2], muP[3], muP[0]);
     double Mmu;
     Mmu = sqrt(Pmu*Pmu);
     if ( RC->getDebug() and i<10) {
       std::cout << "       ----> Mass of muon check: " << Mmu << std::endl;
     }
     hmumass->Fill(Mmu);

     TVector3 b;
     b = -Pmu.BoostVector();
     //b[0] = 0.;
     //b[1] = 0.;
     //b[2] = -0.999844984769905;
     if ( RC->getDebug() and i<10) {
       std::cout << "       ----> Setting up boost parameters check: "
		 << std::endl;
       std::cout << "             Boost vector: ";
       for(int ii=0;ii<3;ii++){
	 std::cout << b[ii] << ", ";
       }
       std::cout << std::endl;
     }
     TLorentzVector PmuRest = Pmu;
     PmuRest.Boost(b);
     if ( RC->getDebug() and i<10) {
       std::cout << "       ----> Muon 4-mmtm in rest frame check: ";
       for(int ii=0;ii<4;ii++){
	 std::cout << PmuRest[ii] << ", ";
       }
       std::cout << std::endl;
     }

     TLorentzVector Pnue, Pnumu;
     Pnue.SetPxPyPzE( nueP[1],  nueP[2],  nueP[3],  nueP[0]);
     Pnumu.SetPxPyPzE(numuP[1], numuP[2], numuP[3], numuP[0]);
     TLorentzVector PnueRest  = Pnue;
     TLorentzVector PnumuRest = Pnumu;
     PnueRest.Boost(b);
     PnumuRest.Boost(b);

     hnueErest->Fill(PnueRest[3]);
     hnumuErest->Fill(PnumuRest[3]);
     
   }

   Double_t dEvt = nEvt;
   Double_t dE   = 0.07/100.;

   TF1 *fnueErest  = new TF1("nueE",   nueErest,  0., 0.07, 2);
   fnueErest->SetParameters(dEvt,dE);
   fnueErest->SetLineColor(kRed); fnueErest->SetLineStyle(2);

   TF1 *fnumuErest = new TF1("nuemu", numuErest, 0., 0.07, 2);
   fnumuErest->SetParameters(dEvt,dE);
   fnumuErest->SetLineColor(kRed); fnumuErest->SetLineStyle(2);

   TCanvas *c = new TCanvas();
   hnueErest->Draw();
   fnueErest->Draw("SAME");
   c->Print("hnueErest.png");
   hnumuErest->Draw();
   fnumuErest->Draw("SAME");
   c->Print("hnumuErest.png");
   
   TFile oRf("histos.root","new");
   hnueE->Write();
   hnueErest->Write();
   hnumuErest->Write();
   hmumass->Write();

}