/*    nuAnalysis class implementation file */

#include <iostream>
#include<dirent.h>

#include <vector>

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"

//#include "TVector.h"
#include "TLorentzVector.h"

#include "RunControl.hpp"
#include "nuAnalysis.hpp"

nuAnalysis::nuAnalysis(bool Dbg) {

  Debug = Dbg;
  if ( Debug ) {
    std::cout << " nuAnalysis: create instance, start of Debug print:" << std::endl;
  }

  RunControl* RC = RunControl::getInstance();
  if ( Debug ) {
    std::cout << "     ----> RunControl: instance with parameters:" << std::endl;
    std::cout << "           Debug         : "
	      << RC->getDebug() << std::endl;
    std::cout << "           ROOT file name: "
	      << RC->getROOTfilename() << std::endl;
    std::cout << "           Directory with nTuples to chain: "
	      << RC->getCHAINdirname() << std::endl;
    std::cout << "           Output directory: "
	      << RC->getOUTPUTdirname() << std::endl;
  }

  // Initialise files for reading:
  eventHist_ch = new TChain("eventHistory");
  //runInfo_ch = new TChain("runInfo");
  //beam_ch    = new TChain("beam");
  //flux_ch    = new TChain("flux");

  // Process files in directory with "chain":
  if ( RC->getChainFlag() ){
    struct dirent *d;
    DIR *dr;
    std::string ChnDir = RC->getCHAINdirname().c_str();
    char *c = const_cast<char*>(ChnDir.c_str());
    dr = opendir(c);
    std::cout << "           Setting up to process files in: "
	      << RC->getCHAINdirname() << std::endl;
    std::cout << "               ----> Directory open: " << c << std::endl;

    std::string FilePath;
    for ( d=readdir(dr) ; d!=NULL ; d=readdir(dr) ){
      if ( strcmp(d->d_name, ".") != 0 && strcmp(d->d_name, "..") != 0 ){
	FilePath = RC->getCHAINdirname() + "/" + d->d_name;
	std::cout << "                     ----> Path: " << FilePath << std::endl;
  eventHist_ch->AddFile(FilePath.c_str());
  //runInfo_ch->AddFile(FilePath.c_str());
	//beam_ch->AddFile(FilePath.c_str());
	//flux_ch->AddFile(FilePath.c_str());
      }
    }
    closedir(dr);
    std::cout << "               <---- Directory closed: " << c << std::endl;
  }
  if ( RC->getFileFlag() ){
    eventHist_ch->AddFile(RC->getROOTfilename().c_str());
    //runInfo_ch->AddFile(RC->getROOTfilename().c_str());
    //beam_ch->AddFile(RC->getROOTfilename().c_str());
    //flux_ch->AddFile(RC->getROOTfilename().c_str());
  }
  if ( Debug ) {
    std::cout << "     ----> nuAnalysis: print TChains: " << std::endl;
    std::cout << "           Title: " << eventHist_ch->GetName() << std::endl;
    //std::cout << "           Title: " << runInfo_ch->GetName() << std::endl;
    //std::cout << "           Title: " << beam_ch->GetName()    << std::endl;
    //std::cout << "           Title: " << flux_ch->GetName()    << std::endl;
  }

}

void nuAnalysis::PreEventLoop( bool Dbg ) {

  if ( Debug ) {
    std::cout << " ----> nuAnalysis: Pre-event loop method entered:"
	      << std::endl;
  }

  if ( Debug ) {
    std::cout << " <---- Leaving Pre-event loop method." << std::endl;
  }

}

void nuAnalysis::EventLoop( bool Dbg ) {

  if ( Debug ) {
    std::cout << " ----> nuAnalysis: Event loop method entered:" << std::endl;
  }

  // Loop over neutrinos in flux ntuple:
  int nEvt = eventHist_ch->GetEntries();
  if ( Debug ) {
    std::cout << "     ----> eventHist ntuple has "
	      << nEvt << " entries" << std::endl;
  }

  for (int i=0 ; i<nEvt ; i++) {
    eventHist_ch->GetEntry(i);
    if ( Debug and i<10) {
      std::cout << "          ----> Event: " << i << std::endl;
    }

  }

  if ( Debug ) {
    std::cout << " <---- Leaving event loop method." << std::endl;
  }

}

void nuAnalysis::PostEventLoop( bool Dbg ) {

  if ( Debug ) {
    std::cout << " ----> nuAnalysis: Post event loop method entered:"
	      << std::endl;
  }
  RunControl* RC = RunControl::getInstance();

  if ( Debug ) {
    std::cout << " <---- Leaving post event loop method." << std::endl;
  }

}

void nuAnalysis::HistFitDo() {

  if ( Debug ) {
    std::cout << " ----> nuAnalysis: HistFitDo method entered:"
	      << std::endl;
  }
  RunControl* RC = RunControl::getInstance();

  std::string OutFile = RC->getOUTPUTdirname() + "nuAnalysis.root";

  TFile oRf(OutFile.c_str(),"RECREATE");
  for (int iH = 0; iH < TH1Flist.size(); iH++)
      TH1Flist[iH]->Write();
  for (int iF = 0; iF < TF1list.size(); iF++)
      TF1list[iF]->Write();

  if ( Debug ) {
    std::cout << " <---- Leaving HistFitDo method." << std::endl;
  }

}

void nuAnalysis::seteventHist_ch( TChain* eHist_ch ) {
  nuAnalysis::eventHist_ch = eHist_ch;
}

//void nuAnalysis::setrunInfo_ch( TChain* rInfo_ch ) {
//  nuAnalysis::runInfo_ch = rInfo_ch;
//}

//void nuAnalysis::setbeam_ch( TChain* bm_ch ) {
//  nuAnalysis::beam_ch = bm_ch;
//}

//void nuAnalysis::setflux_ch( TChain* flx_ch ) {
//  nuAnalysis::flux_ch = flx_ch;
//}

void nuAnalysis::setDebug(bool Dbg) {
  nuAnalysis::Debug = Dbg;
}

Double_t nuAnalysis::nueErest(Double_t *xin, Double_t *par) {
  Double_t mmu  = 0.1056583745;
  Double_t s    = par[0];
  Double_t x    = (2.*xin[0]/mmu);
  Double_t scl  = (12.*s);
  Double_t dx   = (2.*par[1]/mmu);
  Double_t y    = scl*x*x*(1.-x)*dx;
  if ( y < 0. or x > 1.) y = 0.;
  return y;
}

Double_t nuAnalysis::numuErest(Double_t *xin, Double_t *par) {
  Double_t mmu  = 0.1056583745;
  Double_t s    = par[0];
  Double_t x    = (2.*xin[0]/mmu);
  Double_t scl  = (2.*s);
  Double_t dx   = (2.*par[1]/mmu);
  Double_t y    = scl*x*x*(3. - 2.*x)*dx;
  if ( y < 0. or x > 1.) y = 0.;
  return y;
}

TLorentzVector nuAnalysis::DetectorHitPosition(TLorentzVector Xin, TLorentzVector Pin, Double_t t0, Double_t zPosDet) {
  TLorentzVector Xout;
  Double_t deltaZ;
  Double_t ds;
  deltaZ = zPosDet - Xin.Z();
  Xout.SetX(Xin.X() + Pin.Px()*deltaZ/Pin.Pz());
  Xout.SetY(Xin.Y() + Pin.Py()*deltaZ/Pin.Pz());
  Xout.SetZ(zPosDet);
  ds = sqrt(pow((Xin.X()-Xout.X()),2)+pow((Xin.Y()-Xout.Y()),2)+pow((Xin.Z()-Xout.Z()),2));
  Xout.SetT(Xin.T() + ds/0.299792458 + t0);
  return Xout;
}

nuSIMtstRestFrame::nuSIMtstRestFrame(bool Dbg) {

  nuAnalysis::setDebug(Dbg);
  if ( nuAnalysis::getDebug() ) {
    std::cout << " nuAnalysis::nuAnalysis: create instance, start of Debug print:" << std::endl;
  }

}

void nuSIMtstRestFrame::PreEventLoop( bool Dbg ) {

  if ( nuAnalysis::getDebug() ) {
    std::cout << " ----> nuAnalysis: Pre-event loop method entered:"
	      << std::endl;
  }

  // Set up histograms etc. prior to event loop:
  TH1F *hmumass    = new TH1F("hmumass", "muon mass", 150, 0., 0.150);
  nuAnalysis::TH1Flist.push_back(hmumass);

  TH1F *hnueErest  = new TH1F("hnueErest", "nue energy in rest frame", 130, 0., 65.0);
  nuAnalysis::TH1Flist.push_back(hnueErest);

  TH1F *hnumuErest = new TH1F("hnumuErest", "numu energy in rest frame", 130, 0., 65.0);
  nuAnalysis::TH1Flist.push_back(hnumuErest);

  TH1F *hmuE  = new TH1F("hmuE", "muon energy", 50, 1.5, 2.0);
  nuAnalysis::TH1Flist.push_back(hmuE);

  TH1F *hnueE  = new TH1F("hnueE", "nue energy", 50, 0., 2.0);
  nuAnalysis::TH1Flist.push_back(hnueE);

  TH1F *hnumuE = new TH1F("hnumuE", "numu energy", 50, 0., 2.0);
  nuAnalysis::TH1Flist.push_back(hnumuE);

  TH1F *hpiE = new TH1F("hpiE", "pion energy", 50, 2.5, 3.5);
  nuAnalysis::TH1Flist.push_back(hpiE);

  TH1F *hpimass = new TH1F("hpimass", "pion mass", 150, 0., 0.150);
  nuAnalysis::TH1Flist.push_back(hpimass);

  TH2F *hnueDetHit = new TH2F("hnueDetHit", "nue detector hit position", 80, -40., 40., 80, -40., 40.);
  nuAnalysis::TH2Flist.push_back(hnueDetHit);

  TH2F *hnumuDetHit = new TH2F("hnumuDetHit", "numu detector hit position", 80, -40., 40., 80, -40., 40.);
  nuAnalysis::TH2Flist.push_back(hnumuDetHit);

  int nEvt = nuAnalysis::geteventHist_ch()->GetEntries();
  Double_t dEvt = nEvt;
  Double_t dE   = 0.07/100.;

  TF1 *fnueErest  = new TF1("nueE",   nueErest,  0., 0.07, 2);
  fnueErest->SetParameters(dEvt,dE);
  fnueErest->SetLineColor(kRed); fnueErest->SetLineStyle(2);
  nuAnalysis::TF1list.push_back(fnueErest);

  TF1 *fnumuErest = new TF1("nuemu", numuErest, 0., 0.07, 2);
  fnumuErest->SetParameters(dEvt,dE);
  fnumuErest->SetLineColor(kRed); fnumuErest->SetLineStyle(2);
  nuAnalysis::TF1list.push_back(fnumuErest);

  if ( nuAnalysis::getDebug() ) {
    std::cout << "     ----> Booked "
	      << TH1Flist.size() << " histos:" << std::endl;
    for(int iH = 0; iH < TH1Flist.size(); iH++) {
      std::cout << "           Printing histo " << iH
		<< ": Title: " << TH1Flist[iH]->GetTitle() << std::endl;
    }
    std::cout << "     ----> Booked "
	      << TF1list.size() << " fits:" << std::endl;
    for(int iF = 0; iF < TF1list.size(); iF++) {
      std::cout << "           Printing fit " << iF
		<< ": Title: " << TF1list[iF]->GetTitle() << std::endl;
    }
  }

  if ( nuAnalysis::getDebug() ) {
    std::cout << " <---- Leaving Pre-event loop method." << std::endl;
  }

}

void nuSIMtstRestFrame::EventLoop( bool Dbg ) {

  if ( nuAnalysis::getDebug() ) {
    std::cout << " ----> nuAnalysis: Event loop method entered:" << std::endl;
  }

  // Loop over neutrinos in flux ntuple:
  int nEvt = nuAnalysis::geteventHist_ch()->GetEntries();
  if ( nuAnalysis::getDebug() ) {
    std::cout << "     ----> eventHist ntuple has "
	      << nEvt << " entries" << std::endl;
  }

  struct eventHistory trgt[13], prdStrght[13], piDcy[13], muPrdctn[13], piFlshNu[13], muDcy[13], ePrdctn[13], numuPrdctn[13], nuePrdctn[13], numuDtctr[13], nueDtctr[13];
  nuAnalysis::geteventHist_ch()->SetBranchAddress("target",  &trgt);
  nuAnalysis::geteventHist_ch()->SetBranchAddress("productionStraight", &prdStrght);
  nuAnalysis::geteventHist_ch()->SetBranchAddress("pionDecay", &piDcy);
  nuAnalysis::geteventHist_ch()->SetBranchAddress("muonProduction",  &muPrdctn);
  nuAnalysis::geteventHist_ch()->SetBranchAddress("piFlashNu", &piFlshNu);
  nuAnalysis::geteventHist_ch()->SetBranchAddress("muonDecay", &muDcy);
  nuAnalysis::geteventHist_ch()->SetBranchAddress("eProduction",  &ePrdctn);
  nuAnalysis::geteventHist_ch()->SetBranchAddress("numuProduction", &numuPrdctn);
  nuAnalysis::geteventHist_ch()->SetBranchAddress("nueProduction", &nuePrdctn);
  nuAnalysis::geteventHist_ch()->SetBranchAddress("numuDetector",  &numuDtctr);
  nuAnalysis::geteventHist_ch()->SetBranchAddress("nueDetector", &nueDtctr);

  TH1F *hmumass    = nuAnalysis::TH1Flist[0];
  TH1F *hnueErest  = nuAnalysis::TH1Flist[1];
  TH1F *hnumuErest = nuAnalysis::TH1Flist[2];
  TH1F *hmuE       = nuAnalysis::TH1Flist[3];
  TH1F *hnueE      = nuAnalysis::TH1Flist[4];
  TH1F *hnumuE     = nuAnalysis::TH1Flist[5];
  TH1F *hpiE       = nuAnalysis::TH1Flist[6];
  TH1F *hpimass    = nuAnalysis::TH1Flist[7];

  TH2F *hnueDetHit = nuAnalysis::TH2Flist[0];
  TH2F *hnumuDetHit= nuAnalysis::TH2Flist[1];

  TLorentzVector Xpi;
  TLorentzVector Xmu;
  TLorentzVector Xnue, Xnumu;
  TLorentzVector Ppi;
  TLorentzVector Pmu;
  TLorentzVector PmuRest;
  TLorentzVector Pnue, Pnumu;
  TLorentzVector PnueRest;
  TLorentzVector PnumuRest;
  TLorentzVector XnueDet, XnumuDet;
  double Emu;
  double Enue;
  double Enumu;
  double Epi;
  double Mmu;
  TVector3 b;

  int count_nu = 0;

  for (int i=0 ; i<nEvt ; i++) {
    nuAnalysis::geteventHist_ch()->GetEntry(i);

    Epi = sqrt(pow(trgt[12],2)+pow(trgt[7],2)+pow(trgt[8],2)+pow(trgt[9],2));

    Emu = sqrt(pow(muPrdctn[12],2)+pow(muPrdctn[7],2)+pow(muPrdctn[8],2)+pow(muPrdctn[9],2));
    Pmu.SetPxPyPzE(muPrdctn[7], muPrdctn[8], muPrdctn[9], Emu);
    Mmu = sqrt(Pmu*Pmu);

    hpiE->Fill(Epi);
    hpimass->Fill(trgt[12]);

    if (muDcy[11] > 0.){
      hmumass->Fill(Mmu);
      hmuE->Fill(Emu);
    }

    b = -Pmu.BoostVector();

    PmuRest = Pmu;
    PmuRest.Boost(b);

    Enue = sqrt(pow(nuePrdctn[12],2)+pow(nuePrdctn[7],2)+pow(nuePrdctn[8],2)+pow(nuePrdctn[9],2));
    Enumu = sqrt(pow(numuPrdctn[12],2)+pow(numuPrdctn[7],2)+pow(numuPrdctn[8],2)+pow(numuPrdctn[9],2));
    Pnue.SetPxPyPzE(nuePrdctn[7], nuePrdctn[8], nuePrdctn[9], Enue);
    Pnumu.SetPxPyPzE(numuPrdctn[7], numuPrdctn[8], numuPrdctn[9], Enumu);
    PnueRest  = Pnue;
    PnumuRest = Pnumu;
    PnueRest.Boost(b);
    PnumuRest.Boost(b);

    //if ( nuAnalysis::getDebug() and i<10) {
    if (nuePrdctn[11] > 0.) {

      Xnue.SetXYZT(nuePrdctn[3],nuePrdctn[4],nuePrdctn[5],nuePrdctn[10]);
      Xnumu.SetXYZT(numuPrdctn[3],numuPrdctn[4],numuPrdctn[5],numuPrdctn[10]);

      XnueDet = nuAnalysis::DetectorHitPosition(Xnue,Pnue,trgt[10],230.);
      XnumuDet = nuAnalysis::DetectorHitPosition(Xnumu,Pnumu,trgt[10],230.);

      if ( nuAnalysis::getDebug() && (count_nu++ < 10)) {

        std::cout << "          ----> Run:   " << nuePrdctn[0] << std::endl;
        std::cout << "          ----> Event: " << i << std::endl;
        std::cout << "          ----> Weight:   " << nuePrdctn[11] << std::endl;

        std::cout << "              ----> Nue.P: ";
        for(int ii=7;ii<10;ii++){
	         std::cout << nuePrdctn[ii] << ", ";
         }
         std::cout << std::endl;

         std::cout << "              ----> NuMu.P: ";
         for(int ii=7;ii<10;ii++){
	          std::cout << numuPrdctn[ii] << ", ";
         }
         std::cout << std::endl;

         std::cout << "              ----> muon.P: ";
         for(int ii=7;ii<10;ii++){
	          std::cout << muPrdctn[ii] << ", ";
         }
         std::cout << std::endl;

         std::cout << "                  ----> Mass of muon check: "<< Mmu << std::endl;

         std::cout <<"                  ----> Setting up boost parameters check: "<< std::endl;
         std::cout << "                      Boost vector: ";
         for(int ii=0;ii<3;ii++){
	          std::cout << b[ii] << ", ";
         }
         std::cout << std::endl;

         std::cout <<"                  ----> Muon 4-mmtm in rest frame check: ";
         for(int ii=0;ii<4;ii++){
	          std::cout << PmuRest[ii] << ", ";
         }
         std::cout << std::endl;

         std::cout <<"                  ----> NuMu 4-mmtm in rest frame check: ";
         for(int ii=0;ii<4;ii++){
	          std::cout << PnumuRest[ii] << ", ";
         }
         std::cout << std::endl;

         std::cout <<"                  ----> NuE 4-mmtm in rest frame check: ";
         for(int ii=0;ii<4;ii++){
	          std::cout << PnueRest[ii] << ", ";
         }
         std::cout << std::endl;

         // TODO: Also print out XnueDet and XnumuDet and compare with Pauls entries (tree->Show(eventNum))
         //       like I did on Friday with the momentum

         std::cout << "              ----> NuE.X @ DetectorPlane: ";
         for (int iii = 0; iii<4; iii++){
           std::cout<<XnueDet[iii]<<", ";
         }
         std::cout<<std::endl;

         std::cout << "              ----> NuMu.X @ DetectorPlane: ";
         for (int iii = 0; iii<4; iii++){
           std::cout<<XnumuDet[iii]<<", ";
         }
         std::cout<<std::endl;
      }

      hnueDetHit->Fill(XnueDet.X(),XnueDet.Y());
      hnumuDetHit->Fill(XnumuDet.X(),XnumuDet.Y());

      hnueErest->Fill(PnueRest[3]);
      hnumuErest->Fill(PnumuRest[3]);

      hnueE->Fill(Enue);
      hnumuE->Fill(Enumu);
    }

    //Pmu.SetPxPyPzE(muPrdctn[7], muPrdctn[8], muPrdctn[9], sqrt(pow(muPrdctn[12],2)+pow(muPrdctn[7],2)+pow(muPrdctn[8],2)+pow(muPrdctn[9],2)));
    //Mmu = sqrt(Pmu*Pmu);
    //if ( nuAnalysis::getDebug() and i<10) {
    //  std::cout << "                  ----> Mass of muon check: "<< Mmu << std::endl;
    //}
    //hmumass->Fill(Mmu);

    //b = -Pmu.BoostVector();
    //if ( nuAnalysis::getDebug() and i<10) {
    //  std::cout <<"                  ----> Setting up boost parameters check: "<< std::endl;
    //  std::cout << "                      Boost vector: ";
    //  for(int ii=0;ii<3;ii++){
	  //     std::cout << b[ii] << ", ";
    //  }
    //  std::cout << std::endl;
    //}

    //PmuRest = Pmu;
    //PmuRest.Boost(b);
    //if ( nuAnalysis::getDebug() and i<10) {
    //  std::cout <<"                  ----> Muon 4-mmtm in rest frame check: ";
    //  for(int ii=0;ii<4;ii++){
	  //     std::cout << PmuRest[ii] << ", ";
    //  }
    //  std::cout << std::endl;
    //}

    //Pnue.SetPxPyPzE(nuePrdctn[7], nuePrdctn[8], nuePrdctn[9], sqrt(pow(nuePrdctn[12],2)+pow(nuePrdctn[7],2)+pow(nuePrdctn[8],2)+pow(nuePrdctn[9],2)));
    //Pnumu.SetPxPyPzE(numuPrdctn[7], numuPrdctn[8], numuPrdctn[9], sqrt(pow(numuPrdctn[12],2)+pow(numuPrdctn[7],2)+pow(numuPrdctn[8],2)+pow(numuPrdctn[9],2)));
    //PnueRest  = Pnue;
    //PnumuRest = Pnumu;
    //PnueRest.Boost(b);
    //PnumuRest.Boost(b);

    //if ( nuAnalysis::getDebug() and i<10) {
    //  std::cout <<"                  ----> NuMu 4-mmtm in rest frame check: ";
    //  for(int ii=0;ii<4;ii++){
	  //     std::cout << PnumuRest[ii] << ", ";
    //  }
    //  std::cout << std::endl;
    //  std::cout <<"                  ----> NuE 4-mmtm in rest frame check: ";
    //  for(int ii=0;ii<4;ii++){
	  //     std::cout << PnueRest[ii] << ", ";
    //  }
    //  std::cout << std::endl;
    //}

    //hnueErest->Fill(PnueRest[3]);
    //hnumuErest->Fill(PnumuRest[3]);

  }

  if ( nuAnalysis::getDebug() ) {
    std::cout << " <---- Leaving event loop method." << std::endl;
  }

}

void nuSIMtstRestFrame::PostEventLoop( bool Dbg ) {

  if ( nuAnalysis::getDebug() ) {
    std::cout << " ----> nuAnalysis: Post event loop method entered:"
	      << std::endl;
  }
  RunControl* RC = RunControl::getInstance();

  TH1F *hmumass    = nuAnalysis::TH1Flist[0];
  TH1F *hnueErest  = nuAnalysis::TH1Flist[1];
  TH1F *hnumuErest = nuAnalysis::TH1Flist[2];
  TH1F *hmuE       = nuAnalysis::TH1Flist[3];
  TH1F *hnueE      = nuAnalysis::TH1Flist[4];
  TH1F *hnumuE     = nuAnalysis::TH1Flist[5];
  TH1F *hpiE       = nuAnalysis::TH1Flist[6];
  TH1F *hpimass    = nuAnalysis::TH1Flist[7];

  TF1 *fnueErest  = nuAnalysis::TF1list[0];
  TF1 *fnumuErest = nuAnalysis::TF1list[1];

  TH2F *hnueDetHit = nuAnalysis::TH2Flist[0];
  TH2F *hnumuDetHit= nuAnalysis::TH2Flist[1];

  TCanvas *c = new TCanvas();
  gErrorIgnoreLevel = kWarning;

  std::string PltFile = RC->getOUTPUTdirname() + "hmumass.png";
  hmumass->Draw();
  std::cout<<"    ----> Last bin with an entry in hmumass: "<<hmumass->FindLastBinAbove()<<std::endl;
  c->Print(PltFile.c_str());

  PltFile = RC->getOUTPUTdirname() + "hnueErest.png";
  hnueErest->Draw();
  std::cout<<"    ----> Last bin with an entry in hnueErest: "<<hnueErest->FindLastBinAbove()<<std::endl;
  fnueErest->Draw("SAME");
  c->Print(PltFile.c_str());

  PltFile = RC->getOUTPUTdirname() + "hnumuErest.png";
  hnumuErest->Draw();
  std::cout<<"    ----> Last bin with an entry in hnumuErest: "<<hnumuErest->FindLastBinAbove()<<std::endl;
  fnumuErest->Draw("SAME");
  c->Print(PltFile.c_str());

  PltFile = RC->getOUTPUTdirname() + "hmuE.png";
  hmuE->Draw();
  std::cout<<"    ----> Last bin with an entry in hmuE: "<<hmuE->FindLastBinAbove()<<std::endl;
  c->Print(PltFile.c_str());

  PltFile = RC->getOUTPUTdirname() + "hnueE.png";
  hnueE->Draw();
  std::cout<<"    ----> Last bin with an entry in hnueE: "<<hnueE->FindLastBinAbove()<<std::endl;
  c->Print(PltFile.c_str());

  PltFile = RC->getOUTPUTdirname() + "hnumuE.png";
  hnumuE->Draw();
  std::cout<<"    ----> Last bin with an entry in hnumuE: "<<hnumuE->FindLastBinAbove()<<std::endl;
  c->Print(PltFile.c_str());

  PltFile = RC->getOUTPUTdirname() + "hpiE.png";
  hpiE->Draw();
  std::cout<<"    ----> Last bin with an entry in hpiE: "<<hpiE->FindLastBinAbove()<<std::endl;
  c->Print(PltFile.c_str());

  PltFile = RC->getOUTPUTdirname() + "hpimass.png";
  hpimass->Draw();
  std::cout<<"    ----> Last bin with an entry in hpimass: "<<hpimass->FindLastBinAbove()<<std::endl;
  c->Print(PltFile.c_str());

  PltFile = RC->getOUTPUTdirname() + "hnueDetHit.png";
  hnueDetHit->Draw();
  //std::cout<<"    ----> Last bin with an entry in hnueDetHit: "<<hnueDetHit->FindLastBinAbove()<<std::endl;
  c->Print(PltFile.c_str());

  PltFile = RC->getOUTPUTdirname() + "hnumuDetHit.png";
  hnumuDetHit->Draw();
  //std::cout<<"    ----> Last bin with an entry in hnumuDetHit: "<<hnumuDetHit->FindLastBinAbove()<<std::endl;
  c->Print(PltFile.c_str());

  if ( nuAnalysis::getDebug() ) {
    std::cout << " <---- Leaving post event loop method." << std::endl;
  }

}
