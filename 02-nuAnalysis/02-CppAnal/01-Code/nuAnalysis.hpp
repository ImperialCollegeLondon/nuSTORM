/*    nuAnalysis class header file */

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

class nuAnalysis {

private:

  bool        Debug;  // Global debug flag

  TChain* eventHist_ch;// Initalise attribute of type TChain for eventHist ntuple
  //TChain* runInfo_ch;  // Initalise attribute of type TChain for runInfo ntuple
  //TChain*    beam_ch;  // Initalise attribute of type TChain for beam ntuple
  //TChain*    flux_ch;  // Initalise attribute of type TChain for flux ntuple

public:
  std::vector<TH1F*> TH1Flist;
  std::vector<TH2F*> TH2Flist;
  std::vector<TF1*>  TF1list;

  nuAnalysis( bool Dbg=true );

  ~nuAnalysis() {  }

  void PreEventLoop( bool Dbg=true );

  void EventLoop( bool Dbg=true );

  void PostEventLoop( bool Dbg=true );

  void HistFitDo( );

  //--> Setters:
  void                  setDebug( bool Dbg );
  void           seteventHist_ch( TChain* eHist_ch );
  //void             setrunInfo_ch( TChain* rInfo_ch );
  //void                setbeam_ch( TChain* bm_ch );
  //void                setflux_ch( TChain* flx_ch );

  //--> Getters:
  bool                  getDebug(){ return Debug; };
  TChain*        geteventHist_ch(){ return eventHist_ch; };
  //TChain*          getrunInfo_ch(){ return runInfo_ch; };
  //TChain*             getbeam_ch(){ return beam_ch; };
  //TChain*             getflux_ch(){ return flux_ch; };

  //INCLUDE Speed of Light and Detector Distance, etc as constant parameters

  static Double_t  nueErest(Double_t *xin, Double_t *par);
  static Double_t numuErest(Double_t *xin, Double_t *par);
  TLorentzVector DetectorHitPosition(TLorentzVector Xin, TLorentzVector Pin, Double_t zPos);

};

class nuSIMtstRestFrame : public nuAnalysis{

public:
  nuSIMtstRestFrame( bool Dbg=true );
  ~nuSIMtstRestFrame() {  }

  void PreEventLoop( bool Dbg=true );
  void EventLoop( bool Dbg=true );
  void PostEventLoop( bool Dbg=true );

};

class nuSIMUserAnal : public nuAnalysis{

public:
  nuSIMUserAnal( bool Dbg=true );
  ~nuSIMUserAnal() {  }

  void PreEventLoop( bool Dbg=true );
  void EventLoop( bool Dbg=true );
  void PostEventLoop( bool Dbg=true );

};
