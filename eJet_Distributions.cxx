#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include <TVector2.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLatex.h>
#include <TColor.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <vector>
#include <math.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <cmath>
#include <TVectorT.h>
#include <TRandom.h>

float calc_Q_square(float inE, ROOT::Math::PtEtaPhiEVector v)
{
  return 2.*inE*v.E()*(1-TMath::Abs(TMath::Cos(v.Theta())));
}
float constituent_pT_threshold(float eta,float B_field)
{
  // Minimum pT for B = 1.5 T (https://physdiv.jlab.org/DetectorMatrix/):
  float eta_bins[6] = {3.0,2.5,2.0,1.5,1.0,0.0};
  float pT_threshold_array[5] = {0.1,0.13,0.70,0.15,0.2};

  if (B_field == 3.0)
  {
    //    printf("B Field = %f \n",B_field);
    pT_threshold_array[0]=0.15; pT_threshold_array[1]=0.22; 
    pT_threshold_array[2] = 0.16; pT_threshold_array[3]=0.3;
    pT_threshold_array[4]=0.4;
  }
  float pT_threshold = 0;
  for (int i = 0; i < 5; i++)
    if ( (eta < eta_bins[i])&&(eta > eta_bins[i+1]) )
      pT_threshold = pT_threshold_array[i];

  /* std::cout<<"pt threshold = "<<pT_threshold<<std::endl; */
  return pT_threshold;
}

float smear_E(float E)
{
  //smears electron energy according to backward calorimeter in EIC yellow report  
  //http://www.eicug.org/web/sites/default/files/Yellow_Report_v1.1.pdf Page 21.

  float Sigma_EoverE = sqrt(pow(0.02/sqrt(E),2)+pow(0.03,2));
  float Sigma_E = Sigma_EoverE*E;
  float E_smeared =gRandom->Gaus(E,Sigma_E);
  /* fprintf(stderr,"%d: E Truth = %1.3f : E_Sigma = %1.3f : E Smear = %1.3f\n",__LINE__,E,Sigma_E,E_smeared); */

  return E_smeared;
}

bool EoverP_2Sigma(float EoverP)
{
//determines if input is within 3 sigma of electron E/P Curve.
//Curve is a normal distribution obtained from smearing E_Truth
  float Sigma = 0.03;
  float Mean = 1.002; 
  bool within_2sigma = (EoverP >= Mean-2*Sigma && EoverP <= Mean+2*Sigma); 
  return within_2sigma;
}

float dphi_res(float B, float jet_phi, float jet_eta, float e_phi, float e_eta)
{//returns relative error based on 1-sigma of phi resolution


  float mrad_to_rad = 1000;

  //electron Sigma
  float e_sigma = 0;  //electron phi res [mrad]
  if (e_eta < -0.5) e_sigma = 4.;
  else if (e_eta > 0.5) e_sigma = 0.6;
  else e_sigma = 4.;
  e_sigma = e_sigma/mrad_to_rad;
  
  //Jet Sigma
  float eta_bin[] = {-3.0,-1.5,-.5,0.5,1.5,3.0};
  int n_eta = sizeof(eta_bin)/sizeof(*eta_bin);
  float jet_sigma_array[] = {11.187,13.351,13.445,14.039,12.539};//1.4T B Field
  if (B==3.0) {
    jet_sigma_array[0] = 12.294;
    jet_sigma_array[1] = 15.815;
    jet_sigma_array[2] = 16.427;
    jet_sigma_array[3] = 16.237;
    jet_sigma_array[4] = 14.236;
  }

  float jet_sigma=0;
  for (int i_eta = 0; i_eta < n_eta-1; i_eta++)
    if (jet_eta > eta_bin[i_eta] && jet_eta >= eta_bin[i_eta+1])
      jet_sigma = jet_sigma_array[i_eta];

  jet_sigma = jet_sigma/mrad_to_rad;

  float dphi_percent = sqrt( pow(jet_sigma/jet_phi,2) + pow(e_sigma/e_phi,2));
  /* std::cout<<__LINE__<<"Delta Phi (%) Sigma = "<<dphi_percent<<" Jet Eta = "<<jet_eta<<std::endl; */  

  return dphi_percent;

}

using namespace std;
int main(int argc, char *argv[])
{
  if (argc < 3) {
    std::cout<<"Syntax: [Command] [File] [B Field]"<<std::endl;
    exit(EXIT_FAILURE);
  }
  //for (int iarg = 1; iarg < argc; iarg++) {
  int iarg = 1;
  TString root_file = (TString)argv[iarg];

  float B_Field = atof(argv[2]);
  std::cout<< "B = "<<B_Field<<std::endl;
  std::cout << "Opening: " << (TString)argv[iarg] << std::endl;

  float scale_factor = 1.;
  if (B_Field == 3.0)
    scale_factor = 0.4010189811930663;
  else //1.4T
    scale_factor = 0.4165767978773618;
  std::cout<<"Applying Scale Factor of "<<scale_factor<<" to simulate EIC Luminosity \n";

  TFile *file = TFile::Open(root_file);

  if (file == NULL) {
    std::cout << " File Fail" << std::endl;
    exit(EXIT_FAILURE); 
  } 

  file->Print();

  TTree *_tree_event = dynamic_cast<TTree *>(file->Get("T"));

  if (_tree_event == NULL) {
    std::cout << " Tree Fail " << std::endl;
    exit(EXIT_FAILURE);
  }

  // enum {MaxNumJets = 20,kMaxConstituents = 100};
  TString infile = argv[1];
  TFile * F = new TFile(infile);
  if (F == NULL) { std::cout << " File Fail" << std::endl; exit(EXIT_FAILURE); }

  TTree *T = dynamic_cast<TTree *>(F->Get("T"));
  if (_tree_event == NULL) { std::cout << " Tree Fail " << std::endl; exit(EXIT_FAILURE); }
  //TTree * T = (TTree*) F -> Get("T");
  Int_t njets,nAlltruthjets;
  int nEntries = T -> GetEntries();
  nEntries *= scale_factor; //Applies EIC Luminosity Scale Factor

  const int MaxNumJets = 20;
  const int kMaxConstituents = 100;

  array<Float_t, MaxNumJets> E,Eta,Phi,Pt,gE,gEta,gPhi,gPt;
  array<Float_t, MaxNumJets> all_truthE,all_truthP,all_truthPt,all_truthPhi,all_truthEta;
  array<Int_t, MaxNumJets> NComponent,gNComponent, alltruth_NComponent;

  Float_t electron_gE,electron_gEta,electron_gPhi,electron_gPt;
  Float_t electron_E,electron_Eta,electron_Phi,electron_Pt;

  array<array<Float_t, kMaxConstituents >, MaxNumJets > gComponent_Eta,gComponent_PID,
    gComponent_Pt,gComponent_Phi,gComponent_E, gComponent_Charge;

  array<array<Float_t, kMaxConstituents >, MaxNumJets > Component_Eta,Component_PID,
    Component_Pt,Component_Phi,Component_P, Component_Charge;

  array<array<Float_t, kMaxConstituents >, MaxNumJets > all_Component_Eta,all_Component_PID, all_Component_E,
    all_Component_Pt,all_Component_Phi,all_Component_P, all_Component_Charge;

  T -> SetBranchAddress("njets",&njets);
  T -> SetBranchAddress("e",&E);
  T -> SetBranchAddress("eta",&Eta);
  T -> SetBranchAddress("phi",&Phi);
  T -> SetBranchAddress("pt",&Pt);
  T -> SetBranchAddress("nComponent",&NComponent);

  bool reco_branches = true;
  if (reco_branches){
    T -> SetBranchAddress("Constituent_recoEta", Component_Eta.data());
    T -> SetBranchAddress("Constituent_recoPt",Component_Pt.data());
    T -> SetBranchAddress("Constituent_recoP",Component_P.data());
    T -> SetBranchAddress("Constituent_recoPhi",Component_Phi.data());

    T -> SetBranchAddress("electron_recoE",&electron_E);
    T -> SetBranchAddress("electron_recoEta",&electron_Eta);
    T -> SetBranchAddress("electron_recoPhi",&electron_Phi);
    T -> SetBranchAddress("electron_recoPt",&electron_Pt);
  }

  T -> SetBranchAddress("matched_truthE",&gE);
  T -> SetBranchAddress("matched_truthEta",&gEta);
  T -> SetBranchAddress("matched_truthPhi",&gPhi);
  T -> SetBranchAddress("matched_truthPt",&gPt);
  T -> SetBranchAddress("matched_truthNComponent",&gNComponent);

  T -> SetBranchAddress("matched_Constituent_truthEta", gComponent_Eta.data());
  T -> SetBranchAddress("matched_Constituent_truthPID",gComponent_PID.data());
  T -> SetBranchAddress("matched_Constituent_truthPt",gComponent_Pt.data());
  T -> SetBranchAddress("matched_Constituent_truthE",gComponent_E.data());
  T -> SetBranchAddress("matched_Constituent_truthPhi",gComponent_Phi.data());
  T -> SetBranchAddress("matched_Constituent_truthCharge",gComponent_Charge.data());

  T -> SetBranchAddress("nAlltruthjets", &nAlltruthjets);
  T -> SetBranchAddress("all_truthE",&all_truthE);
  T -> SetBranchAddress("all_truthEta",&all_truthEta);
  T -> SetBranchAddress("all_truthPhi",&all_truthPhi);
  T -> SetBranchAddress("all_truthPt",&all_truthPt);
  T -> SetBranchAddress("all_truthNComponent",&alltruth_NComponent);

  T -> SetBranchAddress("all_Constituent_truthEta", all_Component_Eta.data());
  T -> SetBranchAddress("all_Constituent_truthPID",all_Component_PID.data());
  T -> SetBranchAddress("all_Constituent_truthPt",all_Component_Pt.data());
  T -> SetBranchAddress("all_Constituent_truthE",all_Component_E.data());
  T -> SetBranchAddress("all_Constituent_truthPhi",all_Component_Phi.data());
  T -> SetBranchAddress("all_Constituent_truthCharge",all_Component_Charge.data());

  T -> SetBranchAddress("electron_truthE",&electron_gE);
  T -> SetBranchAddress("electron_truthEta",&electron_gEta);
  T -> SetBranchAddress("electron_truthPhi",&electron_gPhi);
  T -> SetBranchAddress("electron_truthPt",&electron_gPt);

  /* T->SetBranchAddress("all_truthE", all_truthE.data()); */

  //gStyle Plotting
  gStyle->SetOptStat("emr");
  gStyle->SetStatY(0.85);
  gStyle->SetStatX(0.87);
  gStyle->SetStatW(0.15);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes  
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  gStyle->SetMarkerColor(4);
  gStyle->SetLineColor(4);

  float central = 0.5;
  //central (abs() < central)
  //forward (> central)
  //backward (< -central)

  //2D Histos
  TH2F * Tjve = new TH2F("ETrueJet_vs_Eelectron", "E^{True}_{Jet} (|#eta^{Jet}|<0.7) vs. E_{e}^{True}",100,0,25,100,0,25);
  TH2F * Rjve = new TH2F("ERecoJet_vs_Eelectron", "E^{Reco}_{Jet} (|#eta^{Jet}|<0.7) vs. E_{e}^{True}",100,0,25,100,0,25);

  //Ratio Histos
  TH1F * eoTj = new TH1F("ETrueJet_over_Eelectron", "E_{Reco}^{True}/E^{True}_{e} (|#eta^{Jet}|<0.7)",80,0,2);
  TH1F * eoRj = new TH1F("Eelectron_over_ERecoJet_over_Eelectron", "E_{Jet}^{Reco}/E^{True}_{e} (|#eta^{Jet}|<0.7)",80,0,2);
  //TH1F * RjoTj = new TH1F("PRecoJet_over_PTrueJet", "P_{Jet}^{True} - P^{Reco}_{Jet} / P^{True}_{Jet}",80,-0.4,0.4);
  TH1F * RjoTj = new TH1F("dP_P", "P_{Jet}^{True} - P^{Reco}_{Jet} / P^{True}_{Jet}",100,-0.25,0.25);
  RjoTj->GetYaxis()->SetTitle("Counts");
  TH1F * justE = new TH1F("Reco Energy", "Energy",100,0,100);

  //Difference Histos
  TH1F * emTj = new TH1F("Eelectron_minus_PTrueJet", "E_{e}^{True} - E^{True}_{Jet} (|#eta^{Jet}|<0.7)",100,-20,30);
  TH1F * emRj = new TH1F("Eelectron_minus_ERecoJet", "E_{e}^{True} - E^{Reco}_{Jet} (|#eta^{Jet}|<0.7)",100,-20,30);
  //Detector Coordinate Histos
  TH1F * dEtaTj = new TH1F("dEta_e_TrueJet", "|#Delta#eta| (#eta_{e} - #eta^{True}_{Jet})", 80,-10,10);
  TH1F * dEtaRj = new TH1F("dEta_e_RecoJet", "|#Delta#eta| (#eta_{e} - #eta^{Reco}_{Jet})", 80,-10,10);
  TH1F * dPhiTj_whole = new TH1F("dPhi_e_TrueJet", "|#Delta#varphi| (#varphi_{e} - #varphi^{True}_{Jet}) ", 128,0,M_PI);
  TH1F * dPhiRj_whole = new TH1F("dPhi_e_RecoJet", "|#Delta#varphi| #varphi_{e} - #varphi(Jet^{Reco}_{Jet}) ", 128,0,M_PI);
  TH1F *aT_edPhi_whole = new TH1F("all_dPhi_e_TrueJet", "|#Delta#varphi| (#varphi_{e} - #varphi^{True}_{Jet}) ", 128,0,M_PI);

  float eta_bin[] = {-3.0,-.5,0.5,3.0};//backward,central,forword
  const int size_eta_bin = sizeof(eta_bin)/sizeof(*eta_bin)-1;

  TH1F ** dPhiTj= new TH1F*[size_eta_bin];
  TH1F ** dPhiRj= new TH1F*[size_eta_bin];
  TH1F ** aT_edPhi= new TH1F*[size_eta_bin];

  for(int et = 0 ; et < size_eta_bin ; et++){
    dPhiTj[et] = new TH1F(Form("dPhi_e_TrueJet_eta_%i",et), "|#Delta#varphi| (#varphi_{e} - #varphi^{True}_{Jet}) ", 128,0,M_PI);
    dPhiRj[et] = new TH1F(Form("dPhi_e_RecoJet_eta_%i",et), "|#Delta#varphi| #varphi_{e} - #varphi(Jet^{Reco}_{Jet}) ", 128,0,M_PI);
    aT_edPhi[et] = new TH1F(Form("all_dPhi_e_TrueJet_eta_%i",et), "|#Delta#varphi| (#varphi_{e} - #varphi^{True}_{Jet}) ", 128,0,M_PI);
  }

  //Simple Distributions
  // gStyle->SetMarkerColor(4);
  // gStyle->SetLineColor(4);
  TH1F *reco_phi = new TH1F("reco_phi","Reconstructed Jet #varphi",16,-M_PI,M_PI); 
  TH1F *reco_eta = new TH1F("reco_eta","Recontsructed Jet #eta",50,-5,5);
  TH1F *reco_E = new TH1F("reco_E","Reconstructed Jet Energy",100,0,50);
  TH1F *reco_nconst = new TH1F("reco_nconst","Reconstructed Jet N Component",20,0,20);
  TH1F *reco_P = new TH1F("reco_P","Reconstructed Jet Momentum",100,0,50);
  TH1F *nconst_diff = new TH1F("nconst_diff","Truth - Reco No. Constituents",20,0,20);
  TH1F *comp_eta = new TH1F("truth_comp_eta","Jet Truth Component #eta",80,-4,4);
  TH1F *comp_pid = new TH1F("truth_comp_pid","Jet Truth Component PID",1000,-500,500);
  TH1F *comp_pt = new TH1F("truth_comp_pt","Jet Truth Component p_{T}",500,0,50);

  TH1F *reco_comp_eta = new TH1F("reco_comp_eta","Jet Reco Component #eta",240,-4,4);
  TH1F *reco_comp_pt = new TH1F("reco_comp_pt","Jet Reco Component p_{T}",500,0,50);
  TH1F* Q2 = new TH1F("Q2","Q^{2}",100,0,500);

  TH2F *P_Component_vs_JetEta = new TH2F("comp_p_vs_jetEta","Component P vs. #eta^{Jet}_{Truth}",25,-5,5,40,0,20);
  TH2F *Pt_Component_vs_JetEta = new TH2F("comp_pT_vs_jetEta","Component p_{T} vs. #eta^{Jet}_{Truth}",25,-5,5,40,0,20);

  TH2F *n_neutrals_vs_dP_P = new TH2F("n_neutrals_vs_dP_P", "No. Neutral Components in Original Truth Jet VS. dP/P",80,-0.4,0.4,9,-1,8);
  TH2F *n_missed_vs_dP_P = new TH2F("n_missed_vs_dP_P", "No. Missed Jet Constituents VS. dP/P",80,-0.4,0.4,16,-8,8);  
  TH1F *n_neutrals = new TH1F("n_neutrals", "No. Original Neutral Components in Truth Jet",10,0,10);
  TH1I *n_charged = new TH1I("n_charged", "No. Original Charged Components in Truth Jet",10,0,10);
  TH1I *n_missed = new TH1I("n_missed","No. Missed Jet Components",10,-5,5);
  TH1I *n_constituents = new TH1I("n_constituents","No. Jet Components",10,0,10);
  TH1I *n_reco_constituents = new TH1I("reco_n_constituents","No. Jet Components",10,0,10);
  TH1I *n_missed_minuseutrals = new TH1I("n_missed_minusneutrals","No. Missed Jet Components (Truth - Reco - Neutral)",10,-5,5);
  TH1I *simple_ndiff = new TH1I("simple_ndiff", "gNComponent - NComponent", 10,-5,5);

  TH1F *central_phi_percent = new TH1F("central_dphi_percent","Sigma Phi (percent)",200,0,100); 
  TH1F *backwards_phi_percent = new TH1F("backwards_dphi_percent","Sigma Phi (percent)",200,0,100); 
  TH1F *forward_phi_percent = new TH1F("forward_dphi_percent","Sigma Phi (percent)",200,0,100); 
  TH1F *whole_phi_percent = new TH1F("dphi_percent","Sigma Phi (percent)",200,0,100); 

  int nEta_bins = 8;
  int root_colors[8] = {6,52,60,70,80,90,94,100};
  int Eta_Bins[9] = {-4,-3,-2,-1,0,1,2,3,4};//Extend past normal range (abs(3)) for viewing
  /* int Eta_Bins[2] = {-1,1}; */
  TH1F ** momentum_in_eta_bins = new TH1F*[nEta_bins];
  for (int ieta=0; ieta < nEta_bins; ++ieta)
  {
    momentum_in_eta_bins[ieta] = new TH1F(Form("momentum_etaBin_%i",ieta),
        Form("Jet p_{T} Distribution, %i <#eta_{jet} < %i",
          Eta_Bins[ieta],Eta_Bins[ieta+1]),50,0,25);
    momentum_in_eta_bins[ieta]->SetMarkerColor(4);
    momentum_in_eta_bins[ieta]->SetLineColor(4);
  }

  TH1F *truth_phi = new TH1F("truth_phi","Truth Jet #varphi",16,-M_PI,M_PI); 
  TH1F *truth_eta = new TH1F("truth_eta","Truth Jet #eta",50,-5,5);
  TH1F *truth_E = new TH1F("truth_E","Truth Jet Energy",100,0,50);
  TH1F *truth_P = new TH1F("truth_P","Truth Jet Momentum",100,0,50);

  gStyle->SetMarkerColor(2);
  gStyle->SetLineColor(2);
  TH1F *reco_phi_anticut = new TH1F("reco_phi(anti-cut)","Reconstructed Jet #varphi(anti-cut)",16,-M_PI,M_PI); 
  TH1F *reco_eta_anticut = new TH1F("reco_eta(anti-cut)","Recontsructed Jet #eta(anti-cut)",50,-5,5);
  TH1F *reco_E_anticut = new TH1F("reco_E(anti-cut)","Reconstructed Jet Energy(anti-cut)",100,0,50);
  TH1F *reco_nconst_anticut = new TH1F("reco_nconst(anti-cut)","Reconstructed Jet N Component(anti-cut)",20,0,20);
  TH1F *reco_P_anticut = new TH1F("reco_P(anti-cut)","Reconstructed Jet Momentum(anti-cut)",100,0,50);
  TH1F *RjoTj_anticut = new TH1F("PRecoJet_over_PTrueJet_anticut", "P_{Jet}^{True} - P_{Jet}^{Reco} / P^{True}_{Jet} (anti-cut) ",100,-0.25,0.25);
  TH1F *nconst_diff_anticut = new TH1F("nconst_diff_anticut","Truth - Reco No. Constituents",20,0,20);
  TH1F *comp_eta_anticut = new TH1F("comp_eta_anticut","Jet Component #eta (Anti-Cut)",40,-4,4);
  TH1F *comp_pt_anticut = new TH1F("comp_pt_anticut","Jet Component p_{T} (Anti-Cut)",500,0,50);
  TH1F *comp_pid_anticut = new TH1F("comp_pid_anticut","Jet Component PID (Anti-Cut)",1000,-500,500);
  TH1F *Q2_anticut = new TH1F("Q2_anticut","Q^{2} (anticut)",100,0,500);
  //TH1I * PID = new TH1I("PID_Histo", "PID",1000,-500,500);
  TH2F *momentum_response = new TH2F("momentum_response","Charged Jet Momentum Response",100,0,50,100,0,50);
  momentum_response->GetXaxis()->SetTitle("p_{Jet}^{Reco} [GeV/c]");
  momentum_response->GetYaxis()->SetTitle("p_{Jet}^{Truth,Ch} [GeV/c]");
  TH2F *P_Component_vs_JetEta_anticut = new TH2F("comp_p_vs_jetEta_anticut","Component P vs. #eta^{Jet}_{Truth} (anti-cut)",25,-5,5,40,0,20);
  TH2F *Pt_Component_vs_JetEta_anticut = new TH2F("comp_pT_vs_jetEta_anticut","Component p_{T} vs. #eta^{Jet}_{Truth} (anticut)",25,-5,5,40,0,20);
  TH1F ** anticut_momentum_in_eta_bins = new TH1F*[nEta_bins];
  for (int ieta=0; ieta < nEta_bins; ++ieta)
  {
    anticut_momentum_in_eta_bins[ieta] = new TH1F(Form("anticut_momentum_etaBin_%i",ieta),
        Form("Jet p_{T} Distribution (anticut), %i <#eta_{jet} < %i",
          Eta_Bins[ieta],Eta_Bins[ieta+1]),50,0,25);
    anticut_momentum_in_eta_bins[ieta]->SetMarkerColor(2);
    anticut_momentum_in_eta_bins[ieta]->SetLineColor(2);
  }

  TH1F *truth_phi_anticut = new TH1F("anticut_truth_phi","Truth Jet #varphi (anticut)",16,-M_PI,M_PI); 
  TH1F *truth_eta_anticut = new TH1F("anticut_truth_eta","Truth Jet #eta (anticut)",50,-5,5);
  TH1F *truth_E_anticut = new TH1F("anticut_truth_E","Truth Jet Energy (anticut)",100,0,50);
  TH1F *truth_P_anticut = new TH1F("anticut_truth_P","Truth Jet Momentum (anticut)",100,0,50);


  //Electron Histograms
  TH1F *truth_electron_E = new TH1F("truth_eletron_E","E_{Truth}^{e}",300,0,60);
  TH1F *smeared_electron_E = new TH1F("smeared_electron_E","E_{smeared}^{e}",300,0,60);
  TH1F *electron_EoverP = new TH1F("electron_EoverP","Electron E/P",400,0.8,1.2);
  TH1F *smeared_EoverP = new TH1F("smeared_electron_EoverP","Smeared Electron E/P",400,0.8,1.2);

  TH1F *truth_FF = new TH1F("truth_fragmentation_fuction","Truth Jet Charged Fragmentation Function",24,0,1);
  TH1F *truth_FF_zT = new TH1F("truth_fragmentation_fuction_zT","Truth Jet Charged Fragmentation Function (zT)",24,0,1);
  TH1F *reco_FF = new TH1F("fragmentation_fuction","Reco Jet Charged Fragmentation Function",24,0,1);
  TH1F *reco_FF_zT = new TH1F("fragmentation_fuction_zT","Reco Jet Charged Fragmentation Function (zT)",24,0,1);

  TH1F *all_truth_FF = new TH1F("all_truth_fragmentation_fuction","Truth Jet Charged Fragmentation Function",24,0,1);
  TH1F *all_truth_FF_zT = new TH1F("all_truth_fragmentation_fuction_zT","Truth Jet Charged Fragmentation Function (zT)",24,0,1);

  TH1F *All_truth_E_histo = new TH1F("all_truth_E","all_truth Jet Energy",100,0,50);


  /* TH1F * dPhiTj_DeltaR = new TH1F("DR_dPhi_e_TrueJet", "|#Delta#varphi| (#varphi_{e} - #varphi^{True}_{Jet}) ", 128,0,M_PI); */
  TH1F * dPhiRj_DeltaR = new TH1F("DR_dPhi_e_RecoJet", "|#Delta#varphi| #varphi_{e} - #varphi(Jet^{Reco}_{Jet}) ", 128,0,M_PI);
  TH1F *reco_FF_DeltaR = new TH1F("DR_fragmentation_fuction","Reco Jet Charged Fragmentation Function",24,0,1);
  /* TH1F *aT_edPhi_DeltaR = new TH1F("DR_all_dPhi_e_TrueJet", "|#Delta#varphi| (#varphi_{e} - #varphi^{True}_{Jet}) ", 128,0,M_PI); */

  //Smearing by DeltaPhi Resolution
  TH1F * dPhiTj_sigmaplus = new TH1F("sigmaplus_dPhi_e_TrueJet", "|#Delta#varphi| (#varphi_{e} - #varphi^{True}_{Jet}) ", 128,0,M_PI);
  TH1F * dPhiTj_sigmamin = new TH1F("sigmamin_dPhi_e_TrueJet", "|#Delta#varphi| (#varphi_{e} - #varphi^{True}_{Jet}) ", 128,0,M_PI);
  TH1F * dPhiRj_sigmaplus = new TH1F("sigmaplus_dPhi_e_RecoJet", "|#Delta#varphi| #varphi_{e} - #varphi(Jet^{Reco}_{Jet}) ", 128,0,M_PI);
  TH1F * dPhiRj_sigmamin= new TH1F("sigmamin_dPhi_e_RecoJet", "|#Delta#varphi| #varphi_{e} - #varphi(Jet^{Reco}_{Jet}) ", 128,0,M_PI);

  float DPhi_res = 0.015059458510564696; //1.4 T field
  if (B_Field == 3.0)
    DPhi_res = 0.017602173880750444;
    

  //--------------------Cut Parameters--------------------//
  float max_DeltaR = 0.1; //reco-truth match
  //float max_dE_E = 0.03;
  int min_comp = 4;
  float minE = 4.0;
  float jet_cut_counter[3] = {0};
  float max_miss_const = 1;
  int n_neutral_max = 1;

  //--------------------Event & Jet Loops--------------------//

  for (Long64_t ev = 0; ev < nEntries; ev++){
    T->GetEntry(ev);
    if (ev%10000==0) fprintf(stderr,"%d: Entry %lli out of %d\n",__LINE__,ev,nEntries);
    //if (ev == 500000) break;

    //ELECTRON SMEARING
    ROOT::Math::PtEtaPhiEVector electron_gLorentz(electron_gPt,electron_gEta,electron_gPhi,electron_gE);//Truth
    ROOT::Math::PtEtaPhiMVector electron_Lorentz(electron_Pt,electron_Eta,electron_Phi,0.0005);//Mass  of e. Reco (electron E is not reco well)

    truth_electron_E->Fill(electron_gE);
    electron_EoverP->Fill(electron_gE/electron_gLorentz.P());
    float_t electron_E_smeared = smear_E(electron_gE);

    smeared_electron_E->Fill(electron_E_smeared);
    float e_EoverP = electron_E_smeared/electron_Lorentz.P();
    smeared_EoverP->Fill(e_EoverP);
    if (not(EoverP_2Sigma(e_EoverP)))
      continue;

    float max_full_E = 0;
    for (int ia = 0; ia < nAlltruthjets; ia++){
      if (isnan(all_truthE[ia])) continue;


      if (alltruth_NComponent[ia] < 4) continue;
      if (all_truthE[ia] < minE) continue;

      if (all_truthE[ia] > max_full_E)
        max_full_E = all_truthE[ia];
    
    }//all_truth_loop

    //There are about 4X more "all-truth" jets than reco, so taking max energy jets rather than background jets 
    for (int ia = 0; ia < nAlltruthjets; ia++){
      if (alltruth_NComponent[ia] < 4) continue;
      if (all_truthE[ia] < minE) continue;

      if (isnan(all_truthE[ia])) continue;
      if (all_truthE[ia] < max_full_E) continue;

      ROOT::Math::PtEtaPhiEVector ALorentz(all_truthPt[ia],all_truthEta[ia],all_truthPhi[ia],all_truthE[ia]);
      All_truth_E_histo->Fill(ALorentz.E());
      if (isnan(ALorentz.Phi())) continue;

      bool pt_const_cut = false;
      //Uncomment to take charged_all truth
      for (int i = 0; i < alltruth_NComponent[ia]; i++)
      {
        pt_const_cut = (all_Component_Pt[ia][i] < .1);
        if (pt_const_cut) continue;
      /*   ROOT::Math::PtEtaPhiEVector AConstLorentz(all_Component_Pt[ia][i],all_Component_Eta[ia][i],all_Component_Phi[ia][i],all_Component_E[ia][i]); */
      /*   /1* if (all_Component_Charge[ia][i] == 0) *1/ */
      /*     /1* ALorentz -= AConstLorentz;//done before jet histograms are filled *1/ */
      /*   //Edited 6/16: Taking a look at full-truth jets */
      }
      if (pt_const_cut) continue;
      Float_t AllTrue_DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(ALorentz.Phi() - electron_gPhi - TMath::Pi()));
      aT_edPhi_whole->Fill(AllTrue_DeltaPhi);
       
      for(int et = 0 ; et < size_eta_bin ; et++)
        if (ALorentz.Eta() >= eta_bin[et] && ALorentz.Eta() < eta_bin[et+1])
          aT_edPhi[et]->Fill(AllTrue_DeltaPhi);
      
      for (int i = 0; i < alltruth_NComponent[i]; i++)
      {
        if (all_Component_Charge[ia][i] != 0)
        {
          ROOT::Math::PtEtaPhiEVector AConstLorentz(all_Component_Pt[ia][i],all_Component_Eta[ia][i],all_Component_Phi[ia][i],all_Component_E[ia][i]);
          all_truth_FF->Fill(ALorentz.P()/AConstLorentz.P());
          all_truth_FF_zT->Fill(ALorentz.Pt()/AConstLorentz.Pt());
        }
      }
    }//all_truth_loop

    //RECO Jet Loop
    for (int n = 0; n < njets; ++n) {

      if (std::isnan(gE[n])) continue; //reco jet must have a truth match
      simple_ndiff->Fill(gNComponent[n] - NComponent[n]);

      ROOT::Math::PtEtaPhiEVector Lorentz(Pt[n],Eta[n],Phi[n],E[n]);
      ROOT::Math::PtEtaPhiEVector gLorentz(gPt[n],gEta[n],gPhi[n],gE[n]);
      ROOT::Math::PtEtaPhiEVector e_vector (electron_gPt,electron_gEta,electron_gPhi,electron_gE);


      //--------------------Reco Jet CUTS--------------------//

      /* if (dR > max_DeltaR) continue; */
      if (E[n] < minE) continue;
      n_reco_constituents->Fill(NComponent[n]);
      if (NComponent[n] < min_comp) continue;

      //--------------------Truth Cuts--------------------//
      //if ((gE[n] < 15) || (gE[n] > 20)) continue;
      // bool gconst_eta = true;
      // for (int i = 0; i < gNComponent[n]; i++){
      //   if (abs(gComponent_Eta[n][i]) > 3.) gconst_eta=false;
      // 	if (!gconst_eta) break;
      // }			  
      // //if !(gconst_eta) continue

      //--------------------Constituent Cuts & Counting-----------------------//
      bool eta_const_cut = false; //cut based on radiation length table, where barrel meets endcap
      bool pt_const_cut = false;//cut away helixes
      int n_neutral = 0;
      int n_ch = 0;
      bool maxed_neutrals = true;
      float constant_p_cut =false;
      for (int i = 0; i < NComponent[n]; i++){

        eta_const_cut = (  ( (abs(Component_Eta[n][i]) > 1.06) && (abs(Component_Eta[n][i]) < 1.13) )
            || (abs(Component_Eta[n][i]) > 3.5)  );

        pt_const_cut = (Component_Pt[n][i] < constituent_pT_threshold(Component_Eta[n][i],B_Field));

        if (B_Field == 3.0)
          constant_p_cut = Component_P[n][i] < 1.0; 
        if (B_Field > 1.399 && B_Field < 1.4001) //no idea why == 1.4 fails
          constant_p_cut = Component_P[n][i] < 0.6;

        // Apply cuts to jet, breaking out of 
        if (eta_const_cut ) break; 
        if (pt_const_cut ) break; 
        if (constant_p_cut) break;
      }

      if (eta_const_cut || pt_const_cut || constant_p_cut) continue;//see break statement above

      if (ROOT::Math::VectorUtil::DeltaR(gLorentz,Lorentz) < max_DeltaR)
        dPhiRj_DeltaR->Fill(TMath::Abs(TVector2::Phi_mpi_pi(Lorentz.Phi() - electron_Phi - TMath::Pi())));

      //Loop to Subtract Neutrals from Truth Jets
      for (int i = 0; i < gNComponent[n]; i++){
        if (gComponent_Charge[n][i] == 0)
          n_neutral++;
        else
          n_ch++;

        ROOT::Math::PtEtaPhiEVector gConstLorentz(gComponent_Pt[n][i],gComponent_Eta[n][i],
            gComponent_Phi[n][i],gComponent_E[n][i]);
        if (gComponent_Charge[n][i] == 0)
          gLorentz -= gConstLorentz;//done before jet histograms are filled
      }

      maxed_neutrals =  n_neutral <  n_neutral_max;

      float dR = ROOT::Math::VectorUtil::DeltaR(Lorentz,gLorentz);
      float dP_P = (gLorentz.P() - Lorentz.P()) / gLorentz.P();
      Float_t Q_square = calc_Q_square(20,e_vector); //electron Beam of 20 GeV/c
      n_neutrals_vs_dP_P->Fill(dP_P,n_neutral);
      n_missed_vs_dP_P->Fill(dP_P,(gNComponent[n] - NComponent[n]) - n_neutral);
      n_missed->Fill(gNComponent[n] - NComponent[n]);
      n_missed_minuseutrals->Fill(gNComponent[n] - NComponent[n] - n_neutral);
      n_neutrals->Fill(n_neutral);
      n_charged->Fill(n_ch);
      n_constituents->Fill(gNComponent[n]);
      momentum_response->Fill(Lorentz.P(),gLorentz.P());

      float dphi_percent= dphi_res(B_Field,Lorentz.Phi(), Lorentz.Eta(), electron_gPhi,electron_gEta);
      whole_phi_percent->Fill(dphi_percent);
      if (Lorentz.Eta() < -0.5)
        backwards_phi_percent->Fill(dphi_percent);

      else if (Lorentz.Eta() > 0.5)
        forward_phi_percent->Fill(dphi_percent);
      else 
        central_phi_percent->Fill(dphi_percent);

      //--------------------Cut To Study--------------------//
      bool cut_to_study = true; //Check for conflicts with continue statements
      //cut_to_study = (abs(dP_P) < max_dP_P);
      //cut_to_study = (NComponent[n] >= min_comp);
      //cut_to_study =( (NComponent[n] > min_comp) && (E[n] > minE) );
      //cut_to_study = ((gNComponent[n] - NComponent[n] - n_neutral) < max_miss_const);
      //cut_to_study = (constant_p_cut);
      //cut_to_study = pt_const_cut;
      //cut_to_study = maxed_neutrals;

      if (cut_to_study)
      {
        reco_E->Fill(Lorentz.E());
        reco_P->Fill(Lorentz.P());
        reco_eta->Fill(Lorentz.Eta());
        reco_phi->Fill(Lorentz.Phi());
        reco_nconst->Fill(NComponent[n]);
        truth_E->Fill(gLorentz.E());
        truth_P->Fill(gLorentz.P());
        truth_eta->Fill(gLorentz.Eta());
        truth_phi->Fill(gLorentz.Phi());

        float reco_z = 0.;
        float reco_zT = 0.;
        if (reco_branches){ 
          for (int i = 0; i < NComponent[n]; i++){
            ROOT::Math::PtEtaPhiEVector ConstLorentz(Component_Pt[n][i],
                Component_Eta[n][i],Component_Phi[n][i],Component_P[n][i]);
            reco_z = ConstLorentz.P()/Lorentz.P();
            reco_FF->Fill(reco_z);
            reco_zT = ConstLorentz.Pt()/Lorentz.Pt();
            reco_FF_zT->Fill(reco_zT);
            reco_comp_eta->Fill(Component_Eta[n][i]);
            reco_comp_pt->Fill(Component_Pt[n][i]);
            /* if (gLorentz.DeltaR(Lorentz) < max_DeltaR) */
            if (ROOT::Math::VectorUtil::DeltaR(gLorentz,Lorentz) < max_DeltaR)
              reco_FF_DeltaR->Fill(reco_z);
          }
        }
        RjoTj->Fill(dP_P);
        nconst_diff->Fill(gNComponent[n] - NComponent[n]);
        Q2->Fill(Q_square);
        float truth_z = 0.;
        float truth_zT = 0.;
        for (int i = 0; i < gNComponent[n]; i++){
          if (gComponent_Charge[n][i] == 0) continue;
          comp_eta->Fill(gComponent_Eta[n][i]);
          comp_pid->Fill(gComponent_PID[n][i]);
          comp_pt->Fill(gComponent_Pt[n][i]);
          ROOT::Math::PtEtaPhiEVector gConstLorentz(gComponent_Pt[n][i],gComponent_Eta[n][i],gComponent_Phi[n][i],gComponent_E[n][i]);
          P_Component_vs_JetEta->Fill(gEta[n],gConstLorentz.P());
          Pt_Component_vs_JetEta->Fill(gEta[n],gComponent_Pt[n][i]);
          truth_z = gConstLorentz.P()/gLorentz.P();
          truth_FF->Fill(truth_z);
          truth_zT = gConstLorentz.Pt()/gLorentz.Pt();
          truth_FF_zT->Fill(truth_zT);
        }
        jet_cut_counter[1]+=1.; //increment
        for (int ieta = 0; ieta < nEta_bins; ieta++)
          if ((Eta[n] >= Eta_Bins[ieta]) && (Eta[n] < Eta_Bins[ieta+1]))
            momentum_in_eta_bins[ieta]->Fill(Pt[n]);

        Float_t True_DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(gLorentz.Phi() - electron_gPhi - TMath::Pi()));
        Float_t Reco_DeltaPhi = 0;
        if (reco_branches)
          Reco_DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(Lorentz.Phi() - electron_Phi - TMath::Pi()));

        dPhiTj_whole->Fill(True_DeltaPhi);
        dPhiTj_sigmaplus->Fill(True_DeltaPhi*(1+DPhi_res));
        dPhiTj_sigmamin->Fill(True_DeltaPhi*(1-DPhi_res));

        dPhiRj_whole->Fill(Reco_DeltaPhi);
        dPhiRj_sigmaplus->Fill(Reco_DeltaPhi*(1+DPhi_res));
        dPhiRj_sigmamin->Fill(Reco_DeltaPhi*(1-DPhi_res));
        /* std::cout<<"DPhi_res = "<<DPhi_res<<std::endl; */
        /* std::cout<<"Delta Thing = "<<Reco_DeltaPhi*(1-DPhi_res)<<std::endl; */
        /* if (gLorentz.DeltaR(Lorentz) < max_DeltaR) */

        for(int et = 0 ; et < size_eta_bin ; et++){
          //electron-jet correlation
          if (gEta[n] >= eta_bin[et] && gEta[n] < eta_bin[et+1])
            dPhiTj[et]->Fill(True_DeltaPhi);
          if (reco_branches)
            if (Eta[n] >= eta_bin[et] && Eta[n] < eta_bin[et+1]){
              dPhiRj[et]->Fill(Reco_DeltaPhi);
            }
        }
        
        dEtaTj->Fill(electron_gEta-gEta[n]);
        dEtaRj->Fill(electron_gEta-Eta[n]);

      }//cut_to_study

      else
      {
        reco_E_anticut->Fill(E[n]);
        reco_P_anticut->Fill(Lorentz.P());
        reco_eta_anticut->Fill(Eta[n]);
        reco_phi_anticut->Fill(Phi[n]);

        truth_E_anticut->Fill(gE[n]);
        truth_P_anticut->Fill(gLorentz.P());
        truth_eta_anticut->Fill(gEta[n]);
        truth_phi_anticut->Fill(gPhi[n]);

        reco_nconst_anticut->Fill(NComponent[n]);
        RjoTj_anticut->Fill(dP_P);
        nconst_diff_anticut->Fill(gNComponent[n] - NComponent[n]);
        Q2_anticut->Fill(Q_square);    
        for (int i = 0; i < gNComponent[n]; i++){
          if (gComponent_Charge[n][i] == 0) continue;
          comp_eta_anticut->Fill(gComponent_Eta[n][i]);
          comp_pid_anticut->Fill(gComponent_PID[n][i]);
          comp_pt_anticut->Fill(gComponent_Pt[n][i]);
          ROOT::Math::PtEtaPhiEVector gConstLorentz(gComponent_Pt[n][i],gComponent_Eta[n][i],gComponent_Phi[n][i],gComponent_E[n][i]);
          P_Component_vs_JetEta_anticut->Fill(gEta[n],gConstLorentz.P());
          Pt_Component_vs_JetEta_anticut->Fill(gEta[n],gComponent_Pt[n][i]);
        }
        for (int ieta = 0; ieta < nEta_bins; ieta++)
          if ((Eta[n] >= Eta_Bins[ieta]) && (Eta[n] < Eta_Bins[ieta+1]))
            anticut_momentum_in_eta_bins[ieta]->Fill(Pt[n]);

      }

      jet_cut_counter[0]+=1.; //Jets that only pass continue statements

      //Inclusive Spectra
      Rjve->Fill(electron_gE,E[n]);
      Tjve->Fill(electron_gE,gE[n]);

      //Kinematic Cuts
      // if (True_DeltaPhi < M_PI/2) continue;
      // if (gE[n] < 3.0) continue;
      // if (TMath::Abs(eta[n]) < 0.7) continue;
      //electron/Jet Comparisons
      float eE = electron_gE;
      eoTj->Fill(gE[n]/eE);
      eoRj->Fill(E[n]/eE);
      emTj->Fill(electron_gE-gE[n]);
      emRj->Fill(electron_gE-E[n]);
    }
  }

  jet_cut_counter[2] = jet_cut_counter[1]/jet_cut_counter[0];
  TVectorT<double> TJet_counter(3); //total jets, jets passed, %passed
  cout<<endl<<"Number of Jets passing Cut:"<<jet_cut_counter[1]<<endl;
  cout<<endl<<"Fraction of Jets passing Cut:"<<jet_cut_counter[2]<<endl;
  for (int i = 0; i < 3; i++) TJet_counter[i] = jet_cut_counter[i];
  //entry loop

  //Write to new root file
  TFile* fout = new TFile(Form("Histograms_Jet_Callibration_%fT.root",B_Field),"RECREATE");
  TJet_counter.Write("TJet_Counter");
  momentum_response->Write();
  justE->Write();
  //PID->Write();

  //reco_E->Scale(1./reco_E->GetEntries());
  //reco_P->Scale(1./reco_P->GetEntries());
  //reco_eta->Scale(1./reco_eta->GetEntries());
  //reco_nconst->Scale(1./reco_nconst->GetEntries());
  reco_E->Write();
  reco_P->Write();
  reco_eta->Write();
  reco_phi->Write();

  truth_E->Write();
  truth_P->Write();
  truth_eta->Write();
  truth_phi->Write();


  truth_electron_E->Write();
  smeared_electron_E->Write();
  electron_EoverP->Write();
  smeared_EoverP->Write();

  reco_nconst->Write();
  RjoTj->GetXaxis()->SetTitle("dp/p");
  RjoTj->Write();
  nconst_diff->Write();
  comp_eta->Write();
  comp_pid->Write();
  comp_pt->Scale(1./comp_pt->GetEntries());
  comp_pt->Write();
  reco_comp_eta->Write();
  reco_comp_pt->Write();
  P_Component_vs_JetEta->Write();
  Pt_Component_vs_JetEta->Write();
  n_neutrals_vs_dP_P->Write();
  n_missed_vs_dP_P->SetYTitle("No. Missed Jet Constituents");
  n_missed_vs_dP_P->SetXTitle("dp/p");
  n_missed_vs_dP_P->Write();
  n_neutrals->Write();
  n_charged->Write();
  n_missed->Write();
  n_missed_minuseutrals->Write();
  simple_ndiff->Write();
  n_constituents->Write();
  n_reco_constituents->Write();
  //reco_P_anticut->Scale(1./reco_P_anticut->GetEntries());
  //reco_E_anticut->Scale(1./reco_E_anticut->GetEntries());
  //reco_eta_anticut->Scale(1./reco_eta_anticut->GetEntries());
  //reco_nconst_anticut->Scale(1./reco_nconst_anticut->GetEntries());
  reco_E_anticut->Write();
  reco_P_anticut->Write();
  reco_eta_anticut->Write();
  reco_phi_anticut->Write();

  truth_E_anticut->Write();
  truth_P_anticut->Write();
  truth_eta_anticut->Write();
  truth_phi_anticut->Write();

  RjoTj_anticut->Write();
  reco_nconst_anticut->Write();
  nconst_diff_anticut->Write();
  comp_eta_anticut->Write();
  comp_pid_anticut->Write();
  comp_pt_anticut->Scale(1./comp_pt_anticut->GetEntries());
  comp_pt_anticut->Write();
  P_Component_vs_JetEta_anticut->Write();
  Pt_Component_vs_JetEta_anticut->Write();  


  Tjve->GetXaxis()->SetTitle("E^{True}_{electron} [GeV]");
  Rjve->GetXaxis()->SetTitle("E^{True}_{electron} [GeV]");
  Tjve->GetYaxis()->SetTitle("E^{True}_{Jet} [GeV]");
  Rjve->GetYaxis()->SetTitle("E^{Reco}_{Jet} [GeV]");
  Tjve->Write();
  Rjve->Write();

  eoTj->GetXaxis()->SetTitle("E_{e}/E_{jet}");
  eoRj->GetXaxis()->SetTitle("E_{e}/E_{jet}");
  eoTj->Write();
  eoRj->Write();

  emTj->GetXaxis()->SetTitle("E_{e}-E_{jet}");
  emRj->GetXaxis()->SetTitle("E_{e}-E_{jet}");
  emTj->Write();
  emRj->Write();

  for(int et = 0 ; et < size_eta_bin ; et++){
    dPhiTj[et]->GetXaxis()->SetTitle("#Delta#varphi");
    dPhiRj[et]->GetXaxis()->SetTitle("#Delta#varphi");
    dPhiTj[et]->Write();
    dPhiRj[et]->Write();
    aT_edPhi[et]->Write();
  }
  dPhiTj_whole->Write();
  dPhiRj_whole->Write();

  dPhiTj_sigmaplus->Write();
  dPhiTj_sigmamin->Write();
  dPhiRj_sigmaplus->Write();
  dPhiRj_sigmamin->Write();

  aT_edPhi_whole->Write();

  dPhiRj_DeltaR->Write();
  reco_FF_DeltaR->Write();

  dEtaTj->GetXaxis()->SetTitle("#Delta#eta");
  dEtaRj->GetXaxis()->SetTitle("#Delta#eta");
  dEtaTj->Write();
  dEtaRj->Write();

  Q2->Write();
  Q2_anticut->Write();
  truth_FF->Write();
  truth_FF_zT->Write();
  reco_FF->Write();
  reco_FF_zT->Write();
  all_truth_FF->Write();
  all_truth_FF_zT->Write();
  All_truth_E_histo->Write();
  // H_dR->Write();
  // H_NExtra_Matches->Write();

  central_phi_percent->Write();
  backwards_phi_percent->Write();
  forward_phi_percent->Write();
  whole_phi_percent->Write();
  for (int ieta = 0; ieta < nEta_bins; ieta++) {
    momentum_in_eta_bins[ieta]->Write();
    anticut_momentum_in_eta_bins[ieta]->Write();
  }  
  fout->Close();

}
