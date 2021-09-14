#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

// Root includes
#include <TROOT.h>
#include "TRint.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TVectorT.h"
#include <TVector2.h>
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

using namespace std;

float constituent_pT_threshold(float eta, float B_Field);
// ============================================================================================================================================
int main(int argc, char ** argv) {

    if(argc < 3){
      cout << "Run this code as:\n\033[1;32m./analysis_momentum_resolution [Cut] [B Field (T)] filename.root [B Field]\033[0m\n";
      cout << "[Cut] 0 -> Don't Apply any N_Missing Cut\n 1 -> Apply N_Missing < 1\n  = 2 -> apply N_Missing >= 1 Selection\n";
      cout << "[B Field] -> either 1.4 or 3.0 T at this time. DONT TRY ANOTHER B FIELD. Make sure file matches B_Field\n";
      cout << "[Root File] -> At this point in time, ../data/3T_combined.root or ../data/1p4T_combined.root\n";
      exit(0);
    }

  bool do_const_cut = true;
  bool N_Missing_Cut = false;
  int N_Missing_Max = 1;
  bool plot_mrad = true;
  string rad_mrad_string = "[rad]";
  if (plot_mrad)
    rad_mrad_string = "[mrad]";

  cout << "\033[1;31m**********\nUSEFUL INFO:\033[0m\nWill be loading data from file: '" << argv[3] << "' assumed to be in directory 'data'" << endl;

  if(atoi(argv[1])==0){do_const_cut = false; cout << " Will NOT apply any N Missing Const. Cut\n";}
  else if(atoi(argv[1])==1){N_Missing_Cut = true ;   cout << "Will apply N Missing Const. < 1 Cut\n" ;}
  else if(atoi(argv[1])==2){N_Missing_Cut = false;   cout << "Will apply N Missing Const. >= 1 Cut\n";}
  else{cout << "Something wrong with your election of input parameter 'Cut'. Bailing out!\n"; exit(0);}

  float B_Field = atof(argv[2]);
  printf("B Field = %f\n",B_Field);
  if ((B_Field != float(3.0)) && (B_Field != float(1.4))) {printf("B Field = %f not valid. Only 1.4 or 3.0 Tesla!\n",B_Field); exit(0);}

  // -------------------------
  // Binning
  float eta_bin[] = {-3.0,-1.5,-.5,0.5,1.5,3.0};
  float mom_bin[] = {4.,6.,8.,10.,12.,15, 20.};

  const int size_eta_bin = sizeof(eta_bin)/sizeof(*eta_bin);
  const int size_mom_bin = sizeof(mom_bin)/sizeof(*mom_bin);

  TVectorT<double> TVT_eta_bin(size_eta_bin);	for(int i = 0 ; i < size_eta_bin ; i++) TVT_eta_bin[i] = eta_bin[i];
  TVectorT<double> TVT_mom_bin(size_mom_bin);	for(int i = 0 ; i < size_mom_bin ; i++) TVT_mom_bin[i] = mom_bin[i];
  // -------------------------
  
  // Filename Strings 
  TString infile = argv[3];

  cout<<endl<<"Input Root File = "<<argv[3]<<endl;

  TString outfile; 

  if (do_const_cut && N_Missing_Cut)
    outfile = "histograms_reco_No_Missing_const_output_mom_res_" + (TString)Form("sigma_eta_%i_p_%i_B_%1.1f",size_eta_bin-1,size_mom_bin-1,B_Field) + ".root";

  if (do_const_cut && !N_Missing_Cut)
    outfile = "histograms_reco_Missing_const_output_mom_res_" + (TString)Form("sigma_eta_%i_p_%i_B_%1.1f",size_eta_bin-1,size_mom_bin-1,B_Field) + ".root";

  if (!do_const_cut)
    outfile = "histograms_reco_NoCuts_output_mom_res_" + (TString)Form("sigma_eta_%i_p_%i_B_%1.1f",size_eta_bin-1,size_mom_bin-1,B_Field) + ".root";

  cout<<endl<<"Input Root File = "<<infile<<endl;
  cout<<endl<<"Output Root File = "<<outfile<<endl;
  TString out_pdf = "output_fits_mom_res_" + (TString)Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".pdf";
  TString out_pdf2 = "results/results_mom_res_" + (TString)Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + ".pdf";
  // -------------------------------------------------------------
  // Some settings
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle -> SetOptStat(0);	
  // -------------------------------------------------------------

  // Loading Data from ROOT file
  TFile * F = new TFile(infile);
  TTree *T = dynamic_cast<TTree *>(F->Get("T"));
  if (T == NULL) { std::cout << " Tree Fail " << std::endl; exit(EXIT_FAILURE); }
  Int_t njets;
  int nEntries = T -> GetEntries();
  const int MaxNumJets = 20;
  const int kMaxConstituents = 100;

  array<Float_t, MaxNumJets> E,Eta,Phi,Pt,gE,gEta,gPhi,gPt;
  array<Int_t, MaxNumJets> NComponent, gNComponent;

  Float_t electron_gE,electron_gEta,electron_gPhi,electron_gPt;
  Float_t electron_E,electron_Eta,electron_Phi,electron_Pt;

  array<array<Float_t, kMaxConstituents >, MaxNumJets > gComponent_Eta,gComponent_PID,
    gComponent_Pt,gComponent_Phi,gComponent_E, gComponent_Charge;

  array<array<Float_t, kMaxConstituents >, MaxNumJets > Component_Eta,Component_Phi,Component_P,Component_Pt;

  T -> SetBranchAddress("njets",&njets);
  T -> SetBranchAddress("e",&E);
  T -> SetBranchAddress("eta",&Eta);
  T -> SetBranchAddress("phi",&Phi);
  T -> SetBranchAddress("pt",&Pt);
  T -> SetBranchAddress("nComponent",&NComponent);

  T -> SetBranchAddress("Constituent_recoEta", Component_Eta.data());
  T -> SetBranchAddress("Constituent_recoPt",Component_Pt.data());
  T -> SetBranchAddress("Constituent_recoP",Component_P.data());
  T -> SetBranchAddress("Constituent_recoPhi",Component_Phi.data());

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

  T -> SetBranchAddress("electron_recoE",&electron_E);
  T -> SetBranchAddress("electron_recoEta",&electron_Eta);
  T -> SetBranchAddress("electron_recoPhi",&electron_Phi);
  T -> SetBranchAddress("electron_recoPt",&electron_Pt);

  T -> SetBranchAddress("electron_truthE",&electron_gE);
  T -> SetBranchAddress("electron_truthEta",&electron_gEta);
  T -> SetBranchAddress("electron_truthPhi",&electron_gPhi);
  T -> SetBranchAddress("electron_truthPt",&electron_gPt);

  // -------------------------------------------------------------

  // -------------------------------------------------------------
  // Defining histograms
  TH1F *** h1_dpp_p_et_bins = new TH1F**[size_eta_bin-1];	// delta p / p vs. p in eta bins
  TH1F *** h1_dth_p_et_bins = new TH1F**[size_eta_bin-1];	// delta theta vs. p in eta bins
  TH1F *** h1_dph_p_et_bins = new TH1F**[size_eta_bin-1];	// delta phi   vs. p in eta bins

  TH1F *** h1_eDelta_dph_p_et_bins = new TH1F**[size_eta_bin-1];	// electron-jet delta phi vs p in eta ibns
  TH1F *** h1_efulljet_Delta_dph_p_et_bins = new TH1F**[size_eta_bin-1];	// electron-jet delta phi vs p in eta ibns

  int n_dpp_bins = 30;
  int n_dth_bins = 30;
  int n_dph_bins = 30;

  for(int et = 0 ; et < size_eta_bin-1 ; et++){
    h1_dpp_p_et_bins[et] = new TH1F*[size_mom_bin-1];
    h1_dth_p_et_bins[et] = new TH1F*[size_mom_bin-1];
    h1_dph_p_et_bins[et] = new TH1F*[size_mom_bin-1];

    h1_eDelta_dph_p_et_bins[et] = new TH1F*[size_mom_bin-1];
    h1_efulljet_Delta_dph_p_et_bins[et] = new TH1F*[size_mom_bin-1];

    for(int p = 0 ; p < size_mom_bin-1 ; p++){

      if (!do_const_cut){
        h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),
            ";dP/P;Counts",300,-0.4,0.8);
        h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),
            ";d#theta [rad];Counts",80,-0.02,0.02);
        h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),
            ";d#phi [rad];Counts",80,-0.02,0.02);
        h1_eDelta_dph_p_et_bins[et][p] = new TH1F(Form("h1_eDelta_dph_p_et_bins_%i_%i",et,p),
            ";e-Jet d#Delta#phi [rad];Counts",45,-0.04,0.04);
        h1_efulljet_Delta_dph_p_et_bins[et][p] = new TH1F(Form("h1_efulljet_Delta_dph_p_et_bins_%i_%i",et,p),
            ";e-Jet_{full} d#Delta#phi [rad];Counts",45,-0.04,0.04);
      }

      else if (N_Missing_Cut){
        h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),
            ";dP/P;Counts",n_dpp_bins,-0.06,0.06);
        h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),
            ";d#theta [rad];Counts",20,-0.005,0.005);
        h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),
            ";d#phi [rad];Counts"  ,40,-0.01,0.01);
        h1_eDelta_dph_p_et_bins[et][p] = new TH1F(Form("h1_eDelta_dph_p_et_bins_%i_%i",et,p),
            ";e-Jet d#Delta#phi [rad];Counts"  ,n_dph_bins,-0.01,0.01);
        h1_efulljet_Delta_dph_p_et_bins[et][p] = new TH1F(Form("h1_efulljet_Delta_dph_p_et_bins_%i_%i",et,p),
            ";e-Jet_{full} d#Delta#phi [rad];Counts"  ,n_dph_bins,-0.01,0.01);
      }

      else if (!N_Missing_Cut) {
        h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),
            ";dP/P;Counts",300,-0.4,0.8);
        h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),
            ";d#theta [rad];Counts",80,-0.02,0.02);
        h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),
            ";d#phi [rad];Counts",80,-0.02,0.02);
        h1_eDelta_dph_p_et_bins[et][p] = new TH1F(Form("h1_eDelta_dph_p_et_bins_%i_%i",et,p),
            ";e-Jet d#Delta#phi [rad];Counts",80,-0.2,0.2);
        h1_efulljet_Delta_dph_p_et_bins[et][p] = new TH1F(Form("h1_efulljet_Delta_dph_p_et_bins_%i_%i",et,p),
            ";e-Jet_{full} d#Delta#phi [rad];Counts",80,-0.2,0.2);
      }

      h1_dpp_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < P < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
      h1_dth_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < P < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
      h1_dph_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < P < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
    } 
  }	
  // -------------------------------------------------------------	

  //Histograms to verify cuts:
  /* truth and constituent pt, p, and eta */ 
  /*   number of missing constituets distribution */
  /*   jet pT,p, and eta distributions */
  TH1F *reco_jet_p = new TH1F("reco_jet_p","Reco Jet Momentum Distribution",200,0,50);
  TH1F *reco_jet_eta = new TH1F("reco_jet_eta","Reco Jet #eta",50,-5,5);
  TH1F *reco_comp_pt = new TH1F("reco_comp_pt","Jet reco Component p_{T}",200,0,20);
  TH1F *reco_comp_eta = new TH1F("reco_comp_eta","Jet reco Component #eta",100,-5,5);
  TH1F *reco_NConstituents = new TH1F("reco_nconstituents","Number of Reconstructed Jet Constituents",20,0,20);

  TH1F *truth_jet_p = new TH1F("truth_jet_p","Truth Jet Momentum Distribution",200,0,50);
  TH1F *truth_jet_eta = new TH1F("truth_jet_eta","Truth Jet #eta",50,-5,5);
  TH1F *truth_comp_pt = new TH1F("truth_comp_pt","Jet Truth Component p_{T}",200,0,20);
  TH1F *truth_comp_eta = new TH1F("truth_comp_eta","Jet Truth Component #eta",100,-5,5);
  TH1F *truth_NConstituents = new TH1F("truth_nconstituents","Number of truthnstructed Jet Constituents",20,0,20);

  // Momentum Response
  TH2F * mom_response = new TH2F("mom_response","P_{Charge}^{Reco} vs. P_{Charge}^{Truth}",50,0,50,50,0,50);
  cout << "\033[1;31m********************************************************************\033[0m\n";
  // -------------------------------------------------------------
  // Loop over entries of the tree	
  for (Long64_t ev = 0; ev < nEntries; ev++){
    T->GetEntry(ev);
    if (ev%10000==0) fprintf(stderr,"\r%d: Entry %lli out of %d",__LINE__,ev,nEntries);
    //if (ev==500000) break;
    for (int n = 0; n < njets; ++n) {

      if (NComponent[n] < 4) continue;
      if (isnan(gE[n])) continue;
      if (E[n] < 4.0) continue;
      ROOT::Math::PtEtaPhiEVector Lorentz(Pt[n],Eta[n],Phi[n],E[n]);
      ROOT::Math::PtEtaPhiEVector gLorentz(gPt[n],gEta[n],gPhi[n],gE[n]);
      ROOT::Math::PtEtaPhiEVector gFullLorentz(gPt[n],gEta[n],gPhi[n],gE[n]);

      bool eta_const_cut = true; //avoid crack where barrel meets disks
      bool pt_const_cut = true;//cut away helixes
      int n_neutral = 0;
      int n_ch = 0;
      int n_neutral_max = 1;
      float max_DeltaR = 0.1;

      //Truth-Level Cuts (Don't use)
      /* for (int i = 0; i < gNComponent[n]; i++) */
      /* eta_const_cut = (  ( (abs(gComponent_Eta[n][i]) > 1.04) && (abs(gComponent_Eta[n][i]) < 1.15) ) */
      /*     || (abs(gComponent_Eta[n][i]) > 3.5)  );//Aluminum Cone between 1.06 and 1.13 eta */
      /* pt_const_cut = (gComponent_Pt[n][i] < constituent_pT_threshold(gComponent_Eta[n][i],B_Field)); */


      //Using Reco-Constituents. Only Cuts that can be used in time of experiment can be used to study performance.
      for (int i = 0; i < NComponent[n]; i++){
        eta_const_cut = (  ( (abs(Component_Eta[n][i]) > 1.04) && (abs(Component_Eta[n][i]) < 1.15) )
            || (abs(Component_Eta[n][i]) > 3.5)  );
        pt_const_cut = (Component_Pt[n][i] < constituent_pT_threshold(Component_Eta[n][i],B_Field));

        if (eta_const_cut || pt_const_cut) break; //skip jets that fail. Break Constituent loop, continue jet loop
      }

      if (eta_const_cut || pt_const_cut) continue;

      for (int t = 0; t < gNComponent[n]; t++){
        if (gComponent_Charge[n][t] == 0)
          n_neutral++;
        else
          n_ch++;

        ROOT::Math::PtEtaPhiEVector gConstLorentz(gComponent_Pt[n][t],
            gComponent_Eta[n][t],
            gComponent_Phi[n][t],
            gComponent_E[n][t]);
        if (gComponent_Charge[n][t] == 0)
          gLorentz -= gConstLorentz;
      }

      //Delta R Check after Neutral Subtraction (Charged Comparison)
      /* float dR = ROOT::Math::VectorUtil::DeltaR(Lorentz,gLorentz); */
      /* if (dR > max_DeltaR) continue; */

      int N_Missing = gNComponent[n] - NComponent[n] - n_neutral;
      if (do_const_cut)
      {
        if ( (N_Missing_Cut) && (N_Missing >= N_Missing_Max) ) {continue;}
        if ( (!N_Missing_Cut) && (N_Missing < N_Missing_Max) ) {continue;} //anticut of above
      }

      //After last continue statement, fill histos that check cuts
      reco_jet_p->Fill(Lorentz.P());
      reco_jet_eta->Fill(Lorentz.Eta());
      reco_NConstituents->Fill(NComponent[n]);
      for (int r = 0; r < NComponent[n]; r++){
        reco_comp_pt->Fill(Component_Pt[n][r]);
        reco_comp_eta->Fill(Component_Eta[n][r]);
      }

      truth_jet_p->Fill(gLorentz.P());
      truth_jet_eta->Fill(gLorentz.Eta());
      truth_NConstituents->Fill(gNComponent[n]);
      for (int t = 0; t < gNComponent[n]; t++){
        truth_comp_pt->Fill(gComponent_Pt[n][t]);
        truth_comp_eta->Fill(gComponent_Eta[n][t]);
      }

      float geta = gLorentz.Eta();
      float P_reco = Lorentz.P();
      float P_truth = gLorentz.P();

      float eTrue_DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(gLorentz.Phi() - electron_gPhi - TMath::Pi()));
      float eReco_DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(Lorentz.Phi() - electron_Phi - TMath::Pi()));
      float eJet_dDeltaPhi = (eTrue_DeltaPhi-eReco_DeltaPhi);

      float eFull_DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(gFullLorentz.Phi() - electron_gPhi - TMath::Pi()));
      float efullJet_dDeltaPhi = (eTrue_DeltaPhi-eFull_DeltaPhi);

      float E_gE = (Lorentz.E()/gLorentz.E());
      //float E_gE = (E[n]/gE[n]);
      mom_response->Fill(gLorentz.P(),Lorentz.P());

      //Most Important Variables
      float dth = Lorentz.Theta() - gLorentz.Theta();
      float dph = Lorentz.Phi() - gLorentz.Phi();
      float dP_P = (P_truth-P_reco)/P_truth;

      //Filling 
      for(int et = 0 ; et < size_eta_bin-1 ; et++){
        if( geta >  eta_bin[et] &&  geta <= eta_bin[et+1] ){
          for(int p = 0 ; p < size_mom_bin-1 ; p++){
            if( P_truth > mom_bin[p] && P_truth <= mom_bin[p+1] ){
              h1_dpp_p_et_bins[et][p] -> Fill( dP_P );
              h1_dth_p_et_bins[et][p] -> Fill( dth  );
              h1_dph_p_et_bins[et][p] -> Fill( dph  );
              h1_eDelta_dph_p_et_bins[et][p] -> Fill(eJet_dDeltaPhi );
              h1_efulljet_Delta_dph_p_et_bins[et][p] -> Fill(efullJet_dDeltaPhi);
              //fprintf(stderr,"%d: Full-Jet Delta Phi = %1.4f, Truth-Reco Delta Phi = %1.4f\n",__LINE__,eFull_DeltaPhi,eJet_dDeltaPhi);
            }	
          }
        }
      }
    }//do_const
  }//jets
  cout << "\033[1;31m********************************************************************\033[0m\n";
  // -------------------------------------------------------------
  // Declaring other useful variables and functions

  // -------------------------------------------------------------
  // -------------------------------------------------------------

  // -------------------------------------------------------------
  // -------------------------------------------------------------
  // Plotting histograms

  // -------------------------------------------------------------
  // Saving histograms
  TFile * Fout = new TFile(outfile,"recreate");
  mom_response->Write();

  reco_jet_p->Write();
  reco_jet_eta->Write();
  reco_NConstituents->Write();
  reco_comp_pt->Write();
  reco_comp_eta->Write();

  truth_jet_p->Write();
  truth_jet_eta->Write();
  truth_NConstituents->Write();
  truth_comp_pt->Write();
  truth_comp_eta->Write();

  TVT_eta_bin.Write("TVT_eta_bin");
  TVT_mom_bin.Write("TVT_mom_bin");

  TDirectory *dpp_dir =Fout->mkdir("dpp_histos");
  dpp_dir->cd();
  for(int p = 0 ; p < size_mom_bin-1 ; p++){
    for(int et = 0 ; et < size_eta_bin-1 ; et++){ 
      h1_dpp_p_et_bins[et][p]->Write(Form("h1_dpp_et_%i_p_%i_bin",et,p));
    }
  } 

  TDirectory *dph_dir =Fout->mkdir("dph_histos");
  dph_dir->cd();
  for(int p = 0 ; p < size_mom_bin-1 ; p++){
    for(int et = 0 ; et < size_eta_bin-1 ; et++){ 
      h1_dph_p_et_bins[et][p]->Write(Form("h1_dph_et_%i_p_%i_bin",et,p));
    }
  } 

  TDirectory *Ddph_dir =Fout->mkdir("eD_dph_histos");
  Ddph_dir->cd();
  for(int p = 0 ; p < size_mom_bin-1 ; p++){
    for(int et = 0 ; et < size_eta_bin-1 ; et++){ 
      h1_eDelta_dph_p_et_bins[et][p]->Write(Form("h1_eDelta_dph_et_%i_p_%i_bin",et,p));
    }
  }
  for(int p = 0 ; p < size_mom_bin-1 ; p++){
    for(int et = 0 ; et < size_eta_bin-1 ; et++){ 
      h1_efulljet_Delta_dph_p_et_bins[et][p]->Write(Form("h1_efulljet_Delta_dph_et_%i_p_%i_bin",et,p));
    }
  } 

  TDirectory *dth_dir =Fout->mkdir("dth_histos");
  dth_dir->cd();
  for(int p = 0 ; p < size_mom_bin-1 ; p++){
    for(int et = 0 ; et < size_eta_bin-1 ; et++){ 
      h1_dth_p_et_bins[et][p]->Write(Form("h1_dth_et_%i_p_%i_bin",et,p));
    }
  } 
  Fout -> Close();
  return 0;
}//end main
// ============================================================================================================================================

//-----------------------------------------------

float constituent_pT_threshold(float eta, float B_Field)
{
  // Minimum pT B = 1.5 T or 3.0 T (https://physdiv.jlab.org/DetectorMatrix/):

  eta = abs(eta);
  float eta_bins[6] = {3.5,2.5,2.0,1.5,1.0,0.0};
  float pT_threshold_array[5] = {0.1,0.13,0.07,0.15,0.2}; 

  if (B_Field == 3.0)
  {
    pT_threshold_array[0]=0.15; pT_threshold_array[1]=0.22;
    pT_threshold_array[2] = 0.16; pT_threshold_array[3]=0.3;
    pT_threshold_array[4]=0.4;
  }
  float pT_threshold = 0;
  for (int i = 0; i < 5; i++)
    if ( (eta < eta_bins[i])&&(eta > eta_bins[i+1]) )
      pT_threshold = pT_threshold_array[i];

  return pT_threshold;
}


