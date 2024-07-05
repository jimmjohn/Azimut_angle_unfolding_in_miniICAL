#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include "TLegend.h"

#include "TVector.h"
#include <vector>
#include <TF1.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TStyle.h>
#include "TPostScript.h"
#include <TPad.h>
#include <TLine.h>
#include <TRandom.h>

#include "TUnfold.h"
#include "TUnfoldDensity.h"


using namespace std;

//bool muMinus = false;
bool mlApplied = false;

bool var_p = true;

//int version = 107;//100 means all models

double chisquare_p = 0;
double chisquare_phi = 0;

//theta
int nthetabins_gen = 40;
int nthetabins_reco = 80;
double theta_low = 0.;
double theta_high = 90.;

//phi
//int nphibins_gen[5] = {50,38,34,34,34};
int nphibins_gen[5] = {60,60,60,60,60};
int nphibins_reco = 70;
double phi_low = -1.*TMath::Pi();
double phi_high = TMath::Pi();

//double phi_low_gen[5] =  {-50.*TMath::Pi()/30., -38.*TMath::Pi()/30., -34.*TMath::Pi()/30., -34.*TMath::Pi()/30., -34.*TMath::Pi()/30.};
//double phi_high_gen[5] = { 50.*TMath::Pi()/30.,  38.*TMath::Pi()/30.,  34.*TMath::Pi()/30.,  34.*TMath::Pi()/30.,  34.*TMath::Pi()/30.};
double phi_low_gen[5] =  {-60.*TMath::Pi()/30., -60.*TMath::Pi()/30., -60.*TMath::Pi()/30., -60.*TMath::Pi()/30., -60.*TMath::Pi()/30.};
double phi_high_gen[5] = { 60.*TMath::Pi()/30.,  60.*TMath::Pi()/30.,  60.*TMath::Pi()/30.,  60.*TMath::Pi()/30.,  60.*TMath::Pi()/30.};

//Momentum
int npbins_gen_init=8;
int npbins_reco = 16;
double p_low = 0.8;
double p_high = 3.0;


double binwidth_gen = (p_high-p_low)/npbins_gen_init;
double binwidth_reco = (p_high-p_low)/npbins_reco;

int extra_bins = (p_low - 0.5)/binwidth_gen;
int npbins_gen = npbins_gen_init;//+2*extra_bins;

void DivideHistogramByBinWidth(TH1D *histogram) {
  for (int i = 1; i <= histogram->GetNbinsX(); ++i) {
    double binContent = histogram->GetBinContent(i);
    double binWidth = histogram->GetBinWidth(i);
    histogram->SetBinContent(i, binContent / binWidth);
    histogram->SetBinError(i, histogram->GetBinError(i) / binWidth);
  }
}

int main(int argc, char** argv){
//void unfolding(){

  bool muMinus = (atoi(argv[3])==1) ? true : false;

  
  TFile* infileRM;//RSA  
  if(!muMinus){
    infileRM = TFile::Open("response_matrix/fileOut/phi_response_matrix_v100_muplus_ndof5.root","READ");
  }
  else{
    infileRM = TFile::Open("response_matrix/fileOut/phi_response_matrix_v100_muminus_ndof5.root","READ");
  }

  int version = atoi(argv[1]);
  string saveDir = "fileOut";
  
//  string filename = "Mc_v107_2018_Data_Magnetic.root";
string filename = argv[2];
char name[100];

double p_low_temp = p_low;

if(muMinus) {p_low = -1*p_high; p_high = -1*p_low_temp;}

  double pival = acos(-1.);
const int numThetaRanges = 5;
double thetaRanges[numThetaRanges][2] = {{0, 17}, {17, 26}, {26, 34}, {34, 44}, {44,90}};
//double thetaRanges[numThetaRanges][2] = {{0, 90}};

  //double reco_bin_sch[23] = {0,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1};
  //double gen_bin_sch[13] = {0,1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3,3.1};

  //double reco_bin_sch[35] = {0, 1, 1.06149, 1.12298, 1.18446, 1.24595, 1.30744, 1.36893, 1.43041, 1.4919, 1.55339, 1.61488, 1.67636, 1.73785, 1.79934, 1.86083, 1.92231, 1.9838, 2.04529, 2.10678, 2.16826, 2.22975, 2.29124, 2.35273, 2.41421, 2.4757, 2.53719, 2.59868, 2.66016, 2.72165, 2.78314, 2.84463, 2.90611, 2.9676, 3.1};
  //double gen_bin_sch[17] = {0, 1, 1.14142, 1.28284, 1.42426, 1.56569, 1.70711, 1.84853, 1.98995, 2.13137, 2.27279, 2.41421, 2.55563, 2.69706, 2.83848, 2.9799, 3.1};

  //double reco_bin_sch[20] = {0, 1, 1.11314, 1.22627, 1.33941, 1.45255, 1.56569, 1.67882, 1.79196, 1.9051, 2.01823, 2.13137, 2.24451, 2.35765, 2.47078, 2.58392, 2.69706, 2.81019, 2.92333, 3.1};
  //double gen_bin_sch[10] = {0, 1, 1.28284, 1.56569, 1.84853, 2.13137, 2.41421, 2.69706, 2.9799, 3.1};

  //Variable bin width based on resolution MuPlus
  //double reco_bin_sch_muplus[18] = {0, 1, 1.11314, 1.22819, 1.3452, 1.4642, 1.58521, 1.70828, 1.83344, 1.96072, 2.09016, 2.2218, 2.35567, 2.49181, 2.63027, 2.77107, 2.91427, 3.1};
  //double gen_bin_sch_muplus[9] = {0, 1, 1.28284, 1.57769, 1.88504, 2.20543, 2.53941, 2.88757, 3.1};
  // double reco_bin_sch_muplus[47] = {0.8, 0.82, 0.85, 0.87, 0.9, 0.92, 0.95, 0.98, 1.01, 1.04, 1.07, 1.1, 1.13, 1.16, 1.2, 1.23, 1.27, 1.3, 1.34, 1.38, 1.42, 1.46, 1.51, 1.55, 1.59, 1.64, 1.69, 1.74, 1.79, 1.84, 1.89, 1.95, 2.01, 2.06, 2.13, 2.19, 2.25, 2.32, 2.38, 2.45, 2.52, 2.6, 2.67, 2.75, 2.83, 2.92, 3};
  // double gen_bin_sch_muplus[23] = {0.8, 0.85, 0.9, 0.96, 1.02, 1.08, 1.15, 1.22, 1.29, 1.37, 1.46, 1.55, 1.65, 1.75, 1.86, 1.97, 2.09, 2.22, 2.36, 2.51, 2.66, 2.83, 3};
  double reco_bin_sch_muplus[17] = {0.8, 0.87, 0.94, 1.02, 1.11, 1.21, 1.31, 1.43, 1.55, 1.68, 1.83, 1.98, 2.16, 2.34, 2.54, 2.76, 3};
  double gen_bin_sch_muplus[9] = {0.8, 0.94, 1.11, 1.31, 1.55, 1.83, 2.16, 2.54, 3};

  //Start from 0.6 - 3.5
  //0.7-0.85, 0.85-1.0

  //suppose 0.8-3.0 its 20 bins total bins =22, number of points =23
  //suppose 0.8-3.0 its 8 bins total bins =10 , number of points =11


  //Variable bin width based on resolution MuMinus
  //double reco_bin_sch_muminus[47] = {-3, -2.92, -2.83, -2.75, -2.67, -2.6, -2.52, -2.45, -2.38, -2.32, -2.25, -2.19, -2.13, -2.06, -2.01, -1.95, -1.89, -1.84, -1.79, -1.74, -1.69, -1.64, -1.59, -1.55, -1.51, -1.46, -1.42, -1.38, -1.34, -1.3, -1.27, -1.23, -1.2, -1.16, -1.13, -1.1, -1.07, -1.04, -1.01, -0.98, -0.95, -0.92, -0.9, -0.87, -0.85, -0.82, -0.8};
  //double gen_bin_sch_muminus[23] = {-3, -2.83, -2.66, -2.51, -2.36, -2.22, -2.09, -1.97, -1.86, -1.75, -1.65, -1.55, -1.46, -1.37, -1.29, -1.22, -1.15, -1.08, -1.02, -0.96, -0.9, -0.85, -0.8,};
  double reco_bin_sch_muminus[17] = {-3, -2.76, -2.54, -2.34, -2.16, -1.98, -1.83, -1.68, -1.55, -1.43, -1.31, -1.21, -1.11, -1.02, -0.94, -0.87, -0.8};
  double gen_bin_sch_muminus[9] = {-3, -2.54, -2.16, -1.83, -1.55, -1.31, -1.11, -0.94, -0.8};

  double reco_bin_sch[17];
  double gen_bin_sch[9];

  if(!muMinus) {
    std::copy(std::begin(reco_bin_sch_muplus), std::end(reco_bin_sch_muplus), std::begin(reco_bin_sch));
    std::copy(std::begin(gen_bin_sch_muplus), std::end(gen_bin_sch_muplus), std::begin(gen_bin_sch));
  } else {
    std::copy(std::begin(reco_bin_sch_muminus), std::end(reco_bin_sch_muminus), std::begin(reco_bin_sch));
    std::copy(std::begin(gen_bin_sch_muminus), std::end(gen_bin_sch_muminus), std::begin(gen_bin_sch));
  }

  TH2D* covariance_matrix_mc[numThetaRanges];
  TH2D* covariance_matrix_data[numThetaRanges];

  
  //Phi
  TH1D *hist_phi_reco[numThetaRanges];
  for (int j = 0; j < numThetaRanges; ++j) {
    if(!muMinus) {sprintf(name,"h_phi_mc_reco_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_phi_reco[j] = new TH1D(name,name,nphibins_reco, phi_low, phi_high); hist_phi_reco[j]->Sumw2();}
    if(muMinus)  {sprintf(name,"h_phi_mc_reco_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_phi_reco[j] = new TH1D(name,name,nphibins_reco, phi_low, phi_high); hist_phi_reco[j]->Sumw2();}
  }
  TH1D *hist_phi_data[numThetaRanges];
  for (int j = 0; j < numThetaRanges; ++j) {
    if(!muMinus) {sprintf(name,"h_phi_data_reco_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_phi_data[j] = new TH1D(name,name,nphibins_reco, phi_low, phi_high); hist_phi_data[j]->Sumw2();}
    if(muMinus)  {sprintf(name,"h_phi_data_reco_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_phi_data[j] = new TH1D(name,name,nphibins_reco, phi_low, phi_high); hist_phi_data[j]->Sumw2();}
  }
  TH1D *hist_phi_gen[numThetaRanges];
  for (int j = 0; j < numThetaRanges; ++j) {
    if(!muMinus) {sprintf(name,"h_phi_gen_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_phi_gen[j] = new TH1D(name,name,nphibins_gen[j], phi_low_gen[j], phi_high_gen[j]); hist_phi_gen[j]->Sumw2();}
    if(muMinus)  {sprintf(name,"h_phi_gen_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_phi_gen[j] = new TH1D(name,name,nphibins_gen[j], phi_low_gen[j], phi_high_gen[j]); hist_phi_gen[j]->Sumw2();}
  }
  if(!muMinus) {sprintf(name,"h_phi_gen_reco_muplus_v%d",version);}
  if(muMinus)  {sprintf(name,"h_phi_gen_reco_muminus_v%d",version);}
  TH2D *mat_phi_rm[numThetaRanges];
  for (int j = 0; j < numThetaRanges; ++j) {
    if(!muMinus) {sprintf(name,"h_phi_gen_reco_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      mat_phi_rm[j] = new TH2D(name,name,nphibins_reco, phi_low, phi_high, nphibins_gen[j], phi_low_gen[j], phi_high_gen[j]); mat_phi_rm[j]->Sumw2();}
    if(muMinus)  {sprintf(name,"h_phi_gen_reco_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      mat_phi_rm[j] = new TH2D(name,name,nphibins_reco, phi_low, phi_high, nphibins_gen[j], phi_low_gen[j], phi_high_gen[j]); mat_phi_rm[j]->Sumw2();}
  }



  //Get the efficiecny and fakes from all models
  TH1D *fake_rate_phi[numThetaRanges];
  TH1D *efficiency_phi[numThetaRanges];

  for (int j = 0; j < numThetaRanges; ++j) {
    if(!muMinus) {
      sprintf(name,"fake_rate_phi_muplus_v100_theta_%d_%d", (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      fake_rate_phi[j] = (TH1D*)infileRM->Get(name);
      sprintf(name,"efficiency_phi_muplus_v100_theta_%d_%d", (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      efficiency_phi[j] = (TH1D*)infileRM->Get(name);
    }
    if(muMinus)  {
      sprintf(name,"fake_rate_phi_muminus_v100_theta_%d_%d", (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      fake_rate_phi[j] = (TH1D*)infileRM->Get(name);
      sprintf(name,"efficiency_phi_muminus_v100_theta_%d_%d", (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      efficiency_phi[j] = (TH1D*)infileRM->Get(name);
    }
  }



  int nbins_gen, nbins_reco;
  double lowedge, upedge;

  double thgen, phgen, pgen;
  double threco, phreco, preco;
  double chi2reco;
  int ndfreco;

  //In tree
  Float_t momin[1], thein[1], phiin[1];
  Float_t momrf[1], therf[1], phirf[1];
  float chisquare[1];
  UInt_t ndof[1];
  UInt_t fitfailed;
  Float_t learnedVal, learnedErr;

  Float_t momrfd[1], therfd[1], phirfd[1];
  float chisquared[1];
  UInt_t ndofd[1];
  UInt_t fitfailedd;
  Float_t learnedVald, learnedErrd;

//sprintf(name,"%s",filename.c_str());
sprintf(name,"/var/nfscondor/rajshah/paper2/mom_unfolding/%s",filename.c_str());
 cout<<name<<endl;
 TFile *file1 = new TFile(name,"read");
 
 file1->cd();
  TTree *T1;
  T1 = (TTree*)file1->Get("T3");
  T1->SetBranchAddress("thein",thein);
  T1->SetBranchAddress("phiin",phiin);
  T1->SetBranchAddress("momin",momin);
  T1->SetBranchAddress("therf",therf);
  T1->SetBranchAddress("phirf",phirf);
  T1->SetBranchAddress("momrf",momrf);
  T1->SetBranchAddress("chisquare",chisquare);
  T1->SetBranchAddress("ndof",ndof);
  T1->SetBranchAddress("fitfailed",&fitfailed);
  T1->SetBranchAddress("learnedVal",&learnedVal);
  T1->SetBranchAddress("learnedErr",&learnedErr);

  TTree *T2;
  T2 = (TTree*)file1->Get("T5");
  T2->SetBranchAddress("therf",therfd);
  T2->SetBranchAddress("phirf",phirfd);
  T2->SetBranchAddress("momrf",momrfd);
  T2->SetBranchAddress("chisquare",chisquared);
  T2->SetBranchAddress("ndof",ndofd);
  T2->SetBranchAddress("fitfailed",&fitfailedd);
  T2->SetBranchAddress("learnedVal",&learnedVald);
  T2->SetBranchAddress("learnedErr",&learnedErrd);
  
  int totalevents= 0;

  //RSA
  //for(int ientry=0; ientry<T1->GetEntries(); ientry++) {
    //    if(ientry%10000==0)cout<<ientry<<endl;
    for(int ientry=0; ientry<2590000; ientry++) {

    T1->GetEntry(ientry);

    //if(abs(momin[0])>3) {continue;}   To check to see it has any effect on assymettry

    if(!mlApplied && !muMinus) {learnedVal=1.;}
    if(!mlApplied && muMinus)  {learnedVal=0.;}

    thgen = thein[0]*180./pival; phgen = phiin[0]*180./180.; pgen = momin[0];
    threco = therf[0]*180./pival; phreco = phirf[0]*180./180.;
    	
    //if(momrf[0]< 0){preco = -0.0398542 + 1.46176*momrf[0];}//pol1 fit in range -2.0 to -0.6   // not used
    //if(momrf[0]< 0){preco = -0.0162093 + 1.59788*momrf[0];}//pol1 fit in range -2.0 to -0.6   // not used
    	
		
    // if(momrf[0]> 0){preco = -0.014765 + 1.60720*momrf[0];}//pol1 fit in range 0.6 to 2.0
    // if(momrf[0]< 0){preco = -0.185643 + 1.21464*momrf[0];}//pol1 fit in range -2.0 to -0.6


    if(momrf[0]>0) {preco = momrf[0];}
    if(momrf[0]<0) {preco = momrf[0];}


    //if(momrf[0]< 0){preco = momrf[0]-0.2;}
    chi2reco = chisquare[0]; ndfreco = ndof[0];
    bool selected =  (fitfailed==1 && ndfreco>=5 && chi2reco/ndfreco<2.0)? true: false;

    if(fitfailed==1 || fitfailed !=1){
      //if(totalevents>1000000) break;

      //////////////////MuPlus//////////////////////////////////////////////////////////////////////////////
     
      if(selected && learnedVal>0.9 && !muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
	    if(preco>=p_low && preco<=p_high){hist_phi_reco[j]->Fill(phreco);}
	    else{hist_phi_reco[j]->Fill(-200.);}
	  }
	}
	totalevents++;
      } else if(!muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
	    hist_phi_reco[j]->Fill(-200.);
	  }
	}
      }

       
      //////////////MuMinus----------------------------------------------------------
   
      if(selected && learnedVal<0.1 && muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
	    if(preco>=p_low && preco<=p_high) {hist_phi_reco[j]->Fill(phreco);}
	    else{hist_phi_reco[j]->Fill(-200.);}
	  }
	}
	totalevents++;
      } else if(muMinus){
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco >= thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
	    hist_phi_reco[j]->Fill(-200.);						
	  }
	}				
      }

    }
    }

   


  cout<<"MC TotalEvents Reco good"<<totalevents<<endl;
  totalevents=0;
  
  for(int ientry=0; ientry<T2->GetEntries(); ientry++){

    T2->GetEntry(ientry);

    if(!mlApplied && !muMinus) {learnedVald=1.;}
    if(!mlApplied && muMinus)  {learnedVald=0.;}

    threco = therfd[0]*180./pival; phreco = phirfd[0]*180./180.;
		
    //if(momrfd[0]< 0){preco = -0.0398542 + 1.46176*momrfd[0];}//pol1 fit in range -2.0 to -0.6 //not used
    //if(momrfd[0]< 0){preco = -0.0162093 + 1.59788*momrfd[0];}//pol1 fit in range -2.0 to -0.6 //not used
    	
    // if(momrfd[0]>0) {preco = -0.014765 + 1.60720*momrfd[0];}//pol1 fit in range 0.6 to 2.0
    // if(momrfd[0]< 0){preco = -0.185643 + 1.21464*momrfd[0];}//pol1 fit in range -2.0 to -0.6
    	
    if(momrfd[0]>0) {preco = momrfd[0];}
    if(momrfd[0]<0) {preco = momrfd[0];}
		
		
    //if(momrfd[0]< 0){preco = momrfd[0]-0.2;}
    chi2reco = chisquared[0]; ndfreco = ndofd[0];

    bool selected = (fitfailedd==1 && ndfreco>=5  && chi2reco/ndfreco<2.0)? true: false;

    if(fitfailedd==1 || fitfailedd!=1){
      //if(totalevents>1000000) break;
      if(selected && learnedVald>0.9 && !muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
	    if(preco>=p_low && preco<=p_high){hist_phi_data[j]->Fill(phreco);}
	    else{hist_phi_data[j]->Fill(-200.);}
	  }
	}
	totalevents++;
      } else if(!muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
	    hist_phi_data[j]->Fill(-200.);
	  }
	}
      }

      if(selected && learnedVald<0.1 && muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
	    if(preco>=p_low && preco<=p_high){hist_phi_data[j]->Fill(phreco);}
	    else{hist_phi_data[j]->Fill(-200.);}
	  }
	}
	totalevents++;
      } else if(muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
	    hist_phi_data[j]->Fill(-200.);
	  }
	}
      }


    }
  }


  for (int j = 0; j < numThetaRanges; ++j) {
    hist_phi_data[j]->Scale(0.8*2590000.0/T2->GetEntries());
  }

  cout<<"Subtract fakes from fraterate estimated from all models"<<endl;

  for (int j = 0; j < numThetaRanges; ++j) {
    for(int ix=0; ix <((hist_phi_reco[j]->GetNbinsX())); ix++){
      //      cout<<"reco "<<ix+1<<" = "<<1-fakerate[ix+1]<<" * "<<hist_phi_reco[j]->GetBinContent(ix+1)<<endl;
      hist_phi_reco[j]->SetBinContent(ix+1,(1-fake_rate_phi[j]->GetBinContent(ix+1))*(hist_phi_reco[j]->GetBinContent(ix+1)));
      hist_phi_reco[j]->SetBinError(ix+1,(1-fake_rate_phi[j]->GetBinContent(ix+1))*(hist_phi_reco[j]->GetBinError(ix+1)));

    }
    hist_phi_reco[j]->SetBinContent(0,0);
    hist_phi_reco[j]->SetBinError(0,0);

    for(int ix=0; ix <((hist_phi_data[j]->GetNbinsX())); ix++){
      //      cout<<"reco "<<ix+1<<" = "<<1-fakerate[ix+1]<<" * "<<hist_phi_data[j]->GetBinContent(ix+1)<<endl;
      hist_phi_data[j]->SetBinContent(ix+1,(1-fake_rate_phi[j]->GetBinContent(ix+1))*(hist_phi_data[j]->GetBinContent(ix+1)));
      hist_phi_data[j]->SetBinError(ix+1,(1-fake_rate_phi[j]->GetBinContent(ix+1))*(hist_phi_data[j]->GetBinError(ix+1)));
    }
    hist_phi_data[j]->SetBinContent(0,0);
    hist_phi_data[j]->SetBinError(0,0);
  
  }

  cout<<"Clone response, gen from root file"<<endl;
  
  for (int j = 0; j < numThetaRanges; ++j) {
      if(!muMinus) {sprintf(name,"h_phi_gen_muplus_v100_theta_%d_%d", (int)thetaRanges[j][0], (int)thetaRanges[j][1]);                        
	hist_phi_gen[j] = (TH1D*)infileRM->Get(name); hist_phi_gen[j]->Sumw2();}
      if(muMinus)  {sprintf(name,"h_phi_gen_muminus_v100_theta_%d_%d", (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
	hist_phi_gen[j] = (TH1D*)infileRM->Get(name); hist_phi_gen[j]->Sumw2();}
      
      if(!muMinus) {sprintf(name,"h_phi_gen_reco_muplus_v100_theta_%d_%d", (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
	mat_phi_rm[j] = (TH2D*)infileRM->Get(name); mat_phi_rm[j]->Sumw2();}
      if(muMinus)  {sprintf(name,"h_phi_gen_reco_muminus_v100_theta_%d_%d", (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
	mat_phi_rm[j] = (TH2D*)infileRM->Get(name); mat_phi_rm[j]->Sumw2();}

      //Scale:
      mat_phi_rm[j]->Scale(1./8.);
      hist_phi_gen[j]->Scale(1./8.);
     
    }

  
  TFile *fileout ;
  if(!muMinus)sprintf(name,"%s/TUnfold_phi_muplus_v%d.root",saveDir.c_str(),version);
  if(muMinus)sprintf(name,"%s/TUnfold_phi_muminus_v%d.root",saveDir.c_str(),version);
  fileout = new TFile(name,"recreate") ;


  //Phi
  TH1D *hist_phi_unf_reco[numThetaRanges];
  for (int j = 0; j < numThetaRanges; ++j) {
    if(!muMinus) {sprintf(name,"h_phi_unfolded_mc_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_phi_unf_reco[j] = new TH1D(name,name,nphibins_gen[j], phi_low_gen[j], phi_high_gen[j]); hist_phi_unf_reco[j]->Sumw2();}
    if(muMinus)  {sprintf(name,"h_phi_unfolded_mc_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_phi_unf_reco[j] = new TH1D(name,name,nphibins_gen[j], phi_low_gen[j], phi_high_gen[j]); hist_phi_unf_reco[j]->Sumw2();}
  }

  TH1D *hist_phi_unf_data[numThetaRanges];
  for (int j = 0; j < numThetaRanges; ++j) {
    if(!muMinus) {sprintf(name,"h_phi_unfolded_data_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_phi_unf_data[j] = new TH1D(name,name,nphibins_gen[j], phi_low_gen[j], phi_high_gen[j]); hist_phi_unf_data[j]->Sumw2();}
    if(muMinus)  {sprintf(name,"h_phi_unfolded_data_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_phi_unf_data[j] = new TH1D(name,name,nphibins_gen[j], phi_low_gen[j], phi_high_gen[j]); hist_phi_unf_data[j]->Sumw2();}
  }

    ///////////Phi//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    cout<<"unfolding phi"<<endl;

    if(!muMinus){sprintf(name,"Lcurve_mc_reco_v%d_muplus",version);}
    else{sprintf(name,"Lcurve_mc_reco_v%d_muminus",version);}
    TCanvas* canvasm = new TCanvas(name,name,800,1200);
    canvasm->Divide(3,5,1.e-6,1.e-6);
    gStyle->SetLabelSize(0.06,"XYZ");
    gStyle->SetTitleSize(0.06,"XYZ");

    if(!muMinus){sprintf(name,"Lcurve_data_reco_v%d_muplus",version);}
    else{sprintf(name,"Lcurve_data_reco_v%d_muminus",version);}
    
    TCanvas* canvasd = new TCanvas(name,name,800,1200);
    canvasd->Divide(3,5,1.e-6,1.e-6);
    gStyle->SetLabelSize(0.06,"XYZ");
    gStyle->SetTitleSize(0.06,"XYZ");
    
    for (int j = 0; j < numThetaRanges; ++j) {

     TUnfold::ERegMode regMode=TUnfold::kRegModeCurvature;
     //TUnfold::ERegMode regMode=TUnfold::kRegModeSize;

    // preserve the area
    TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;

    TUnfold unfoldBbBreco_phi(mat_phi_rm[j],TUnfold::kHistMapOutputVert,regMode, constraintMode);
    TUnfold unfoldBbBdata_phi(mat_phi_rm[j],TUnfold::kHistMapOutputVert, regMode, constraintMode);
    
    
    unfoldBbBreco_phi.SetInput(hist_phi_reco[j],1.0);
    unfoldBbBdata_phi.SetInput(hist_phi_data[j],1.0);
    //Regulatisation
    unfoldBbBreco_phi.RegularizeBins(3,1,30, regMode);
    unfoldBbBdata_phi.RegularizeBins(3,1,30, regMode);

    Double_t tauMin_phi=pow(10,-1);
    Double_t tauMax_phi=pow(10,-0.2);
    Int_t nScan_phi=50;
    Int_t iBest_phi;
    TSpline *logTauX_phi,*logTauY_phi;
    TSpline *logTauCurvature;
    TGraph *lCurvem_phi;
    TGraph *lCurved_phi;
    cout<<"check2"<<endl;

    Int_t *binMap=new Int_t[nphibins_reco+2];
    for(Int_t i=1;i<=nphibins_reco;i++) binMap[i]=i;
    binMap[0]=-1;
    binMap[nphibins_reco+1]=-1;

    iBest_phi = unfoldBbBreco_phi.ScanLcurve(nScan_phi,tauMin_phi,tauMax_phi,&lCurvem_phi,&logTauX_phi,&logTauY_phi, &logTauCurvature);
    std::cout<<"chi**2_reco="<<unfoldBbBreco_phi.GetChi2A()<<"+"<<unfoldBbBreco_phi.GetChi2L()<<" / "<<unfoldBbBreco_phi.GetNdf()<<"\n";
    Double_t t_phi[1],x_phi[1],y_phi[1];
    logTauX_phi->GetKnot(iBest_phi,t_phi[0],x_phi[0]);
    logTauY_phi->GetKnot(iBest_phi,t_phi[0],y_phi[0]);
    TGraph *bestLcurvem=new TGraph(1,x_phi,y_phi);
    TGraph *bestLogTauLogChi2m=new TGraph(1,t_phi,x_phi);
    
    cout<<"tau = "<<t_phi[0] <<"\tx = "<<x_phi[0]<<"\ty = "<<y_phi[0]<<endl;
    //hist_phi_unf_reco[j] = (TH1D*)unfoldBbBreco_phi.GetOutput("Unfolded");
    unfoldBbBreco_phi.GetOutput(hist_phi_unf_reco[j],binMap);	
    cout<<"check3"<<endl;

    canvasm->cd(3*j+1);
    logTauX_phi->Draw();
    bestLogTauLogChi2m->SetMarkerColor(kRed);
    bestLogTauLogChi2m->Draw("*");
    // show the L curve
    canvasm->cd(3*j+2);
    lCurvem_phi->Draw("AL");
    bestLcurvem->SetMarkerColor(kRed);
    bestLcurvem->Draw("*");
    //logtaucurvature
    canvasm->cd(3*j+3);
    logTauCurvature->Draw();
    if(j==4){
      sprintf(name,"%s/%s.png",saveDir.c_str(),canvasm->GetName());
      canvasm->SaveAs(name);
    }
            
    iBest_phi = unfoldBbBdata_phi.ScanLcurve(nScan_phi,tauMin_phi,tauMax_phi,&lCurved_phi,&logTauX_phi,&logTauY_phi, &logTauCurvature);
    std::cout<<"chi**2_reco="<<unfoldBbBdata_phi.GetChi2A()<<"+"<<unfoldBbBdata_phi.GetChi2L()<<" / "<<unfoldBbBdata_phi.GetNdf()<<"\n";
    logTauX_phi->GetKnot(iBest_phi,t_phi[0],x_phi[0]);
    logTauY_phi->GetKnot(iBest_phi,t_phi[0],y_phi[0]);
    cout<<"tau = "<<t_phi[0] <<"\tx = "<<x_phi[0]<<"\ty = "<<y_phi[0]<<endl;
    TGraph *bestLcurved=new TGraph(1,x_phi,y_phi);
    TGraph *bestLogTauLogChi2d=new TGraph(1,t_phi,x_phi);

    //hist_phi_unf_data[j] = (TH1D*)unfoldBbBdata_phi.GetOutput("Unfolded");
    unfoldBbBdata_phi.GetOutput(hist_phi_unf_data[j],binMap);

    
    
    if(!muMinus) {sprintf(name,"h_phi_unfolded_data_muplus_covariance_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);}
    if(muMinus)  {sprintf(name,"h_phi_unfolded_data_muminus_covariance_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);}
    //unfoldBbBdata_phi.GetEmatrix(covariance_matrix_data[j]);
    //covariance_matrix_data[j]->SetName(name);
    
    if(!muMinus) {sprintf(name,"h_phi_unfolded_mc_muplus_covariance_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);}
    if(muMinus)  {sprintf(name,"h_phi_unfolded_mc_muminus_covariance_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);}
    //unfoldBbBreco_phi.GetEmatrix(covariance_matrix_mc[j]);
    //covariance_matrix_mc[j]->SetName(name);

    
    
    canvasd->cd(3*j+1);
    logTauX_phi->Draw();
    bestLogTauLogChi2d->SetMarkerColor(kRed);
    bestLogTauLogChi2d->Draw("*");
    // show the L curve                                                                                                                                                                                            
    canvasd->cd(3*j+2);
    lCurved_phi->Draw("AL");
    bestLcurved->SetMarkerColor(kRed);
    bestLcurved->Draw("*");
    canvasd->cd(3*j+3);
    logTauCurvature->Draw();

    if(j==4){sprintf(name,"%s/%s.png",saveDir.c_str(),canvasd->GetName());
      canvasd->SaveAs(name);
    }

    //Correct for inefficiecny:   
    for(int ix=0; ix<(hist_phi_unf_reco[j]->GetNbinsX()); ix++){
      cout<<"reco "<<ix+1<<" = "<<efficiency_phi[j]->GetBinContent(ix+1)<<" * "<<hist_phi_unf_reco[j]->GetBinContent(ix+1)<<endl;
      double efficiency;
      efficiency = efficiency_phi[j]->GetBinContent(ix+1);
      if(efficiency<0.0001) {efficiency=0.001;}
      if(hist_phi_unf_reco[j]->GetBinCenter(ix+1)<-1*TMath::Pi() || hist_phi_unf_reco[j]->GetBinCenter(ix+1)>TMath::Pi()) {
	//efficiency = efficiency/hist_phi_unf_reco[j]->GetEntries();
      }      
      //hist_phi_unf_reco[j]->SetBinContent(ix+1,(hist_phi_unf_reco[j]->GetBinContent(ix+1))/efficiency_phi[j]->GetBinContent(ix+1));
      hist_phi_unf_reco[j]->SetBinContent(ix+1,(hist_phi_unf_reco[j]->GetBinContent(ix+1))/efficiency);
      hist_phi_unf_reco[j]->SetBinError(ix+1,(hist_phi_unf_reco[j]->GetBinError(ix+1))/efficiency);
    }
    hist_phi_unf_reco[j]->SetBinContent(0,0);
    hist_phi_unf_reco[j]->SetBinError(0,0);

    for(int ix=0; ix<(hist_phi_unf_data[j]->GetNbinsX()); ix++){
      cout<<"reco "<<ix+1<<" = "<<efficiency_phi[j]->GetBinContent(ix+1)<<" * "<<hist_phi_unf_data[j]->GetBinContent(ix+1)<<endl;

      double efficiency;
      efficiency = efficiency_phi[j]->GetBinContent(ix+1);
      if(efficiency<0.0001) {efficiency=0.001;}
      if(hist_phi_unf_data[j]->GetBinCenter(ix+1)<-1*TMath::Pi() || hist_phi_unf_data[j]->GetBinCenter(ix+1)>TMath::Pi()) {
	//efficiency = efficiency/hist_phi_unf_data[j]->GetEntries();
      }      
      //hist_phi_unf_data[j]->SetBinContent(ix+1,(hist_phi_unf_data[j]->GetBinContent(ix+1))/efficiency_phi[j]->GetBinContent(ix+1));
      hist_phi_unf_data[j]->SetBinContent(ix+1,(hist_phi_unf_data[j]->GetBinContent(ix+1))/efficiency);
      hist_phi_unf_data[j]->SetBinError(ix+1,(hist_phi_unf_data[j]->GetBinError(ix+1))/efficiency);

      
    }
    hist_phi_unf_data[j]->SetBinContent(0,0);
    hist_phi_unf_data[j]->SetBinError(0,0);  
    
   

    hist_phi_data[j]->SetLineColor(kBlack);
    hist_phi_reco[j]->SetLineColor(kBlue);
    hist_phi_gen[j]->SetLineColor(kRed);
    hist_phi_unf_reco[j]->SetLineColor(kMagenta);
    hist_phi_unf_data[j]->SetLineColor(kGreen);

    hist_phi_reco[j]->SetStats(0);
    hist_phi_reco[j]->SetStats(0);
    hist_phi_gen[j]->SetStats(1);
    hist_phi_unf_reco[j]->SetStats(0);
    hist_phi_unf_data[j]->SetStats(0);
    hist_phi_gen[j]->SetLineStyle(kDashed);

    // unfoldBbBreco_phi.Delete();
    // unfoldBbBdata_phi.Delete();
   
    //JMJ correction for phi

    for(int kl=1; kl<=hist_phi_unf_reco[j]->GetNbinsX();kl++) {
      if(hist_phi_unf_reco[j]->GetBinCenter(kl)<-1*TMath::Pi()) {
	double binerr=hist_phi_unf_reco[j]->GetBinError(hist_phi_unf_reco[j]->FindBin(2*TMath::Pi()+hist_phi_unf_reco[j]->GetBinCenter(kl)));
	hist_phi_unf_reco[j]->Fill(2*TMath::Pi()+hist_phi_unf_reco[j]->GetBinCenter(kl), hist_phi_unf_reco[j]->GetBinContent(kl));
	hist_phi_unf_reco[j]->SetBinError(hist_phi_unf_reco[j]->FindBin(2*TMath::Pi()+hist_phi_unf_reco[j]->GetBinCenter(kl)),binerr);//Setting bin error
	hist_phi_unf_reco[j]->SetBinContent(kl,0.);
      }
      else if(hist_phi_unf_reco[j]->GetBinCenter(kl)>TMath::Pi()) {
	double binerr=hist_phi_unf_reco[j]->GetBinError(hist_phi_unf_reco[j]->FindBin(-2*TMath::Pi()+hist_phi_unf_reco[j]->GetBinCenter(kl)));
	hist_phi_unf_reco[j]->Fill(-2*TMath::Pi()+hist_phi_unf_reco[j]->GetBinCenter(kl), hist_phi_unf_reco[j]->GetBinContent(kl));
	hist_phi_unf_reco[j]->SetBinError(hist_phi_unf_reco[j]->FindBin(-2*TMath::Pi()+hist_phi_unf_reco[j]->GetBinCenter(kl)),binerr);//Setting bin error
	hist_phi_unf_reco[j]->SetBinContent(kl,0.);
      }	
    }

    for(int kl=1; kl<=hist_phi_unf_data[j]->GetNbinsX();kl++) {
      if(hist_phi_unf_data[j]->GetBinCenter(kl)<-1*TMath::Pi()) {
	double binerr=hist_phi_unf_data[j]->GetBinError(hist_phi_unf_data[j]->FindBin(2*TMath::Pi()+hist_phi_unf_data[j]->GetBinCenter(kl)));
	hist_phi_unf_data[j]->Fill(2*TMath::Pi()+hist_phi_unf_data[j]->GetBinCenter(kl), hist_phi_unf_data[j]->GetBinContent(kl));
	hist_phi_unf_data[j]->SetBinError(hist_phi_unf_data[j]->FindBin(2*TMath::Pi()+hist_phi_unf_data[j]->GetBinCenter(kl)),binerr);//Setting bin error	
	hist_phi_unf_data[j]->SetBinContent(kl,0.);
      }
      else if(hist_phi_unf_data[j]->GetBinCenter(kl)>TMath::Pi()) {
	double binerr=hist_phi_unf_data[j]->GetBinError(hist_phi_unf_data[j]->FindBin(-2*TMath::Pi()+hist_phi_unf_data[j]->GetBinCenter(kl)));
	hist_phi_unf_data[j]->Fill(-2*TMath::Pi()+hist_phi_unf_data[j]->GetBinCenter(kl), hist_phi_unf_data[j]->GetBinContent(kl));
	hist_phi_unf_data[j]->SetBinError(hist_phi_unf_data[j]->FindBin(-2*TMath::Pi()+hist_phi_unf_data[j]->GetBinCenter(kl)),binerr);//Setting bin error	
	hist_phi_unf_data[j]->SetBinContent(kl,0.);
      }	
    }

   //Gen Correction                                                                                                                                                                                               
    for(int kl=1; kl<=hist_phi_gen[j]->GetNbinsX();kl++) {
      if(hist_phi_gen[j]->GetBinCenter(kl)<-1*TMath::Pi()) {
        double binerr=hist_phi_gen[j]->GetBinError(hist_phi_gen[j]->FindBin(2*TMath::Pi()+hist_phi_gen[j]->GetBinCenter(kl)));
        hist_phi_gen[j]->Fill(2*TMath::Pi()+hist_phi_gen[j]->GetBinCenter(kl), hist_phi_gen[j]->GetBinContent(kl));
	hist_phi_gen[j]->SetBinError(hist_phi_gen[j]->FindBin(2*TMath::Pi()+hist_phi_gen[j]->GetBinCenter(kl)),binerr);//Setting bin error                   
        hist_phi_gen[j]->SetBinContent(kl,0.);
      }
      else if(hist_phi_gen[j]->GetBinCenter(kl)>TMath::Pi()) {
        double binerr=hist_phi_gen[j]->GetBinError(hist_phi_gen[j]->FindBin(-2*TMath::Pi()+hist_phi_gen[j]->GetBinCenter(kl)));
	hist_phi_gen[j]->Fill(-2*TMath::Pi()+hist_phi_gen[j]->GetBinCenter(kl), hist_phi_gen[j]->GetBinContent(kl));
        hist_phi_gen[j]->SetBinError(hist_phi_gen[j]->FindBin(-2*TMath::Pi()+hist_phi_gen[j]->GetBinCenter(kl)),binerr);//Setting bin error                  
	hist_phi_gen[j]->SetBinContent(kl,0.);
      }
    }

    

    }//j loop


  if(muMinus){
    gStyle->SetStatX(0.3);
    gStyle->SetStatY(0.3);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.25);
  }

  //Phi

  if(!muMinus){
    sprintf(name,"Tunfold_output_mc_reco_phi_muplus_v%d",version);
  }
  else{
    sprintf(name,"Tunfold_output_mc_reco_phi_muminus_v%d",version);
  }



  TCanvas *c3 = new TCanvas(name, name, 800, 1200);
  c3->Divide(2,3);

  for (int j = 0; j < numThetaRanges; ++j) {
    TLegend *leg1 = new TLegend(0.4,0.2,0.6,0.4);
    leg1->SetBorderSize(0);
    c3->cd(j+1);

    //gPad->SetLogx(1);
    gPad->SetLogy(1);

    hist_phi_gen[j]->GetYaxis()->SetRangeUser(100,300000);
    //hist_phi_gen[j]->SetTitle("MC Reco Unfolding Without ML");

    hist_phi_gen[j]->Draw("hist");  leg1->AddEntry(hist_phi_gen[j],"GEN","l");
    hist_phi_data[j]->Draw("hist:same");  leg1->AddEntry(hist_phi_data[j],"Data Reco","l");
    hist_phi_reco[j]->Draw("hist:sames"); leg1->AddEntry(hist_phi_reco[j],"MC Reco","l");
    hist_phi_unf_reco[j]->Draw("sames"); leg1->AddEntry(hist_phi_unf_reco[j],"Unfolded","l");

    leg1->Draw();

  }

  sprintf(name,"%s/%s.png",saveDir.c_str(), c3->GetName());
  c3->SaveAs(name);
  fileout->cd();
  c3->Write();


  if(!muMinus){
    sprintf(name, "Tunfold_output_data_reco_phi_muplus_v%d", version); 
  }
  else{
    sprintf(name, "Tunfold_output_data_reco_phi_muminus_v%d", version); 
  }


  TCanvas *c4 = new TCanvas(name, name, 800, 1200);
  c4->Divide(2,3);
  for (int j = 0; j < numThetaRanges; ++j) {
    TLegend *leg2 = new TLegend(0.4,0.2,0.6,0.4);
    leg2->SetBorderSize(0);
    c4->cd(j+1);
    gPad->SetLogy(1);

    hist_phi_gen[j]->GetYaxis()->SetRangeUser(100,300000);
    //hist_phi_gen[j]->SetTitle("Data Reco Unfolding Without ML");

    hist_phi_gen[j]->Draw("hist");  leg2->AddEntry(hist_phi_gen[j],"GEN","l");
    hist_phi_data[j]->Draw("hist:same");  leg2->AddEntry(hist_phi_data[j],"Data Reco","l");
    hist_phi_reco[j]->Draw("hist:sames"); leg2->AddEntry(hist_phi_reco[j],"MC Reco","l");
    hist_phi_unf_data[j]->Draw("sames"); leg2->AddEntry(hist_phi_unf_data[j],"Unfolded","l");

    leg2->Draw();

        
  }

  sprintf(name,"%s/%s.png",saveDir.c_str(), c4->GetName());
  c4->SaveAs(name);
  fileout->cd();
  c4->Write();


  fileout->cd();

  for (int j = 0; j < numThetaRanges; ++j) {
    hist_phi_gen[j]->Write();
    hist_phi_reco[j]->Write();
    hist_phi_data[j]->Write();
    hist_phi_unf_reco[j]->Write();
    hist_phi_unf_data[j]->Write();
    mat_phi_rm[j]->Write();
    //covariance_matrix_mc[j]->Write();
    //covariance_matrix_data[j]->Write();
  }
  fileout->Write();
  fileout->Close();

  infileRM->cd();
  infileRM->Close();

  return 0;
}
