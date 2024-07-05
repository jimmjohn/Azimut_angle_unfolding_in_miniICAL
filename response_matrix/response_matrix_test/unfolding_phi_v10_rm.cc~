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
//int nphibins_gen[5] = {30,30,30,30,30};
int nphibins_gen[5] = {60,60,60,60,60};
int nphibins_reco = 70;
double phi_low = -1.*TMath::Pi();
double phi_high = TMath::Pi();

//double phi_low_gen[5] =  {-50.*TMath::Pi()/30., -38.*TMath::Pi()/30., -34.*TMath::Pi()/30., -34.*TMath::Pi()/30., -34.*TMath::Pi()/30.};
//double phi_high_gen[5] = { 50.*TMath::Pi()/30.,  38.*TMath::Pi()/30.,  34.*TMath::Pi()/30.,  34.*TMath::Pi()/30.,  34.*TMath::Pi()/30.};
//double phi_low_gen[5] =  {-30.*TMath::Pi()/30., -30.*TMath::Pi()/30., -30.*TMath::Pi()/30., -30.*TMath::Pi()/30., -30.*TMath::Pi()/30.};
//double phi_high_gen[5] = { 30.*TMath::Pi()/30.,  30.*TMath::Pi()/30.,  30.*TMath::Pi()/30.,  30.*TMath::Pi()/30.,  30.*TMath::Pi()/30.};
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


  int version = atoi(argv[1]);
  string saveDir = "fileOut";
  
//  string filename = "Mc_v107_2018_Data_Magnetic.root";
string filename = argv[2];
 int myndof = atoi(argv[4]);
 
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


  TH1D *fake_rate_phi[numThetaRanges];
  TH1D *efficiency_phi[numThetaRanges];

  for (int j = 0; j < numThetaRanges; ++j) {
    if(!muMinus) {
      sprintf(name,"fake_rate_phi_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      fake_rate_phi[j] = new TH1D(name,name,nphibins_reco, phi_low, phi_high);
      sprintf(name,"efficiency_phi_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      efficiency_phi[j] = new TH1D(name,name,nphibins_gen[j], phi_low_gen[j], phi_high_gen[j]);
    }
    if(muMinus)  {
      sprintf(name,"fake_rate_phi_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      fake_rate_phi[j] = new TH1D(name,name,nphibins_reco, phi_low, phi_high);
      sprintf(name,"efficiency_phi_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      efficiency_phi[j] = new TH1D(name,name,nphibins_gen[j], phi_low_gen[j], phi_high_gen[j]);
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


  //sprintf(name,"%s",filename.c_str());
  sprintf(name,"/var/nfscondor/rajshah/paper2/mom_unfolding_fiducial/%s",filename.c_str());
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


  //RSA
  for(int ientry=0; ientry<T1->GetEntries(); ientry++) {
    //    if(ientry%10000==0)cout<<ientry<<endl;
  //  for(int ientry=0; ientry<2590000; ientry++) {

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
    bool selected =  (fitfailed==1 && ndfreco>=myndof && chi2reco/ndfreco<2.0)? true: false;


    if(selected==1 && abs(phgen-phreco)>TMath::Pi() && preco>p_low && preco<p_high) { //JMJ
      if(phgen<0) {phgen=phgen+2.*TMath::Pi();}
      else {phgen=phgen-2.*TMath::Pi();}
    }

    bool skipOutRange = false;
    for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1] && selected==1) {
	    if(phgen>hist_phi_gen[j]->GetXaxis()->GetXmax() || phgen<hist_phi_gen[j]->GetXaxis()->GetXmin()) {skipOutRange=true;break;}
	  }
    }
    //    if(skipOutRange) {continue;}

    if(fitfailed==1 || fitfailed !=1){
      //if(totalevents>1000000) break;

      //////////////////MuPlus//////////////////////////////////////////////////////////////////////////////
      //if(pgen>p_low && pgen<p_high) {   //pgen should be outside
      if(!muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
            if(pgen>=p_low && pgen<=p_high ){hist_phi_gen[j]->Fill(phgen);}
	  }
	}
      }

      if(selected && learnedVal>0.9 && !muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
            if(preco>=p_low && preco<=p_high){hist_phi_reco[j]->Fill(phreco);}
	    else{hist_phi_reco[j]->Fill(-200.);}
	  }
	}
      } else if(!muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
            hist_phi_reco[j]->Fill(-200.);
	  }
	}
      }

      // if(selected && preco>p_low && preco<p_high && ientry<T1->GetEntries()/2.) {
      // 		hist_p_data->Fill(preco);
      // } else if(ientry<T1->GetEntries()/2.){hist_p_data->Fill(-100);}

      // if(selected && learnedVal>0.9 && !muMinus) {
      // 	mat_th_rm->Fill(threco,thgen);
      // 	for (int j = 0; j < numThetaRanges; ++j) {
      // 		if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
      // 			mat_p_rm[j]->Fill(preco,pgen);                                     //Selected events
      // 			if(preco>=p_low && preco<=p_high && pgen>=p_low && pgen<=p_high)mat_phi_rm[j]->Fill(phreco,phgen);
      // 			if(preco>=p_low && preco<=p_high && (pgen<p_low || pgen>p_high)) mat_phi_rm[j]->Fill(phreco,-200.); //Fakes
      // 			if((preco<p_low || preco>p_high) && pgen>=p_low && pgen<=p_high )mat_phi_rm[j]->Fill(-200., phgen); //Misses
      // 			if((preco<p_low || preco>p_high) && (pgen<p_low || pgen>p_high)) mat_phi_rm[j]->Fill(-200.,-200.); //out of range
      // 		}
      // 	}

      // } else if(!muMinus) {
      // 	mat_th_rm->Fill(-100., thgen);
      // 	for (int j = 0; j < numThetaRanges; ++j) {
      // 		if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
      // 			mat_p_rm[j]->Fill(-100.,pgen);                                    //Selected events
      // 			mat_phi_rm[j]->Fill(-200.,phgen);                                 //Misses
      // 		}
      // 	}
      // }

      if(selected && learnedVal>0.9 && !muMinus && preco>=p_low && preco<=p_high && pgen>=p_low && pgen<=p_high) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
	    mat_phi_rm[j]->Fill(phreco,phgen);
	  }
	}
      }

      if(selected && learnedVal>0.9 && !muMinus && preco>=p_low && preco<=p_high && (pgen<p_low || pgen>p_high)) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
	    //Fakes
	    mat_phi_rm[j]->Fill(phreco,-200.);
	  }
	}
      }

      if((!selected || learnedVal<0.9 || preco<p_low || preco>p_high) && !muMinus && pgen>=p_low && pgen<=p_high) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
	    //Misses
	    mat_phi_rm[j]->Fill(-200,phgen);
	  }
	}
      }


      //////////////MuMinus----------------------------------------------------------
      //Gen
      if(muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
            if(pgen<0 && pgen>=p_low && pgen<=p_high){hist_phi_gen[j]->Fill(phgen);}
	  }
	}
      }
      //Reco
      // if(selected && learnedVal<0.1 && muMinus) {
      // 	hist_th_reco->Fill(threco);
      // 	for (int j = 0; j < numThetaRanges; ++j) {
      // 		if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
      // 			hist_p_reco[j]->Fill(preco);
      // 			if(preco>=p_low && preco<=p_high){hist_phi_reco[j]->Fill(phreco);}
      // 			else{hist_phi_reco[j]->Fill(-200.);}
      // 		}
      // 	}
      // 		totalevents++;
      // } else if(muMinus) {
      // 	hist_th_reco->Fill(-100.);
      // 	for (int j = 0; j < numThetaRanges; ++j) {
      // 		if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
      // 			hist_p_reco[j]->Fill(100);
      // 			hist_phi_reco[j]->Fill(-200.);
      // 		}
      // 	}
      // }

      
      if(selected && learnedVal<0.1 && muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
	    if(preco>=p_low && preco<=p_high) {hist_phi_reco[j]->Fill(phreco);}
            else {hist_phi_reco[j]->Fill(-200.);}
	  }
	}
      } else if(muMinus){
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco >= thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
	    hist_phi_reco[j]->Fill(-200.);						
	  }
	}				
      }


      //Response Matrix
      // if(selected && learnedVal<0.1 && muMinus) {
      // 	mat_th_rm->Fill(threco, thgen);
      // 	for (int j = 0; j < numThetaRanges; ++j) {
      // 		if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
      // 			mat_p_rm[j]->Fill(preco,pgen);                                    //Selected events
      // 			if(preco>=p_low && preco<=p_high && pgen>=p_low && pgen<=p_high) {mat_phi_rm[j]->Fill(phreco,phgen);}
      // 			if(preco>=p_low && preco<=p_high && (pgen<p_low || pgen>p_high)) {mat_phi_rm[j]->Fill(phreco,-200.);} //Fakes
      // 			if((preco<p_low || preco>p_high) && pgen>=p_low && pgen<=p_high){mat_phi_rm[j]->Fill(-200., phgen);} //Misses
      // 			if((preco<p_low || preco>p_high) && (pgen<p_low || pgen>p_high)) {mat_phi_rm[j]->Fill(-200.,-200.);} //out of range
      // 		}
      // 	}

      // } else if(muMinus) {
      // 	mat_th_rm->Fill(-100.,thgen);
      // 	for (int j = 0; j < numThetaRanges; ++j) {
      // 		if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
      // 			mat_p_rm[j]->Fill(100.,pgen);                                    //Selected events
      // 			mat_phi_rm[j]->Fill(-200.,phgen);                                //Misses
      // 		}
      // 	}
      // }

      if(selected && learnedVal<0.1 && muMinus && preco>=p_low && preco<=p_high && pgen>=p_low && pgen<=p_high ) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
            mat_phi_rm[j]->Fill(phreco,phgen);
	  }
	}
      }

      if(selected && learnedVal<0.1 && muMinus && preco>=p_low && preco<=p_high && (pgen<p_low || pgen>p_high)) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
	    //Fakes
	    mat_phi_rm[j]->Fill(phreco,-200.);
	  }
	}
      }

      if((!selected || learnedVal>0.1 || preco<p_low || preco>p_high) && muMinus && pgen>=p_low && pgen<=p_high) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco <= thetaRanges[j][1]) {
	    //Misses
	    mat_phi_rm[j]->Fill(-200,phgen);
	  }
	}
      }



      /*
	if(pgen>p_low && pgen<p_high && preco>p_low && preco<p_high && selected) {
	mat_p_rm->Fill(preco,pgen);
	}
	if((pgen<p_low || pgen>p_high) && preco>p_low && preco<p_high && selected) {   //Fakes
	mat_p_rm->Fill(preco,-100);
	}
	// if(pgen>p_high  && preco>p_low && preco<p_high) {
	// 	mat_p_rm->Fill(preco,900);
	// }
	if(pgen>p_low && pgen<p_high && (preco<p_low || preco>p_high || (!selected))) {   //Misses
	mat_p_rm->Fill(-100,pgen);
	}
      */
    }
    }

    

  /*
    for(int ix=0; ix<(hist_p_gen->GetNbinsX()); ix++){
    cout<<"gen "<<ix+1<<" = "<<effi[ix+1]<<" * "<<hist_p_gen->GetBinContent(ix+1)<<endl;
    hist_p_gen->SetBinContent(ix+1,(hist_p_gen->GetBinContent(ix+1))*effi[ix+1]);
    hist_p_gen->SetBinError(ix+1,(hist_p_gen->GetBinError(ix+1))*effi[ix+1]);
    }

    hist_p_gen->SetBinContent(0,0);
    hist_p_gen->SetBinError(0,0);
  */


  //Fake Rate Calcualtion for phi
  int nbinx = mat_phi_rm[0]->GetNbinsX();
  int nbiny = mat_phi_rm[0]->GetNbinsY();

  static const int nbinmx = 120; //max(hResp->GetNbinsY(),hResp->GetNbinsX()) + 5;
  double totalgen[nbinmx]={0.};
  double totalreco[nbinmx]={0.};

  double fakerate[nbinmx];
  double effi[nbinmx];


  //TH2D* mat_p_rm_initial = (TH2D*)mat_p_rm->Clone();

  for (int j = 0; j < numThetaRanges; ++j) {
    cout<<"j: "<<j<<endl;
    for(int ib=0; ib< nbinmx; ib++){
      totalgen[ib] = 0;
      totalreco[ib] = 0;
      fakerate[ib] = 0;
      effi[ib] = 0;
    }

    for (int ix=0; ix<nbinx+1; ix++) {
      for (int iy=0; iy<nbiny+1; iy++) {
	if(ix==0&&iy==0) continue;
	if(mat_phi_rm[j]->GetBinContent(ix,iy) < -0.1)	mat_phi_rm[j]->SetBinContent(ix,iy,0);
	totalreco[ix] += mat_phi_rm[j]->GetBinContent(ix, iy);          // Total number of reco events in bin ix
	if (iy==0) fakerate[ix] = mat_phi_rm[j]->GetBinContent(ix,iy);  //Number of muMinus events accepted in bin ix
	totalgen[iy] +=mat_phi_rm[j]->GetBinContent(ix, iy);            //Total number of gen events in bin iy
	if (ix==0) effi[iy] = mat_phi_rm[j]->GetBinContent(ix, iy) ;     //Number events not reconstructed but generated in bin iy   // Actually it should be ineffi
	if (ix==0 || iy==0) {
	  mat_phi_rm[j]->SetBinContent(ix, iy, 0.0);
	  mat_phi_rm[j]->SetBinError(ix, iy, 0.0);
	}
      }//iy
    }//ix

    for (int iy=1; iy<nbiny+1; iy++) {
      effi[iy] = (totalgen[iy] - effi[iy])/max(1.e-10, totalgen[iy]);      //Now inefficiency changed to efficiency
      if (abs(effi[iy]) > 1) effi[iy] = 1;
      else if ( effi[iy] < 0) effi[iy] = 0.0000001;
      efficiency_phi[j]->SetBinContent(iy, effi[iy]);
    } //iy

    for (int ix=1; ix<nbinx+1; ix++) {
      //      cout<<"jim, ix="<<ix<<"\tfake="<<fakerate[ix]<<"\ttotal="<<totalreco[ix]<<endl;//Actually fakes are too much for the same width bins
      fakerate[ix] = fakerate[ix] / max(1.e-10, totalreco[ix]);
      if(abs(fakerate[ix]) > 1) fakerate[ix] = 0.99999999;
      else if ( fakerate[ix] < 0) fakerate[ix] = 0.0000001;
      fake_rate_phi[j]->SetBinContent(ix, fakerate[ix]);
    }//ixfakefake

    // for(int ix=0; ix <((hist_phi_reco[j]->GetNbinsX())); ix++){
    //   //      cout<<"reco "<<ix+1<<" = "<<1-fakerate[ix+1]<<" * "<<hist_phi_reco[j]->GetBinContent(ix+1)<<endl;
    //   hist_phi_reco[j]->SetBinContent(ix+1,(1-fakerate[ix+1])*(hist_phi_reco[j]->GetBinContent(ix+1)));
    //   //hist_phi_reco[j]->SetBinError(ix+1,(1-fakerate[ix+1])*(hist_phi_reco[j]->GetBinError(ix+1)));
    // }
    // hist_phi_reco[j]->SetBinContent(0,0);
    // hist_phi_reco[j]->SetBinError(0,0);

    // for(int ix=0; ix <((hist_phi_data[j]->GetNbinsX())); ix++){
    //   //      cout<<"reco "<<ix+1<<" = "<<1-fakerate[ix+1]<<" * "<<hist_phi_data[j]->GetBinContent(ix+1)<<endl;
    //   hist_phi_data[j]->SetBinContent(ix+1,(1-fakerate[ix+1])*(hist_phi_data[j]->GetBinContent(ix+1)));
    //   //hist_phi_data[j]->SetBinError(ix+1,(1-fakerate[ix+1])*(hist_phi_data[j]->GetBinError(ix+1)));
    // }
    // hist_phi_data[j]->SetBinContent(0,0);
    // hist_phi_data[j]->SetBinError(0,0);
  
  } 
 

  
  ////RSA construct response matrix from all models

  TFile* outfileforRM;
  if(!muMinus) {
    outfileforRM = new TFile(TString::Format("fileOut/phi_response_matrix_v100_muplus_ndof%d.root",myndof),"RECREATE");
  }
  else{
    outfileforRM = new TFile(TString::Format("fileOut/phi_response_matrix_v100_muminus_ndof%d.root",myndof),"RECREATE");

  }

  outfileforRM->cd();
  for (int j = 0; j < numThetaRanges; ++j) {	  
    mat_phi_rm[j]->Write();  
    hist_phi_gen[j]->Write();
    hist_phi_reco[j]->Write();
    fake_rate_phi[j]->Write();
    efficiency_phi[j]->Write();
  }

  outfileforRM->Write();       
  outfileforRM->Close();

  
  return 0;
}
