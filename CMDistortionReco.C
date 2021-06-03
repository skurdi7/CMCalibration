// step 2
#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TTree.h"
#include "/sphenix/u/skurdi/CMCalibration/PHG4TpcCentralMembrane.h"
R__LOAD_LIBRARY(build/.libs/libg4tpccentralmembrane)

//from phg4tpcsteppingaction.cc
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
R__LOAD_LIBRARY(libphg4hit.so)


using namespace std;

class Shifter {
public:
Shifter(TString sourcefilename);
  TFile *forward, *average;
  TH3F *hX, *hY, *hZ, *hR, *hPhi, *hXave, *hYave, *hZave, *hRave, *hPhiave;  
};

Shifter::Shifter(TString sourcefilename){
  //single event distortion file
  forward=TFile::Open(sourcefilename,"READ"); 

  hX=(TH3F*)forward->Get("hIntDistortionPosX");
  hY=(TH3F*)forward->Get("hIntDistortionPosY");
  hZ=(TH3F*)forward->Get("hIntDistortionPosZ");

  hR=(TH3F*)forward->Get("hIntDistortionPosR");
  hPhi=(TH3F*)forward->Get("hIntDistortionPosP");

  //average distortion file
  average=TFile::Open("/sphenix/user/rcorliss/distortion_maps/2021.04/apr07.average.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ"); 
  
  hXave=(TH3F*)average->Get("hIntDistortionX");
  hYave=(TH3F*)average->Get("hIntDistortionY");
  hZave=(TH3F*)average->Get("hIntDistortionZ");
  
  hRave=(TH3F*)average->Get("hIntDistortionR");
  hPhiave=(TH3F*)average->Get("hIntDistortionP");

  //subtract average from total distortions to study fluctuations
  hX->Add(hXave,-1);
  hY->Add(hYave,-1);
  hZ->Add(hZave,-1);
  
  hR->Add(hRave,-1);
  hPhi->Add(hPhiave,-1);
}

int CMDistortionReco() {
  Shifter *shifter;
  int nbins = 35; 
  double low = -80.0;
  double high = 80.0;
  double deltaX, deltaY, deltaZ, deltaR, deltaPhi;

  int nEvents = 3; //change based on number of tree files available in source directory
    
  //take in events
  const char * inputpattern="/sphenix/user/rcorliss/distortion_maps/2021.04/*h_Charge_*.root"; 
  
  //find all files that match the input string (includes wildcards)
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  TString sourcefilename;
  
  TCanvas *canvas=new TCanvas("canvas","CMDistortionReco",2000,3000);

  int nsumbins = 20;
  int minsum = -10;
  int maxsum = 10;
  
  //set up summary plots
  TH1F *hDifferenceMeanR = new TH1F("hDifferenceMeanR", "Average Difference between R Model and True of All Events (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    TH1F *hDifferenceStdDevR = new TH1F("hDifferenceStdDevR", "Std Dev of Difference between R Model and True of All Events (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    
    TH1F *hTrueMeanR = new TH1F("hTrueMeanR", "Mean True R Distortion Model of All Events (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    TH1F *hTrueStdDevR = new TH1F("hTrueStdDevR", "Std Dev of True R Distortion Model of All Events (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    
    TH1F *hDifferenceMeanPhi = new TH1F("hDifferenceMeanPhi", "Average Difference between Phi Model and True of All Events (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
    TH1F *hDifferenceStdDevPhi = new TH1F("hDifferenceStdDevPhi", "Std Dev of Difference between Phi Model and True of All Events (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
    
    TH1F *hTrueMeanPhi = new TH1F("hTrueMeanPhi", "Mean True Phi Distortion Model of All Events (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
    TH1F *hTrueStdDevPhi = new TH1F("hTrueStdDevPhi", "Std Dev of True Phi Distortion Model of All Events (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);

    TVector3 *position, *newposition;
     position = new TVector3(1.,1.,1.);
     newposition = new TVector3(1.,1.,1.);
 
  for (int ifile=0;ifile < nEvents;ifile++){
    //for each file, find all histograms in that file
    sourcefilename=((TFileInfo*)(filelist->GetList()->At(ifile)))->GetCurrentUrl()->GetFile();

    //create shifter
    shifter = new Shifter(sourcefilename);
    
    //get data from ttree
    char const *treename="cmDistHitsTree";
    TFile *input=TFile::Open(Form("cmDistHitsTree_Event%d.root", ifile));
    TTree *inTree=(TTree*)input->Get("tree");
    
    inTree->SetBranchAddress("position",&position);
    inTree->SetBranchAddress("newposition",&newposition);
    

    // 0 to 2pi for phi, 0 to 90 for r
    
    //for forward only
   
    TH2F *hStripesPerBin = new TH2F("hStripesPerBin","CM Stripes Per Bin (z in stripes); x (cm); y (cm)",nbins,low,high,nbins,low,high); 

    //TH2F *hStripesPerBinRPhi = new TH2F("hStripesPerBinRPhi","CM Stripes Per Bin (z in stripes); phi (rad); r (cm)",nbins,low,high,nbins,low,high); 
    
    TH2F *hCartesianForward[3];
    hCartesianForward[0] = new TH2F("hForwardX","X Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianForward[1] = new TH2F("hForwardY","Y Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianForward[2] = new TH2F("hForwardZ","Z Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);

    TH2F *hCylindricalForward[2];
    hCylindricalForward[0] = new TH2F("hForwardR","Radial Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalForward[1] = new TH2F("hForwardPhi","Phi Shift Forward of Stripe Centers (rad); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  
  
    for (int i=0;i<inTree->GetEntries();i++){
      inTree->GetEntry(i);
      
      double r = position->Perp();
    
      double phi = position->Phi();
      if(position->Phi() < 0.0){
	phi = position->Phi() + TMath::TwoPi(); 
      }

      hStripesPerBin->Fill(position->X(),position->X(),1);
      
      deltaX = (newposition->X() - position->X())*(1e4); //convert from cm to micron 
      deltaY = (newposition->Y() - position->Y())*(1e4);
      deltaZ = (newposition->Z() - position->Z())*(1e4);

      deltaR = (newposition->Perp() - position->Perp())*(1e4);
      deltaPhi = newposition->DeltaPhi(*position);

      hCartesianForward[0]->Fill(position->X(),position->Y(),deltaX);
      hCartesianForward[1]->Fill(position->X(),position->Y(),deltaY);
      hCartesianForward[2]->Fill(position->X(),position->Y(),deltaZ);

      hCylindricalForward[0]->Fill(position->X(),position->Y(),deltaR);
      hCylindricalForward[1]->Fill(position->X(),position->Y(),deltaPhi);
    }

    TH2F *hCartesianAveShift[3];
    hCartesianAveShift[0] = new TH2F("AveShiftX","Average of CM Model X over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCartesianAveShift[1] = new TH2F("AveShiftY","Average of CM Model Y over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCartesianAveShift[2] = new TH2F("AveShiftZ","Average of CM Model Z over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 

    TH2F *hCylindricalAveShift[4];
    hCylindricalAveShift[0] = new TH2F("AveShiftR","Average of CM Model R over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCylindricalAveShift[1] = new TH2F("AveShiftPhi","Average of CM Model Phi over Stripes per Bin (rad); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCylindricalAveShift[2] = new TH2F("AveShiftRCart","Average of CM Model R over Stripes per Bin from Cartesian (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCylindricalAveShift[3] = new TH2F("AveShiftPhiCart","Average of CM Model Phi over Stripes per Bin from Cartesian (rad); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    

    for (int i = 0; i < 3; i ++){
      hCartesianAveShift[i]->Divide(hCartesianForward[i],hStripesPerBin);
    }
  
    hCylindricalAveShift[0]->Divide(hCylindricalForward[0],hStripesPerBin);
    hCylindricalAveShift[1]->Divide(hCylindricalForward[1],hStripesPerBin);

    
    //r from cart in loop over bins below
    for(int i = 0; i < nbins; i++){
      double x = low + ((high - low)/(1.0*nbins))*(i+0.5); //center of bin
      for(int j = 0; j < nbins; j++){
	double y = low + ((high - low)/(1.0*nbins))*(j+0.5); //center of bin
      
	int xbin = hCartesianAveShift[0]->FindBin(x,y);
	int ybin = hCartesianAveShift[1]->FindBin(x,y);
	double xaveshift = (hCartesianAveShift[0]->GetBinContent(xbin))*(1e-4); // converts  microns to cm 
	double yaveshift = (hCartesianAveShift[1]->GetBinContent(ybin))*(1e-4);
	
	TVector3 shifted, original;
	original.SetX(x);
	original.SetY(y);
	shifted.SetX(x+xaveshift);
	shifted.SetY(y+yaveshift);
	
	double raveshift = (shifted.Perp() - original.Perp())*(1e4);
       	double phiaveshift = shifted.DeltaPhi(original);
	
	hCylindricalAveShift[2]->Fill(x,y,raveshift);
	hCylindricalAveShift[3]->Fill(x,y,phiaveshift);
      } 
    } 

    //same range and bins for each coordinate, can use hR for all, binned in cm
    int nphi = shifter->hR->GetXaxis()->GetNbins();
    int nr = shifter->hR->GetYaxis()->GetNbins();
    int nz = (shifter->hR->GetZaxis()->GetNbins())/2;
    
    double minphi = shifter->hR->GetXaxis()->GetXmin();
    double minr = shifter->hR->GetYaxis()->GetXmin();
    double minz = 5.0;
    
    double maxphi = shifter->hR->GetXaxis()->GetXmax();
    double maxr = shifter->hR->GetYaxis()->GetXmax();
    double maxz = shifter->hR->GetZaxis()->GetXmax();
    
    TH3F *hCartesianCMModel[3];
    hCartesianCMModel[0]=new TH3F("hCMModelX", "CM Model: X Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz); //rad, cm, cm
    hCartesianCMModel[1]=new TH3F("hCMModelY", "CM Model: Y Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
    hCartesianCMModel[2]=new TH3F("hCMModelZ", "CM Model: Z Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);

    TH3F *hCylindricalCMModel[4];
    hCylindricalCMModel[0]=new TH3F("hCMModelR", "CM Model: Radial Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
    hCylindricalCMModel[1]=new TH3F("hCMModelPhi", "CM Model: Phi Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
    hCylindricalCMModel[2]=new TH3F("hCMModelRCart", "CM Model: Radial Shift Forward of Stripe Centers from Cartesian", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
    hCylindricalCMModel[3]=new TH3F("hCMModelPhiCart", "CM Model: Phi Shift Forward of Stripe Centers from Cartesian", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);  
    
    double xshift, yshift, zshift, rshift, phishift, rshiftcart, phishiftcart;
  
    for(int i = 0; i < nphi; i++){
      double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
      for(int j = 0; j < nr; j++){
	double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin

	double x = r*cos(phi); //cm
	double y = r*sin(phi);

	for(int k = 0; k < nz; k++){
	  double z = minz + ((maxz - minz)/(1.0*nz))*(k+0.5); //center of bin

	  xshift=(hCartesianAveShift[0]->Interpolate(x,y))*(1e-4);//coordinate of your stripe
	  yshift=(hCartesianAveShift[1]->Interpolate(x,y))*(1e-4);//convert micron to cm
	  zshift=(hCartesianAveShift[2]->Interpolate(x,y))*(1e-4);
 
	  rshift=(hCylindricalAveShift[0]->Interpolate(x,y))*(1e-4);
	  phishift=(hCylindricalAveShift[1]->Interpolate(x,y));

	  //rshift calculated from xshift n yshift
	  rshiftcart=(hCylindricalAveShift[2]->Interpolate(x,y))*(1e-4);
	  phishiftcart=hCylindricalAveShift[3]->Interpolate(x,y);

	  hCartesianCMModel[0]->Fill(phi,r,z,xshift*(1-z/105.5));
	  hCartesianCMModel[1]->Fill(phi,r,z,yshift*(1-z/105.5));
	  hCartesianCMModel[2]->Fill(phi,r,z,zshift*(1-z/105.5));
	
	  hCylindricalCMModel[0]->Fill(phi,r,z,rshift*(1-z/105.5));
	  hCylindricalCMModel[1]->Fill(phi,r,z,phishift*(1-z/105.5));

	  //radial and phi models from cartesian models 0,1
	  hCylindricalCMModel[2]->Fill(phi,r,z,rshiftcart*(1-z/105.5));
	  hCylindricalCMModel[3]->Fill(phi,r,z,phishiftcart*(1-z/105.5));
	
	}
      }
    }

    int ndiff = 300;
    int mindiff = -20;
    int maxdiff = 20;
  
    TH1F *hCartesianShiftDifference[3];
    hCartesianShiftDifference[0] = new TH1F("hShiftDifferenceX", "Difference between CM Model X and True (R > 30); #Delta X (#mum)", ndiff, mindiff, maxdiff);
    hCartesianShiftDifference[1] = new TH1F("hShiftDifferenceY", "Difference between CM Model Y and True (R > 30); #Delta Y (#mum)", ndiff, mindiff, maxdiff);
    hCartesianShiftDifference[2] = new TH1F("hShiftDifferenceZ", "Difference between CM Model Z and True (R > 30); #Delta Z (#mum)", ndiff, mindiff, maxdiff);
  
    TH1F *hCylindricalShiftDifference[4];
    hCylindricalShiftDifference[0] = new TH1F("hShiftDifferenceR", "Difference between CM Model R and True (R > 30); #Delta R (#mum)", ndiff, mindiff, maxdiff);
    hCylindricalShiftDifference[1] = new TH1F("hShiftDifferencePhi", "Difference between CM Model Phi and True (R > 30); #Delta Phi (#mum)", ndiff, mindiff, maxdiff);
    hCylindricalShiftDifference[2] = new TH1F("hShiftDifferenceRCart", "Difference between CM Model R from Cartesian and True (R > 30); #Delta R (#mum)", ndiff, mindiff, maxdiff);
    hCylindricalShiftDifference[3] = new TH1F("hShiftDifferencePhiCart", "Difference between CM Model Phi from Cartesian and True (R > 30); #Delta Phi (#mum)", ndiff, mindiff, maxdiff);

    TH1F *hRShiftDifference = new TH1F("hRShiftDifference", "Difference between CM Model R from Cartesian and CM Model R from R data (R > 30); #Delta R (#mum)", ndiff, mindiff, maxdiff);
    TH1F *hPhiShiftDifference = new TH1F("hPhiShiftDifference", "Difference between CM Model Phi from Cartesian and CM Model Phi from Phi data (R > 30); #Delta Phi (#mum)", ndiff, mindiff, maxdiff);

    TH1F *hRShiftTrue = new TH1F("hRShiftTrue", "True R Distortion Model (R > 30); #Delta R (#mum)", ndiff, mindiff, maxdiff);
    TH1F *hPhiShiftTrue = new TH1F("hPhiShiftTrue", "True Phi Distortion Model (R > 30); #Delta Phi (#mum)", ndiff, mindiff, maxdiff);
  
    TH2F *hCartesianDiff[6];
    hCartesianDiff[0] = new TH2F("hDiffXYX", "Difference in XY for CM Model X; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianDiff[1] = new TH2F("hDiffRZX", "Difference in RZ for CM Model X; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCartesianDiff[2] = new TH2F("hDiffXYY", "Difference in XY for CM Model Y; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianDiff[3] = new TH2F("hDiffRZY", "Difference in RZ for CM Model Y; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCartesianDiff[4] = new TH2F("hDiffXYZ", "Difference in XY for CM Model Z; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianDiff[5] = new TH2F("hDiffRZZ", "Difference in RZ for CM Model Z; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  
    TH2F *hCylindricalDiff[8];
    hCylindricalDiff[0] = new TH2F("hDiffXYR", "Difference in XY for CM Model R; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalDiff[1] = new TH2F("hDiffRZR", "Difference in RZ for CM Model R; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCylindricalDiff[2] = new TH2F("hDiffXYPhi", "Difference in XY for CM Model Phi; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalDiff[3] = new TH2F("hDiffRZPhi", "Difference in RZ for CM Model Phi; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);  
    hCylindricalDiff[4] = new TH2F("hDiffXYRCart", "Difference in XY for CM Model R from Cartesian; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalDiff[5] = new TH2F("hDiffRZRCart", "Difference in RZ for CM Model R from Cartesian; z (cm); r (cm)",nz,minz,maxz,nr,minr,maxr);
    hCylindricalDiff[6] = new TH2F("hDiffXYPhiCart", "Difference in XY for CM Model Phi from Cartesian; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalDiff[7] = new TH2F("hDiffRZPhiCart", "Difference in RZ for CM Model Phi from Cartesian; z (cm); r (cm)",nz,minz,maxz,nr,minr,maxr);

    TH2F *hRDiff[2];
    hRDiff[0] = new TH2F("hRDiffXY", "Difference between R Models in XY; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hRDiff[1] = new TH2F("hRDiffRZ", "Difference between R Models in RZ; z (cm); r (cm)",nz,minz,maxz,nr,minr,maxr);
    TH2F *hPhiDiff[2];
    hPhiDiff[0] = new TH2F("hPhiDiffXY", "Difference between Phi Models in XY; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hPhiDiff[1] = new TH2F("hPhiDiffRZ", "Difference between Phi Models in RZ; z (cm); r (cm)",nz,minz,maxz,nr,minr,maxr);
  
    TH2F *hCartesianAveDiff[6];
    hCartesianAveDiff[0] = new TH2F("hAveDiffXYX", "X Model - Truth Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianAveDiff[1] = new TH2F("hAveDiffRZX", "X Model - Truth Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCartesianAveDiff[2] = new TH2F("hAveDiffXYY", "Y Model - Truth Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianAveDiff[3] = new TH2F("hAveDiffRZY", "Y Model - Truth Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCartesianAveDiff[4] = new TH2F("hAveDiffXYZ", "Z Model - Truth Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianAveDiff[5] = new TH2F("hAveDiffRZZ", "Z Model - Truth Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  
    TH2F *hCylindricalAveDiff[8];
    hCylindricalAveDiff[0] = new TH2F("hAveDiffXYR", "R Model - Truth Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalAveDiff[1] = new TH2F("hAveDiffRZR", "R Model - Truth Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCylindricalAveDiff[2] = new TH2F("hAveDiffXYPHi", "Phi Model - Truth Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalAveDiff[3] = new TH2F("hAveDiffRZPhi", "Phi Model - Truth Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCylindricalAveDiff[4] = new TH2F("hAveDiffXYRCart", "R Model from Cartesian - Truth Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalAveDiff[5] = new TH2F("hAveDiffRZRCart", "R Model from Cartesian - Truth Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCylindricalAveDiff[6] = new TH2F("hAveDiffXYPhiCart", "Phi Model from Cartesian - Truth Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalAveDiff[7] = new TH2F("hAveDiffRZPhiCart", "Phi Model from Cartesian - Truth Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);

    TH2F *hRAveDiff[2];
    hRAveDiff[0] = new TH2F("hRAveDiffXY", "R Model from Cartesian - Original R Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hRAveDiff[1] = new TH2F("hRAveDiffRZ", "R Model from Cartesian - Original R Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);

    TH2F *hPhiAveDiff[2];
    hPhiAveDiff[0] = new TH2F("hPhiAveDiffXY", "Phi Model from Cartesian - Original Phi Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hPhiAveDiff[1] = new TH2F("hPhiAveDiffRZ", "Phi Model from Cartesian - Original Phi Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  
    TH2F *hSamplePerBinXY = new TH2F("hSamplePerBinXY", "Filling each xy bin; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    TH2F *hSamplePerBinRZ = new TH2F("hSamplePerBinRZ", "Filling each rz bin; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);

    TH2F *hCompareXY = new TH2F("hCompareXY", "Compare Difference in X and Y Models; x diff (#mum); y diff (#mum)",nbins,low,high,nbins,low,high);

    TH2F *hCompareRTrue = new TH2F("hCompareRTrue", "Compare Difference from R Model and True (R > 30, 10 < z < 90); reco shift (#mum); true shift (#mum)",nbins,-550,550,nbins,-550,550);
    TH2F *hComparePhiTrue = new TH2F("hComparePhiTrue", "Compare Difference from Phi Model and True (R > 30, 10 < z < 90); reco shift (#mum); true shift (#mum)",nbins,-550,550,nbins,-550,550);

    TH2F *hRDiffvR = new TH2F("hRDiffvR", "Difference between R Model and True vs. r (R > 30, 10 < z < 90); r (cm); shift difference (#mum)",nr,minr,maxr,ndiff,mindiff,maxdiff);
    TH2F *hRDiffvZ = new TH2F("hRDiffvZ", "Difference between R Model and True vs. z (R > 30); z (cm); shift difference (#mum)",nz,minz,maxz,ndiff,mindiff,maxdiff);
    TH2F *hRDiffvPhi = new TH2F("hRDiffvPhi", "Difference between R Model and True vs. phi (R > 30, 10 < z < 90); phi (rad); shift difference (#mum)",nphi,minphi,maxphi,ndiff,mindiff,maxdiff);

    TH2F *hPhiDiffvR = new TH2F("hPhiDiffvR", "Difference between Phi Model and True vs. r (R > 30, 10 < z < 90); r (cm); shift difference (#mum)",nr,minr,maxr,ndiff,mindiff,maxdiff);
    TH2F *hPhiDiffvZ = new TH2F("hPhiDiffvZ", "Difference between Phi Model and True vs. z (R > 30); z (cm); shift difference (#mum)",nz,minz,maxz,ndiff,mindiff,maxdiff);
    TH2F *hPhiDiffvPhi = new TH2F("hPhiDiffvPhi", "Difference between Phi Model and True vs. phi (R > 30, 10 < z < 90); phi (rad); shift difference (#mum)",nphi,minphi,maxphi,ndiff,mindiff,maxdiff);

    TH2F *hCMmodelSliceRvTrue = new TH2F("hCMmodelSliceRvTrue", "Difference between R Model and True as a function of R and Phi; r (cm); phi (rad)",nr,minr,maxr,nphi,minphi,maxphi);

    for(int i = 1; i < nphi - 1; i++){
      double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
      for(int j = 1; j < nr - 1; j++){
	double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin
	for(int k = 1; k < nz - 1; k++){
	  double z = minz + ((maxz - minz)/(1.0*nz))*(k+0.5); //center of bin

	  double shiftrecoCart[3];
	  double shifttrueCart[3];
	  double differenceCart[3];
	
	  double shiftrecoCyl[4];
	  double shifttrueCyl[4];
	  double differenceCyl[4];

	  double differenceR, differencePhi;

	  int bin = hCartesianCMModel[0]->FindBin(phi,r,z); //same for all

	  if((r > 30.0) && (r < 76.0)){
	  //if ((z > 20) && (z < 90)){
	    shifttrueCart[0] = (shifter->hX->Interpolate(phi,r,z))*(1e4); //convert from cm to micron
	    shifttrueCart[1] = (shifter->hY->Interpolate(phi,r,z))*(1e4); //convert from cm to micron 
	    shifttrueCart[2] = (shifter->hZ->Interpolate(phi,r,z))*(1e4); //convert from cm to micron 

	    //x y and z
	    for(int l = 0; l < 3; l ++){
	      shiftrecoCart[l] =  (hCartesianCMModel[l]->GetBinContent(bin))*(1e4);
	  
	      differenceCart[l] = shiftrecoCart[l] - shifttrueCart[l]; 

	      hCartesianShiftDifference[l]->Fill(differenceCart[l]);
	    }

	    //r
	    for(int l = 0; l < 3; l = l + 2){  
	      shiftrecoCyl[l] =  (hCylindricalCMModel[l]->GetBinContent(bin))*(1e4);
	      shifttrueCyl[l] = (shifter->hR->Interpolate(phi,r,z))*(1e4); //convert from cm to micron 
	      differenceCyl[l] = shiftrecoCyl[l] - shifttrueCyl[l]; 

	      //if ((z > 20) && (z < 90)){
		hCylindricalShiftDifference[l]->Fill(differenceCyl[l]);
		//}
	    }

	    //phi 
	    for(int l = 1; l < 4; l = l + 2){  
	      shiftrecoCyl[l] = r*(1e4)*(hCylindricalCMModel[l]->GetBinContent(bin));
	      shifttrueCyl[l] = (shifter->hPhi->Interpolate(phi,r,z))*(1e4); 
	      differenceCyl[l] = (shiftrecoCyl[l] - shifttrueCyl[l]); 

	      //if ((z > 20) && (z < 90)){
		hCylindricalShiftDifference[l]->Fill(differenceCyl[l]);
		//}
	    }
	
	    differenceR = differenceCyl[2]-differenceCyl[0];
	    hRShiftDifference->Fill(differenceR);

	    differencePhi = differenceCyl[3]-differenceCyl[1];
	    hPhiShiftDifference->Fill(differencePhi);

	    hRShiftTrue->Fill(shifttrueCyl[2]);
	    hPhiShiftTrue->Fill(shifttrueCyl[3]);

	
	    if (k == 1){
	      hCMmodelSliceRvTrue->Fill(r,phi,differenceCyl[2]);
	    }
	
	    double x = r*cos(phi);
	    double y = r*sin(phi);
       	
	    hCompareXY->Fill(differenceCart[0],differenceCart[1],1); 

	    //x
	    hCartesianDiff[0]->Fill(x,y, differenceCart[0]);
	    hCartesianDiff[1]->Fill(z,r, differenceCart[0]);
	    //y
	    hCartesianDiff[2]->Fill(x,y, differenceCart[1]);	  
	    hCartesianDiff[3]->Fill(z,r, differenceCart[1]);
	    //z
	    hCartesianDiff[4]->Fill(x,y, differenceCart[2]);
	    hCartesianDiff[5]->Fill(z,r, differenceCart[2]);

	    //r
	    hCylindricalDiff[0]->Fill(x,y, differenceCyl[0]);
	    hCylindricalDiff[1]->Fill(z,r, differenceCyl[0]);
	    //phi
	    hCylindricalDiff[2]->Fill(x,y, differenceCyl[1]);
	    hCylindricalDiff[3]->Fill(z,r, differenceCyl[1]);

	    //r cart
	    hCylindricalDiff[4]->Fill(x,y, differenceCyl[2]);
	    hCylindricalDiff[5]->Fill(z,r, differenceCyl[2]);
	    //phi cart
	    hCylindricalDiff[6]->Fill(x,y, differenceCyl[3]);
	    hCylindricalDiff[7]->Fill(z,r, differenceCyl[3]);
	  
	    //compare r and phi models from cartesian to originals
	    hRDiff[0]->Fill(x,y,differenceR);
	    hPhiDiff[0]->Fill(x,y,differencePhi);

	    //if ((z > 20) && (z < 90)){
	      hRDiff[1]->Fill(z,r,differenceR);
		hPhiDiff[1]->Fill(z,r,differencePhi);
		hSamplePerBinRZ->Fill(z,r,1);
		// }
	    //exclude ends
	    //if ((z > 10) && (z < 90)){
	      hCompareRTrue->Fill(shiftrecoCyl[2],shifttrueCyl[2]);
	      hComparePhiTrue->Fill(shiftrecoCyl[3],shifttrueCyl[3]);

	      hRDiffvR->Fill(r,differenceCyl[2],1);
	      hRDiffvPhi->Fill(phi,differenceCyl[2],1);

	      hPhiDiffvR->Fill(r,differenceCyl[3],1);
	      hPhiDiffvPhi->Fill(phi,differenceCyl[3],1);
	      //}
	    
	    hRDiffvZ->Fill(z,differenceCyl[2],1);
	    
	    hPhiDiffvZ->Fill(z,differenceCyl[3],1);
	    
	    hSamplePerBinXY->Fill(x,y,1);
	    // }
	  }
	}
      }
    }

    for (int m = 0; m < 6; m = m+2){
      hCartesianAveDiff[m]->Divide(hCartesianDiff[m],hSamplePerBinXY);
      hCylindricalAveDiff[m]->Divide(hCylindricalDiff[m],hSamplePerBinXY);
    }
    for (int m = 1; m < 6; m = m+2){
      hCartesianAveDiff[m]->Divide(hCartesianDiff[m],hSamplePerBinRZ);
      hCylindricalAveDiff[m]->Divide(hCylindricalDiff[m],hSamplePerBinRZ);
    }

    hCylindricalAveDiff[6]->Divide(hCylindricalDiff[6],hSamplePerBinXY);
    hCylindricalAveDiff[7]->Divide(hCylindricalDiff[7],hSamplePerBinRZ);
  
    hRAveDiff[0]->Divide(hRDiff[0],hSamplePerBinXY);
    hRAveDiff[1]->Divide(hRDiff[1],hSamplePerBinRZ);

    hPhiAveDiff[0]->Divide(hPhiDiff[0],hSamplePerBinXY);
    hPhiAveDiff[1]->Divide(hPhiDiff[1],hSamplePerBinRZ);

    //summary plots
    hDifferenceMeanR->Fill(hCylindricalShiftDifference[2]->GetMean(1));
    hDifferenceStdDevR->Fill(hCylindricalShiftDifference[2]->GetStdDev(1));

    hTrueMeanR->Fill(hRShiftTrue->GetMean(1));
    hTrueStdDevR->Fill(hRShiftTrue->GetStdDev(1));
    
    hDifferenceMeanPhi->Fill(hCylindricalShiftDifference[3]->GetMean(1));
    hDifferenceStdDevPhi->Fill(hCylindricalShiftDifference[3]->GetStdDev(1));

    hTrueMeanPhi->Fill(hPhiShiftTrue->GetMean(1));
    hTrueStdDevPhi->Fill(hPhiShiftTrue->GetStdDev(1));
   
    //  TFile *plots;

  
    // plots=TFile::Open("shift_plots.root","RECREATE");
    /* hXplots.Write();
       hYplots.Write();
       hZplots.Write();
       hRplots.Write();
       hPhiplots.Write(); */

    /* for(int i = 0; i < 3; i++){
      hCartesianCMModel[i]->Write();
    }
    hCylindricalCMModel[2]->Write();
    hCylindricalCMModel[3]->Write();
    plots->Close();
    */
  
    for (int i = 0; i < 3; i++){
      hCartesianForward[i]->SetStats(0);
      hCartesianAveShift[i]->SetStats(0);
      hCylindricalForward[i]->SetStats(0);
      hCylindricalAveShift[i]->SetStats(0);
    }

    hStripesPerBin->SetStats(0);

    for (int m = 0; m < 6; m = m+2){
      hCartesianAveDiff[m]->SetStats(0);
      hCylindricalAveDiff[m]->SetStats(0);
    }
    for (int m = 1; m < 6; m = m+2){
      hCartesianAveDiff[m]->SetStats(0);
      hCylindricalAveDiff[m]->SetStats(0);
    }

    hRAveDiff[0]->SetStats(0);
    hRAveDiff[1]->SetStats(0);

    hCylindricalAveShift[3]->SetStats(0);
    hCylindricalAveDiff[6]->SetStats(0);
    hCylindricalAveDiff[7]->SetStats(0);
    hPhiAveDiff[0]->SetStats(0);
    hPhiAveDiff[1]->SetStats(0);
    hCompareXY->SetStats(0);
    hCompareRTrue->SetStats(0);
    hComparePhiTrue->SetStats(0);

    hRDiffvR->SetStats(0);
    hRDiffvZ->SetStats(0);
    hRDiffvPhi->SetStats(0);
  
    hPhiDiffvR->SetStats(0);
    hPhiDiffvZ->SetStats(0);
    hPhiDiffvPhi->SetStats(0);

    /*double tsize = 0.07;
    double stsize = 0.03;
    double plotsize;

    plotsize = 1 - tsize - (6.0*stsize); */
    
    
    TPad *c1=new TPad("c1","",0.0,0.8,1.0,0.93); //can i do an array of pads?
    TPad *c2=new TPad("c2","",0.0,0.64,1.0,0.77);
    TPad *c3=new TPad("c3","",0.0,0.48,1.0,0.61);
    TPad *c4=new TPad("c4","",0.0,0.32,1.0,0.45);
    TPad *c5=new TPad("c5","",0.0,0.16,1.0,0.29);
    TPad *c6=new TPad("c6","",0.0,0.0,1.0,0.13);
    
    TPad *titlepad=new TPad("titlepad","",0.0,0.96,1.0,1.0);

    TPad *stitlepad1=new TPad("stitlepad1","",0.0,0.93,1.0,0.96);
    TPad *stitlepad2=new TPad("stitlepad2","",0.0,0.77,1.0,0.8);
    TPad *stitlepad3=new TPad("stitlepad3","",0.0,0.61,1.0,0.64);
    TPad *stitlepad4=new TPad("stitlepad4","",0.0,0.45,1.0,0.48);
    TPad *stitlepad5=new TPad("stitlepad5","",0.0,0.29,1.0,0.32);
    TPad *stitlepad6=new TPad("stitlepad6","",0.0,0.13,1.0,0.16);
    
    TLatex * title = new TLatex(0.0,0.0,"");

    TLatex * stitle1 = new TLatex(0.0,0.0,""); //array?
    TLatex * stitle2 = new TLatex(0.0,0.0,"");
    TLatex * stitle3 = new TLatex(0.0,0.0,"");
    TLatex * stitle4 = new TLatex(0.0,0.0,"");
    TLatex * stitle5 = new TLatex(0.0,0.0,"");
    TLatex * stitle6 = new TLatex(0.0,0.0,"");
    
    title->SetNDC();
    stitle1->SetNDC();
    stitle2->SetNDC();
    stitle3->SetNDC();
    stitle4->SetNDC();
    stitle5->SetNDC();
    stitle6->SetNDC();
    
    title->SetTextSize(0.32);
    stitle1->SetTextSize(0.35);
    stitle2->SetTextSize(0.35);
    stitle3->SetTextSize(0.35);
    stitle4->SetTextSize(0.35);
    stitle5->SetTextSize(0.35);
    stitle6->SetTextSize(0.35);
    
    canvas->cd();
    c1->Draw();
    stitlepad1->Draw();
    c2->Draw();
    stitlepad2->Draw();
    c3->Draw();
    stitlepad3->Draw();
    c4->Draw();
    stitlepad4->Draw();
    c5->Draw();
    stitlepad5->Draw();
    c6->Draw();
    stitlepad6->Draw();
    titlepad->Draw();

    //x plots
    c1->Divide(4,1);
    c1->cd(1);
    hCartesianAveDiff[0]->Draw("colz");
    c1->cd(2);
    hCartesianAveDiff[1]->Draw("colz");
    c1->cd(3);
    hCartesianShiftDifference[0]->Draw();
    //c1->cd(4)->Clear();  
    c1->cd(4);
    //hCMmodelSliceRvTrue->Draw("colz");
    hSamplePerBinRZ->Draw("colz");
    //y plots
    c2->Divide(4,1);
    c2->cd(1);
    hCartesianAveDiff[2]->Draw("colz");
    c2->cd(2);
    hCartesianAveDiff[3]->Draw("colz");
    c2->cd(3);
    hCartesianShiftDifference[1]->Draw();
    //c2->cd(4)->Clear();
    c2->cd(4);
    //hStripesPerBin->Draw("colz");
    hSamplePerBinXY->Draw("colz");
    
    //r cart
    c3->Divide(4,1);
    c3->cd(1);
    hCylindricalAveDiff[4]->Draw("colz");
    c3->cd(2);
    hCylindricalAveDiff[5]->Draw("colz");
    c3->cd(3);
    hCylindricalShiftDifference[2]->Draw();
    c3->cd(4);
    hRShiftTrue->Draw();
    
    //phi cart
    c4->Divide(4,1);
    c4->cd(1);
    hCylindricalAveDiff[6]->Draw("colz");
    c4->cd(2);
    hCylindricalAveDiff[7]->Draw("colz");
    c4->cd(3);
    hCylindricalShiftDifference[3]->Draw();
    c4->cd(4);
    hPhiShiftTrue->Draw();

    //r to true comparison
    c5->Divide(4,1);
    c5->cd(1);
    hCompareRTrue->Draw("colz");
    c5->cd(2);
    hRDiffvR->Draw("colz");
    c5->cd(3);
    hRDiffvZ->Draw("colz");
    c5->cd(4);
    hRDiffvPhi->Draw("colz");

    //phi to true comparison
    c6->Divide(4,1);
    c6->cd(1);
    hComparePhiTrue->Draw("colz");
    c6->cd(2);
    hPhiDiffvR->Draw("colz");
    c6->cd(3);
    hPhiDiffvZ->Draw("colz");
    c6->cd(4);
    hPhiDiffvPhi->Draw("colz");

    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.01,0.4,Form("Event %d; %s", ifile, sourcefilename.Data())); 
    title->Draw();
    
    stitlepad1->cd();
    stitlepad1->Clear();
    stitle1->DrawLatex(0.45,0.2,"X Model"); 
    stitle1->Draw();
     
    stitlepad2->cd();
    stitlepad2->Clear();
    stitle2->DrawLatex(0.45,0.2,"Y Model"); 
    stitle2->Draw();

    stitlepad3->cd();
    stitlepad3->Clear();
    stitle3->DrawLatex(0.45,0.2,"R Model"); 
    stitle3->Draw();

    stitlepad4->cd();
    stitlepad4->Clear();
    stitle4->DrawLatex(0.45,0.2,"Phi Model"); 
    stitle4->Draw();

    stitlepad5->cd();
    stitlepad5->Clear();
    stitle5->DrawLatex(0.4,0.2,"Comparing R Model to True"); 
    stitle5->Draw();

    stitlepad6->cd();
    stitlepad6->Clear();
    stitle6->DrawLatex(0.4,0.2,"Comparing Phi Model to True"); 
    stitle6->Draw();

    if(ifile == 0){ 
      //if(ifile == 1){
      canvas->Print("CMDistortionReco.pdf(","pdf");
    }
    else{
      canvas->Print("CMDistortionReco.pdf","pdf");
    }
  }
  
  TCanvas *summary = new TCanvas("summary","ShiftPlotsSummary",2000,3000);
  
  TPad *sumtitlepad = new TPad("sumtitlepad","",0.0,0.96,1.0,1.0);
  TPad *sumplots = new TPad("sumplotspad","",0.0,0.0,1.0,0.96);
  
  TLatex *sumtitle = new TLatex(0.0,0.0,"");
  
  sumtitle->SetNDC();
  sumtitle->SetTextSize(0.4);
  
  summary->cd();
  sumplots->Draw();
  sumtitlepad->Draw();

  sumplots->Divide(4,6);
  sumplots->cd(1);
  hDifferenceMeanR->Draw();
  sumplots->cd(2);
  hDifferenceStdDevR->Draw();
  sumplots->cd(3);
  hTrueMeanR->Draw();
  sumplots->cd(4);
  hTrueStdDevR->Draw();
  sumplots->cd(5);
  hDifferenceMeanPhi->Draw();
  sumplots->cd(6);
  hDifferenceStdDevPhi->Draw();
  sumplots->cd(7);
  hTrueMeanPhi->Draw();
  sumplots->cd(8);
  hTrueStdDevPhi->Draw();
  sumplots->cd(9);
  sumplots->cd(10)->Clear();
  sumplots->cd(11)->Clear();
  sumplots->cd(12)->Clear();
  sumplots->cd(13)->Clear();
  sumplots->cd(14)->Clear();
  sumplots->cd(15)->Clear();
  sumplots->cd(16)->Clear();
  sumplots->cd(17)->Clear();
  sumplots->cd(18)->Clear();
  sumplots->cd(19)->Clear();
  sumplots->cd(20)->Clear();
  sumplots->cd(21)->Clear();
  sumplots->cd(22)->Clear();
  sumplots->cd(23)->Clear();
  sumplots->cd(24)->Clear();

  sumtitlepad->cd();
  sumtitlepad->Clear();
  sumtitle->DrawLatex(0.4,0.4,"Summary of Events"); 
  summary->Print("CMDistortionReco.pdf)","pdf");
  
  return 0;
}
