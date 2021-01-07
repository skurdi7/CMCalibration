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


// all distances in mm, all angles in rad
// class that generates stripes and dummy hit coordinates
// stripes have width of one mm, length of one pad width, and are centered in middle of sector gaps

using namespace std;

class Shifter {
public:
  Shifter();
  TVector3 Shift(TVector3 position);
  TVector3 ShiftForward(TVector3 position); //only shift with forward histogram
  TVector3 ShiftBack(TVector3 position); //
  TFile *forward, *back;
  TH3F *hX, *hY, *hZ, *hR, *hXBack, *hYBack, *hZBack;  
};

Shifter::Shifter(){
  //forward=TFile::Open("/sphenix/user/rcorliss/distortion_maps/res_scan/Summary_bX1508071_0_10_events.root.h_Charge_evt_0.real_B1.5_E-400.0.ross_phi1_sphenix_phislice_lookup_r23xp23xz35.distortion_map.hist.root","READ"); //using temporary histogram for testing
  forward=TFile::Open("/gpfs/mnt/gpfs02/sphenix/user/rcorliss/distortion_maps/elevatorpitch/fluct_single.1side.3d.file0.h_Charge_0.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ"); 
  
  back=TFile::Open("/sphenix/user/rcorliss/distortion_maps/averages/empty.2sides.3d.file0.h_Charge_0.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ"); //tells it only to read, not to write anything you make there.

  hX=(TH3F*)forward->Get("hIntDistortionX");
  hY=(TH3F*)forward->Get("hIntDistortionY");
  hZ=(TH3F*)forward->Get("hIntDistortionZ");

  hR=(TH3F*)forward->Get("hIntDistortionR");
   
  hXBack=(TH3F*)back->Get("hIntDistortionX");
  hYBack=(TH3F*)back->Get("hIntDistortionY");
  hZBack=(TH3F*)back->Get("hIntDistortionZ");
}

TVector3 Shifter::ShiftForward(TVector3 position){
double x, y, z, xshift, yshift, zshift;
  const double mm = 1.0;
  const double cm = 10.0;
  TVector3 shiftposition;

  x= position.X();
  y= position.Y();
  z= position.Z();

  double r=position.Perp();
  double phi=position.Phi();
  if(position.Phi() < 0.0){
    phi = position.Phi() + 2.0*TMath::Pi(); 
  }
  
  xshift=hX->Interpolate(phi,r,z);//coordinate of your stripe
  yshift=hY->Interpolate(phi,r,z);
  zshift=hZ->Interpolate(phi,r,z);

  TVector3 forwardshift(x+xshift,y+yshift,z);

  return forwardshift;
}

TVector3 Shifter::ShiftBack(TVector3 forwardshift){
double x, y, z, xshift, yshift, zshift;
  const double mm = 1.0;
  const double cm = 10.0;
  TVector3 shiftposition;

  x= forwardshift.X();
  y= forwardshift.Y();
  z= forwardshift.Z();

  double rforward=forwardshift.Perp();
  double phiforward=forwardshift.Phi();
  if(forwardshift.Phi() < 0.0){
    phiforward += 2.0*TMath::Pi();
  }
  
  double xshiftback=-1*hXBack->Interpolate(phiforward,rforward,z);
  double yshiftback=-1*hYBack->Interpolate(phiforward,rforward,z);
  double zshiftback=-1*hZBack->Interpolate(phiforward,rforward,z);
    
  shiftposition.SetXYZ(x+xshiftback,y+yshiftback,z);

  return shiftposition;
}

TVector3 Shifter::Shift(TVector3 position){
  
  return ShiftBack(ShiftForward(position));
}

void ScanHist(int nbins, double low, double high, double x, double y);
void IDLabels();

int cmShiftPlots() {
  Shifter shifter;
  StripesClass stripes;
  vector<PHG4Hitv1*> Hits = stripes.PHG4Hits;
  int nbins; 
  double x, y, z;
  TVector3 position, newposition;
  double low = -80.0;
  double high = 80.0;
  double deltaX, deltaY, deltaZ, deltaR, deltaPhi;
  
  nbins = 40;
  /*rsteps = 100;
  phisteps = 100;
  
  rstepsize = (stripes.end_CM - stripes.begin_CM)/rsteps;
  phistepsize = 2*TMath::Pi()/phisteps; */

  //ScanHist(nbins, low, high, x, y);
  //IDLabels();
   
  TH2F *RShift = new TH2F("RShift","Radial shift of stripe centers (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see

  TH2F *hStripesPerBin = new TH2F("hStripesPerBin","CM Stripes Per Bin (z in stripes); x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see

  //TH2F *AveShift = new TH2F("AveShift","Average of CM Model over Stripes per Bin; x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see
  TH2F *hPhiCheck2d = new TH2F("hPhiCheck2d","what phi am i using; x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see

  TH1F *PhiCheck = new TH1F("PhiCheck","what phi am i using; phi (radians)",200,-10.0,10.0);
  
  for (int i = 0; i < Hits.size(); i++){ 
    x = (Hits[i]->get_x(0) + Hits[i]->get_x(1))/2; //stripe center
    y = (Hits[i]->get_y(0) + Hits[i]->get_y(1))/2;
    z = 0.5;
    
    position.SetXYZ(x,y,z);

    double phi=position.Phi();
    if(position.Phi() < 0.0){
      phi = position.Phi() + 2.0*TMath::Pi(); 
    }
  
    PhiCheck->Fill(phi);
    hPhiCheck2d->Fill(x,y,phi);
    
    newposition = shifter.Shift(position);

    deltaR = newposition.Perp() - position.Perp();
    RShift->Fill(x,y,deltaR);
    hStripesPerBin->Fill(x,y,1);
	  
    // cout << i << endl;
  }

  //AveShift->Divide(RShift,hStripesPerBin);
  //hPhiCheck2d->Divide(hStripesPerBin);


  //repeat for forward only
  //TH2F *hForwardR = new TH2F("hForwardR","Radial Shift Forward of Stripe Centers; x (cm); y (cm)",nbins,low,high,nbins,low,high);

  TH2F *hCartesianForward[3];
  hCartesianForward[0] = new TH2F("hForwardX","X Shift Forward of Stripe Centers (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCartesianForward[1] = new TH2F("hForwardY","Y Shift Forward of Stripe Centers (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCartesianForward[2] = new TH2F("hForwardZ","Z Shift Forward of Stripe Centers (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);

  TH2F *hCylindricalForward[2];
  hCylindricalForward[0] = new TH2F("hForwardR","Radial Shift Forward of Stripe Centers (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCylindricalForward[1] = new TH2F("hForwardPhi","Phi Shift Forward of Stripe Centers (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  //hCylindricalForward[2] = new TH2F("hForwardRCart","R Shift Forward of Stripe Centers from Cartesian; x (cm); y (cm)",nbins,low,high,nbins,low,high);
  
  for (int i = 0; i < Hits.size(); i++){
    x = (Hits[i]->get_x(0) + Hits[i]->get_x(1))/2; //stripe center
    y = (Hits[i]->get_y(0) + Hits[i]->get_y(1))/2;
    z = 0.5;
    
    position.SetXYZ(x,y,z);
    
    double phi=position.Phi();
    if(position.Phi() < 0.0){
      phi = position.Phi() + 2.0*TMath::Pi(); 
    }
  
    PhiCheck->Fill(phi);
    hPhiCheck2d->Fill(x,y,phi);
      
    newposition = shifter.ShiftForward(position);

    deltaX = (newposition.X() - position.X())*(1e4);
    deltaY = (newposition.Y() - position.Y())*(1e4);
    deltaZ = (newposition.Z() - position.Z())*(1e4);

    deltaR = (newposition.Perp() - position.Perp())*(1e4);
    deltaPhi = newposition.Phi() - position.Phi();
    //deltaRCart = sqrt(deltaX*deltaX + deltaY*deltaY); // wrong calculation

    hCartesianForward[0]->Fill(x,y,deltaX);
    hCartesianForward[1]->Fill(x,y,deltaY);
    hCartesianForward[2]->Fill(x,y,deltaZ);

    hCylindricalForward[0]->Fill(x,y,deltaR);
    hCylindricalForward[1]->Fill(x,y,deltaPhi);
    //hCylindricalForward[2]->Fill(x,y,deltaRCart);
    //hForwardR->Fill(x,y,deltaR);
  
  }

  TH2F *hCartesianAveShift[3];
  hCartesianAveShift[0] = new TH2F("AveShiftX","Average of CM Model X over Stripes per Bin (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see
  hCartesianAveShift[1] = new TH2F("AveShiftY","Average of CM Model Y over Stripes per Bin (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see
  hCartesianAveShift[2] = new TH2F("AveShiftZ","Average of CM Model Z over Stripes per Bin (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see

  TH2F *hCylindricalAveShift[3];
  hCylindricalAveShift[0] = new TH2F("AveShiftR","Average of CM Model R over Stripes per Bin (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see
  hCylindricalAveShift[1] = new TH2F("AveShiftPhi","Average of CM Model Phi over Stripes per Bin (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see
  hCylindricalAveShift[2] = new TH2F("AveShiftRCart","Average of CM Model R over Stripes per Bin from Cartesian (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see
  
  //AveShift->Divide(hForwardR,hStripesPerBin);

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
      // try interpolate
      double xaveshift = hCartesianAveShift[0]->Interpolate(x,y);
      double yaveshift = hCartesianAveShift[1]->Interpolate(x,y);
      //fill with r from x n y
      double raveshift = sqrt(xaveshift*xaveshift + yaveshift*yaveshift);
      hCylindricalAveShift[2]->Fill(x,y,raveshift);
    }
  }
  hPhiCheck2d->Divide(hStripesPerBin);

  //same range and bins for each coordinate, can use hR for all
  int nphi = shifter.hR->GetXaxis()->GetNbins();
  int nr = shifter.hR->GetYaxis()->GetNbins();
  int nz = shifter.hR->GetZaxis()->GetNbins();
  
  double minphi = shifter.hR->GetXaxis()->GetXmin();
  double minr = shifter.hR->GetYaxis()->GetXmin();
  double minz = shifter.hR->GetZaxis()->GetXmin();
  
  double maxphi = shifter.hR->GetXaxis()->GetXmax();
  double maxr = shifter.hR->GetYaxis()->GetXmax();
  double maxz = shifter.hR->GetZaxis()->GetXmax();


  TH3F *hCartesianCMModel[3];
  hCartesianCMModel[0]=new TH3F("hCMModelX", "CM Model: X Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
  hCartesianCMModel[1]=new TH3F("hCMModelY", "CM Model: Y Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
  hCartesianCMModel[2]=new TH3F("hCMModelZ", "CM Model: Z Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);

  TH3F *hCylindricalCMModel[3];
  hCylindricalCMModel[0]=new TH3F("hCMModelR", "CM Model: Radial Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
  hCylindricalCMModel[1]=new TH3F("hCMModelPhi", "CM Model: Phi Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
  hCylindricalCMModel[2]=new TH3F("hCMModelRCart", "CM Model: Radial Shift Forward of Stripe Centers from Cartesian", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);  
  //TH3F *hCMModel = new TH3F("hCMModel", "CM Model: Radial Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);

  double xshift, yshift, zshift, rshift, phishift, rshiftcart;
  
  for(int i = 0; i < nphi; i++){
    double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
    for(int j = 0; j < nr; j++){
      double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin

      double x = r*cos(phi);
      double y = r*sin(phi);

      //cout << "x" << x << endl;
      //cout << "y" << y << endl;
      
      for(int k = 0; k < nz; k++){
	double z = minz + ((maxz - minz)/(1.0*nz))*(k+0.5); //center of bin
	
	xshift=hCartesianAveShift[0]->Interpolate(x,y);//coordinate of your stripe
	yshift=hCartesianAveShift[1]->Interpolate(x,y);
	zshift=hCartesianAveShift[2]->Interpolate(x,y);
 
	rshift=hCylindricalAveShift[0]->Interpolate(x,y);
	phishift=hCylindricalAveShift[1]->Interpolate(x,y);

	//rshift calculated from xshift n yshift
	rshiftcart=sqrt(xshift*xshift + yshift*yshift);

	hCartesianCMModel[0]->Fill(phi,r,z,xshift*(1-z/105.5));
	hCartesianCMModel[1]->Fill(phi,r,z,yshift*(1-z/105.5));
	hCartesianCMModel[2]->Fill(phi,r,z,zshift*(1-z/105.5));
	
	hCylindricalCMModel[0]->Fill(phi,r,z,rshift*(1-z/105.5));
	hCylindricalCMModel[1]->Fill(phi,r,z,phishift*(1-z/105.5));

	//radial model from cartesian models 0,1
	hCylindricalCMModel[2]->Fill(phi,r,z,rshiftcart*(1-z/105.5));
	
      }
    }
  }

  TH1F *hCartesianShiftDifference[3];
  hCartesianShiftDifference[0] = new TH1F("hShiftDifferenceX", "Difference between CM Model X and True (R > 30); #Delta X (#mu m)", 300, -200, 200);
  hCartesianShiftDifference[1] = new TH1F("hShiftDifferenceY", "Difference between CM Model Y and True (R > 30); #Delta Y (#mu m)", 300, -200, 200);
  hCartesianShiftDifference[2] = new TH1F("hShiftDifferenceZ", "Difference between CM Model Z and True (R > 30); #Delta Z (#mu m)", 300, -200, 200);
  
  TH1F *hCylindricalShiftDifference[3];
  hCylindricalShiftDifference[0] = new TH1F("hShiftDifferenceR", "Difference between CM Model R and True (R > 30); #Delta R (#mu m)", 300, -200, 200);
  hCylindricalShiftDifference[1] = new TH1F("hShiftDifferencePhi", "Difference between CM Model Phi and True (R > 30); #Delta Phi (#mu m)", 300, -200, 200);
  hCylindricalShiftDifference[2] = new TH1F("hShiftDifferenceRCart", "Difference between CM Model R from Cartesian and True (R > 30); #Delta R (#mu m)", 300, -200, 200);

  TH1F *hRShiftDifference = new TH1F("hRShiftDifference", "Difference between CM Model R from Cartesian and CM Model R from R data (R > 30); #Delta R (#mu m)", 300, -200, 200);
  
TH2F *hCartesianDiff[6];
  hCartesianDiff[0] = new TH2F("hDiffXYX", "Difference in XY for CM Model X; x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCartesianDiff[1] = new TH2F("hDiffRZX", "Difference in RZ for CM Model X; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  hCartesianDiff[2] = new TH2F("hDiffXYY", "Difference in XY for CM Model Y; x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCartesianDiff[3] = new TH2F("hDiffRZY", "Difference in RZ for CM Model Y; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  hCartesianDiff[4] = new TH2F("hDiffXYZ", "Difference in XY for CM Model Z; x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCartesianDiff[5] = new TH2F("hDiffRZZ", "Difference in RZ for CM Model Z; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  
  TH2F *hCylindricalDiff[6];
  hCylindricalDiff[0] = new TH2F("hDiffXYR", "Difference in XY for CM Model R; x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCylindricalDiff[1] = new TH2F("hDiffRZR", "Difference in RZ for CM Model R; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  hCylindricalDiff[2] = new TH2F("hDiffXYPhi", "Difference in XY for CM Model Phi; x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCylindricalDiff[3] = new TH2F("hDiffRZPhi", "Difference in RZ for CM Model Phi; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);  
  hCylindricalDiff[4] = new TH2F("hDiffXYRCart", "Difference in XY for CM Model R from Cartesian; x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCylindricalDiff[5] = new TH2F("hDiffRZRCart", "Difference in RZ for CM Model R from Cartesian; z (cm); r (cm)",nz,minz,maxz,nr,minr,maxr);

  TH2F *hRDiff[2];
  hRDiff[0] = new TH2F("hRDiffXY", "Difference between R Models in XY; x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hRDiff[1] = new TH2F("hRDiffRZ", "Difference between R Models in RZ; z (cm); r (cm)",nz,minz,maxz,nr,minr,maxr);
  
  TH2F *hCartesianAveDiff[6];
  hCartesianAveDiff[0] = new TH2F("hAveDiffXYX", "X Model - Truth Averaged Over z (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCartesianAveDiff[1] = new TH2F("hAveDiffRZX", "X Model - Truth Averaged Over phi (z in cm); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  hCartesianAveDiff[2] = new TH2F("hAveDiffXYY", "Y Model - Truth Averaged Over z (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCartesianAveDiff[3] = new TH2F("hAveDiffRZY", "Y Model - Truth Averaged Over phi (z in cm); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  hCartesianAveDiff[4] = new TH2F("hAveDiffXYZ", "Z Model - Truth Averaged Over z (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCartesianAveDiff[5] = new TH2F("hAveDiffRZZ", "Z Model - Truth Averaged Over phi (z in cm); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  
  TH2F *hCylindricalAveDiff[6];
  hCylindricalAveDiff[0] = new TH2F("hAveDiffXYR", "R Model - Truth Averaged Over z (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCylindricalAveDiff[1] = new TH2F("hAveDiffRZR", "R Model - Truth Averaged Over phi (z in cm); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  hCylindricalAveDiff[2] = new TH2F("hAveDiffXYPHi", "Phi Model - Truth Averaged Over z (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCylindricalAveDiff[3] = new TH2F("hAveDiffRZPhi", "Phi Model - Truth Averaged Over phi (z in cm); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  hCylindricalAveDiff[4] = new TH2F("hAveDiffXYRCart", "R Model from Cartesian - Truth Averaged Over z (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCylindricalAveDiff[5] = new TH2F("hAveDiffRZRCart", "R Model from Cartesian - Truth Averaged Over phi (z in cm); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);

  TH2F *hRAveDiff[2];
  hRAveDiff[0] = new TH2F("hRAveDiffXY", "R Model from Cartesian - Original R Averaged Over z (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hRAveDiff[1] = new TH2F("hRAveDiffRZ", "R Model from Cartesian - Original R Averaged Over phi (z in cm); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  
  TH2F *hSamplePerBinXY = new TH2F("hSamplePerBinXY", "Filling each xy bin; x (cm); y (cm)",nbins,low,high,nbins,low,high);
  TH2F *hSamplePerBinRZ = new TH2F("hSamplePerBinRZ", "Filling each rz bin; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);

  for(int i = 0; i < nphi; i++){
    double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
    for(int j = 0; j < nr; j++){
      double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin
      for(int k = 0; k < nz; k++){
	double z = minz + ((maxz - minz)/(1.0*nz))*(k+0.5); //center of bin

	double shiftrecoCart[3];
	double shifttrueCart[3];
	double differenceCart[3];
	
	double shiftrecoCyl[3];
	double shifttrueCyl[3];
	double differenceCyl[3];

	double differenceR;
	
	/*cout << "phi: " << phi;
	cout << "r: " << r;
	cout << "z: " << z << endl;
	*/

	int bin = shifter.hR->FindBin(phi,r,z);

	for(int l = 0; l < 3; l ++){
	  shiftrecoCart[l] =  hCartesianCMModel[l]->GetBinContent(bin);
	  shifttrueCart[l] = (shifter.hR->GetBinContent(bin))*(1e4);
	  differenceCart[l] = shiftrecoCart[l] - shifttrueCart[l]; 

	  hCartesianShiftDifference[l]->Fill(differenceCart[l]);
	  
	  shiftrecoCyl[l] =  hCylindricalCMModel[l]->GetBinContent(bin);
	  shifttrueCyl[l] = (shifter.hR->GetBinContent(bin))*(1e4);
	  differenceCyl[l] = shiftrecoCyl[l] - shifttrueCyl[l]; 

	  hCylindricalShiftDifference[l]->Fill(differenceCyl[l]);
	}

	differenceR = differenceCyl[2]-differenceCyl[0];
	hRShiftDifference->Fill(differenceR);
	  
	if (k == nz/2){
	  //cmmodelslice -> fill(shift reco)
	  //trueslice -> fill(shift true)
	}
	
	double x = r*cos(phi);
	double y = r*sin(phi);

	//if difference < -0.8
	if(r > 30.0){
	  for(int l = 0; l < 3; l ++){
	    for (int m = 0; m < 6; m = m+2){
	      hCartesianDiff[m]->Fill(x,y, differenceCart[l]);
	      hCylindricalDiff[m]->Fill(x,y, differenceCyl[l]);
	    }
	    for (int m = 1; m < 6; m = m+2){
	      hCartesianDiff[m]->Fill(z,r, differenceCart[l]);
	      hCylindricalDiff[m]->Fill(z,r, differenceCyl[l]);
	    }

	    hRDiff[0]->Fill(x,y,differenceR);
	    hRDiff[1]->Fill(z,r,differenceR);
	  }
	  
	  hSamplePerBinXY->Fill(x,y,1);
	  hSamplePerBinRZ->Fill(z,r,1);
	}

	  //also compare r model to r model from cartesian
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

  hRAveDiff[0]->Divide(hRDiff[0],hSamplePerBinXY);
  hRAveDiff[1]->Divide(hRDiff[1],hSamplePerBinRZ);
  
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
  
  // gStyle->SetOptStat(0);
  
  TCanvas *c=new TCanvas("c","ShiftPlots",1500,1000);
  // x plots
  c->Divide(3,2);
  c->cd(1);
  hCartesianForward[0]->Draw("colz");
  c->cd(2);
  hStripesPerBin->Draw("colz");
  c->cd(3);
  hCartesianAveShift[0]->Draw("colz");
  c->cd(4);
  hCartesianAveDiff[0]->Draw("colz");
  c->cd(5);
  hCartesianAveDiff[1]->Draw("colz");
  c->cd(6);
  hCartesianShiftDifference[0]->Draw();
  c->Print("ShiftPlots.pdf(","pdf");
  
  // y plots
  //c->Divide(3,2);
  c->cd(1);
  hCartesianForward[1]->Draw("colz");
  c->cd(2);
  hStripesPerBin->Draw("colz");
  c->cd(3);
  hCartesianAveShift[1]->Draw("colz");
  c->cd(4);
  hCartesianAveDiff[2]->Draw("colz");
  c->cd(5);
  hCartesianAveDiff[3]->Draw("colz");
  c->cd(6);
  hCartesianShiftDifference[1]->Draw();
  c->Print("ShiftPlots.pdf","pdf");
  
  // z plots
  // c->Divide(3,2);
  c->cd(1);
  hCartesianForward[2]->Draw("colz");
  c->cd(2);
  hStripesPerBin->Draw("colz");
  c->cd(3);
  hCartesianAveShift[2]->Draw("colz");
  c->cd(4);
  hCartesianAveDiff[4]->Draw("colz");
  c->cd(5);
  hCartesianAveDiff[5]->Draw("colz");
  c->cd(6);
  hCartesianShiftDifference[2]->Draw();
  c->Print("ShiftPlots.pdf","pdf");
  
  // r plots
  //c->Divide(3,2);
  c->cd(1);
  hCylindricalForward[0]->Draw("colz");
  c->cd(2);
  hStripesPerBin->Draw("colz");
  c->cd(3);
  hCylindricalAveShift[0]->Draw("colz");
  c->cd(4);
  //PhiCheck->Draw();
  hCylindricalAveDiff[0]->Draw("colz");
  c->cd(5);
  // hPhiCheck2d->Draw("colz");
  hCylindricalAveDiff[1]->Draw("colz");
  c->cd(6);
  hCylindricalShiftDifference[0]->Draw();
  c->Print("ShiftPlots.pdf","pdf");

  // r plots from cart
  //c->Divide(3,2);
  c->cd(1);
  
  c->cd(2);
  hStripesPerBin->Draw("colz");
  c->cd(3);
  hCylindricalAveShift[2]->Draw("colz");
  c->cd(4);
  hCylindricalAveDiff[4]->Draw("colz");
  c->cd(5);
  hCylindricalAveDiff[5]->Draw("colz");
  c->cd(6);
  hCylindricalShiftDifference[2]->Draw();
  c->Print("ShiftPlots.pdf","pdf");

  // compare the two R models
  c->cd(1);
  hCylindricalAveShift[0]->Draw("colz");
  c->cd(2);
  hCylindricalAveShift[2]->Draw("colz");
  c->cd(3);
  
  c->cd(4);
  hRAveDiff[0]->Draw("colz");
  c->cd(5);
  hRAveDiff[1]->Draw("colz");
  c->cd(6);
  hRShiftDifference->Draw();
  c->Print("ShiftPlots.pdf","pdf");
  
  // phi plots
  //c->Divide(3,2);
  c->cd(1);
  hCylindricalForward[1]->Draw("colz");
  c->cd(2);
  hStripesPerBin->Draw("colz");
  c->cd(3);
  hCylindricalAveShift[1]->Draw("colz");
  c->cd(4);
  hCylindricalAveDiff[2]->Draw("colz");
  c->cd(5);
  hCylindricalAveDiff[3]->Draw("colz");
  c->cd(6);
  hCylindricalShiftDifference[1]->Draw();
  c->Print("ShiftPlots.pdf)","pdf");
  
  // c->SaveAs("RShift.pdf"); // replace w print
  
  return 0;
}

void ScanHist(int nbins, double low, double high, double x, double y){
  StripesClass stripes;
  int stripeID; 
  //histogram from search
  TH2F *Pattern1 = new TH2F("Pattern1","X,Y Scan if in Stripe;X (mm);Y (mm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see
  
  //TLatex *tex=new TLatex(x,y,"Stripe");
  //tex->SetTextSize(0.005);
  // for (r = stripes.begin_CM; r < stripes.end_CM; r = r + rstepsize){ // radii spanning full CM
  // for (phi = 0.0; phi < 2*TMath::Pi(); phi = phi + phistepsize){ // angles spanning full CM
  for (int i=0;i<nbins;i++){
    x=low+(high-low)/nbins*(i+0.5);
    for (int j=0;j<nbins;j++){
      y=low+(high-low)/nbins*(j+0.5);
 
    // x = r*cos(phi);
    // y = r*sin(phi);
      // cout << x << endl;
      //cout << y << endl; 
      
      stripeID = stripes.getStripeID(x, y);
      /* TLatex *tex=new TLatex(x,y,"StripeID");
	 tex->SetTextSize(0.005);
	 tex->DrawLatex(x,y,Form("%d",stripeID));
      */
      
      if(stripeID == -1){
        Pattern1->Fill(x,y,0);
      } else{
        Pattern1->Fill(x,y,1);
      }
      
    }
  }

  TCanvas *c=new TCanvas("a","CheckStripeID.cpp",500,500); 
  Pattern1->Draw("colz"); //check if theres still white spaces
  c->SaveAs("cmScan.pdf");
  //
}

void IDLabels(){
  StripesClass stripes;
  int stripeID;
  vector<PHG4Hitv1*> Hits = stripes.PHG4Hits;
  const double mm = 1.0;
  const double cm = 10.0;
  
  vector<double> xhit;
  vector<double> yhit;
  
  //build tgraph from dummy hits
  for (int i = 0; i < Hits.size(); i++){
    xhit.push_back(Hits[i]->get_x(0)*cm/mm); 
    yhit.push_back(Hits[i]->get_y(0)*cm/mm);
    xhit.push_back(Hits[i]->get_x(1)*cm/mm);
    yhit.push_back(Hits[i]->get_y(1)*cm/mm);
  }
  
  int npts = 2*Hits.size();
  TGraph *gDummyHits = new TGraph(npts, &xhit[0], &yhit[0]);
  gDummyHits->SetMarkerColor(2);
  gDummyHits->SetMarkerSize(0.5);
  
  gStyle->SetOptStat(0);
  TCanvas *c=new TCanvas("a","CheckStripeID.cpp",500,500);
  gDummyHits->Draw("AP");

  //loop thru hits again
  double xav, yav, xa, xb, ya, yb;

  for (int i = 0; i < Hits.size(); i++){
    //avg x0 n x1, y0 n y1 and use to draw stripeID
    xav = (Hits[i]->get_x(0)*cm/mm + Hits[i]->get_x(1)*cm/mm)/2;
    yav = (Hits[i]->get_y(0)*cm/mm + Hits[i]->get_y(1)*cm/mm)/2;
    //cout << "i: " << i << endl;
    xa = Hits[i]->get_x(0)*cm/mm ;
    xb= Hits[i]->get_x(1)*cm/mm ;
    ya= Hits[i]->get_y(0)*cm/mm ;
    yb = Hits[i]->get_y(1)*cm/mm ;
    //cout << "xav: " << xav << endl;
    //cout << "yav: " << yav << endl; 
    stripeID = stripes.getStripeID(xav, yav);
    TLatex *tex=new TLatex(xav,yav,"StripeID");
    tex->SetTextSize(0.005);
    tex->DrawLatex(xav,yav,Form("%d",stripeID));
  }

  c->SaveAs("cmStripeID.pdf");
}  
  //loop thru hits again
  /*
  double xav, yav, xa, xb, ya, yb;
 for (int i = 0; i < Hits.size(); i++){
  
    //avg x0 n x1, y0 n y1 and use to draw stripeID
    xav = (Hits[i]->get_x(0)*cm/mm + Hits[i]->get_x(1)*cm/mm)/2;
    yav = (Hits[i]->get_y(0)*cm/mm + Hits[i]->get_y(1)*cm/mm)/2;
    //cout << "i: " << i << endl;
    xa = Hits[i]->get_x(0)*cm/mm ;
    xb= Hits[i]->get_x(1)*cm/mm ;
    ya= Hits[i]->get_y(0)*cm/mm ;
    yb = Hits[i]->get_y(1)*cm/mm ;
    //cout << "xav: " << xav << endl;
    //cout << "yav: " << yav << endl; 
    stripeID = stripes.getStripeID(xav, yav);
    TLatex *tex=new TLatex(xav,yav,"StripeID");
    tex->SetTextSize(0.005);
    tex->DrawLatex(xav,yav,Form("%d",stripeID));
    

    TLatex *tex=new TLatex(xav,yav,"Stripe");
    tex->SetTextSize(0.005);
    if(stripeID == -1){
      tex->DrawLatex(xav,yav,Form("%d",0));
    } else{
      tex->DrawLatex(xav,yav,Form("%d",1));
    }
}
  */
      
    //TLine *line=new TLine;
    //line->DrawLine(Hits[i]->get_x(0)*cm/mm ,Hits[i]->get_y(0)*cm/mm,Hits[i]->get_x(1)*cm/mm, Hits[i]->get_y(1)*cm/mm);
    //line->DrawLine(0, 0, 200, 200*((yb- ya)/(xb-xa)));
    
  
  
  /*for (r = stripes.begin_CM; r < stripes.end_CM; r = r + rstepsize){ // radii spanning full CM
    for (phi = 0.0; phi < 2.0*TMath::Pi()/9.0; phi = phi + phistepsize){ // angles spanning full CM
      
      x = r*cos(phi);
      y = r*sin(phi);
      
      
      
      //TLatex tex;
      
      //if(stripeID == 1)
      //Pattern1->Fill(x,y);
      cout << stripeID << endl;
    }
  }
  */
  
 
