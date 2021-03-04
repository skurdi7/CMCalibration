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
  TVector3 Shift(TVector3 position);
  TVector3 ShiftForward(TVector3 position); //only shift with forward histogram
  TVector3 ShiftBack(TVector3 position); //
  TFile *forward, *back, *average;
  TH3F *hX, *hY, *hZ, *hR, *hPhi, *hXave, *hYave, *hZave, *hRave, *hPhiave, *hXBack, *hYBack, *hZBack;  
};

Shifter::Shifter(TString sourcefilename){
  //forward=TFile::Open("/sphenix/user/rcorliss/distortion_maps/res_scan/Summary_bX1508071_0_10_events.root.h_Charge_evt_0.real_B1.5_E-400.0.ross_phi1_sphenix_phislice_lookup_r23xp23xz35.distortion_map.hist.root","READ"); //using temporary histogram for testing, try running
  //forward=TFile::Open("/gpfs/mnt/gpfs02/sphenix/user/rcorliss/distortion_maps/elevatorpitch/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ"); //distortions due to the average amount of space charge in the TPC with effect of distortions due to the external electric and magnetic fields removed
  //forward=TFile::Open("/gpfs/mnt/gpfs02/sphenix/user/rcorliss/distortion_maps/elevatorpitch/fluct_single.1side.3d.file0.h_Charge_0.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ"); //distortions due to single-event differences from the average
  forward=TFile::Open(sourcefilename,"READ"); //single event distortion

  hX=(TH3F*)forward->Get("hIntDistortionX");
  hY=(TH3F*)forward->Get("hIntDistortionY");
  hZ=(TH3F*)forward->Get("hIntDistortionZ");

  hR=(TH3F*)forward->Get("hIntDistortionR");
  hPhi=(TH3F*)forward->Get("hIntDistortionP");

  average=TFile::Open("/gpfs/mnt/gpfs02/sphenix/user/rcorliss/distortion_maps/averages/empty.2sides.3d.file0.h_Charge_0.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ"); 
  hXave=(TH3F*)average->Get("hIntDistortionX");
  hYave=(TH3F*)average->Get("hIntDistortionY");
  hZave=(TH3F*)average->Get("hIntDistortionZ");

  hRave=(TH3F*)average->Get("hIntDistortionR");
  hPhiave=(TH3F*)average->Get("hIntDistortionP");

  hX->Add(hXave,-1);
  hY->Add(hYave,-1);
  hZ->Add(hZave,-1);

  hR->Add(hRave,-1);
  hPhi->Add(hPhiave,-1);
  
  back=TFile::Open("/sphenix/user/rcorliss/distortion_maps/averages/empty.2sides.3d.file0.h_Charge_0.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ"); //tells it only to read, not to write anything you make there.
   
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

  TVector3 forwardshift(x+xshift,y+yshift,z+zshift);

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
    
  shiftposition.SetXYZ(x+xshiftback,y+yshiftback,z+zshiftback);

  return shiftposition;
}

TVector3 Shifter::Shift(TVector3 position){
  
  return ShiftBack(ShiftForward(position));
}

void ScanHist(int nbins, double low, double high, double x, double y);
void IDLabels();

int cmShiftPlots() {
  Shifter *shifter;
  PHG4TpcCentralMembrane stripes;
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
  TCanvas *canvas=new TCanvas("canvas","ShiftPlotsAllEvents",2500,2000);
  
  const char * inputpattern="/gpfs/mnt/gpfs02/sphenix/user/rcorliss/distortion_maps/Oct20/full_maps/*.root";
  //find all files that match the input string (includes wildcards)
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  //TFile *infile;
  TString sourcefilename;
  //int nEvents = filelist->GetNFiles();
  int nEvents = 10;

  for (int ifile=0;ifile < nEvents;ifile++){
    //for each file, find all histograms in that file.
    sourcefilename=((TFileInfo*)(filelist->GetList()->At(ifile)))->GetCurrentUrl()->GetFile();//gross
    //infile=TFile::Open(sourcefilename.Data(),"READ");

    shifter = new Shifter(sourcefilename);
    
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
    
      newposition = shifter->Shift(position);

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
    hCartesianForward[0] = new TH2F("hForwardX","X Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianForward[1] = new TH2F("hForwardY","Y Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianForward[2] = new TH2F("hForwardZ","Z Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);

    TH2F *hCylindricalForward[2];
    hCylindricalForward[0] = new TH2F("hForwardR","Radial Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalForward[1] = new TH2F("hForwardPhi","Phi Shift Forward of Stripe Centers (rad); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    //hCylindricalForward[2] = new TH2F("hForwardRCart","R Shift Forward of Stripe Centers from Cartesian (z in cm); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  
    for (int i = 0; i < Hits.size(); i++){
      x = (Hits[i]->get_x(0) + Hits[i]->get_x(1))/2; //stripe center
      y = (Hits[i]->get_y(0) + Hits[i]->get_y(1))/2;
      z = 0.5;
    
      position.SetXYZ(x,y,z);
    
      double phi=position.Phi();
      if(position.Phi() < 0.0){
	phi = position.Phi() + TMath::TwoPi(); 
      }
  
      PhiCheck->Fill(phi);
      hPhiCheck2d->Fill(x,y,phi);
      
      newposition = shifter->ShiftForward(position);

      deltaX = (newposition.X() - position.X())*(1e4); //convert from cm to micron 
      deltaY = (newposition.Y() - position.Y())*(1e4);
      deltaZ = (newposition.Z() - position.Z())*(1e4);

      deltaR = (newposition.Perp() - position.Perp())*(1e4);
      deltaPhi = newposition.DeltaPhi(position);
      //double newR = sqrt(newposition.X()*newposition.X() + newposition.Y()*newposition.Y());
    
      // deltaRCart = (sqrt(newposition.X()*newposition.X() + newposition.Y()*newposition.Y()) - sqrt(position.X()*position.X() + position.Y()*position.Y()))*(1e4);

      hCartesianForward[0]->Fill(x,y,deltaX);
      hCartesianForward[1]->Fill(x,y,deltaY);
      hCartesianForward[2]->Fill(x,y,deltaZ);

      hCylindricalForward[0]->Fill(x,y,deltaR);
      hCylindricalForward[1]->Fill(x,y,deltaPhi);
      // hCylindricalForward[2]->Fill(x,y,deltaRCart);
      //hForwardR->Fill(x,y,deltaR);
  
    }

    TH2F *hCartesianAveShift[3];
    hCartesianAveShift[0] = new TH2F("AveShiftX","Average of CM Model X over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see
    hCartesianAveShift[1] = new TH2F("AveShiftY","Average of CM Model Y over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCartesianAveShift[2] = new TH2F("AveShiftZ","Average of CM Model Z over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 

    TH2F *hCylindricalAveShift[4];
    hCylindricalAveShift[0] = new TH2F("AveShiftR","Average of CM Model R over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see
    hCylindricalAveShift[1] = new TH2F("AveShiftPhi","Average of CM Model Phi over Stripes per Bin (rad); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCylindricalAveShift[2] = new TH2F("AveShiftRCart","Average of CM Model R over Stripes per Bin from Cartesian (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCylindricalAveShift[3] = new TH2F("AveShiftPhiCart","Average of CM Model Phi over Stripes per Bin from Cartesian (rad); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
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
      
	int xbin = hCartesianAveShift[0]->FindBin(x,y);
	int ybin = hCartesianAveShift[1]->FindBin(x,y);
	double xaveshift = (hCartesianAveShift[0]->GetBinContent(xbin))*(1e-4); // converts  microns to cm 
	double yaveshift = (hCartesianAveShift[1]->GetBinContent(ybin))*(1e-4);
	
	TVector3 shifted, original;
	original.SetX(x);
	original.SetY(y);
	shifted.SetX(x+xaveshift);
	shifted.SetY(y+yaveshift);
	// have x n y above for orig
	//shifted is orig + ave shift
	
	double raveshift = (shifted.Perp() - original.Perp())*(1e4);
	//double phiaveshift = shifted.Phi() - original.Phi();
	double phiaveshift = shifted.DeltaPhi(original);

	//if(phiaveshift > TMath::TwoPi()){
	// phiaveshift = phiaveshift - TMath::TwoPi(); 
	//}
	
	//fill with r from x n y
	// double raveshift = sqrt(xaveshift*xaveshift + yaveshift*yaveshift);
	hCylindricalAveShift[2]->Fill(x,y,raveshift);
	hCylindricalAveShift[3]->Fill(x,y,phiaveshift);
      } 
    } 
  
    hPhiCheck2d->Divide(hStripesPerBin);

    //same range and bins for each coordinate, can use hR for all, binned in cm
    int nphi = shifter->hR->GetXaxis()->GetNbins();
    int nr = shifter->hR->GetYaxis()->GetNbins();
    int nz = shifter->hR->GetZaxis()->GetNbins();
  
    double minphi = shifter->hR->GetXaxis()->GetXmin();
    double minr = shifter->hR->GetYaxis()->GetXmin();
    double minz = shifter->hR->GetZaxis()->GetXmin();

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
    //TH3F *hCMModel = new TH3F("hCMModel", "CM Model: Radial Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);

    double xshift, yshift, zshift, rshift, phishift, rshiftcart, phishiftcart;
  
    for(int i = 0; i < nphi; i++){
      double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
      for(int j = 0; j < nr; j++){
	double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin

	double x = r*cos(phi); //cm
	double y = r*sin(phi);

	//cout << "x" << x << endl;
	//cout << "y" << y << endl;
      
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

    TH2F *hCompareRTrue = new TH2F("hCompareRTrue", "Compare Difference from R Model and True (R > 30); reco shift (#mum); true shift (#mum)",nbins,low,high,nbins,low,high);
    TH2F *hComparePhiTrue = new TH2F("hComparePhiTrue", "Compare Difference from Phi Model and True (R > 30); reco shift (#mum); true shift (#mum)",nbins,low,high,nbins,low,high);

    TH2F *hRDiffvR = new TH2F("hRDiffvR", "Difference between R Model and True vs. r (R > 30); r (cm); shift difference (#mum)",nr,minr,maxr,ndiff,mindiff,maxdiff);
    TH2F *hRDiffvZ = new TH2F("hRDiffvZ", "Difference between R Model and True vs. z (R > 30); z (cm); shift difference (#mum)",nz,minz,maxz,ndiff,mindiff,maxdiff);
    TH2F *hRDiffvPhi = new TH2F("hRDiffvPhi", "Difference between R Model and True vs. phi (R > 30); phi (rad); shift difference (#mum)",nphi,minphi,maxphi,ndiff,mindiff,maxdiff);

    TH2F *hPhiDiffvR = new TH2F("hPhiDiffvR", "Difference between Phi Model and True vs. r (R > 30); r (cm); shift difference (#mum)",nr,minr,maxr,ndiff,mindiff,maxdiff);
    TH2F *hPhiDiffvZ = new TH2F("hPhiDiffvZ", "Difference between Phi Model and True vs. z (R > 30); z (cm); shift difference (#mum)",nz,minz,maxz,ndiff,mindiff,maxdiff);
    TH2F *hPhiDiffvPhi = new TH2F("hPhiDiffvPhi", "Difference between Phi Model and True vs. phi (R > 30); phi (rad); shift difference (#mum)",nphi,minphi,maxphi,ndiff,mindiff,maxdiff);

    //TH1F *hCMmodelslicePhi = 
  
    for(int i = 0; i < nphi; i++){
      double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
      for(int j = 0; j < nr; j++){
	double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin
	for(int k = 0; k < nz; k++){
	  double z = minz + ((maxz - minz)/(1.0*nz))*(k+0.5); //center of bin

	  double shiftrecoCart[3];
	  double shifttrueCart[3];
	  double differenceCart[3];
	
	  double shiftrecoCyl[4];
	  double shifttrueCyl[4];
	  double differenceCyl[4];

	  double differenceR, differencePhi;
	
	  /*cout << "phi: " << phi;
	    cout << "r: " << r;
	    cout << "z: " << z << endl;
	  */

	  int bin = shifter->hR->FindBin(phi,r,z);

	  if(r > 30.0){
	  
	    shifttrueCart[0] = (shifter->hX->GetBinContent(bin))*(1e4); //convert from cm to micron
	    shifttrueCart[1] = (shifter->hY->GetBinContent(bin))*(1e4); //convert from cm to micron 
	    shifttrueCart[2] = (shifter->hZ->GetBinContent(bin))*(1e4); //convert from cm to micron 

	    for(int l = 0; l < 3; l ++){
	      shiftrecoCart[l] =  (hCartesianCMModel[l]->GetBinContent(bin))*(1e4);
	  
	      differenceCart[l] = shiftrecoCart[l] - shifttrueCart[l]; 

	      hCartesianShiftDifference[l]->Fill(differenceCart[l]);
	    }
	
	    for(int l = 0; l < 3; l = l + 2){  
	      shiftrecoCyl[l] =  (hCylindricalCMModel[l]->GetBinContent(bin))*(1e4);
	      shifttrueCyl[l] = (shifter->hR->GetBinContent(bin))*(1e4); //convert from cm to micron 
	      differenceCyl[l] = shiftrecoCyl[l] - shifttrueCyl[l]; 
	    
	      hCylindricalShiftDifference[l]->Fill(differenceCyl[l]);
	    }
	  
	    for(int l = 1; l < 4; l = l + 2){  
	      shiftrecoCyl[l] = r*(1e4)*(hCylindricalCMModel[l]->GetBinContent(bin));
	      shifttrueCyl[l] = (shifter->hPhi->GetBinContent(bin))*(1e4); 
	      differenceCyl[l] = (shiftrecoCyl[l] - shifttrueCyl[l]); 

	      hCylindricalShiftDifference[l]->Fill(differenceCyl[l]);
	    }
	
	    differenceR = differenceCyl[2]-differenceCyl[0];
	    hRShiftDifference->Fill(differenceR);

	    differencePhi = differenceCyl[3]-differenceCyl[1];
	    hPhiShiftDifference->Fill(differencePhi);

	    hRShiftTrue->Fill(shifttrueCyl[2]);
	    hPhiShiftTrue->Fill(shifttrueCyl[3]);

	
	    if (k == nz/2){
	      //hCMmodelslicePhi->Fill(shiftrecoCyl[3]);
	      //hTrueslicePhi->Fill(shifttrueCyl[3]);
	      //
	    }
	
	    double x = r*cos(phi);
	    double y = r*sin(phi);
       	
	    hCompareXY->Fill(differenceCart[0],differenceCart[1],1); 

	    hCompareRTrue->Fill(shiftrecoCyl[2],shifttrueCyl[2]);
	    hComparePhiTrue->Fill(shiftrecoCyl[3],shifttrueCyl[3]);

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
	    hRDiff[1]->Fill(z,r,differenceR);

	    hPhiDiff[0]->Fill(x,y,differencePhi);
	    hPhiDiff[1]->Fill(z,r,differencePhi);

	    hRDiffvR->Fill(r,differenceCyl[2],1);
	    hRDiffvZ->Fill(z,differenceCyl[2],1);
	    hRDiffvPhi->Fill(phi,differenceCyl[2],1);
	  
	    hPhiDiffvR->Fill(r,differenceCyl[3],1);
	    hPhiDiffvZ->Fill(z,differenceCyl[3],1);
	    hPhiDiffvPhi->Fill(phi,differenceCyl[3],1);
	  
	    hSamplePerBinXY->Fill(x,y,1);
	    hSamplePerBinRZ->Fill(z,r,1);
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

    
    TPad *c=new TPad("c","",0.0,0.0,1.0,0.9);
    TPad *titlepad=new TPad("titlepad","",0.0,0.9,1.0,1.0);
    TLatex * title = new TLatex(0.0,0.0,"");
    title->SetNDC();
    title->SetTextSize(0.5);
    canvas->cd();
    c->Draw();
    titlepad->Draw();

    c->Divide(5,4);
    //x plots
    c->cd(1);
    hCartesianAveDiff[0]->Draw("colz");
    c->cd(2);
    hCartesianAveDiff[1]->Draw("colz");
    c->cd(3);
    hCartesianShiftDifference[0]->Draw();
    //y plots
    c->cd(4);
    hCartesianAveDiff[2]->Draw("colz");
    c->cd(5);
    hCartesianAveDiff[3]->Draw("colz");
    c->cd(6);
    hCartesianShiftDifference[1]->Draw();
    //r cart
    c->cd(7);
    hCylindricalAveDiff[4]->Draw("colz");
    c->cd(8);
    hCylindricalAveDiff[5]->Draw("colz");
    c->cd(9);
    hCylindricalShiftDifference[2]->Draw();
    //phi cart
    c->cd(10);
    hCylindricalAveDiff[6]->Draw("colz");
    c->cd(11);
    hCylindricalAveDiff[7]->Draw("colz");
    c->cd(12);
    hCylindricalShiftDifference[3]->Draw();
    //r to true comparison
    c->cd(13);
    hCompareRTrue->Draw("colz");
    c->cd(14);
    hRDiffvR->Draw("colz");
    c->cd(15);
    hRDiffvZ->Draw("colz");
    c->cd(16);
    hRDiffvPhi->Draw("colz");
    //phi to true comparison
    c->cd(17);
    hComparePhiTrue->Draw("colz");
    c->cd(18);
    hPhiDiffvR->Draw("colz");
    c->cd(19);
    hPhiDiffvZ->Draw("colz");
    c->cd(20);
    hPhiDiffvPhi->Draw("colz");
    
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,Form("Event %d", ifile)); //how do i change the number
    title->Draw();
    if(ifile == 0){
      canvas->Print("ShiftPlotsAllEvents.pdf(","pdf");
    }
    else if (ifile == nEvents - 1){
      canvas->Print("ShiftPlotsAllEvents.pdf)","pdf"); 
    }
    else{
      canvas->Print("ShiftPlotsAllEvents.pdf","pdf");
    }

    canvas->Print(Form("ShiftPlotsEvent%d.gif", ifile),"gif");
    
    /*   TCanvas *canvas=new TCanvas("canvas","ShiftPlots",1500,1000);
    TPad *c=new TPad("c","",0.0,0.0,1.0,0.9);
    TPad *titlepad=new TPad("titlepad","",0.0,0.9,1.0,1.0);
    TLatex * title = new TLatex(0.0,0.0,"");
    title->SetNDC();
    title->SetTextSize(0.5);
    canvas->cd();
    c->Draw();
    titlepad->Draw();
  
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
    //c->cd();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"X Shift Model");
    title->Draw();
    canvas->Print("ShiftPlots.pdf(","pdf");
  
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
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"Y Shift Model");
    canvas->Print("ShiftPlots.pdf","pdf");
  
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
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"Z Shift Model");
    canvas->Print("ShiftPlots.pdf","pdf");
  
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
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"R Shift Model");
    canvas->Print("ShiftPlots.pdf","pdf");

    // r plots from cart
    //c->Divide(3,2);
    c->cd(1)->Clear();
    //hCylindricalForward[2]->Draw("colz");
    c->cd(2)->Clear();
    //hStripesPerBin->Draw("colz");
    c->cd(3);
    hCylindricalAveShift[2]->Draw("colz");
    c->cd(4);
    hCylindricalAveDiff[4]->Draw("colz");
    c->cd(5);
    hCylindricalAveDiff[5]->Draw("colz");
    c->cd(6);
    hCylindricalShiftDifference[2]->Draw();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.3,0.2,"R Shift from Cartesian Model");
    canvas->Print("ShiftPlots.pdf","pdf");

    // compare the two R models
    c->cd(1);
    hCylindricalAveShift[0]->Draw("colz");
    c->cd(2);
    hCylindricalAveShift[2]->Draw("colz");
    c->cd(3)->Clear();
    c->cd(4);
    hRAveDiff[0]->Draw("colz");
    c->cd(5);
    hRAveDiff[1]->Draw("colz");
    c->cd(6);
    hRShiftDifference->Draw();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"Comparing R Models");
    canvas->Print("ShiftPlots.pdf","pdf");

    //compare R cart to true
    c->cd(1);
    hCylindricalShiftDifference[2]->Draw();
    c->cd(2);
    hRShiftTrue->Draw();
    c->cd(3);
    hCompareRTrue->Draw("colz");
    c->cd(4);
    hRDiffvR->Draw("colz");
    c->cd(5);
    hRDiffvZ->Draw("colz");
    c->cd(6);
    hRDiffvPhi->Draw("colz");
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.3,0.2,"Comparing R from Cart Model to True");
    canvas->Print("ShiftPlots.pdf","pdf");
  
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
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"Phi Shift Model");
    canvas->Print("ShiftPlots.pdf","pdf");

    // phi plots from cart
    c->cd(1)->Clear();
    c->cd(2)->Clear();
    c->cd(3);
    hCylindricalAveShift[3]->Draw("colz");
    c->cd(4);
    hCylindricalAveDiff[6]->Draw("colz");
    c->cd(5);
    hCylindricalAveDiff[7]->Draw("colz");
    c->cd(6);
    hCylindricalShiftDifference[3]->Draw();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.3,0.2,"Phi Shift from Cartesian Model");
    canvas->Print("ShiftPlots.pdf","pdf");

    // compare the two Phi models
    c->cd(1);
    hCylindricalAveShift[1]->Draw("colz");
    c->cd(2);
    hCylindricalAveShift[3]->Draw("colz");
    c->cd(3)->Clear();
    c->cd(4);
    hPhiAveDiff[0]->Draw("colz");
    c->cd(5);
    hPhiAveDiff[1]->Draw("colz");
    c->cd(6);
    hPhiShiftDifference->Draw();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"Comparing Phi Models");
    canvas->Print("ShiftPlots.pdf","pdf");
  
    //compare Phi cart to true
    c->cd(1);
    hCylindricalShiftDifference[3]->Draw();
    c->cd(2);
    hPhiShiftTrue->Draw();
    c->cd(3);
    hComparePhiTrue->Draw("colz");
    c->cd(4);
    hPhiDiffvR->Draw("colz");
    c->cd(5);
    hPhiDiffvZ->Draw("colz");
    c->cd(6);
    hPhiDiffvPhi->Draw("colz");
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.3,0.2,"Comparing Phi from Cart Model to True");
    canvas->Print("ShiftPlots.pdf","pdf");
  
    //Compare X and Y diff
    c->cd(1);
    hCompareXY->Draw("colz");
    c->cd(2)->Clear();
    c->cd(3)->Clear();
    c->cd(4)->Clear();
    c->cd(5)->Clear();
    c->cd(6)->Clear();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"Comparing X and Y Models");
    canvas->Print("ShiftPlots.pdf)","pdf"); */
  
    // c->SaveAs("RShift.pdf"); // replace w print
  }
  return 0;
}

void ScanHist(int nbins, double low, double high, double x, double y){
  PHG4TpcCentralMembrane stripes;
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
  PHG4TpcCentralMembrane stripes;
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
  
 
