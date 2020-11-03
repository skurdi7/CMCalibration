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
  forward=TFile::Open("/sphenix/user/rcorliss/distortion_maps/res_scan/Summary_bX1508071_0_10_events.root.h_Charge_evt_0.real_B1.5_E-400.0.ross_phi1_sphenix_phislice_lookup_r23xp23xz35.distortion_map.hist.root","READ"); //using temporary histogram for testing
  
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
  double deltaR;
  
  nbins = 40;
  /*rsteps = 100;
  phisteps = 100;
  
  rstepsize = (stripes.end_CM - stripes.begin_CM)/rsteps;
  phistepsize = 2*TMath::Pi()/phisteps; */

  //ScanHist(nbins, low, high, x, y);
  //IDLabels();
   
  TH2F *RShift = new TH2F("RShift","Radial shift of stripe centers; x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see

  TH2F *hStripesPerBin = new TH2F("hStripesPerBin","Stripes Per Bin; x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see

  TH2F *AveShift = new TH2F("AveShift","Divide RShift by ShiftCheck; x (cm); y (cm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see
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
	  
    cout << i << endl;
  }

  //AveShift->Divide(RShift,hStripesPerBin);
  //hPhiCheck2d->Divide(hStripesPerBin);


  //repeat for forward only
  TH2F *hForward = new TH2F("hForward","Radial Shift Forward of Stripe Centers",nbins,low,high,nbins,low,high); 
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
    
    deltaR = newposition.Perp() - position.Perp();
    hForward->Fill(x,y,deltaR);
  
  }

  AveShift->Divide(hForward,hStripesPerBin);
  hPhiCheck2d->Divide(hStripesPerBin);
  	  
  int nphi = shifter.hX->GetXaxis()->GetNbins();
  int nr = shifter.hX->GetYaxis()->GetNbins();
  int nz = shifter.hX->GetZaxis()->GetNbins();
  
  double minphi = shifter.hX->GetXaxis()->GetMinimum();
  double minr = shifter.hX->GetYaxis()->GetMinimum();
  double minz = shifter.hX->GetZaxis()->GetMinimum();
  
  double maxphi = shifter.hX->GetXaxis()->GetMaximum();
  double maxr = shifter.hX->GetYaxis()->GetMaximum();
  double maxz = shifter.hX->GetZaxis()->GetMaximum();

  TH3F *hCMModel = new TH3F("hCMModel", "Radial Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);

  double rshift;
  
  for(int i = 0; i < nphi; i++){
    double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
    for(int j = 0; j < nr; j++){
      double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin

      double x = r*cos(phi);
      double y = r*sin(phi);
      
      for(int k = 0; k < nz; k++){
	double z = minz + ((maxz - minz)/(1.0*nz))*(k+0.5); //center of bin

	
	rshift=AveShift->Interpolate(x,y,z);//coordinate of your stripe
	
	hCMModel->Fill(phi,r,z,rshift*z/105.5);
      }
    }
  }

  TH1F *hShiftDifference = new TH1F("hShiftDifference", "Difference between Radial Shift Reco and True", 300, -0.5, 0.5);

  for(int i = 0; i < nphi; i++){
    double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
    for(int j = 0; j < nr; j++){
      double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin
      for(int k = 0; k < nz; k++){
	double z = minz + ((maxz - minz)/(1.0*nz))*(k+0.5); //center of bin

	double difference = hCMModel->Interpolate(phi,r,z) - shifter.hR->Interpolate(phi,r,z);
	hShiftDifference->Fill(difference);
      }
    }
  }
  
  
  TCanvas *c=new TCanvas("c","RShift",1500,1000);
  c->Divide(3,2);
  c->cd(1);
  hForward->Draw("colz");
  c->cd(2);
  hStripesPerBin->Draw("colz");
  c->cd(3);
  AveShift->Draw("colz");
  c->cd(4);
  PhiCheck->Draw();
  c->cd(5);
  hPhiCheck2d->Draw("colz");
  c->cd(6);
  hShiftDifference->Draw();
  
  c->SaveAs("RShift.pdf");
  
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
  
 
