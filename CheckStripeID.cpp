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

void ScanHist(int nbins, double low, double high, double x, double y);
void IDLabels();

int CheckStripeID() {
  int nbins, rsteps, phisteps, result; 
  double r, phi, x, y, xmod, ymod, phimod, rstepsize, phistepsize;
  double low = 250.0;
  double high = 290.0;
  
  nbins = 1000;
  /*rsteps = 100;
  phisteps = 100;
  
  rstepsize = (stripes.end_CM - stripes.begin_CM)/rsteps;
  phistepsize = 2*TMath::Pi()/phisteps; */

  //ScanHist(nbins, low, high, x, y);
  IDLabels();
  
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
  Pattern1->Draw();
  c->SaveAs("cmScan.pdf");
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
  
 
