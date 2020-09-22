#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TTree.h"
#include "/sphenix/u/skurdi/CMCalibration/PHG4TpcCentralMembrane.h"
R__LOAD_LIBRARY(build/.libs/libg4tpccentralmembrane)

// generate xml code for all stripes to put into eagle board 

int fullCMstripes() {
  StripesClass stripes;
  ofstream eagle;
  eagle.open ("eaglestripes.txt");
  
  vector<PHG4Hitv1*> BV = stripes.BotVertices;
  vector<PHG4Hitv1*> TV = stripes.TopVertices;
  double xbl, xbr, xtl, xtr;
  double ybl, ybr, ytl, ytr;
  //vector<double> xbotleft, xbotright, xtopleft, xtopright;
  //vector<double> ybotleft, ybotright, ytopleft, ytopright;

  for (int i = 0; i < BV.size(); i++){
    xbl = BV[i]->get_x(0);
    xbr = BV[i]->get_x(1);
    xtl = TV[i]->get_x(0);
    xtr = TV[i]->get_x(1);
    ybl = BV[i]->get_y(0);
    ybr = BV[i]->get_y(1);
    ytl = TV[i]->get_y(0);
    ytr = TV[i]->get_y(1);

  
    eagle << "wire x1=\"" << xbl << "\" y1=\"" << ybl << "\" x2=\""<< xtl << "\" y2=\"" << ytl << "\" width=\"0.1524\" layer=\"46\" curve=\"-180\"/> " << endl;
    eagle << "<wire x1=\"" << xtl << "\" y1=\"" << ytl << "\" x2=\""<< xtr << "\" y2=\"" << ytr << "\" width=\"0.1524\" layer=\"46\"/> " << endl;
    eagle << "<wire x1=\"" << xtr << "\" y1=\"" << ytr << "\" x2=\""<< xbr << "\" y2=\"" << ybr << "\" width=\"0.1524\" layer=\"46\" curve=\"-180\"/> " << endl;  
    eagle << "<wire x1=\"" << xbr << "\" y1=\"" << ybr << "\" x2=\""<< xbl << "\" y2=\"" << ybl << "\" width=\"0.1524\" layer=\"46\"/> " << endl;
    eagle << endl;
 
  }
  eagle.close();
  return 0;
}
