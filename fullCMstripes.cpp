#include <iostream>
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

  vector<PHG4Hitv1*> BotVertices = stripes.BotVertices;
  vector<PHG4Hitv1*> TopVertices = stripes.TopVertices;
  double xbl, xbr, xtl, xtr;
  double ybl, ybr, ytl, ytr;
  //vector<double> xbotleft, xbotright, xtopleft, xtopright;
  //vector<double> ybotleft, ybotright, ytopleft, ytopright;

  for (int i = 0; i < BotVertices.size(); i++){
    xbl = BotVertices[i]->get_x(0);
    xbr = BotVertices[i]->get_x(1);
    xtl = TopVertices[i]->get_x(0);
    xtr = TopVertices[i]->get_x(1);
    ybl = BotVertices[i]->get_y(0);
    ybr = BotVertices[i]->get_y(1);
    ytl = TopVertices[i]->get_y(0);
    ytr = TopVertices[i]->get_y(1);
    
    cout << "wire x1=\"" << xbl << "\" y1=\"" << ybl << "\" x2=\""<< xtl << "\" y2=\"" << ytl << "\" width=\"0.1524\" layer=\"46\" curve=\"-180\"/> " << endl;
    cout << "<wire x1=\"" << xtl << "\" y1=\"" << ytl << "\" x2=\""<< xtr << "\" y2=\"" << ytr << "\" width=\"0.1524\" layer=\"46\"/> " << endl;
    cout << "<wire x1=\"" << xtr << "\" y1=\"" << ytr << "\" x2=\""<< xbr << "\" y2=\"" << ybr << "\" width=\"0.1524\" layer=\"46\" curve=\"-180\"/> " << endl;  
    cout << "<wire x1=\"" << xbr << "\" y1=\"" << ybr << "\" x2=\""<< xbl << "\" y2=\"" << ybl << "\" width=\"0.1524\" layer=\"46\"/> " << endl;
    cout << endl;
  }
  
  return 0;
}
