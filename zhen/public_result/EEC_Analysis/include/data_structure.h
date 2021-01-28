//
//  Ltunfold.h
//  
//
//  Created by 林桢 on 2020/10/30.
//


#include <iostream>

#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TF2.h>
#include "TUnfoldDensity.h"
#ifndef data_structure_h
#define data_structure_h

using namespace std;
struct binninginfo{
    int numberofbins;
    int numberofhalfbins;
    double bin[100];
    double halfbin[100];
};
struct objrequirement{
    vector<float> *obj_parameter;
    double xmin;
    double xmax;
    vector<float> *level;
};
struct objweight{
    vector<float> *obj_parameter;
    vector<float> *level;
};

#endif


