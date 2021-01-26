#include "../include/Lbinning.h"
#include <iostream>
#include "../include/Ltunfold.h"
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
#include "../include/TUnfoldDensity.h"
#include "../include/LgetNPcorrelation.h"
#include "../include/LsetJetPtWeight.h"


#include<functional>


int GetNpointpairinfo(TString inputfile,TString outputfile)
{
   // return 0;

    //TFile *f_tqH = new TFile("../../sample_root/gg_HToZZTo4L_M125_13TeV.root");
    /*input file*/
    
    TFile *f_tqH = new TFile(inputfile);
   /*define output location and name*/
    
//TFile *f_tqH = new TFile("../../sample_root/yyl/gg_HToZZTo4L_M125_13TeV.root","UPDATE");
    //TFile *f_tqH = new TFile("../../sample_root/yyl/BACKUP.root","UPDATE");
    
    TTree *t_tqH = (TTree *)f_tqH->Get("jetInfos/JetsAndDaughters");
    
    
  // LsetJetPtWeight(t_tqH,"Reco",tqH_gen_jetPt,30,30,400,2,tempRequirement1);
    LsetJetPtWeight(t_tqH,"Reco","JetPt",100,30,400);
    
    f_tqH->Close();
   TFile *f3= new TFile("Reco.root");
   t_tqH = (TTree *)f3->Get("JetsAndDaughters");
    
    //{tqH_reco_jetpartonflavour,-9,10,tqH_reco_Jettag}
    /*part2*/
    
     LsetJetPtWeight(t_tqH,"Gen","JetPt",100,30,400);
    f3->Close();
    //return 0;
     TFile *f2= new TFile("Gen.root");
     t_tqH = (TTree *)f2->Get("JetsAndDaughters");
    //return 0;
    /*part3*/
    
    cout<<1111<<endl;
    LgetNPcorrelation(t_tqH,"Reco","DaughterEta","DaughterPhi","DaughterEnergy","DaughterMatching","DaughterJetId","JetMatching",2);
    //cout<<i<<endl;
    cout<<1<<endl;
    f2->Close();
    TFile *f4= new TFile("Reco.root");
    t_tqH = (TTree *)f4->Get("JetsAndDaughters");
    /*part4*/
    LgetNPcorrelation(t_tqH,"Gen","DaughterEta","DaughterPhi","DaughterEnergy","DaughterMatching","DaughterJetId","JetMatching",2);
    f4->Close();
    
    //t_tqH->AddFriend(tSetJetPtWeight);
    
    gSystem->Exec("rm Reco.root");
    gSystem->Exec("mv Gen.root "+outputfile);

    
    return 0;
}
/*
int test(){
    
    
    struct binninginfo binningweight=Lbinning("../../chunk3/ZZ4lAnalysis.root","dr2",0,0,0.8,0,1.0/3);
   // struct binninginfo binningweight=Lreadoutbinfile("objbins.txt");
    for (int i=0;i<=binningweight.numberofbins;i++){
        cout<<binningweight.bin[i]<<endl;
    }
    
    double a[2]={1,2};
    double *b;
    b=a;
    cout<<b[1]<<endl;
    
    
    return 0;
}
*/

