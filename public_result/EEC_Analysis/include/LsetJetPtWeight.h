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
#include "data_structure.h"


#ifndef LsetJetPtWeight_h
#define LsetJetPtWeight_h
using namespace std;


TFile  * LsetJetPtWeight(TTree *t,TString level,TString bJetPt1,int nbins,double min,double max){
    //1.Jettag
    //2.Particle info Daughter.(Eta,Phi,Energy,Matchtag,Jettag)
    vector<float> *JetPt1 = new vector<float>;
    vector<float> *tqH_gen_jetpartonflavour = new vector<float>;
    cout<<"t "<<t->GetEntries()<<endl;
    TFile *f1=new TFile(level+".root","RECREATE");
    auto t1 = t->CloneTree(0);
    t1->CopyEntries(t);
    bJetPt1=level+bJetPt1;
    TString bPartonFlavour=level+"JetPartonFlavourHard";
    
    t1->SetBranchAddress(bJetPt1, &JetPt1);
    t1->SetBranchAddress(bPartonFlavour, &tqH_gen_jetpartonflavour);
    
    struct objrequirement particle_selected[2]{
        {tqH_gen_jetpartonflavour,20,22,NULL},
        {tqH_gen_jetpartonflavour,-9,9,NULL}
    };
    int numberofrequirements=2;
    

    
    
    vector<float> *Weight = new vector<float>;
    vector<float> *TempWeight = new vector<float>;
    struct objrequirement tempRequirement[2];
    double number[numberofrequirements+1];
    double density[numberofrequirements+1];
    
    //t1->Branch(level+"JetPtWeight",Weight);
    //Weight=t->AddBranch(level+"JetPtWeight");
    TBranch *b=t1->Branch(level+"JetPtWeight",Weight);
    
   // t1->SetBranchStatus("*",0);
   // t1->SetBranchStatus(level+"JetPtWeight",1);
    
    for(int r=0;r<numberofrequirements+1;r++){
        number[r]=0;
    }
    int c1=0;
    for(int k=0;k<t1->GetEntries();k++){
        t1->GetEntry(k);
        double weight[JetPt1->size()];
        bool condition[numberofrequirements];
        for(int j=0;j<JetPt1->size();j++){
            bool tempbool;
            bool tempbool1=true;
            for(int r=0;r<numberofrequirements;r++){
                tempRequirement[0]=particle_selected[r];
                condition[r]=Lselectedevent_multi_requirements(1,tempRequirement,j);
            }
            for(int r=0;r<numberofrequirements-1;r++){
                if(condition[r]){
                    number[r]=number[r]+1;
                    tempbool=true;
                    break;
                }
                tempbool1=tempbool1&&(!condition[r]);
                condition[r+1]=tempbool1&&condition[r+1];
            }
            if(condition[numberofrequirements-1]){
                number[numberofrequirements-1]=number[numberofrequirements-1]+1;
                tempbool=true;
                
            }
            if(tempbool){
                continue;
            }
                number[numberofrequirements]=number[numberofrequirements]+1;
            }
    }
    
    double ptrange[nbins+1];
    double Bweight[numberofrequirements+1][nbins];
    ptrange[0]=min;
    ptrange[nbins]=max;
    for(int i=1;i<nbins;i++){
        ptrange[i]=min+(max-min)/nbins*i;
    }
    
    for(int r=0;r<numberofrequirements+1;r++){
        cout<<number[r]<<endl;
        density[r]=number[r]*1.0/(max-min);
        //cout<<density[r]<<endl;
        number[r]=0;
    }
    
    for(int i=0;i<nbins;i++){
        for(int k=0;k<t1->GetEntries();k++){
            t1->GetEntry(k);
            double weight[JetPt1->size()];
            bool condition[numberofrequirements];
            bool condition1[numberofrequirements];
            for(int j=0;j<JetPt1->size();j++){
                for(int r=0;r<numberofrequirements;r++){
                    tempRequirement[0]=particle_selected[r];
                    tempRequirement[1]={JetPt1,ptrange[i],ptrange[i+1],NULL};
                    condition[r]=Lselectedevent_multi_requirements(2,tempRequirement,j);
                    condition1[r]=Lselectedevent_multi_requirements(1,tempRequirement,j);
                }
                bool tempbool1=true;
                for(int r=0;r<numberofrequirements-1;r++){
                    if(condition[r]&&condition1[r]){
                        number[r]=number[r]+1;
                    }
                    tempbool1=tempbool1&&(!condition1[r]);
                    condition1[r+1]=tempbool1&&condition1[r+1];
                }
                
                if(condition[numberofrequirements-1]&&condition1[numberofrequirements-1]){
                    number[numberofrequirements-1]=number[numberofrequirements-1]+1;
                }
                
                tempbool1=tempbool1&&(!condition1[numberofrequirements-1]);
                condition1[numberofrequirements]=tempbool1;
    
                tempRequirement[0]={JetPt1,ptrange[i],ptrange[i+1],NULL};
                condition[numberofrequirements]=(condition1[numberofrequirements-1])&&Lselectedevent_multi_requirements(1,tempRequirement,j);
                if(condition[numberofrequirements]){
                    number[numberofrequirements]=number[numberofrequirements]+1;
                }
            }
            
        }
        
        for(int r=0;r<numberofrequirements+1;r++){
            //cout<<number[r]<<" "<<r<<endl;
            if(number[r]==0){
                Bweight[r][i]=1;
            }
            else
            Bweight[r][i]=density[r]/(number[r]/(ptrange[i+1]-ptrange[i]));
            number[r]=0;
            //cout<<Bweight[r][i]<<endl;
        }
    }
    for(int i=0;i<nbins;i++){
        for(int r=0;r<numberofrequirements+1;r++){
            //cout<<"r="<<r<<"i="<<i<<"Bweight"<<Bweight[r][i]<<endl;
        }
    }
    
    
    
    double aaa=0;
    for(int k=0;k<t1->GetEntries();k++){
        
        t1->GetEntry(k);
        double weight[JetPt1->size()];
        bool condition[numberofrequirements+1];
        bool condition1[numberofrequirements+1];
        for(int j=0;j<JetPt1->size();j++){
            int flavourid=-1;
            int ptid=-1;
            bool getit=false;
            for(int i=0;i<nbins;i++){
                
                for(int r=0;r<numberofrequirements;r++){
                    tempRequirement[0]=particle_selected[r];
                    tempRequirement[1]={JetPt1,ptrange[i],ptrange[i+1],NULL};
                    condition[r]=Lselectedevent_multi_requirements(2,tempRequirement,j);
                    condition1[r]=Lselectedevent_multi_requirements(1,tempRequirement,j);
                }
                
                bool tempbool1=true;
                for(int r=0;r<numberofrequirements-1;r++){
                    if(condition[r]&&condition1[r]){
                        flavourid=r;
                        ptid=i;
                        //cout<<ptid<<" "<<flavourid<<endl;
                        getit=true;
                    }
                    tempbool1=tempbool1&&(!condition1[r]);
                    condition1[r+1]=tempbool1&&condition1[r+1];
                }
                if(condition[numberofrequirements-1]&&condition1[numberofrequirements-1]){
                    flavourid=numberofrequirements-1;
                    ptid=i;
                    getit=true;
                }
                tempbool1=tempbool1&&(!condition1[numberofrequirements-1]);
                condition1[numberofrequirements]=tempbool1;
                
                tempRequirement[0]={JetPt1,ptrange[i],ptrange[i+1],NULL};
                condition[numberofrequirements]=(condition1[numberofrequirements])&&Lselectedevent_multi_requirements(1,tempRequirement,j);
                if(condition[numberofrequirements]){
                    flavourid=numberofrequirements;
                    ptid=i;
                    getit=true;
                }
                if(getit){
                    break;
                }
            }
            if(flavourid==0&&ptid!=0)
            //cout<<flavourid<<"flavourid"<<ptid<<"ptid"<<Bweight[flavourid][ptid]<<endl;
            aaa=aaa+Bweight[flavourid][ptid];
            //cout<<Bweight[flavourid][ptid]<<endl;
            if(ptid==-1||flavourid==-1){
                Weight->push_back(1);
            }
            else{
            Weight->push_back(Bweight[flavourid][ptid]);
            }
        }
        //cout<<aaa<<endl;
        //t1->Fill();
       
        b->Fill();
        //t1->Fill();
        Weight->clear();
    }

    cout<<"t1  "<<t1->GetEntries()<<endl;
    //t1->Write();
    //t1->Write();
    f1->Write();
    f1->Close();
    
    return f1;
}


#endif /* Ltunfold_h */
