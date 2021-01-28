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


using namespace std;

double LDeltaR(double Eta1,double Phi1,double Eta2,double Phi2){
    double dEta=Eta1-Eta2;
    double pi=3.1415926;
    
    
    while(Phi1>2*pi){
        Phi1=Phi1-2*pi;
    }
    while(Phi2>2*pi){
        Phi2=Phi2-2*pi;
    }
    
    while(Phi1<0){
        Phi1=Phi1+2*pi;
    }
    while(Phi2<0){
        Phi2=Phi2+2*pi;
    }
    
    
    double dPhi=fabs(Phi1-Phi2);
    
    if(dPhi>pi){
        dPhi=2*pi-dPhi;
    }
    else{
        dPhi=dPhi;
    }
    
    return TMath::Sqrt(dEta*dEta+dPhi*dPhi);
}


void LgetNPcorrelation(TTree *t, TString name_of_root_file,TString bDaughter_Eta,TString bDaughter_Phi,TString bDaughter_Energy,TString bDaughter_Matchtag,TString bDaughter_Jettag,TString bJet_Matchtag,int NPoint=2){
    //1.Jettag
    //2.Particle info Daughter.(Eta,Phi,Energy,Matchtag,Jettag)
    int nDaughter=100;
    int maxpairnum=1;
    for(int j=0;j<NPoint;j++){
        maxpairnum=nDaughter*maxpairnum;
    }
    
    
    TString branch_part1_name=name_of_root_file+TString::Format("N%d_",NPoint);
    
    int paircount=0;
    int matchpaircount=0;
    int paircount1=0;
     cout<<"t "<<t->GetEntries()<<endl;
    TFile *f1=new TFile(name_of_root_file+".root","RECREATE");
    auto t1 = t->CloneTree(0);
    t1->CopyEntries(t);
    
    vector<float> *pair_Jettag = new vector<float>;
    vector<float> *pair_Matchtag = new vector<float>;
    vector<float> *pair_Weight = new vector<float>;
    vector<float> *pair_DeltaR = new vector<float>;
    
    vector<float> *Daughter_Eta= new vector<float>;
    vector<float> *Daughter_Phi= new vector<float>;
    vector<float> *Daughter_Energy= new vector<float>;
    vector<float> *Daughter_Matchtag= new vector<float>;
    vector<float> *Daughter_Jettag= new vector<float>;
    vector<float> *Jet_Matchtag= new vector<float>;
    //return 0;
    
    
    
    
    //t1->Branch("paircount1",&paircount1,"paircount1/I");
    t1->SetBranchAddress(name_of_root_file+bDaughter_Eta,&Daughter_Eta);
    t1->SetBranchAddress(name_of_root_file+bDaughter_Phi,&Daughter_Phi);
    t1->SetBranchAddress(name_of_root_file+bDaughter_Energy,&Daughter_Energy);
    t1->SetBranchAddress(name_of_root_file+bDaughter_Matchtag,&Daughter_Matchtag);
    t1->SetBranchAddress(name_of_root_file+bDaughter_Jettag,&Daughter_Jettag);
    t1->SetBranchAddress(name_of_root_file+bJet_Matchtag,&Jet_Matchtag);
    
    
    
    TBranch *b1=t1->Branch(branch_part1_name+"pair_Jettag",pair_Jettag);
    TBranch *b2=t1->Branch(branch_part1_name+"pair_Matchtag",pair_Matchtag);
    TBranch *b3=t1->Branch(branch_part1_name+"pair_Weight",pair_Weight);
    TBranch *b4=t1->Branch(branch_part1_name+"pair_DeltaR",pair_DeltaR);
    
    
    
    
    
    
    cout<<1<<endl;
    // return 0;
    double step[11];
    step[0]=0;
    int s=1;
    for(Int_t i=0; i<t1->GetEntries();i++){
        t1->GetEntry(i);
        step[s]=i*1.0/t1->GetEntries();
        if(step[s]-step[s-1]>0.1){
            cout << setprecision(1)<<i*1.0/t1->GetEntries()<<endl;
            s=s+1;
            //break;
        }
        int jetnum=Jet_Matchtag->size();
        double Q[jetnum];
        
        for(int k=0;k<jetnum;k++){
            Q[k]=0;
        }
        for (Int_t k=0;k<Daughter_Matchtag->size();k++){
            int temp=Daughter_Jettag->at(k);
            Q[temp]=Daughter_Energy->at(k)+Q[temp];
            if(temp>10)
                cout<<temp<<endl;
            //cout<<temp<<endl;
        }
        
        int k[NPoint];
        int pairnum=1;
        int Daughter_num=Daughter_Matchtag->size();
        for(int j=0;j<NPoint;j++){
            pairnum=pairnum*Daughter_num;
        }
        
        for(int kn=0;kn<pairnum;kn++){
            int kn_temp=kn;
            for(int ki=0;ki<NPoint;ki++){
                k[ki]= kn_temp%Daughter_num;
                kn_temp=kn_temp/Daughter_num;
            }
            bool is_same_jet=true;
            for(int ki=1;ki<NPoint;ki++){    is_same_jet=is_same_jet&&(Daughter_Jettag->at(k[0])==Daughter_Jettag->at(k[ki]));
            }
            if(is_same_jet){
                int jetid=Daughter_Jettag->at(k[0]);
                pair_Jettag->push_back(Daughter_Jettag->at(k[0]));
                
                //if(pair_Jettag[paircount]>10)
                //cout<<pair_Jettag[paircount]<<endl;
                bool ismatch_particle=true;
                int dismatchpaircount=0;
                for(int ki=0;ki<NPoint;ki++){
                    ismatch_particle=ismatch_particle&&(Daughter_Matchtag->at(k[ki])>=0);
                }
                if(ismatch_particle){
                    
                    pair_Matchtag->push_back(matchpaircount);
                    matchpaircount++;
                }
                else{
                    pair_Matchtag->push_back(-1);
                    dismatchpaircount++;
                }
                double Weight_temp=1;
                for(int ki=0;ki<NPoint;ki++){
                    Weight_temp=Weight_temp*Daughter_Energy->at(k[ki])/Q[jetid];
                }
                pair_Weight->push_back(Weight_temp);
                if(Weight_temp>1){
                    cout<<Weight_temp<<pairnum<<" "<<Q[0]<<" "<<jetid<<endl;
                }
                double DeltaR_temp=0;
                double Max_DeltaR_temp=0;
                for(int ki=0;ki<NPoint;ki++){
                    for(int ki2=ki+1;ki2<NPoint;ki2++){
                        DeltaR_temp=LDeltaR(Daughter_Eta->at(k[ki]),Daughter_Phi->at(k[ki]),Daughter_Eta->at(k[ki2]),Daughter_Phi->at(k[ki2]));
                        if(DeltaR_temp>Max_DeltaR_temp){
                            Max_DeltaR_temp=DeltaR_temp;
                        }
                    }
                }
                pair_DeltaR->push_back(Max_DeltaR_temp);
                // cout<<Max_DeltaR_temp<<endl;
                paircount++;
            }
            
        }
        paircount1=paircount;
        //cout<<paircount1<<endl;
        b1->Fill();
        b2->Fill();
        b3->Fill();
        b4->Fill();
        //t1->Fill();
        paircount=0;
        matchpaircount=0;
        pair_Jettag->clear();
        pair_Matchtag->clear();
        pair_Weight->clear();
        pair_DeltaR->clear();
        
    }
    f1->cd();
    cout<<1<<endl;
    
    t1->Write("",TObject::kOverwrite);
    f1->Close();

    return 0;

    //return ;
    cout<<2222<<endl;
}
