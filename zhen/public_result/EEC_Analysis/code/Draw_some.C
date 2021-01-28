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


#include<functional>


int Draw_some()
{
    // return 0;
    
    //TFile *f_tqH = new TFile("../../sample_root/yyl/tqH_HToZZTo4L_M125.root");
    TFile *f_tqH = new TFile("../../sample_root/pair/tqH_HToZZTo4L_M125_pair.root");
    TFile *f_ggH = new TFile("../../sample_root/pair/gg_HToZZTo4L_M125_13TeV_pair.root");
    
    //TFile *f_tqH = new TFile("../../sample_root/pair/Gen.root");
    
    TTree *t_tqH = (TTree *)f_tqH->Get("JetsAndDaughters");
    TTree *t_ggH = (TTree *)f_ggH->Get("JetsAndDaughters");
    
    vector<float> *tqH_reco_Phi = new vector<float>;
    vector<float> *tqH_reco_Eta = new vector<float>;
    vector<float> *tqH_reco_Energy = new vector<float>;
    vector<float> *tqH_reco_Jettag = new vector<float>;
    vector<float> *tqH_reco_JetMtag = new vector<float>;
    vector<float> *tqH_reco_jetpartonflavour = new vector<float>;
    vector<float> *tqH_reco_jetPt = new vector<float>;
    vector<float> *tqH_reco_JetPtWeight = new vector<float>;
    
    vector<float> *tqH_reco_N2_pair_Jettag = new vector<float>;
    vector<float> *tqH_reco_N2_pair_Matchtag = new vector<float>;
    vector<float> *tqH_reco_N2_pair_Weight = new vector<float>;
    vector<float> *tqH_reco_N2_pair_DeltaR = new vector<float>;
    vector<float> *tqH_reco_Daughter_Matchtag = new vector<float>;
    
    
    
    vector<float> *ggH_gen_Phi = new vector<float>;
    vector<float> *ggH_gen_Eta = new vector<float>;
    vector<float> *ggH_gen_Energy = new vector<float>;
    vector<float> *ggH_gen_Jettag = new vector<float>;
    vector<float> *ggH_gen_JetMtag = new vector<float>;
    vector<float> *ggH_gen_jetpartonflavour = new vector<float>;
    vector<float> *ggH_gen_jetPt = new vector<float>;
    vector<float> *ggH_gen_JetPtWeight = new vector<float>;
    vector<float> *ggH_gen_N2_pair_Jettag = new vector<float>;
    vector<float> *ggH_gen_N2_pair_Matchtag = new vector<float>;
    vector<float> *ggH_gen_N2_pair_Weight = new vector<float>;
    vector<float> *ggH_gen_N2_pair_DeltaR = new vector<float>;
    vector<float> *ggH_gen_Daughter_Matchtag = new vector<float>;
    
    vector<float> *tqH_gen_Phi = new vector<float>;
    vector<float> *tqH_gen_Eta = new vector<float>;
    vector<float> *tqH_gen_Energy = new vector<float>;
    vector<float> *tqH_gen_Jettag = new vector<float>;
    vector<float> *tqH_gen_JetMtag = new vector<float>;
    vector<float> *tqH_gen_jetpartonflavour = new vector<float>;
    vector<float> *tqH_gen_jetPt = new vector<float>;
    vector<float> *tqH_gen_JetPtWeight = new vector<float>;
    vector<float> *tqH_gen_N2_pair_Jettag = new vector<float>;
    vector<float> *tqH_gen_N2_pair_Matchtag = new vector<float>;
    vector<float> *tqH_gen_N2_pair_Weight = new vector<float>;
    vector<float> *tqH_gen_N2_pair_DeltaR = new vector<float>;
    vector<float> *tqH_gen_Daughter_Matchtag = new vector<float>;
    vector<float> *tqH_gen_Pdgid = new vector<float>;
    vector<float> *ggH_gen_Pdgid = new vector<float>;

    
    
    
    
    
    //t_tqH->SetBranchAddress("RecoDaughterEnergy", &tqH_reco_Energy);
    t_tqH->SetBranchAddress("RecoJetPt", &tqH_reco_Energy);
    t_tqH->SetBranchAddress("RecoDaughterEta", &tqH_reco_Eta);
    t_tqH->SetBranchAddress("RecoDaughterJetId", &tqH_reco_Jettag);
    t_tqH->SetBranchAddress("RecoDaughterPhi", &tqH_reco_Phi);
    t_tqH->SetBranchAddress("RecoDaughterEta", &tqH_reco_Eta);
    //t_tqH->SetBranchAddress("RecoDaughterMatching", &tqH_reco_JetMtag);
    t_tqH->SetBranchAddress("RecoJetMatching", &tqH_reco_JetMtag);
    t_tqH->SetBranchAddress("RecoJetPtWeight", &tqH_reco_JetPtWeight);
    //t_tqH->SetBranchAddress("detjetPartonflavour", &tqH_reco_jetpartonflavour);
    //t_tqH->SetBranchAddress("RecoJetPt", &tqH_reco_jetPt);
    t_tqH->SetBranchAddress("RecoPartonFlavor", &tqH_reco_jetpartonflavour);
    t_tqH->SetBranchAddress("RecoN2_pair_Jettag", &tqH_reco_N2_pair_Jettag);
    t_tqH->SetBranchAddress("RecoN2_pair_Matchtag", &tqH_reco_N2_pair_Matchtag);
    t_tqH->SetBranchAddress("RecoN2_pair_Weight", &tqH_reco_N2_pair_Weight);
    t_tqH->SetBranchAddress("RecoN2_pair_DeltaR", &tqH_reco_N2_pair_DeltaR);
    t_tqH->SetBranchAddress("RecoDaughterMatching", &tqH_reco_Daughter_Matchtag);
    
    
    t_tqH->SetBranchAddress("GenDaughterEta", &tqH_gen_Eta);
    //t_tqH->SetBranchAddress("GenDaughterEnergy", &tqH_gen_Energy);
    t_tqH->SetBranchAddress("GenJetPt", &tqH_gen_Energy);
    t_tqH->SetBranchAddress("GenDaughterJetId", &tqH_gen_Jettag);
    t_tqH->SetBranchAddress("GenDaughterPdgId", &tqH_gen_Pdgid);
    t_tqH->SetBranchAddress("GenDaughterPhi", &tqH_gen_Phi);
    t_tqH->SetBranchAddress("GenDaughterEta", &tqH_gen_Eta);
    //t_tqH->SetBranchAddress("GenDaughterMatching", &tqH_gen_JetMtag);
    t_tqH->SetBranchAddress("GenJetMatching", &tqH_gen_JetMtag);
    t_tqH->SetBranchAddress("GenJetPtWeight", &tqH_gen_JetPtWeight);
    
    t_tqH->SetBranchAddress("GenDaughterMatching", &tqH_gen_Daughter_Matchtag);
    
    
    t_tqH->SetBranchAddress("GenPartonFlavor", &tqH_gen_jetpartonflavour);
    //t_tqH->SetBranchAddress("GenJetPt", &tqH_gen_jetPt);
   // t_tqH->SetBranchAddress("GenPartonFlavor", &tqH_gen_jetPt);
    t_tqH->SetBranchAddress("GenN2_pair_Jettag", &tqH_gen_N2_pair_Jettag);
    t_tqH->SetBranchAddress("GenN2_pair_Matchtag", &tqH_gen_N2_pair_Matchtag);
    t_tqH->SetBranchAddress("GenN2_pair_Weight", &tqH_gen_N2_pair_Weight);
    t_tqH->SetBranchAddress("GenN2_pair_DeltaR", &tqH_gen_N2_pair_DeltaR);
    
    t_ggH->SetBranchAddress("GenDaughterEta", &ggH_gen_Eta);
     //t_ggH->SetBranchAddress("GenDaughterEnergy", &ggH_gen_Energy);
     t_ggH->SetBranchAddress("GenJetPt", &ggH_gen_Energy);
     t_ggH->SetBranchAddress("GenDaughterJetId", &ggH_gen_Jettag);
    t_ggH->SetBranchAddress("GenDaughterPdgId", &ggH_gen_Pdgid);
     t_ggH->SetBranchAddress("GenDaughterPhi", &ggH_gen_Phi);
     t_ggH->SetBranchAddress("GenDaughterEta", &ggH_gen_Eta);
     //t_ggH->SetBranchAddress("GenDaughterMatching", &ggH_gen_JetMtag);
     t_ggH->SetBranchAddress("GenJetMatching", &ggH_gen_JetMtag);
     t_ggH->SetBranchAddress("GenJetPtWeight", &ggH_gen_JetPtWeight);
     
     t_ggH->SetBranchAddress("GenDaughterMatching", &ggH_gen_Daughter_Matchtag);
     
     
     t_ggH->SetBranchAddress("GenPartonFlavor", &ggH_gen_jetpartonflavour);
     //t_ggH->SetBranchAddress("GenJetPt", &ggH_gen_jetPt);
    // t_ggH->SetBranchAddress("GenPartonFlavor", &ggH_gen_jetPt);
     t_ggH->SetBranchAddress("GenN2_pair_Jettag", &ggH_gen_N2_pair_Jettag);
     t_ggH->SetBranchAddress("GenN2_pair_Matchtag", &ggH_gen_N2_pair_Matchtag);
     t_ggH->SetBranchAddress("GenN2_pair_Weight", &ggH_gen_N2_pair_Weight);
     t_ggH->SetBranchAddress("GenN2_pair_DeltaR", &ggH_gen_N2_pair_DeltaR);
    
    /*test for quark&&gluon particle number in different JetPt  */
    /*
    struct objrequirement t_quark_only_select_pjet_gen1[10]={
        {tqH_gen_JetMtag,0,10000,tqH_gen_Jettag},
        {tqH_gen_Daughter_Matchtag,0,10000,NULL},
        {tqH_gen_jetpartonflavour,-9,9,tqH_gen_Jettag}
    };
    struct objrequirement t_gluon_only_select_pjet_gen1[10]={
        {tqH_gen_JetMtag,0,10000,tqH_gen_Jettag},
        {tqH_gen_Daughter_Matchtag,0,10000,NULL},
        {tqH_gen_jetpartonflavour,21,22,tqH_gen_Jettag}
    };
    struct objweight t_gen_weight1[10]={
        {tqH_gen_JetPtWeight,tqH_gen_N2_pair_Jettag}
    };
    
    
   */
    
    
    
    /*test for quark&&gluon EEC distribution at different JetPt  */
    
    struct objrequirement t_quark_only_select_pjet_gen[10]={
        {tqH_gen_JetMtag,-1,10000,tqH_gen_N2_pair_Jettag},
        {tqH_gen_N2_pair_Matchtag,-1,10000,NULL},
        {tqH_gen_jetpartonflavour,-9,9,tqH_gen_N2_pair_Jettag}
    };
    struct objrequirement t_gluon_only_select_pjet_gen[10]={
        {ggH_gen_JetMtag,-1,10000,ggH_gen_N2_pair_Jettag},
        {ggH_gen_N2_pair_Matchtag,-1,10000,NULL},
        {ggH_gen_jetpartonflavour,20,22,ggH_gen_N2_pair_Jettag}
    };
    
    
    
    struct objrequirement t_quark_only_select_pjet_reco[10]={
        {tqH_gen_JetMtag,-1,10000,tqH_gen_N2_pair_Jettag},
        {tqH_gen_N2_pair_Matchtag,-1,10000,NULL},
        {tqH_gen_jetpartonflavour,20,22,tqH_gen_N2_pair_Jettag}
    };
    
    
    
    struct objrequirement t_quark_only_select_pjet_gen1[10]={
        {tqH_gen_JetMtag,-1,10000,tqH_gen_Jettag},
        {tqH_gen_Daughter_Matchtag,-1,10000,NULL},
        {tqH_gen_jetpartonflavour,-9,9,tqH_gen_Jettag}
    };
    struct objrequirement t_quark_only_select_pjet_reco1[10]={
        {ggH_gen_JetMtag,-1,10000,ggH_gen_Jettag},
        {ggH_gen_Daughter_Matchtag,-1,10000,NULL},
        {ggH_gen_jetpartonflavour,20,22,ggH_gen_Jettag}
    };
    
    struct objrequirement t_quark_only_select_pjet_gen2[10]={
        {tqH_gen_JetMtag,-1,10000,NULL},
        {tqH_gen_jetpartonflavour,-9,9,NULL}
    };
    struct objrequirement t_quark_only_select_pjet_reco2[10]={
        {ggH_gen_JetMtag,-1,10000,NULL},
        {ggH_gen_jetpartonflavour,20,22,NULL}
    };
    
    
    
    struct objweight tqH_gen_weight[10]={
        {tqH_gen_N2_pair_Weight,NULL},
        {tqH_gen_JetPtWeight,tqH_gen_N2_pair_Jettag}
    };
    struct objweight ggH_gen_weight[10]={
        {ggH_gen_N2_pair_Weight,NULL},
        {ggH_gen_JetPtWeight,ggH_gen_N2_pair_Jettag}
    };
    
    struct objweight t_reco_weight[10]={
        {tqH_gen_N2_pair_Weight,NULL},
        {tqH_gen_JetPtWeight,tqH_gen_N2_pair_Jettag}
    };
    
    //{tqH_reco_jetpartonflavour,-9,10,tqH_reco_Jettag}
    struct binninginfo binningEnergy=Lbinning(t_tqH,"jetPt",tqH_gen_jetPt,tqH_gen_jetPt,2,30,400,0,1.0/3,0.1,16);
    struct binninginfo binningdR=Lbinning(t_tqH,"pairdR",tqH_gen_jetPt,tqH_gen_jetPt,3,0,0.8,0,1,0.1,30);
    
    TH1D *tqH_perjet_pairs= new TH1D("tqH_perjet_pairs","",binningEnergy.numberofbins,binningEnergy.bin);
    TH1D *ggH_perjet_pairs= new TH1D("ggH_perjet_pairs","",binningEnergy.numberofbins,binningEnergy.bin);
    TH1D *tqH_perjet_cpairs= new TH1D("tqH_perjet_cpairs","",binningEnergy.numberofbins,binningEnergy.bin);
    TH1D *ggH_perjet_cpairs= new TH1D("ggH_perjet_cpairs","",binningEnergy.numberofbins,binningEnergy.bin);
    TH1D *tqH_perjet_npairs= new TH1D("tqH_perjet_npairs","",binningEnergy.numberofbins,binningEnergy.bin);
    TH1D *ggH_perjet_npairs= new TH1D("ggH_perjet_npairs","",binningEnergy.numberofbins,binningEnergy.bin);
    TH1D *tqH_perjet_cpairs_ratio= new TH1D("tqH_perjet_cpairs_ratio","",binningEnergy.numberofbins,binningEnergy.bin);
    TH1D *ggH_perjet_cpairs_ratio= new TH1D("ggH_perjet_cpairs_ratio","",binningEnergy.numberofbins,binningEnergy.bin);
    TH1D *tqH_perjet_npairs_ratio= new TH1D("tqH_perjet_npairs_ratio","",binningEnergy.numberofbins,binningEnergy.bin);
    TH1D *ggH_perjet_npairs_ratio= new TH1D("ggH_perjet_npairs_ratio","",binningEnergy.numberofbins,binningEnergy.bin);
    
    TH1D *ratio_perjet_pairs1= new TH1D("ratio_perjet_pairs1","",binningEnergy.numberofbins,binningEnergy.bin);
    TH1D *ratio_perjet_pairs2= new TH1D("ratio_perjet_pairs2","",binningEnergy.numberofbins,binningEnergy.bin);
    
    TH1D *hist_tqH_EEC_Energy_Resolution_gen[binningEnergy.numberofbins];
    TH1D *hist_ggH_EEC_Energy_Resolution_gen[binningEnergy.numberofbins];
    
    TH1D *hist_tqH_EEC_Energy_Resolution_reco[binningEnergy.numberofbins];
    TH1D *result= new TH1D("JetPtDifference","",binningdR.numberofbins,binningdR.bin);
    //hist_tqH_EEC_Energy_Resolution_gen[0]->DrawNormalized("same");
    //return 0;
    TH1F* a=new TH1F("a","a",100,0,10);
    t_tqH->Draw("GenDaughterPt>>a");
    return 0;
    
    for(int ptrange=1;ptrange<=binningEnergy.numberofbins;ptrange++){
        cout<<ptrange<<endl;
        int particlenum1=0;
        int cparticlenum1=0;
       int jetnum1=0;
        int particlenum2=0;
        int cparticlenum2=0;
        int jetnum2=0;
        t_quark_only_select_pjet_gen[3]={tqH_gen_Energy,binningEnergy.bin[ptrange-1],binningEnergy.bin[ptrange],tqH_gen_N2_pair_Jettag};
        t_gluon_only_select_pjet_gen[3]={ggH_gen_Energy,binningEnergy.bin[ptrange-1],binningEnergy.bin[ptrange],ggH_gen_N2_pair_Jettag};
        
        t_quark_only_select_pjet_reco[3]={tqH_gen_Energy,binningEnergy.bin[ptrange-1],binningEnergy.bin[ptrange],tqH_gen_N2_pair_Jettag};
        
        t_quark_only_select_pjet_gen1[3]={tqH_gen_Energy,binningEnergy.bin[ptrange-1],binningEnergy.bin[ptrange],tqH_gen_Jettag};
        t_quark_only_select_pjet_reco1[3]={ggH_gen_Energy,binningEnergy.bin[ptrange-1],binningEnergy.bin[ptrange],ggH_gen_Jettag};
        t_quark_only_select_pjet_gen2[3]={tqH_gen_Energy,binningEnergy.bin[ptrange-1],binningEnergy.bin[ptrange],NULL};
        t_quark_only_select_pjet_reco2[3]={ggH_gen_Energy,binningEnergy.bin[ptrange-1],binningEnergy.bin[ptrange],NULL};
        for(int k=0;k<t_tqH->GetEntries();k++){
            t_tqH->GetEntry(k);
            for(int i=0;i<tqH_gen_Daughter_Matchtag->size();i++){
                if(Lselectedevent_multi_requirements(4,t_quark_only_select_pjet_gen1,i)){
                    int pdgid=tqH_gen_Pdgid->at(i);
                    pdgid=abs(pdgid);
                    pdgid=pdgid/10;
                    int a[5];
                     a[2]=pdgid%10;
                    pdgid=pdgid/10;
                     a[3]=pdgid%10;
                    pdgid=pdgid/10;
                     a[4]=pdgid%10;
                    bool isodd=((a[2]+a[3])%2==0);
                    bool two_odd_one_even=false;
                    if((a[2]+a[3]+a[4])%2==0){
                        for(int iin=2;iin<=4;iin++){
                            if(a[iin]%2!=0){
                                two_odd_one_even=true;
                                break;
                            }
                        }
                    }
                    
                    bool ischarged=false;
                    
                    if(pdgid>0){
                        if(!two_odd_one_even){
                            ischarged=true;
                        }
                    }
                    if(pdgid<1){
                        if(!isodd){
                            ischarged=true;
                        }
                    }
                    if(ischarged){
                        cparticlenum1++;
                    }
                    
                    particlenum1++;
                }
            }

            for(int i=0;i<tqH_gen_Energy->size();i++){
                if(Lselectedevent_multi_requirements(4,t_quark_only_select_pjet_gen2,i)){
                    jetnum1++;
                }
            }
        }
        for(int k=0;k<t_ggH->GetEntries();k++){
            t_ggH->GetEntry(k);
            for(int i=0;i<ggH_gen_Daughter_Matchtag->size();i++){
                if(Lselectedevent_multi_requirements(4,t_quark_only_select_pjet_reco1,i)){
                    int pdgid=ggH_gen_Pdgid->at(i);
                    pdgid=abs(pdgid);
                    pdgid=pdgid/10;
                    int a[5];
                     a[2]=pdgid%10;
                    pdgid=pdgid/10;
                     a[3]=pdgid%10;
                    pdgid=pdgid/10;
                     a[4]=pdgid%10;
                    bool isodd=((a[2]+a[3])%2==0);
                    bool two_odd_one_even=false;
                    if((a[2]+a[3]+a[4])%2==0){
                        for(int iin=2;iin<=4;iin++){
                            if(a[iin]%2!=0){
                                two_odd_one_even=true;
                                break;
                            }
                        }
                    }
                    
                    bool ischarged=false;
                    
                    if(pdgid>0){
                        if(two_odd_one_even){
                            ischarged=true;
                        }
                    }
                    if(pdgid<1){
                        if(!isodd){
                            ischarged=true;
                        }
                    }
                    if(ischarged){
                        cparticlenum2++;
                    }
                    
                    particlenum2++;
                }
            }
            for(int i=0;i<ggH_gen_Energy->size();i++){
                if(Lselectedevent_multi_requirements(4,t_quark_only_select_pjet_reco2,i)){
                    jetnum2++;
                }
            }
        }
        if(jetnum1==0){
            jetnum1=1;
        }
        if(jetnum2==0){
            jetnum2=1;
        }
        tqH_perjet_pairs->SetBinContent(ptrange,particlenum1*1.0/jetnum1);
        ggH_perjet_pairs->SetBinContent(ptrange,particlenum2*1.0/jetnum2);
        tqH_perjet_cpairs->SetBinContent(ptrange,cparticlenum1*1.0/jetnum1);
        ggH_perjet_cpairs->SetBinContent(ptrange,cparticlenum2*1.0/jetnum2);
        tqH_perjet_npairs->SetBinContent(ptrange,(particlenum1-cparticlenum1)*1.0/jetnum1);
        ggH_perjet_npairs->SetBinContent(ptrange,(particlenum2-cparticlenum2)*1.0/jetnum2);
        tqH_perjet_cpairs_ratio->SetBinContent(ptrange,cparticlenum1*1.0/particlenum1);
        ggH_perjet_cpairs_ratio->SetBinContent(ptrange,cparticlenum2*1.0/particlenum2);
        tqH_perjet_npairs_ratio->SetBinContent(ptrange,(particlenum1-cparticlenum1)*1.0/particlenum1);
        ggH_perjet_npairs_ratio->SetBinContent(ptrange,(particlenum2-cparticlenum2)*1.0/particlenum2);
         double ratio1=particlenum2*1.0/jetnum2/(particlenum1*1.0/jetnum1);
        ratio_perjet_pairs1->SetBinContent(ptrange,ratio1);
        //hist_tqH_EEC_Energy_Resolution_reco[ptrange-1]=Fillhistwithmultiweight_tag(t_tqH,tqH_gen_N2_pair_DeltaR,2,t_reco_weight,TString::Format("hist_ggH_EEC_Energy_Resolution_gen_gluon@%d",ptrange),binningdR.numberofbins,binningdR.bin,0,1.0,4,t_quark_only_select_pjet_reco,true);
        hist_tqH_EEC_Energy_Resolution_gen[ptrange-1]=Fillhistwithmultiweight_tag(t_tqH,tqH_gen_N2_pair_DeltaR,2,tqH_gen_weight,TString::Format("hist_tqH_EEC_Energy_Resolution_gen_quark@%d",ptrange),binningdR.numberofbins,binningdR.bin,0,1.0/20000,4,t_quark_only_select_pjet_gen,true);
        hist_ggH_EEC_Energy_Resolution_gen[ptrange-1]=Fillhistwithmultiweight_tag(t_ggH,ggH_gen_N2_pair_DeltaR,2,ggH_gen_weight,TString::Format("hist_ggH_EEC_Energy_Resolution_gen_Gluon@%d",ptrange),binningdR.numberofbins,binningdR.bin,0,1.0/10000,4,t_gluon_only_select_pjet_gen,true);
        
        cout<<"xx"<<hist_tqH_EEC_Energy_Resolution_gen[ptrange-1]->GetBinContent(1)<<endl;
        hist_tqH_EEC_Energy_Resolution_gen[ptrange-1]=GetNormalized(hist_tqH_EEC_Energy_Resolution_gen[ptrange-1]);
        hist_ggH_EEC_Energy_Resolution_gen[ptrange-1]=GetNormalized(hist_ggH_EEC_Energy_Resolution_gen[ptrange-1]);
       hist_tqH_EEC_Energy_Resolution_reco[ptrange-1]=GetNormalized(hist_tqH_EEC_Energy_Resolution_gen[ptrange-1]);
         
        double ratio2=hist_tqH_EEC_Energy_Resolution_gen[ptrange-1]->GetBinContent(1)*1.0/hist_ggH_EEC_Energy_Resolution_gen[ptrange-1]->GetBinContent(1);
        cout<<"ratio1="<<ratio1<<" "<<"numberofP per Qjet="<<particlenum1*1.0/jetnum1<<" "<<"numberofP per Gjet="<<particlenum2*1.0/jetnum2<<endl;
        cout<<"ratio2="<<ratio2<<" "<<"tqH Bin1="<<hist_tqH_EEC_Energy_Resolution_gen[ptrange-1]->GetBinContent(1)*1.0<<" "<<"ggH Bin1="<<hist_ggH_EEC_Energy_Resolution_gen[ptrange-1]->GetBinContent(1)*1.0<<endl;
        ratio_perjet_pairs2->SetBinContent(ptrange,ratio2);
       //hist_tqH_EEC_Energy_Resolution_reco[ptrange-1]=Fillhistwithoutweight(t_tqH,tqH_reco_Energy,TString::Format("hist_tqH_EEC_Energy_Resolution_reco@%d",ptrange),binningEnergy1.numberofbins,binningEnergy1.bin,0,1,2,t_quark_only_select_pjet_reco);
        //LgetNPcorrelation(t_tqH,tqH_gen_Eta,tqH_gen_Phi,tqH_gen_Energy,tqH_gen_Jettag,tqH_gen_Jettag,tqH_gen_JetMtag,2);
        
    }
    
    TCanvas *c1particleratio=new TCanvas();
    gStyle->SetOptStat(0);
    ratio_perjet_pairs2->GetYaxis()->SetRangeUser(0,2);
    ratio_perjet_pairs2->SetLineColor(1);
    ratio_perjet_pairs1->SetLineColor(2);
    ratio_perjet_pairs1->GetYaxis()->SetRangeUser(1,2);
    ratio_perjet_pairs1->Draw();
    ratio_perjet_pairs2->Draw("same");
    TLegend *L1particleratio = new TLegend(.50, .15, .70, .475);
    L1particleratio->AddEntry(ratio_perjet_pairs1,"EEC distribution Bin1 ratio(Q/G)");
    L1particleratio->AddEntry(ratio_perjet_pairs2,"Number of particles per Jet ratio(G/Q)");
    L1particleratio->Draw();
    
    
    
    TCanvas *c1particle=new TCanvas();
    gStyle->SetOptStat(0);
    ggH_perjet_pairs->GetYaxis()->SetRangeUser(0,45);
    ggH_perjet_pairs->SetLineColor(1);
    tqH_perjet_pairs->SetLineColor(2);
    ggH_perjet_cpairs->SetLineColor(3);
    tqH_perjet_cpairs->SetLineColor(4);
    ggH_perjet_npairs->SetLineColor(5);
    tqH_perjet_npairs->SetLineColor(6);
    ggH_perjet_pairs->Draw();
    tqH_perjet_pairs->Draw("same");
    ggH_perjet_cpairs->Draw("same");
    tqH_perjet_cpairs->Draw("same");
    ggH_perjet_npairs->Draw("same");
    tqH_perjet_npairs->Draw("same");
    
    TLegend *L1particle = new TLegend(.50, .15, .70, .375);
    L1particle->AddEntry(tqH_perjet_pairs,"Quark jet");
    L1particle->AddEntry(ggH_perjet_pairs,"Gluon jet");
    L1particle->AddEntry(tqH_perjet_cpairs,"Quark jet charged particle");
    L1particle->AddEntry(ggH_perjet_cpairs,"Gluon jet charged particle");
    L1particle->AddEntry(tqH_perjet_npairs,"Quark jet neutral particle");
    L1particle->AddEntry(ggH_perjet_npairs,"Gluon jet neutral particle");
    L1particle->Draw();
    
    TCanvas *c2particle=new TCanvas();
    gStyle->SetOptStat(0);
    ggH_perjet_cpairs_ratio->GetYaxis()->SetRangeUser(0,1);
    ggH_perjet_cpairs_ratio->SetLineColor(3);
    tqH_perjet_cpairs_ratio->SetLineColor(4);
    ggH_perjet_npairs_ratio->SetLineColor(5);
    tqH_perjet_npairs_ratio->SetLineColor(6);
    
    ggH_perjet_cpairs_ratio->Draw();
    tqH_perjet_cpairs_ratio->Draw("same");
    ggH_perjet_npairs_ratio->Draw("same");
    tqH_perjet_npairs_ratio->Draw("same");
    
    TLegend *L2particle = new TLegend(.50, .15, .70, .375);
    L2particle->AddEntry(tqH_perjet_cpairs_ratio,"Quark jet charged particle");
    L2particle->AddEntry(ggH_perjet_cpairs_ratio,"Gluon jet charged particle");
    L2particle->AddEntry(tqH_perjet_npairs_ratio,"Quark jet neutral particle");
    L2particle->AddEntry(ggH_perjet_npairs_ratio,"Gluon jet neutral particle");
    L2particle->Draw();
    
    
    
    return 0;
     
    
    TCanvas *c1p=new TCanvas();
    gStyle->SetOptStat(0);
    gPad->SetLogx();
    gPad->SetLogy();
    hist_tqH_EEC_Energy_Resolution_gen[0]->GetXaxis()->SetRangeUser(0.0001,0.8);
    hist_tqH_EEC_Energy_Resolution_reco[0]->GetXaxis()->SetRangeUser(0.0001,0.8);
    hist_tqH_EEC_Energy_Resolution_gen[0]->GetYaxis()->SetRangeUser(0.000001,0.01);
    hist_tqH_EEC_Energy_Resolution_reco[0]->GetYaxis()->SetRangeUser(0.000001,0.01);
    hist_tqH_EEC_Energy_Resolution_gen[0]->SetLineColor(1);
    hist_tqH_EEC_Energy_Resolution_reco[0]->SetLineColor(2);
    hist_ggH_EEC_Energy_Resolution_gen[0]->SetLineColor(2);
    hist_tqH_EEC_Energy_Resolution_gen[0]->Draw();
    hist_ggH_EEC_Energy_Resolution_gen[0]->Draw("same");
    //hist_tqH_EEC_Energy_Resolution_reco[0]->Draw("same");
    cout<<222<<endl;
    TLegend *Legend1p = new TLegend(.50, .15, .70, .375);
    Legend1p->AddEntry(hist_tqH_EEC_Energy_Resolution_gen[0],TString::Format("Quark_JetPt@[%.0f,%.0f]",binningEnergy.bin[0],binningEnergy.bin[1]));
    Legend1p->AddEntry(hist_ggH_EEC_Energy_Resolution_gen[0],TString::Format("Gluon_JetPt@[%.0f,%.0f]",binningEnergy.bin[0],binningEnergy.bin[1]));
    //Legend1p->AddEntry(hist_tqH_EEC_Energy_Resolution_reco[0],TString::Format("Gluon_JetPt@[%.0f,%.0f]",binningEnergy.bin[0],binningEnergy.bin[1]));
    for(int i=1;i<binningEnergy.numberofbins;i++){
        hist_tqH_EEC_Energy_Resolution_gen[i]->SetLineColor(2*i+1);
        hist_ggH_EEC_Energy_Resolution_gen[i]->SetLineColor(2*i+2);
        hist_tqH_EEC_Energy_Resolution_reco[i]->SetLineColor(2*i+2);
        
        hist_tqH_EEC_Energy_Resolution_gen[i]->Draw("same");
        hist_ggH_EEC_Energy_Resolution_gen[i]->Draw("same");
        //hist_tqH_EEC_Energy_Resolution_reco[i]->Draw("same");
        //Legend1p->AddEntry(hist_tqH_EEC_Energy_Resolution_gen[i],TString::Format("JetPt@[%.0f,%.0f]",binningEnergy.bin[i],binningEnergy.bin[i+1]));
        Legend1p->AddEntry(hist_tqH_EEC_Energy_Resolution_gen[i],TString::Format("Quark_JetPt@[%.0f,%.0f]",binningEnergy.bin[i],binningEnergy.bin[i+1]));
        Legend1p->AddEntry(hist_ggH_EEC_Energy_Resolution_gen[i],TString::Format("Gluon_JetPt@[%.0f,%.0f]",binningEnergy.bin[i],binningEnergy.bin[i+1]));
        //Legend1p->AddEntry(hist_tqH_EEC_Energy_Resolution_reco[i],TString::Format("Gluon_JetPt@[%.0f,%.0f]",binningEnergy.bin[i],binningEnergy.bin[i+1]));
        
    }
    Legend1p->SetBorderSize(0);
    Legend1p->Draw();
    c1p->SaveAs("a.eps");
    
    /*                           */
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

