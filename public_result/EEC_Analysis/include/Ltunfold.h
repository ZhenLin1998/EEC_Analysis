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


#ifndef Ltunfold_h
#define Ltunfold_h
using namespace std;

bool Lselectedevent(vector<float> *obj_parameter,int i,double xmin,double xmax){
    if(obj_parameter==NULL){
        return true;
    }
    else{
        if(obj_parameter->size()==0){
            return true;
        }
    }
    
    if(obj_parameter->at(i)>xmin&&obj_parameter->at(i)<xmax){
        //cout<<xmin<< " "<<xmax<<" "<<obj_parameter->at(i)<<" "<<endl;
        return true;
    }
    else {
        return false;
    }
}
/*
 bool Lselectedevent_multi_requirements(int n,struct objrequirement *input,int i){
 bool result=true;
 int i1=i;
 for(int j=0;j<n;j++){
 if((input+j)->level!=NULL)
 i=(input+j)->level->at(i1)-1;
 else{
 //cout<<"y"<<(input+j)->obj_parameter->size()<<" "<<i<<endl;
 //cout<<(input+j)->obj_parameter->at(i)<<" "<<(input+j)->obj_parameter->size()<<" "<<i<<endl;
 }
 //cout<<"a="<<i<<" "<<i1<<" "<<(input+j)->obj_parameter->size()<<"  "<<(input+j)->level->size() <<" "<<j<<endl;
 
 if((input+j)->obj_parameter!=NULL){
 //Fillhistwithweightcout<<(input+j)->obj_parameter->size()<<endl;
 }
 //cout<<"a1="<<i<<" "<<(input+j)->obj_parameter->size()<<endl;
 result=result&&Lselectedevent((input+j)->obj_parameter, i,(input+j)->xmin,(input+j)->xmax);
 //cout<<"a1="<<i<<endl;
 //cout<<"a="<<i<<" "<<i1<<endl;
 i=i1;
 }
 return result;
 }
 */
bool Lselectedevent_multi_requirements(int n,struct objrequirement *input,int i){
    bool result=true;
    int i1=i;
    for(int j=0;j<n;j++){
        if((input+j)->level!=NULL)
            i=(input+j)->level->at(i1);
        else{
            //cout<<"y"<<(input+j)->obj_parameter->size()<<" "<<i<<endl;
            //cout<<(input+j)->obj_parameter->at(i)<<" "<<(input+j)->obj_parameter->size()<<" "<<i<<endl;
        }
        //cout<<"a="<<i<<" "<<i1<<" "<<(input+j)->obj_parameter->size()<<"  "<<(input+j)->level->size() <<" "<<j<<endl;
        
        if((input+j)->obj_parameter!=NULL){
            //Fillhistwithweightcout<<(input+j)->obj_parameter->size()<<endl;
        }
        //cout<<"a1="<<i<<" "<<(input+j)->obj_parameter->size()<<endl;
        result=result&&Lselectedevent((input+j)->obj_parameter, i,(input+j)->xmin,(input+j)->xmax);
        //cout<<"a1="<<i<<endl;
        //cout<<"a="<<i<<" "<<i1<<endl;
        i=i1;
        if(result){
            
        }
        else{
            break;
        }
    }
    return result;
}


TH1 *histogram1Dcreation(TTree *t,vector<float> *obj1,vector<float> *obj2,TString name_of_histogram,TUnfoldBinning *obj_binning,TUnfoldBinning *obj_distribution,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected){
    
    TH1 *result= obj_binning->CreateHistogram(name_of_histogram);
    for (Int_t k = t->GetEntries()*entrymin; k < t->GetEntries()*entrymax ; k++)
    {
        t->GetEntry(k);
        for (Int_t i = 0; i < obj1->size(); i++)
        {
            
            if(Lselectedevent_multi_requirements(numberofrequirements,particle_selected,i)){
                Int_t binNumber = obj_distribution->GetGlobalBinNumber(obj1->at(i), obj2->at(i));
                result->Fill(binNumber);
            }
        }
    }
    return result;
}

TH2 *histogram2Dcreation(TTree *t,vector<float> *obj1,vector<float> *obj2,vector<float> *obj3,vector<float> *obj4,TString name_of_histogram,TUnfoldBinning *obj_binning,TUnfoldBinning *obj_distribution,TUnfoldBinning *obj_binning1,TUnfoldBinning *obj_distribution1,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected){
    
    TH2 *result= TUnfoldBinning::CreateHistogramOfMigrations(obj_binning, obj_binning1, name_of_histogram);
    for (Int_t k = t->GetEntries()*entrymin; k < t->GetEntries()*entrymax ; k++)
    {
        t->GetEntry(k);
        for (Int_t i = 0; i < obj1->size(); i++)
        {
            if(Lselectedevent_multi_requirements(numberofrequirements,particle_selected,i)){
                Int_t binNumber1 = obj_distribution->GetGlobalBinNumber(obj1->at(i), obj2->at(i));
                Int_t binNumber2 = obj_distribution1->GetGlobalBinNumber(obj3->at(i), obj4->at(i));
                result->Fill(binNumber1,binNumber2);
            }
        }
    }
    return result;
}


TUnfoldDensity Ltunfold(TH2 *histMCGenReco,TH1 *histMCRecoBackground,TH1 *histDataReco,TUnfoldBinning * generatorBinning,TUnfoldBinning * detectorBinning){
    TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;
    
    // basic choice of regularisation scheme:
    //    curvature (second derivative)
    TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
    
    // density flags
    TUnfoldDensity::EDensityMode densityFlags =
    TUnfoldDensity::kDensityModeBinWidth;
    
    // detailed steering for regularisation
    const char *REGULARISATION_DISTRIBUTION = 0;
    const char *REGULARISATION_AXISSTEERING = "*[B]";
    
    // set up matrix of migrations
    TUnfoldDensity unfold(histMCGenReco, TUnfold::kHistMapOutputHoriz,
                          regMode, constraintMode, densityFlags,
                          generatorBinning, detectorBinning,
                          REGULARISATION_DISTRIBUTION,
                          REGULARISATION_AXISSTEERING);
    //test
    TH2D *inputEmatrix=
    detectorBinning->CreateErrorMatrixHistogram("input_covar",true);
    for(int i=1;i<=inputEmatrix->GetNbinsX();i++) {
        Double_t e=histDataReco->GetBinError(i);
        inputEmatrix->SetBinContent(i,i,e*e);
    }
    unfold.SubtractBackground(histMCRecoBackground,"bgr");
    unfold.SetInput(histDataReco,0.0,0.0,inputEmatrix);
    
    Int_t nScan = 30;
    Double_t tauMin = 0.0;
    Double_t tauMax = 0.0;
    Int_t iBest;
    TSpline *logTauX, *logTauY;
    TGraph *lCurve=0;
    
    
    TSpline *rhoLogTau=0;
    
    // for determining tau, scan the correlation coefficients
    // correlation coefficients may be probed for all distributions
    // or only for selected distributions
    // underflow/overflow bins may be included/excluded
    //
    const char *SCAN_DISTRIBUTION="signal";
    const char *SCAN_AXISSTEERING=0;
    
    iBest=unfold.ScanTau(50,tauMin,tauMax,&rhoLogTau,TUnfoldDensity::kEScanTauRhoAvg);
    
    
    return unfold;
    
}


TH1D *Ltransfer_2Daxisto1D(TH1 *input,int ndrbins,int nweightbins,double *drbins,double *weightbins){
    const char *hist_name=input->GetName();
    TH1D *output=new TH1D(hist_name,"",ndrbins,drbins);
    double Weight=0;
    double WeightError=0;
    double intweight=0;
    double interror=0;
    for(int i=1;i<=ndrbins;i++){
        intweight=0;
        interror=0;
        for(int j=1;j<=nweightbins;j++){
            Weight=input->GetBinContent((j-1)*(ndrbins)+i)*(weightbins[j-1]+weightbins[j])/2/(drbins[i]-1.0*drbins[i-1]);
            WeightError=input->GetBinError((j-1)*(ndrbins)+i)*(weightbins[j-1]+weightbins[j])/2/(drbins[i]-1.0*drbins[i-1])+input->GetBinContent((j-1)*(ndrbins)+i)*(-1*weightbins[j-1]+weightbins[j]);
            intweight=intweight+Weight;
            interror=interror+WeightError;
            //cout<<double(WeightError/Weight)<<endl;
        }
        output->SetBinContent(i,intweight);
        output->SetBinError(i,interror);
    }
    return output;
}

TH1 *Lsub_histogram(TH1 *input1,TH1 *input2,TString name_of_histogram,double weight=-1){
    
    int nxbins=input1->GetNbinsX();
    TH1 *output= (TH1 *)input1->Clone();
    double a;
    for(int i=1;i<=nxbins;i++){
        a=input1->GetBinContent(i)+weight*input2->GetBinContent(i);
        output->SetBinContent(i,a);
    }
    return output;
    
}
TH1D *Lsub_histogram_TH1D(TH1D *input1,TH1D *input2,TString name_of_histogram,double weight=-1){
    
    int nxbins=input2->GetNbinsX();
    TH1D *output= (TH1D *)input1->Clone();
    double a;
    for(int i=1;i<=nxbins;i++){
        a=input1->GetBinContent(i)+weight*input2->GetBinContent(i);
        output->SetBinContent(i,a);
    }
    return output;
    
}

TH1 *D1histogram1Dcreation(TTree *t,vector<float> *obj1,TString name_of_histogram,TUnfoldBinning *obj_binning,TUnfoldBinning *obj_distribution,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected){
    
    TH1 *result= obj_binning->CreateHistogram(name_of_histogram);
    for (Int_t k = t->GetEntries()*entrymin; k < t->GetEntries()*entrymax ; k++)
    {
        t->GetEntry(k);
        for (Int_t i = 0; i < obj1->size(); i++)
        {
            if(Lselectedevent_multi_requirements(numberofrequirements,particle_selected,i)){
                
                int num=obj_distribution->GetGlobalBinNumber(obj1->at(i));
                result->Fill(num);
            }
        }
    }
    return result;
}

TH2 *D1histogram2Dcreation(TTree *t,vector<float> *obj1,vector<float> *obj2,TString name_of_histogram,TUnfoldBinning *obj_binning,TUnfoldBinning *obj_distribution,TUnfoldBinning *obj_binning1,TUnfoldBinning *obj_distribution1,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected){
    
    TH2 *result= TUnfoldBinning::CreateHistogramOfMigrations(obj_binning, obj_binning1, name_of_histogram);
    for (Int_t k = t->GetEntries()*entrymin; k < t->GetEntries()*entrymax ; k++)
    {
        t->GetEntry(k);
        for (Int_t i = 0; i < obj1->size(); i++)
        {
            if(Lselectedevent_multi_requirements(numberofrequirements,particle_selected,i)){
                int num1=obj_distribution->GetGlobalBinNumber(obj1->at(i));
                int num2=obj_distribution1->GetGlobalBinNumber(obj2->at(i));
                result->Fill(num1,num2);
            }
        }
    }
    return result;
}
TH1 *Ltransfer1DTunfold(TH1 *input1,TH1 *input2,int ibin){
    
    double weight=0;;
    double weighterror=0;
    int nweightbins=input2->GetNbinsX();
    for(int i=1;i<=nweightbins;i++){
        weight=weight+input1->GetBinContent(i)*input1->GetBinCenter(i);
        weighterror=weighterror+input1->GetBinContent(i)*input1->GetBinWidth(i);
    }
    
    input1->SetBinContent(ibin,weight);
    return input1;
}

double Lintweight(TH1 *input,double *xbins1){
    int nweightbins=input->GetNbinsX();
    double weight=0;
    for(int i=1;i<=nweightbins;i++){
        weight=weight+input->GetBinContent(i)*(xbins1[i]+xbins1[i-1])/2;
    }
    cout<<weight<<endl;
    return weight;
}
TH1D *Lc_1Dhist(const char *hist_name,int ndrbins,double *drbins,double *intweight){
    
    TH1D *output=new TH1D(hist_name,"",ndrbins,drbins);
    for(int i=1;i<=ndrbins;i++){
        output->SetBinContent(i,intweight[i]/(drbins[i]-drbins[i-1]));
        cout<<intweight[i]/(drbins[i]-drbins[i-1])<<endl;
    }
    return output;
}

TH1D *Lcombine_hist(TH1D *h1,TH1D *h2){
    
    int nbins1D=h1->GetNbinsX();
    int nbins2D=h2->GetNbinsX();
    int nbins=nbins1D+nbins2D;
    double bin[nbins+1];
    for(int i=0;i<nbins;i++){
        if(i<nbins1D){
            bin[i]=h1->GetBinLowEdge(i+1);
        }
        else{
            cout<<i+1-nbins1D<<" "<<nbins2D<<" "<<h2->GetBinLowEdge(i+1-nbins1D)<<endl;
            
            bin[i]=h2->GetBinLowEdge(i+1-nbins1D);
        }
    }
    
    bin[nbins]=h2->GetBinLowEdge(nbins2D)+h2->GetBinWidth(nbins2D);
    
    TH1D *EEC_hist=new TH1D("EEC_distribution","",nbins,bin);
    for(int i=1;i<=nbins;i++){
        if(i<=nbins1D){
            EEC_hist->SetBinContent(i,h1->GetBinContent(i));
            EEC_hist->SetBinError(i,h1->GetBinError(i));
        }
        
        else{
            EEC_hist->SetBinContent(i,h2->GetBinContent(i-nbins1D));
            EEC_hist->SetBinError(i,h2->GetBinError(i-nbins1D));
        }
    }
    return EEC_hist;
    cout<<"Compelete combination"<<endl;
}


struct binninginfo LCloneBin(TH1D *input){
    
    int nbins=input->GetNbinsX();
    struct binninginfo output;
    cout<<1<<endl;
    
    for(int i=0;i<nbins;i++){
        
        output.bin[i]=input->GetBinLowEdge(i+1);
        cout<<1<<endl;
        
    }
    output.bin[nbins]=input->GetBinLowEdge(nbins)+input->GetBinWidth(nbins);
    
    cout<<2<<endl;
    output.numberofbins=nbins;
    return output;
}
TH1D *GetNormalized(TH1D *input){
    int nbins=input->GetNbinsX();
    TH1D *output= (TH1D *)input->Clone();
    double sum=0;
    for(int i=1;i<=nbins;i++){
        sum=input->GetBinContent(i)+sum;
        cout<<input->GetBinContent(i)<<endl;
    }
    for(int i=1;i<=nbins;i++){
        output->SetBinContent(i,input->GetBinContent(i)*1.0/sum);
        output->SetBinError(i,input->GetBinError(i)*1.0/sum);
        cout<<i<<" "<<output->GetBinContent(i)<<" "<<sum<<endl;
        
    }
    return output;
}



TH1D *Fillhistwithweight(TTree *t,vector<float> *obj1,vector<float> *obj2,TString name_of_histogram,int nbins,double *bin,double entrymin,double entrymax,int numberofrequirements=0,struct objrequirement *particle_selected=NULL,bool ifdensity=true){
    TH1D *result= new TH1D(name_of_histogram,"",nbins,bin);
    //int pairnum=0;
    //int jetnum=0;
    for (Int_t k = t->GetEntries()*entrymin; k < t->GetEntries()*entrymax ; k++)
    {
        t->GetEntry(k);
        for (Int_t i = 0; i < obj1->size(); i++)
        {
            //cout<<"i="<<i<<endl;
            if(Lselectedevent_multi_requirements(numberofrequirements,particle_selected,i)){
                result->Fill(obj1->at(i),obj2->at(i));
            }
        }
        //pairnum=pairnum+obj1->size();
        //if(particle_selected!=NULL)
        //jetnum=jetnum+particle_selected->obj_parameter->size();
    }
    if(ifdensity){
        for(int i=1;i<=nbins;i++){
            result->SetBinContent(i,result->GetBinContent(i)/(bin[i]-bin[i-1]));
            result->SetBinError(i,result->GetBinError(i)/(bin[i]-bin[i-1]));
        }
    }
    //cout<<"jetnum="<<jetnum<<"  "<<"pairnum="<<pairnum<<endl;
    cout<<"complete fill hist with weight"<<endl;
    return result;
    
}

TH1D *Fillhistwithweight_tag(TTree *t,vector<float> *obj1,vector<float> *obj2,vector<float> *obj_tag,TString name_of_histogram,int nbins,double *bin,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected,bool ifdensity=true){
    TH1D *result= new TH1D(name_of_histogram,"",nbins,bin);
    for (Int_t k = t->GetEntries()*entrymin; k < t->GetEntries()*entrymax ; k++)
    {
        t->GetEntry(k);
        for (Int_t i = 0; i < obj1->size(); i++)
        {
            if(obj1->size()*obj2->size()==0)
                continue;
            int jetid=obj_tag->at(i);
            //cout<<"i="<<i<<endl;
            if(Lselectedevent_multi_requirements(numberofrequirements,particle_selected,i)){
                result->Fill(obj1->at(jetid),obj2->at(jetid));
            }
        }
    }
    if(ifdensity){
        for(int i=1;i<=nbins;i++){
            result->SetBinContent(i,result->GetBinContent(i)/(bin[i]-bin[i-1]));
            result->SetBinError(i,result->GetBinError(i)/(bin[i]-bin[i-1]));
        }
    }
    cout<<"complete fill hist with weight"<<endl;
    return result;
    
}
TH1D *Fillhistwithmultiweight_tag(TTree *t,vector<float> *obj1,int numberofweights,struct objweight *particle_weight,TString name_of_histogram,int nbins,double *bin,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected,bool ifdensity=true){
    TH1D *result= new TH1D(name_of_histogram,"",nbins,bin);
    int round=0;
    for (Int_t k = t->GetEntries()*entrymin; k < t->GetEntries()*entrymax ; k++)
    {
        
        int process_percent=k*100/t->GetEntries();
        if(process_percent==round*10){
            cout<<process_percent<<"%"<<endl;
            round++;
        }
            
        
        
        t->GetEntry(k);
        for (Int_t i = 0; i < obj1->size(); i++)
        {
            bool ifvector0=true;
            for(Int_t j = 0; j < numberofweights; j++){
                ifvector0=ifvector0&&(((particle_weight+j)->obj_parameter->size())==0);
            }
            ifvector0=ifvector0&&(obj1->size()==0);
            if(ifvector0)
            continue;
            int i1;
            double weight=1;
            for(Int_t j = 0; j < numberofweights; j++){
                i1=i;
                if((particle_weight+j)->level!=NULL){
                i1=(particle_weight+j)->level->at(i);
                    //cout<<i1<<endl;
                    //cout<<(particle_weight+j)->obj_parameter->at(i1)<<endl;
                }
                //cout<<j<<"j"<<(particle_weight+j)->obj_parameter->at(i1)<<endl;
                //cout<<i<<" "<<obj1->size()<<" "<<j<<" "<<i1<<" "<<(particle_weight+j)->obj_parameter->size();
                weight=(particle_weight+j)->obj_parameter->at(i1)*weight;
            }
            
            //cout<<"i="<<i<<endl;
            if(Lselectedevent_multi_requirements(numberofrequirements,particle_selected,i)){
                result->Fill(obj1->at(i),weight);
                //cout<<obj1->at(i)<<" "<<weight<<endl;
            }
        }
    }
    //cout<<result->GetBinContent(5)<<"ss"<<endl;
    if(ifdensity){
        for(int i=1;i<=nbins;i++){
            result->SetBinContent(i,result->GetBinContent(i)/(bin[i]-bin[i-1]));
            result->SetBinError(i,result->GetBinError(i)/(bin[i]-bin[i-1]));
        }
    }
    cout<<"complete fill hist with weight"<<endl;
    //cout<<result->GetBinContent(1)<<"ss"<<endl;
    
    return result;
    
}



TH1D *Fillhistwithoutweight(TTree *t,vector<float> *obj1,TString name_of_histogram,int nbins,double *bin,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected,bool ifdensity=true){
    cout<<"Begin filling hist without weight"<<endl;
    TH1D *result= new TH1D(name_of_histogram,"",nbins,bin);
    for (Int_t k = t->GetEntries()*entrymin; k < t->GetEntries()*entrymax ; k++)
    {
        t->GetEntry(k);
        for (Int_t i = 0; i < obj1->size(); i++)
        {
            if(Lselectedevent_multi_requirements(numberofrequirements,particle_selected,i)){
                result->Fill(obj1->at(i));
            }
        }
    }
    if(ifdensity){
        for(int i=1;i<=nbins;i++){
            result->SetBinContent(i,result->GetBinContent(i)/(bin[i]-bin[i-1]));
            result->SetBinError(i,result->GetBinError(i)/(bin[i]-bin[i-1]));
        }
    }
    return result;
}


TH1D *Fillhistwithoutweight_tag(TTree *t,vector<float> *obj1,vector<float> *obj2,TString name_of_histogram,int nbins,double *bin,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected,bool ifdensity=true){
    TH1D *result= new TH1D(name_of_histogram,"",nbins,bin);
    int count2=0;
    int count1=0;
    for (Int_t k = t->GetEntries()*entrymin; k < t->GetEntries()*entrymax ; k++)
    {
        t->GetEntry(k);
        for (Int_t i = 0; i < obj1->size(); i++)
        {
            count1++;
            if(obj2->size()==0){
                count2++;
                
                cout<<count1<<" "<<count2<<endl;
                continue;
            }
            int jetid=obj1->at(i);
            // cout<<obj2->size()<<" "<<obj1->at(i)-1<<"x"<<particle_selected->obj_parameter->size()<<endl;
            if(Lselectedevent_multi_requirements(numberofrequirements,particle_selected,i)){
                //cout<<obj2->size()<<" "<<obj1->at(i)-1<<endl;
                result->Fill(obj2->at(jetid));
            }
        }
    }
    if(ifdensity){
        for(int i=1;i<=nbins;i++){
            result->SetBinContent(i,result->GetBinContent(i)/(bin[i]-bin[i-1]));
            result->SetBinError(i,result->GetBinError(i)/(bin[i]-bin[i-1]));
        }
    }
    return result;
}



TH1D *Lreweighting(TTree *t,TTree *t1,vector<float> *obj,vector<float> *obj1,vector<float> *obj_weight_pair,TString name_of_histogram,vector<float> *obj_weight,vector<float> *obj_weight1,int nbins,double *bin,int nbins_weight,double *bin_weight,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected,struct objrequirement *particle_selected_jet,struct objrequirement *particle_selected1,struct objrequirement *particle_selected1_jet){
    
    cout<<1<<endl;
    TH1D *temp;
    double weight_num;
    cout<<"begin reweighting"<<endl;
    TH1D *result= new TH1D(name_of_histogram,"",nbins,bin);
    TH1D *weight=Fillhistwithoutweight(t,obj_weight,"",nbins_weight,bin_weight,entrymin,entrymax, numberofrequirements,particle_selected_jet);
    cout<<"1"<<endl;
    
    TH1D *weight1=Fillhistwithoutweight(t1,obj_weight1,"",nbins_weight,bin_weight,entrymin,entrymax, numberofrequirements,particle_selected1_jet);
    cout<<"1"<<endl;
    for(int i=1;i<=nbins_weight;i++){
        particle_selected[numberofrequirements]={obj_weight_pair,weight->GetBinLowEdge(i),weight->GetBinLowEdge(i)+weight->GetBinWidth(i),NULL};
        cout<<"i="<<i<<endl;
        weight_num=weight1->GetBinContent(i)*1.0/weight->GetBinContent(i);
        if(obj1==NULL)
            temp=Fillhistwithoutweight(t,obj,"",nbins,bin,entrymin,entrymax, numberofrequirements+1,particle_selected);
        else
            temp=Fillhistwithweight(t,obj,obj1,"",nbins,bin,entrymin,entrymax, numberofrequirements+1,particle_selected);
        result=Lsub_histogram_TH1D(result,temp,"",weight_num);
        
    }
    cout<<"complete reweighting"<<endl;
    return result;
    
    
    
}


TH1D *Lreweighting_tag(TTree *t,TTree *t1,vector<float> *obj,vector<float> *obj1,vector<float> *obj_weight_pair,vector<float> *obj_weight_pair_tag,TString name_of_histogram,vector<float> *obj_weight,vector<float> *obj_weight_tag,vector<float> *obj_weight1,vector<float> *obj_weight1_tag,int nbins,double *bin,int nbins_weight,double *bin_weight,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected,struct objrequirement *particle_selected_jet,struct objrequirement *particle_selected1,struct objrequirement *particle_selected1_jet,TH1D *Extra_weight=NULL){
    
    cout<<1<<endl;
    TH1D *temp;
    double weight_num;
    cout<<"begin reweighting"<<endl;
    TH1D *result= new TH1D(name_of_histogram,"",nbins,bin);
    TH1D *weight=Fillhistwithoutweight_tag(t,obj_weight_tag,obj_weight,"",nbins_weight,bin_weight,entrymin,entrymax, numberofrequirements,particle_selected_jet);
    cout<<"1"<<endl;
    
    TH1D *weight1=Fillhistwithoutweight_tag(t1,obj_weight1_tag,obj_weight1,"",nbins_weight,bin_weight,entrymin,entrymax, numberofrequirements,particle_selected1_jet);
    cout<<"1"<<endl;
    for(int i=1;i<=nbins_weight;i++){
        particle_selected[numberofrequirements]={obj_weight_pair,weight->GetBinLowEdge(i),weight->GetBinLowEdge(i)+weight->GetBinWidth(i),obj_weight_pair_tag};
        cout<<"i="<<i<<endl;
        if(Extra_weight!=NULL)
            weight_num=weight1->GetBinContent(i)*1.0/weight->GetBinContent(i)*Extra_weight->GetBinContent(i);
        else
            weight_num=weight1->GetBinContent(i)*1.0/weight->GetBinContent(i);
        if(obj1==NULL)
            temp=Fillhistwithoutweight_tag(t,obj_weight_pair_tag,obj,"",nbins,bin,entrymin,entrymax, numberofrequirements+1,particle_selected);
        else
            temp=Fillhistwithweight(t,obj,obj1,"",nbins,bin,entrymin,entrymax, numberofrequirements+1,particle_selected);
        result=Lsub_histogram_TH1D(result,temp,"",weight_num);
        
    }
    cout<<"complete reweighting"<<endl;
    return result;
}



#endif /* Ltunfold_h */
