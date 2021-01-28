#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <TError.h>
#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TFitter.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TSystem.h>
#include "TUnfoldDensity.h"
#include "Ltunfold.h"
#include "data_structure.h"
#ifndef Lbinning_h
#define Lbinning_h

using namespace std;


struct binninginfo Lbinning(TTree *t,TString object,vector<float> *recoobj,vector<float> *genobj,int binningtool,double xmin,double xmax,double entrymin,double entrymax,double minpurity,int bins=30,vector<float> *obj_parameter=NULL,double xmin1=0, double xmax1=0,vector<float> *obj_parameter1=NULL,double xmin2=0, double xmax2=0)
{

    
    struct binninginfo output;
    cout<<"start binning with function:"<<binningtool<<endl;

    if(binningtool==0){
        double objmax=xmax,objmin=xmin;
        double minstability;
        minstability=minpurity;
        int Nbins_maximum=200;
        
        
        double genobjbins[Nbins_maximum],recoobjbins[Nbins_maximum],stability[Nbins_maximum],purity[Nbins_maximum];
        genobjbins[0]=objmin,recoobjbins[0]=objmin;
        int Nobjbins;
        int nflag=8;

        int pgencount=0,precocount=0,sgencount=0,srecocount=0,gennbin=1,reconbin=1;
        double binmax=(xmax-xmin)/300*(nflag+2);
        double stepobjmax=binmax/nflag;//1/250
        double stepobj=(xmax-xmin)/1000;
        
        double binwidth=0;
        genobjbins[gennbin]=genobjbins[0]+stepobj;
        recoobjbins[reconbin]=recoobjbins[0]+stepobj;
        int nstep=0;

        while(genobjbins[gennbin]<objmax&&recoobjbins[reconbin]<objmax){
            precocount=0;
            pgencount=0;
            srecocount=0;
            sgencount=0;
            nstep++;
            for (int i=entrymin*t->GetEntries();i<entrymax*t->GetEntries();i++){
                t->GetEntry(i);
                

                for (Int_t j = 0; j < genobj->size(); j++)
                {
                    if(Lselectedevent(obj_parameter,j,xmin1,xmax1)&&Lselectedevent(obj_parameter1,j,xmin2,xmax2)){
                    if( genobj->at(j)<genobjbins[gennbin]&&genobj->at(j)>genobjbins[gennbin-1]){
                        pgencount++;
                        if(recoobj->at(j)<recoobjbins[reconbin]&&recoobj->at(j)>recoobjbins[reconbin-1]){
                            precocount++;
                            
                        }
                    }
                    
                    
                    if( recoobj->at(j)<recoobjbins[reconbin]&&recoobj->at(j)>recoobjbins[reconbin-1]){
                        srecocount++;
                        if(genobj->at(j)<genobjbins[gennbin]&&genobj->at(j)>genobjbins[gennbin-1]){
                            sgencount++;
                        }
                    }
                    }
                }
            }
            
            purity[gennbin-1]=1.0*precocount/pgencount;
            stability[reconbin-1]=1.0*sgencount/srecocount;
            /*
             if(purity[gennbin-1]<minpurity){
             if(stability[reconbin-1]>=purity[gennbin-1]){
             recoobjbins[reconbin]+=stepobj;
             }
             }
             if(stability[reconbin-1]<minstability){
             if(stability[reconbin-1]<purity[gennbin-1]){
             genobjbins[gennbin]+=stepobj;
             }
             }
             */
            genobjbins[gennbin]+=stepobj;
            recoobjbins[reconbin]+=stepobj;
            binwidth=recoobjbins[reconbin]-recoobjbins[reconbin-1];
            /*
            if(binwidth>binmax){
               cout<<precocount<<"    "<<pgencount<<endl;
               cout<<srecocount<<"    "<<sgencount<<endl;
               cout<<"new bin defined by biggest binlength"<<genobjbins[gennbin]<<" "<<recoobjbins[reconbin]<<" "<<purity[gennbin-1]<<" "<<stability[reconbin-1]<<endl;
                gennbin++;
                reconbin++;
                
                stepobj=binmax;
                if(stepobj>stepobjmax) stepobj=stepobjmax;
                recoobjbins[reconbin]=recoobjbins[reconbin-1]+stepobj;
                genobjbins[gennbin]=genobjbins[gennbin-1]+stepobj;
                nstep=1;
                continue;
            }*/
            cout<<nstep<<endl;
            cout<<genobjbins[gennbin]<<" "<<recoobjbins[reconbin]<<"  "<<purity[gennbin-1]<<" "<<stability[reconbin-1]<<endl;
            if(purity[gennbin-1]>=minpurity&&stability[reconbin-1]>=minstability){
                
                if(nstep>=nflag){
                    cout<<precocount<<"    "<<pgencount<<endl;
                    cout<<srecocount<<"    "<<sgencount<<endl;
                    cout<<"new bin defined by purity and stability"<<genobjbins[gennbin]<<" "<<recoobjbins[reconbin]<<"  "<<purity[gennbin-1]<<" "<<stability[reconbin-1]<<endl;
                    gennbin++;
                    reconbin++;
                    cout<<gennbin<<endl;
                    
                    stepobj=(genobjbins[gennbin-1]-genobjbins[gennbin-2])/nflag;
                    if(stepobj>stepobjmax) stepobj=stepobjmax;
                    recoobjbins[reconbin]=recoobjbins[reconbin-1]+stepobj;
                    genobjbins[gennbin]=genobjbins[gennbin-1]+stepobj;
                    nstep=1;
                }
                else{
                    stepobj=stepobj/(nflag-nstep)/2;
                    if(stepobj>stepobjmax) stepobj=stepobjmax;
                    genobjbins[gennbin]=genobjbins[gennbin-1]+stepobj;
                    recoobjbins[reconbin]=recoobjbins[reconbin-1]+stepobj;
                //    cout<<gennbin<<"    "<<reconbin<<"  "<<stepobj<<endl;
                    nstep=1;
                }
                
            }
            
        }
        genobjbins[gennbin]=objmax;
        recoobjbins[reconbin]=objmax;
        ofstream outobjbins(object+TString::Format("bins_%.2f-%.2f_%.2f.txt",minpurity,xmin,xmax));
        output.numberofbins=gennbin;
        output.numberofhalfbins=(reconbin-reconbin%2)/2;
    
        for(Int_t i=0;i<=gennbin;i++){
            outobjbins<<gennbin<<"  "<<output.numberofhalfbins<<"   "<<i<<" "<<genobjbins[i]<<"  "<<recoobjbins[i-i%2]<<endl;
            output.bin[i]=genobjbins[i];
            output.halfbin[(i-i%2)/2]=recoobjbins[i-i%2];
        }
        outobjbins.close();
        cout<<"complete binning with function:"<<binningtool<<endl;
        return output;
    }
    if(binningtool==1){
        int binnumx=1000000;
        
        //int nxbins = a1[round];
        Double_t xbins[bins + 1];
        xbins[0] = xmin;
        xbins[bins] = xmax;

        TH1D *hx1 = new TH1D("", "", binnumx, xmin, xmax);
        for (Int_t k = 0; k < t->GetEntries() / 3; k++)
        {
            t->GetEntry(k);
            for (Int_t i = 0; i < recoobj->size(); i++)
            {
                 if(Lselectedevent(obj_parameter,i,xmin1,xmax1)&&Lselectedevent(obj_parameter1,i,xmin2,xmax2)){
                hx1->Fill(recoobj->at(i));
                 }
            }
        }
        double scale = 1 / hx1->Integral();
        hx1->Scale(scale);
        TH1 *hx1c = hx1->GetCumulative(kTRUE, "_cumulative");
        
        double x,y;
        int num;
        TGraph *gx1 = new TGraph();
        for (Int_t k = 1; k <= binnumx; k++)
        {
            x = hx1c->GetBinContent(k);
            y = hx1c->GetBinCenter(k);
            num++;
            gx1->SetPoint(num, x, y);
        }
        for (Int_t i = 1; i <= bins-1; i++)
        {
            xbins[i] = gx1->Eval(1.0 / bins * i);
        }
        ofstream outobjbins(object+TString::Format("bins_%d.txt",bins));
        output.numberofbins=bins;
        output.numberofhalfbins=(bins-bins%2)/2;
        
        for(Int_t i=0;i<=bins;i++){
            outobjbins<<bins<<"  "<<output.numberofhalfbins<<"   "<<i<<" "<<xbins[i]<<"  "<<xbins[i-i%2]<<endl;
            output.bin[i]=xbins[i];
            output.halfbin[(i-i%2)/2]=xbins[i-i%2];
        }
        outobjbins.close();
        cout<<"complete binning with function:"<<binningtool<<endl;
        return output;
    }
     if(binningtool==2){
        int binnumx=1000000;

        //int nxbins = a1[round];
        Double_t xbins[bins + 1];
        xbins[0] = xmin;
        xbins[bins] = xmax;

    
        for (Int_t i = 0; i <= bins; i++)
        {
            xbins[i] = xmin+(xmax-xmin)/bins*i;
        }
        ofstream outobjbins(object+TString::Format("2bins_%d_%.2f_%.2f.txt",bins,xmin,xmax));
        output.numberofbins=bins;
        output.numberofhalfbins=(bins-bins%2)/2;
        
        for(Int_t i=0;i<=bins;i++){
            outobjbins<<bins<<"  "<<output.numberofhalfbins<<"   "<<i<<" "<<xbins[i]<<"  "<<xbins[i-i%2]<<endl;
            
            output.bin[i]=xbins[i];
            output.halfbin[(i-i%2)/2]=xbins[i-i%2];
        }
        outobjbins.close();
        cout<<"complete binning with function:"<<binningtool<<endl;
        return output;

    }
    if(binningtool==3){
        int binnumx=1000000;
int nbins=bins;
        //int nxbins = a1[round];
        Double_t xbins[bins + 1];
        
        xbins[0] = xmin;
        xbins[bins] = xmax;
        int init=1;
        if(xmin==0){
            xbins[1]=0.0001;
            xmin=0.0001;
            init=2;
        }
        Double_t xlogmin=TMath::Log10(xmin);
        Double_t xlogmax=TMath::Log10(xmax);
        Double_t dlogx=(xlogmax-xlogmin)/(nbins);
        
        for(int i=init;i<=nbins;i++){
            Double_t xlog=xlogmin+i*dlogx;
            xbins[i]=TMath::Exp(TMath::Log(10) * xlog);
        }
        ofstream outobjbins(object+TString::Format("2bins_%d_%.2f_%.2f.txt",bins,xmin,xmax));
        output.numberofbins=bins;
        output.numberofhalfbins=(bins-bins%2)/2;
        
        for(Int_t i=0;i<=bins;i++){
            outobjbins<<bins<<"  "<<output.numberofhalfbins<<"   "<<i<<" "<<xbins[i]<<"  "<<xbins[i-i%2]<<endl;
            
            output.bin[i]=xbins[i];
            output.halfbin[(i-i%2)/2]=xbins[i-i%2];
        }
        outobjbins.close();
        cout<<"complete binning with function:"<<binningtool<<endl;
        return output;

    }
    return output;
}
struct binninginfo Lgethalfbinning(struct binninginfo input){
    
    
    struct binninginfo output;
    int bins=input.numberofbins;
    int halfbins=input.numberofhalfbins;
    output.numberofbins=(bins-bins%2)/2;
    output.numberofhalfbins=(halfbins-halfbins%2)/2;
    for(int i=0;i<=bins;i++){
        output.bin[(i-i%2)/2]=input.bin[i-i%2];
    }
    for(int i=0;i<=halfbins;i++){
        output.halfbin[(i-i%2)/2]=input.halfbin[i-i%2];
    }
    
    
    return output;
    
}




struct binninginfo Lreadoutbinfile(TString path){
    
    ifstream in(path);
    string line;
    struct binninginfo output;
    int i;

    while(getline(in,line)){
        istringstream inl(line);
        inl>>output.numberofbins>>output.numberofhalfbins>>i>>output.bin[i]>>output.halfbin[(i-i%2)/2];
    }
    
    return output;
    
}
#endif
