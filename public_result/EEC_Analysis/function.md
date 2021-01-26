Updated at 11.19


Part   ①: binning

struct binninginfo{
int numberofbins;
double bin[200];
double halfbin[100];
};

function    1: get out of the binning of the measurement by several binning tools
struct binninginfo  Lbinning(TString "path of root file " ,TString "binningobject" , int binningtool , double xmin, double x max,double entrymin, double entrymax, double numberofbins,vector<float> *obj_parameter=NULL,double xmin1=0, double xmax1=0,vector<float> *obj_parameter1=NULL,double xmin2=0, double xmax2=0);


binningtool =
0: refers to  binning by purity and stability (both of them > 0.2), can't determin a specific value of numberofbins
1: refers to binning  by getting equal events in each bin
2:refers to binning by equal  numbers of bin's width
xmin & xmax: determin the boundaries of bins
entrymin & entrymax: determin the entries range in root file
numberofbins=
default=30;
this parameter can only effect the binningtool=1;

obj_parameter & xmin1,xmax1: determined the Event selected (at most two selected conditions)
obj_parameter                         : refers to the parameter that need to be selected

We can get the output: binningobject
1. binningobject.numberofbins    =  bin's number of the binning result
2. binningobject.bin[i]          =  The i'th bin value (   0<  i   <=numberofbins)
3. binningobject.halfbin[i]      =  The i'th bin value (   0<  i   <=numberofbins/2)


fuction     2: readout the setbin file
struct binninginfo  Lreadoutbinfile(Tstring "path of setbins file");


Setbins file has a format with:
number of bins          i           genbin_i        recobin_i
number of bins           i+1       genbin_i+1    recobin_i+1

fuction     3:  Get out half number of bins
struct binninginfo Lgethalfbinning(struct binninginfo input)


Part   ②: unfolding


struct objrequirement{
vector<float> *obj_parameter;
double xmin;
double xmax;
};

function    4： judge the event if it is fit the condition
bool Lselectedevent(vector<float> *obj_parameter,int i,double xmin, double xmax);
default=true;

obj_parameter:  refers to the objected parameter like JetPt, dR etc..
i            :  refers to the event i;
xmin,xmax    :  refers to the seleted range of the obj_parameter.

function    5: Selected event in  multi conditions
bool Lselectedevent_multi_requirements(int n,struct objrequirement *input,int i)

n                 :  refers to the number of conditions or requirements
input           : refers to the multi requirements 
i                  : refers to the event i


function    6： creat histogram through existed binning and distribution ( whose bin number consistes with two different objects )
TH1 *histogram1Dcreation(TTree *t,vector<float> *obj1,vector<float> *obj2,TString name_of_histogram,TUnfoldBinning *obj_binning,TUnfoldBinning *obj_distribution,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected);
t: refers to the tree that object stay in
obj1,obj2: refers to two objects
name_of_histogram:refers to the name of the created histogram
obj_binning,obj_distribution: refers to the binning & distribution consisted with two objects
entrymin,entrymax:refers to the specified percentage range of the entries 
numberofrequirements: refers to the number of requirements (help to selecte events)
particle_selected:refers to the multi requirements 


function    7:      creat histogram through existed 2 binnings and distributions ( whose global bin number consists with two different objects )           

TH2 *histogram2Dcreation(TTree *t,vector<float> *obj1,vector<float> *obj2,vector<float> *obj3,vector<float> *obj4,TString name_of_histogram,TUnfoldBinning *obj_binning,TUnfoldBinning *obj_distribution,TUnfoldBinning *obj_binning1,TUnfoldBinning *obj_distribution1,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected);

the parameter is almost same as function6
obj1, obj2, obj_binning, obj_distribution: refers the x axis
obj3, obj4, obj_binning1, obj_distribution1: refers the y axis

function    8: get Tunfold result
                    TUnfoldDensity Ltunfold(TH2 *histMCGenReco,TH1 *histMCRecoBackground,TH1 *histDataReco,TUnfoldBinning * generatorBinning,TUnfoldBinning * detectorBinning)

using scan tau and nscan=30;

histMCGenReco:refers to the migration matrix
histMCRecoBackground:refers to the background on reconstructed level
histDataReco:refers to the Data(input) on the reconstructed level
generatorBinning: refers to the binning on generated level
detectorBinning: refers to the binning on reconstructed level

function 9: transfer the histogram with a binning( consisted with two object) 
TH1D *Ltransfer_2Daxisto1D(TH1 *input,int ndrbins,int nweightbins,double *drbins,double *weightbins)

input:refers to the input histogram
ndrbins,nweightbins:refers to the number of bins of two objected parameter
drbins,weightbins:refers to the bin of two objected parameter

function 10:get a new histogram with: hist1-hist2
TH1 *Lsub_histogram(TH1 *input1,TH1 *input2,TString name_of_histogram)

input1:refers to hist1
input2:refers to hist2
name_of_histogram:refers to the name of new histogram

function 11: creat histogram through existed binning and distribution ( whose bin number consistes with one object )
TH1 *D1histogram1Dcreation(TTree *t,vector<float> *obj1,TString name_of_histogram,TUnfoldBinning *obj_binning,TUnfoldBinning *obj_distribution,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected)

almost the same as function 6

function 12: creat histogram through existed 2 binnings and distributions ( whose global bin number consists with one object )   
TH2 *D1histogram2Dcreation(TTree *t,vector<float> *obj1,vector<float> *obj2,TString name_of_histogram,TUnfoldBinning *obj_binning,TUnfoldBinning *obj_distribution,TUnfoldBinning *obj_binning1,TUnfoldBinning *obj_distribution1,double entrymin,double entrymax,int numberofrequirements,struct objrequirement *particle_selected)

almost the same as function 7

function 13:get sum(bincontent(i)*(bincenter(i))) 
                    Lintweight(TH1 *input,double *xbins1)
                    input : refers to the input histogram
                    xbins1:refer to the bin of this histogram
                    
function 14:  creat a new histogram with determined name ,bin,and content
TH1D *Lc_1Dhist(const char *hist_name,int ndrbins,double *drbins,double *intweight)

hist_name:refer to the name of histogram
ndrbins:refers to the number of drbins
drbins:refers  to the bins of dr
intweight: refers to the weight of bin (we will get bincontent i with intweitht/binwidth)


