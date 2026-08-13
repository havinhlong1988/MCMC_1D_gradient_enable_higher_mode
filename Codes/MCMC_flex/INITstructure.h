#include<iostream>
#include<vector>
#include<cstdlib>
using namespace std;

struct invcov {
  double cov0;
};

struct Qmdef{
  int Qn;  // Qn stores the layerized Qs model, which is defined as 1km/layer
  vector<double> Qmod;
  vector<double> Qlayer;
};

struct groupdef {
int np,nlay;
int flagBs; // flag for Bspline, it will change to 1 once the Bspling fuctions are stored in
int flag,r_flag,q_flag,p_flag;   //1--layered model 2--Bspline 3--gradient 5-- water layer
double thick,vpvs,fvpvs,vpvs1; //thickness of each group,, vpvs ratio, flag_vpvs, vpvs1
double drho,frho,drho1;
double ddep; //thickness step of each group that separate the group into smaller ones (should be less than thick)
double dvs,fdvs,dvs1;
double fvpvs1,fdvs1,frho1;
double Qs;  // Qs value for the group; 
vector<double> ratio,value,value1,value1p,value1r,value1q,thick1,Bsplines;

//ratio, from readin file of layered model;value from readin file,vel;
//value1, layered vel; thick1, layered h;
//Bsplines, Bspline function the this group;
};

struct layermoddef {
int nlayer;
vector<double> vs,vpvs,vp,rho,qs,qp,thick,q_flag,cdep;
//parameters for layered model
};

struct rfdef {
int nrfo,rt;
double L,misfit,gaus;
vector<double> to,rfo,unrfo,tn,rfn,tnn,rfnn;
invcov** covariance;
int nlen;
double nn,nc;
int nflag;
//parameters for receiver function
//to , rfo and unrfo are readin file
//tn, rfn are computed ones.
};

struct dispdef{
int npper,ngper,neper;
int fphase,fgroup,fellip;

double pL,pmisfit,gL,gmisfit,eL,emisfit,L,misfit;
//period, liability,misfit
vector<int> flags;
vector<double> pper,pvelo,pvel,unpvelo;
vector<double> gper,gvelo,gvel,ungvelo;
vector<double> eper,evelo,evel,unevelo;
vector<double> period1;
//readin period, readin vel, model-derived vel, readin unc;
//period1 stores the period that will be computed 
};

struct datadef{
rfdef rf;
dispdef disp;
double p,L,misfit;
//weight of rf;
};

// -------------------- HVL: Data avaiable flag ----------------------------------------------
struct indatadef{
  // Flag for group - phase velocities, hv data, receiver function data
  // =1: data does exist, = 0 data does not exist
  // weighting data flag =1: using weight; =0 equal weight use based on the data
  // weight for each data type
int gvflag,phflag,hvflag,rfflag;
int isequalweight;
double gvw,phw,hvw,rfw;
};

///////////////////////////////////////////////////////////////////////////////////////////////
struct paradef{
int npara,flag;
int nparav; //number of pertubation vel
int npara_tt; //number of pertubation total thickness
int npara_tr; //number of pertubation thickness ratio
//number of parameter, vel and h; onece updated para, falg =1, after 1st model_derived para;
double L,misfit;
vector<double> parameter;
vector< vector<string> > para0;
vector<vector<double> > space1;
//??,??,range of para
};

struct modeldef{
int ngroup,flag,cc;
double tthick;
double mohodepth;//Added by EMB

// initial spline values used in calculating in b/w points of crust
//vector<double> pdBs0;

/*double percdiffpara1vs0;
double percdiffpara1vs2;
double percdiffpara3vs2;
double percdiffpara3vs4;
double percdiffpara5vs4;
double percdiffpara5vs6;
double percdiffpara7vs6;
double percdiffpara7vs8;
double percdiffpara9vs8;
*///end of adds

//#of readin groups;if layered falg=1;
vector<groupdef> groups; // a sequence of groups
layermoddef laym0;
datadef data;
Qmdef Qmodel;
indatadef indata;
};

struct ave_mod_def {
double max_thick;
vector<double> v_mins;
vector<double> v_maxs;
vector<double> v_means;
vector <double> v_stds;
};

