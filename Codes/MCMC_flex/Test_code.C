#include<stdio.h>
#include<iostream>
#include<vector>
#include<algorithm>
#include<vector>
#include<cmath>
#include<fstream>
#include"/home/jiayi/progs/jy/HEAD/head_c++/string_split.C"
#include"/home/jiayi/progs/jy/HEAD/head_c++/generate_Bs.C"
#include"/home/jiayi/progs/jy/HEAD/head_c++/gen_random.C"
#include"INITstructure.h"
#include"CALgroup.C"
#include"CALmodel.C"
#include"CALforward.C"
#include"CALpara.C"
using namespace std;
int main()
{
vector<string> inf;
inf.push_back("/home/jiayi/progs/cv/CODE_oct5/TEST/Q22A.com.txt");
modeldef model;
initmodel(model);
fprintf(stderr,"ngper=%d npper=%d fphase=%d fgroup=%d  \n",model.data.disp.ngper,model.data.disp.npper,model.data.disp.npper,model.data.disp.fphase,model.data.disp.fgroup);
//======read disp
readdisp(model,inf,1);
fprintf(stderr,"ngper=%d npper=%d fphase=%d fphase=%d fgroup=%d , \n",model.data.disp.ngper,model.data.disp.npper,model.data.disp.fphase,model.data.disp.fgroup);
for(int i=0;i<model.data.disp.ngper;i++){cout<<model.data.disp.gper[i]<<" "<<model.data.disp.gvelo[i]<<" "<<model.data.disp.ungvelo[i]<<endl;}
//======read rf
char rfnm[200];
sprintf(rfnm,"/home/jiayi/progs/cv/CODE_oct5/TEST/in.rf");
readrf(model,rfnm);
cout<<model.data.rf.nrfo<<endl;
for(int i=0;i<10;i++){cout<<model.data.rf.to[i]<<" "<<model.data.rf.rfo[i]<<" "<<model.data.rf.unrfo[i]<<endl;}
//------read model
char mdnm[200];
sprintf(mdnm,"/home/jiayi/progs/cv/CODE_oct5/TEST/Q22A.mod1");
readmod(model,mdnm);
for(int i=0;i<model.ngroup;i++)for(int j=0;j<model.groups[i].np;j++){cout<<model.groups[i].flag<<" "<<model.groups[i].value[j]<<" "<<model.groups[i].ratio[j]<<" "<<model.groups[i].vpvs<<endl;}
//---init and -read para
paradef para;
initpara(para);
cout<<para.npara<<" "<<para.L<<" "<<para.misfit<<endl;
char prnm[200];
sprintf(prnm,"/home/jiayi/progs/cv/CODE_oct5/TEST/in.para");
readpara(para,prnm);
for(int i=0;i<para.npara;i++){for(int j=0;j<para.para0[i].size();j++){cout<<para.para0[i][j]<<"\t";}cout<<endl;}
//====update model
updatemodel(model);
for(int i=0;i<model.ngroup;i++){for(int j=0;j<model.groups[i].value1.size();j++){cout<<"\t"<<model.groups[i].value1[j]<<" "<<model.groups[i].thick1[j]<<endl;}cout<<endl;}
cout<<"============\n";
for(int i=0;i<model.laym0.nlayer;i++)cout<<model.laym0.vpvs[i]<<" "<<model.laym0.vs[i]<<" "<<model.laym0.vp[i]<<" "<<model.laym0.rho[i]<<" "<<model.laym0.qs[i]<<" "<<model.laym0.qp[i]<<" "<<model.laym0.thick[i]<<endl;
//====compute rf
return 1;
}
