//read and write binary files
/*#include<iostream>
#include<fstream>
#include<vector>
*/
using namespace std;

#include <stdlib.h>
#include <cmath>
#include <stdio.h>

//int write_bin1(ofstream &out, paradef &para, int sign, int iiter, int iaccp, int icore)
int write_bin1(ofstream &out, paradef &para, boost::array<double, 10> &tinfo)
{

  // fprintf(stderr,"> BIN_rw: write_bin1\n");
  // sleep(1);

int *ids, *n_len, *size, i;
double *Lmis,*param,*para_laym0;
ids=new int[4];
n_len=new int[2];
Lmis=new double[4];
param=new double[para.npara];
para_laym0=new double[4];
size = new int[4];

ids[0]=int(tinfo[0]);
ids[1]=int(tinfo[1]);
ids[2]=int(tinfo[2]);
ids[3]=int(tinfo[3]); 
size[0]=sizeof(int)*4;

Lmis[0]=tinfo[4];
Lmis[1]=tinfo[5];
Lmis[2]=tinfo[6];
Lmis[3]=tinfo[7];
size[1]=sizeof(double)*4;

n_len[0]=int(tinfo[8]);
n_len[1]=int(tinfo[9]);
size[2]=sizeof(int)*2;

for(i=0;i<para.npara;i++)param[i]=para.parameter[i];
size[3]=sizeof(double)*para.npara;

out.write((char *)ids,size[0]);
out.write((char *)Lmis,size[1]);
out.write((char *)n_len,size[2]);
out.write((char *)param,size[3]);

delete [] ids;
delete [] Lmis;
delete [] param;
delete [] n_len;
delete [] size;
return 1;
}



//int write_bin(modeldef &model,char *fbname,paradef &para,int sign, int iiter, int iaccp)
int write_bin(modeldef &model,ofstream &out,paradef &para,int sign, int iiter, int iaccp, int icore)
{

  // fprintf(stderr,"> BIN_rw: write_bin\n");
  // sleep(1);

//cout<<"now beging to write"<<endl;
//ofstream out(fbname,ios::in|ios::binary);
//if(!out)
//{cout<<"####cannot open file to write binary!!\n";exit (0);}
//-----parameter---------
int *ids,*n_len,*size,i;
double *Lmis,*param,*gthick,*para_laym0;
ids=new int[4];
n_len=new int[3];
Lmis=new double[4];
param=new double[para.npara];
gthick=new double[model.ngroup];
para_laym0=new double[4];
size = new int[6];

//cout<<"ok here"<<endl;
//------value passing---------
//--ids
ids[0]=sign; 
ids[1]=iiter; 
ids[2]=iaccp;
ids[3]=icore;  // thread number
size[0]=sizeof(int)*4;
//--Lmis
if (std::isnan(model.data.disp.misfit) || std::isnan(model.data.rf.misfit) || std::isnan(model.data.misfit)) {
  fprintf(stderr,"misfit is nan function! be careful!!!\n");
  }

Lmis[0]=model.data.L; 
Lmis[1]=model.data.misfit;
Lmis[2]=model.data.rf.misfit; 
Lmis[3]=model.data.disp.misfit;
size[1]=sizeof(double)*4;
//--n_len
n_len[0]=para.npara;
n_len[1]=model.ngroup;
n_len[2]=model.laym0.nlayer;

//cout<<"ok here1"<<endl;

size[2]=sizeof(int)*3;
//--param
for(i=0;i<para.npara;i++)param[i]=para.parameter[i];
size[3]=sizeof(double)*para.npara;
//--gthick
for(i=0;i<model.ngroup;i++)gthick[i]=model.groups[i].thick;
size[4]=sizeof(double)*model.ngroup;
//--para_laym0
size[5]=sizeof(double)*4;
//------write binary file----------------
/*
sign iiter iaccp
model.data.L   model.data.misfit    ...rf.misfit    ...disp.misfit
para.npara   model.ngroup    model.laym0.nlayer
para.parameter[]
model.groups[].thick
model.laym0.thick[i] ...vs[i] ...vp[i] ...rho[i]
*/
//cout<<"ok here2"<<endl;
out.write((char *)ids,size[0]);
out.write((char *)Lmis,size[1]);
out.write((char *)n_len,size[2]);
out.write((char *)param,size[3]);
out.write((char *)gthick,size[4]);
//cout<<"ok here3"<<endl;
//cout<<size[0]<<" "<<size[1]<<" "<<size[2]<<endl;
//cout<<Lmis[0]<<" "<<Lmis[1]<<" "<<Lmis[2]<<" "<<Lmis[3]<<endl;
//cout<<ids<<" "<<Lmis<<" "<<n_len<<endl;

/*
for(i=0;i<model.laym0.nlayer;i++)
{
   para_laym0[0]=model.laym0.thick[i];
   para_laym0[1]=model.laym0.vs[i];
   para_laym0[2]=model.laym0.vp[i];
   para_laym0[3]=model.laym0.rho[i];
   out.write((char *)para_laym0,size[5]);
}
*/
//-----END of writing------------------
//out.close();
delete [] ids;
delete [] n_len;
delete [] Lmis;
delete [] param;
delete [] gthick;
delete [] para_laym0;
delete [] size;
return 1;
}

// 
int read_bin1(ifstream &fin , paradef &para, int &iaccp0, boost::array<double, 10> &tinfo)
{

  // fprintf(stderr,"> BIN_rw: read_bin1\n");
  // sleep(1);

int sign, iiter, iaccp, icore;

int *size,*ids,*n_len,i;
double *misL,*gthick,*param;
double *para_laym0;
size = new int[4];
ids = new int[4];
n_len = new int[2];
misL=new double[4];
size[0]=sizeof(int)*4;
size[1]=sizeof(double)*4;
size[2]=sizeof(int)*2;

fin.read((char *)ids,size[0]);//sign iiter iaccp
fin.read((char *)misL,size[1]);//model.data.L   model.data.misfit    ...rf.misfit    ...disp.misfit
fin.read((char *)n_len,size[2]);//para.npara   model.ngroup    model.laym0.nlayer

size[3]=sizeof(double)*n_len[0];
param=new double[size[3]];
//cout<<size[3]<<endl;
fin.read((char *)param,size[3]);//para.parameter[]
double tvalue;
for(i=0;i<n_len[0];i++) { tvalue = param[i]; para.parameter[i] = tvalue; }

sign=ids[0];
iiter=ids[1];
iaccp=ids[2];
icore=ids[3];
//cout<<sign<<" "<<iaccp<<" "<<n_len[0]<<endl;
tinfo[0] = double(sign);
tinfo[1] = double(iiter);
tinfo[2] = double(iaccp);
tinfo[3] = double(icore);

iaccp0 = iaccp;

para.L=misL[0];
para.misfit=misL[1];

tinfo[4]=misL[0];
tinfo[5]=misL[1];
tinfo[6]=misL[2];
tinfo[7]=misL[3];
tinfo[8]=double(n_len[0]);
tinfo[9]=double(n_len[1]);

return 1;
}


//int read_bin(modeldef &model, paradef &para,char *fbname, int &sign, int &iiter, int &iaccp)
int read_bin(ifstream &fin, modeldef &model , paradef &para, int &iaccp0)
{

  // fprintf(stderr,"> BIN_rw: read_bin\n");
  // sleep(1);

//modeldef model;
//paradef para;
int sign, iiter, iaccp, icore;

int *size,*ids,*n_len,i;
double *misL,*gthick,*param;
double *para_laym0;
size = new int[6];
ids = new int[4];
n_len = new int[3];
misL=new double[4];
size[0]=sizeof(int)*4;
size[1]=sizeof(double)*4;
size[2]=sizeof(int)*3;
size[5]=sizeof(double)*4;


fin.read((char *)ids,size[0]);//sign iiter iaccp
fin.read((char *)misL,size[1]);//model.data.L   model.data.misfit    ...rf.misfit    ...disp.misfit
fin.read((char *)n_len,size[2]);//para.npara   model.ngroup    model.laym0.nlayer

//cout<<"read ok!!"<<endl;
//cout<<misL[0]<<" "<<misL[1]<<endl;

size[3]=sizeof(double)*n_len[0];//double*para.npara
size[4]=sizeof(double)*n_len[1];//double*model.ngroup
//cout<<n_len[0]<<endl;
param=new double[size[3]];
gthick=new double[size[4]];
//cout<<gthick<<endl;
//para_laym0=new double(size[5]);
fin.read((char *)param,size[3]);//para.parameter[]
fin.read((char *)gthick,size[4]);//model.groups[].thick

//cout<<"read all info ok!!"<<endl;
//cout<<misL[0]<<" "<<misL[1]<<endl;

sign=ids[0];
iiter=ids[1];
iaccp=ids[2];
icore=ids[3];
iaccp0 = iaccp;

para.L=misL[0];
para.misfit=misL[1];

if (std::isnan(misL[3]) || std::isnan(misL[2]) || std::isnan(misL[1])) {fprintf(stderr,"it is nan somthing wrong here!!! misfit is nan!\n");}

if (std::isnan(misL[1])) misL[1]=100.;
if (std::isnan(misL[2])) misL[2]=100.;
if (std::isnan(misL[3])) misL[3]=100.;

model.data.L=misL[0];
model.data.misfit=misL[1];
model.data.rf.misfit=misL[2];
model.data.disp.misfit=misL[3];
para.npara=n_len[0];
model.ngroup=n_len[1];
model.laym0.nlayer=n_len[2];
//cout<<para.npara<<" "<<model.ngroup<<endl;
//para.parameter.clear();
double tvalue = 0.;
//cout<<"ok here~"<<endl;
//for (i=0;i<para.npara;i++) cout<<param[i]<<" "<<i<<endl;
for(i=0;i<para.npara;i++) { /*cout<<"ok here~ "<<param[i]<<" "<<para.parameter.size()<<endl;*/ tvalue = param[i]; para.parameter[i] = tvalue; }
//cout<<"ok here~"<<endl;
//model.groups.clear();

//cout<<model.ngroup<<endl;
//cout<<gthick<<endl;
for(i=0;i<model.ngroup;i++){ tvalue = gthick[i];model.groups[i].thick=tvalue;
  //cout<<"i "<<i<<" "<<gthick[i]<<endl;
  }
//if (gthick != NULL){
//cout<<gthick<<endl;
//delete gthick;
//}
//----------------------------------------------------------------
//cout<<"ok here~"<<endl;
//fprintf(stderr,"sign:%d iiter:%d iaccp:%d\n",sign,iiter,iaccp);
//fprintf(stderr,"model.data.L: %f misfit: %f rf.misfit: %f disp.misfit: %f\n",model.data.L,model.data.misfit, model.data.rf.misfit, model.data.disp.misfit);
//in.close();

delete [] size;
delete [] ids;
delete [] n_len;
delete [] misL;
delete [] gthick;
delete [] param;

//delete [] para_laym0;
//cout<<"now returen!"<<endl;
return 1;
}
