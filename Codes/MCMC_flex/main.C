#include<iostream>
#include<vector>
#include<algorithm>
#include<vector>
#include<cmath>
#include<fstream>
#include"/home/jiayi/progs/jy/HEAD/head_c++/string_split.C"// " " can not represent "\t"
#include"/home/jiayi/progs/jy/HEAD/head_c++/generate_Bs.C"
#include"/home/jiayi/progs/jy/HEAD/head_c++/gen_random.C"
#include"INITstructure.h"
#include"CALgroup.C"
#include"CALmodel.C"
#include"CALforward.C"
#include"CALpara.C"
int main(){
modeldef model;
vector<string> vnm;
//modelcal mc;
char name[100],name2[100];
sprintf(name,"disp.txt");
sprintf(name2,"mod.txt");
initmodel(model);
vnm.push_back(name2);
readdisp(model,vnm,1);
cout<<"begin read mod\n";
readmod(model,name2);
/*for(int i=0;i<model.data.disp.nper;i++)
	{
	cout<<model.data.disp.per[i]<<endl;
	cout<<model.data.disp.velo[i]<<endl;
	cout<<model.data.disp.unvelo[i]<<endl;
	}
*/
cout<<"model.ngroup "<<model.ngroup<<endl;
for(int i = 0;i<model.ngroup;i++)
	{
	 cout<< "groups i:"<<i<<endl;
	 cout<<"model.falg: "<<model.groups[i].flag<<endl;
	 cout<<"model.group.value: "<<"np="<<model.groups[i].np<<"; size="<<model.groups[i].value.size()<<endl;
	 for(int j=0;j<model.groups[i].np;j++){cout<<model.groups[i].value[j]<<"  "<<endl;}
	 }
int nBs=5,degBs=4,npts=60;
double disfacBs=3.,zminBs=0.,zmaxBs=160.;
vector<float> Bs_basis;
gen_B_spline(nBs,degBs,zminBs,zmaxBs,disfacBs,npts,Bs_basis);
return 1;
}
