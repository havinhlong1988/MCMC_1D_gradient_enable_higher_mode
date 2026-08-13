#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

//#define MAIN
//#include <iostream>
//#include<stdlib.h>
//#include<stdio.h>
//#include <math.h>
//using namespace std;


//void lf_B_spline(int nBs, int degBs, double zmin_Bs, double zmax_Bs, double disfacBs, int npts, double *depth, double *Bs_basis)
int main(int na, char *arg[])
{ 
  if(na!=7)
    {
      //  cout<<"usage:lf_B_spline zmin_Bs zmax_Bs distribution_factor"<<endl;
      fprintf(stderr,"usage:lf_B_spline zmin_Bs zmax_Bs distribution_factor nBs degBs npts\n");
      return 0;
    }
   int m,i,j,pp,n_temp;
   //nBs=4; //5 curves
   int nBs=atoi(arg[4]);
   int degBs=atoi(arg[5]); // cubic B-spline
   m=nBs-1+degBs;
   double t[m+1];
   double step,temp;
   double zmax_Bs=atof(arg[2]);
   double zmin_Bs=atof(arg[1]);
   double disfacBs=atof(arg[3]);
   for(i=0;i<degBs;i++)
     {
       //t[i]=zmin_Bs+i*(zmax_Bs-zmin_Bs)/10000;//Weisen
       t[i]=zmin_Bs+i*(zmax_Bs-zmin_Bs)/10000000;//FC
    }
  for(i=degBs;i<m+1-degBs;i++)
    {
      n_temp=m+1-degBs-degBs+1;
      if(disfacBs!=1)
	{
	  temp=(zmax_Bs-zmin_Bs)*(disfacBs-1)/(pow(disfacBs,n_temp)-1);
	  t[i]=temp/(disfacBs-1)*(pow(disfacBs,(i+1-degBs))-1)+zmin_Bs;
	}
      else
	{
	  temp=(zmax_Bs-zmin_Bs)/n_temp*(i+1-degBs);
	  t[i]=temp+zmin_Bs;
	}
	  //fprintf(stderr,"test %d %lf\n",i,t[i]);
    }
  for(i=m+1-degBs;i<m+1;i++)
    {
      t[i]=zmax_Bs-(zmax_Bs-zmin_Bs)/10000000*(m-i);
    }
  //  for(i=0;i<m+1;i++)
  //t[i]=200.0/12.0*i;
  for(i=0;i<m+1;i++)
    fprintf(stderr,"test %d %lf\n",i,t[i]);
  //cout<<i<<" "<<t[i]<<endl;
  int npts=atoi(arg[6]);
  //int npts=60;
  step=(zmax_Bs-zmin_Bs)/(npts-1);
  double basis_old[m][npts];
  double basis_new[m][npts];
  double depth[npts];
  for(j=0;j<npts;j++)
    depth[j]=j*step+zmin_Bs;
  for(i=0;i<m;i++)
    {
      for(j=0;j<npts;j++)
	{
	  //temp=j*step+zmin_Bs;
	  if(depth[j]>=t[i]&&depth[j]<t[i+1])
	    basis_old[i][j]=1;
	  else
	    basis_old[i][j]=0;
	}
    }
  for(pp=1;pp<degBs;pp++)
    {
      for(i=0;i<m-pp;i++)
	{
	  for(j=0;j<npts;j++)
	    {
	      //temp=j*step+zmin_Bs;
	      basis_new[i][j]=(depth[j]-t[i])/(t[i+pp]-t[i])*basis_old[i][j]+(t[i+pp+1]-depth[j])/(t[i+pp+1]-t[i+1])*basis_old[i+1][j];
	      // cout<<pp<<" "<<i<<" "<<j<<" "<<basis_new[i][j]<<" "<<(t[i+pp+1]-t[i+1])<<endl;
	    }
	}
      for(i=0;i<m-pp;i++)
	{
	  for(j=0;j<npts;j++)
	    {
	      basis_old[i][j]=basis_new[i][j];
	    }
	}
    }
  basis_new[0][0]=1;
  basis_new[m-pp][npts-1]=1;
  FILE *f1;
  char filename[300];
  int k=0;
  double depth_old;
  for(i=0;i<=m-pp;i++)
    {
      sprintf(filename,"B_spline_%d.txt",i);
      f1=fopen(filename,"w");
      depth_old=0;
      for(j=0;j<npts;j++)
	{
	  //	  if(depth[j]-depth_old<depth_old*0.05&&j!=npts-1)
	  //continue;

	  fprintf(f1,"%lf %lf\n",depth[j],basis_new[i][j]);
	  depth_old=depth[j];
	      //cout<<j*step+zmin_Bs<<" "<<basis_new[i][j]<<endl;
	  //Bs_basis[k]=basis_new[i][j];
	  k++;
	}
      fclose(f1);
      //cout<<endl;
    }
  //  for(j=0;j<npts;j++)
  //{
  //  temp=0;
  //  for(i=0;i<=m-pp;i++)
  //{
  //  temp+=basis_new[i][j];
  //}
  //  cout<<depth[j]<<" "<<temp<<endl;
  //}
  //  return 1;
  return 0;
}
