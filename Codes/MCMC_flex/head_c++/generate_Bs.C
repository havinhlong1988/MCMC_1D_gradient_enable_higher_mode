#include<iostream>
#include<math.h>
using namespace std;

//int gen_B_spline(int nBs, int degBs,double zmin_Bs, double zmax_Bs, double disfacBs, int npts)

template <class A, class B, class C,class D>
A gen_B_spline(A nBs, A degBs, B zmin_Bs, B zmax_Bs, C disfacBs, A npts,D &Bs_basis)
{
  //nBs=5 	degBs =4 for cubic Bspine
  int i,j,k,m=nBs-1+degBs,n_temp;
  double temp,step;
  vector<double> t(m+1,0.),depth;
  vector< vector<double> > obasis(m,vector<double>(npts,0.)),nbasis(m,vector<double>(npts,0.));
  FILE *ff;
  char fnm[100];
  for(i=0;i<degBs;i++)
	 { t[i]=zmin_Bs+i*(zmax_Bs-zmin_Bs)/10000.;
        }
  for(i=degBs;i<m+1-degBs;i++)
	{
	  n_temp=m+1-degBs-degBs+1;
	  //n_temp=m+1-degBs-degBs;
	  if(disfacBs!=1)temp=(zmax_Bs-zmin_Bs)*(disfacBs-1)/(pow(disfacBs,n_temp)-1);
	  else temp = (i+1-degBs)*(zmax_Bs-zmin_Bs)/n_temp;

          if (i!=degBs) t[i]=t[i-1]+temp*pow(disfacBs,(i-degBs));
          else t[i]=temp*pow(disfacBs,(i-degBs))+zmin_Bs;

          fprintf (stderr,"t{i}: %g temp: %g i: %d\n",t[i],temp,i);
	}
  for (i=m+1-degBs;i<m+1;i++)
	{t[i]=zmax_Bs-(zmax_Bs-zmin_Bs)/10000.*(m-i);}

  step=(zmax_Bs-zmin_Bs)/(npts-1);
  for (i=0;i<npts;i++)
	{depth.push_back(i*step+zmin_Bs);}
  for(i=0;i<m;i++)
	{
	  for(j=0;j<npts;j++)
		{
		  if(depth[j]>=t[i] and depth[j]<t[i+1])obasis[i][j]=1.;
		  else obasis[i][j]=0.;
		}//forj
	}
	  
  for(k=1;k<degBs;k++)
	{
	  for(i=0;i<m-k;i++)
		{
		  for(j=0;j<npts;j++)
			{
			  nbasis[i][j]=(depth[j]-t[i])/(t[i+k]-t[i])*obasis[i][j]+(t[i+k+1]-depth[j])/(t[i+k+1]-t[i+1])*obasis[i+1][j];
			}//forj
		}//fori
	  for(i=0;i<m-k;i++){for(j=0;j<npts;j++){obasis[i][j]=nbasis[i][j];}}
	}//fork 
  nbasis[0][0]=1;
  nbasis[nBs-1][npts-1]=1;

  for(i=0;i<nBs;i++)
	{
	  sprintf(fnm,"Bs.%d.dat",i);
          ff=fopen(fnm,"w");
	  for(j=0;j<npts;j++)
		{
		  fprintf(ff,"%g %g\n",depth[j],nbasis[i][j]);
		  Bs_basis.push_back(nbasis[i][j]);
		}
  	  fclose(ff);
	}
   return 0;
}//gen Bspline


