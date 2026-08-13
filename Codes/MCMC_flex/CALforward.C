#include<iostream>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
using namespace std;

extern"C" 
{ 
//void theo_(int *n,double *fbeta,double *h,double *vps,double *qa,double *fs,double *din,double *a0,double *c0,double *t0,int *nd,double *rx);
void theo_ (int *n,float *fbeta,float *h,float *vps,float *qa,float *qb,float *fs,float *din,float *a0,float *c0,float *t0,int *nd,float *rx);
//void fast_surf__(int *n_layer0,int *kind0,double *a_ref0,double *b_ref0,double *rho_ref0,double *d_ref0,double *qs_ref0,double *cvper,int *ncvper,double *uR0,double *uL0,double *cR0,double *cL0);
void fast_surf_(int *n_layer0,int *kind0,float *a_ref0,float *b_ref0,float *rho_ref0,float *d_ref0,float *qs_ref0,float *cvper,int *ncvper,float *uR0,float *uL0,float *cR0,float *cL0,float *rR0, float *rL0);
}
/*========CONTENT==========
int compute_disp(modeldef &model)
int compute_rf(modeldef &model)
//=========================
*/

int compute_disp(modeldef &model,int invtype)// WILL BE CHANGED LATER
{
  // fprintf(stderr,"> CALforward: compute_disp\n");
  // sleep(1);

  int newnlayer,nn,i,k,nper,nper1,N=1000;
  float *tvs, *tvp, *tqs, *trho, *tthick; // can use a 2d array for simplicity. but to make code clear, use 1d here.
  float *uR0,*uL0,*cR0,*cL0,*rR0,*rL0,*period;
  double ttt;
  vector<double>::iterator id;
  vector<double> period1,tvel;
  tvs = new float[N]; tvp = new float[N]; tqs= new float[N]; trho= new float[N]; tthick= new float[N];
  uR0 = new float[N]; uL0 = new float[N]; cR0= new float[N]; cL0= new float[N]; period = new float[N];
  rR0 = new float[N]; rL0 = new float[N];
  for(i=0;i<N;i++) // initial setting
	{
	 *(tvs+i)=*(tvp+i)=*(tqs+i)=*(trho+i)=*(tthick+i)=0.;
	 *(uR0+i)=*(uL0+i)=*(cR0+i)=*(cL0+i)=*(period+i)=0.;
         *(rR0+i)=*(rL0+i)=0.;
	}
  if(model.flag!=1)
	{
	  cout<<"#####model hasn't been updated/layered!!\n"; 
	  updatemodel(model,invtype);
	}
  newnlayer=model.laym0.nlayer+1;
  nn=0;
  for(i=0;i<newnlayer-1;i++)
   {
    ttt=0.;
    if(model.laym0.thick[i]<0.005 and i>1)
      {
        ttt=model.laym0.thick[i];
        //cout<<"####too thin layer, thick="<<ttt<<" at layer "<<i<<endl;
        continue;
      }//if
    tvs[nn]=model.laym0.vs[i];
    tvp[nn]=model.laym0.vp[i];
    trho[nn]=model.laym0.rho[i];
    tthick[nn]=model.laym0.thick[i]+ttt;//depth
          //cout<<"tvs_nn "<<tvs[nn]<<"nn "<<nn<<endl;
    nn=nn+1; //nn equivalent to number of layers+1
   }//for i 
  //cout<<"nn "<<nn<<endl;
  k=0;
  //cout<<" check everything: "<<endl;
  for(i=0;i<nn;i++)
   {
    if(model.laym0.qs[i]>0.)
      tqs[k]=1./model.laym0.qs[i];
    else
      tqs[k]=1./(model.laym0.qs[i]+10000.);
      //cout<<k<<" "<<tvs[k]<<" "<<tvp[k]<<" "<<trho[k]<<" "<<tthick[k]<<" "<<tqs[k]<<endl;
      // fprintf(stderr,"check everything: %d %g %g %g %g\n",k,tvs[k],tvp[k],trho[k],tthick[k],tqs[k]);
    k=k+1;
   }//for i
  //the last element in arrays
  k=model.laym0.nlayer-1; // HVL: why loop for k then fixed k here ?
  tvs[nn]=model.laym0.vs[k]+0.01;
  tvp[nn]=model.laym0.vp[k]+0.01;
  tqs[nn]=tqs[nn-2];
  trho[nn]=model.laym0.rho[k]+0.01;
  tthick[nn]=0.;
  //cout<<"last element check.. nn "<<nn<<" k "<<k<<" "<<tvs[nn]<<endl;
  nn=nn+1;//#of element inside the arrays

  nper=0;  
  if(model.data.disp.fphase+model.data.disp.fgroup+model.data.disp.fellip>0)
	{
	  period1=model.data.disp.period1;//temo
	}
  else if(model.data.disp.fphase>0)
	  period1=model.data.disp.pper;
  else if(model.data.disp.fgroup>0)
	  period1=model.data.disp.gper;
  else if (model.data.disp.fellip>0)
    period1=model.data.disp.eper;
  else
	  {cout<<"#####wrong model period!!!\n";
	  exit (0);}
  //
  sort(period1.begin(),period1.end());
  nper1=period1.size();
  for (i=0;i<nper1;i++) {
	period[i]=period1[i];
      //  fprintf(stderr,"%g %g %d %d %d\n",period[i],period1[i],nper,i,nper1);
        }
  nper=nper1;  
/*subroutine FAST_SURF(n_layer0,kind0,
     &          a_ref0,b_ref0,rho_ref0,d_ref0,qs_ref0,
     &          cvper, ncvper,
     &          uR0,uL0,cR0,cL0)
void FAST_SURF_(int *n_layer0,int *kind0,double *a_ref0,double *b_ref0,double *rho_ref0,double *d_ref0,double *qs_ref0,double *cvper,int *ncvper,double *uR0,double *uL0,double *cR0,double *cL0);
*/
  int cflag = 2;
  // fprintf(stderr,"now we'll do the fast_surf!!!! with %d elements and %d period\n",nn, nper);
  fast_surf_(&nn,&cflag,tvp,tvs,trho,tthick, tqs,period,&nper,uR0,uL0,cR0,cL0,rR0,rL0);
  // fprintf(stderr,"now we finish fast_surf!!!! %d %d %d\n",model.data.disp.npper,model.data.disp.ngper,model.data.disp.neper);

  // fprintf(stderr,"now we push phase velocity!!!! %d\n",model.data.disp.npper);
  if(model.data.disp.npper>0)
 	{
	  model.data.disp.pvel.clear();
	  for(i=0;i<model.data.disp.npper;i++)
		{
		  id=find(period1.begin(),period1.end(),model.data.disp.pper[i]);
      // fprintf(stderr,"pper now: %g !\n",model.data.disp.pper[i]);
		  if(id==period1.end())
			{        
        cout<<"#######!!!canot find pper "<<model.data.disp.pper[i]<<" in perio1\n";
        fprintf(stderr,"cannot find pper: %g in perio1!\n exit!!\n",model.data.disp.pper[i]);
        sleep(1);
        exit(0);
        }
		  model.data.disp.pvel.push_back(cR0[id-period1.begin()]);
            //cout<<period1[id-period1.begin()]<<" "<<cR0[id-period1.begin()]<<endl;
		}//for i
	}//if 

  
  if(model.data.disp.ngper>0)
  // fprintf(stderr,"now push the group velocity!\n");
	{
	  model.data.disp.gvel.clear();
	  for(i=0;i<model.data.disp.ngper;i++)
		{
		  id=find(period1.begin(),period1.end(),model.data.disp.gper[i]);
		  if(id==period1.end())
			{
        cout<<"######!!!canot find gper"<<model.data.disp.gper[i]<<"in perio1 i "<<i<<"\n";
        // fprintf(stderr,"cannot find gper: %g in perio1!\n exit!!\n",model.data.disp.gper[i]);
        exit(0);
        }
		  model.data.disp.gvel.push_back(uR0[id-period1.begin()]);
		}
	}//if

  if (model.data.disp.neper>0)
    // fprintf(stderr,"now push the ellipticity!\n");
       {
         model.data.disp.evel.clear();
          for(i=0;i<model.data.disp.neper;i++)
                {
		  //cout<<model.data.disp.eper[i]<<endl;
                  id=find(period1.begin(),period1.end(),model.data.disp.eper[i]);
                  if(id==period1.end())
                        {cout<<"######!!!canot find eper"<<model.data.disp.eper[i]<<"in perio1 i "<<i<<"\n";exit(0);}
                  // fprintf(stderr,"cannot find eper: %g in perio1!\n exit!!\n",model.data.disp.eper[i]);      
                  model.data.disp.evel.push_back(rR0[id-period1.begin()]);
                }
       }
  // fprintf(stderr,"Push finish everything seem ok!\n");
  delete[] tvs;
  delete[] tvp;
  delete[] tqs;
  delete[] trho;
  delete[] tthick;
  delete[] uR0;
  delete[] uL0;
  delete[] cR0;
  delete[] rR0; //added by EMB
  delete[] cL0;
  delete[] rL0; //added by EMB
  delete[] period;   
  //cout<<"finish compute DISP!"<<endl;
  // fprintf(stderr,"finish compute DISP! \n");
  return 1;
}//compute disp
//--------------------------------

int compute_rf(modeldef &model,int invtype)
{

  // fprintf(stderr,"> CALforward: compute_rf\n");
  // sleep(1);

  int i,newnlayer,nl,nn,N=100,NN=1000; //original
  float slow,pi,din,rt;
  float *tvs,*tvpvs,*tqs,*tqp,*trho,*tthick,*rx;

  tvs = new float[N]; 
  tvpvs=new float[N]; 
  tqs=new float[N]; 
  tqp=new float[N]; 
  trho=new float[N]; 
  tthick=new float[N]; 
  rx=new float[NN];

  for(i=0;i<N;i++)
	{tvs[i]=tvpvs[i]=tqs[i]=tqp[i]=trho[i]=tthick[i]=0.;}
  for (i=0;i<NN;i++) rx[i]=0.;

  if(model.flag!=1)
	{
	  cout<<"######model hasn't been layered/updated!\n";
	  updatemodel(model,invtype);  
	} 
  newnlayer=model.laym0.nlayer<100?(model.laym0.nlayer+1):100;
  nl=0;
  for(i=0;i<newnlayer-1;i++)
	{
	  if(model.laym0.vs[i]>0.)
		{
		  tvs[nl]=model.laym0.vs[i];
	  	  tvpvs[nl]=model.laym0.vpvs[i];
		  tqs[nl]=model.laym0.qs[i];
		  tqp[nl]=model.laym0.qp[i];
		  tthick[nl]=model.laym0.thick[nl];
		  nl=nl+1;
		}//if
	}//for i
  tvs[nl]=tvs[nl-1];
  tvpvs[nl]=tvpvs[nl-1];
  tqs[nl]=tqs[nl-1];
  tqp[nl]=tqp[nl-1];
  tthick[nl]=0.; 
  nl=nl+1;

  nn=1000;
  slow=0.06;//ray parameter
  pi=atan(1)*4;
//  cout<<"===== pi="<<pi<<endl;
  din=180.*asin(tvs[nl-1]*tvpvs[nl-1]*slow)/pi;
  rt=model.data.rf.rt;
/*subroutine theo(
     &     n, fbeta, h, vps, qa, qb, fs, din, a0, c0, t0,
     &     nd, rx )
void theo_(int *n,double *fbeta,double *h,double *vps,double *qa,double *fs,double *din,double *a0,double *c0,double *t0,int *nd,double *rx);

*/
  //float gau = 10;//gaussian of 10, just in case
  float gau = 2.5;//gaussian width
  if (model.data.rf.gaus > 0.) gau = (float)(model.data.rf.gaus);
  //fprintf(stderr," GAUSSIAN %g \n\n",gau);
  //float dt = 0.02;//parameter of minimum amplitude level
  float dt = 0.005; //original; parameter of minimum amp level
  float t0 = 0.;
  int nn1 = (int)newnlayer;
  

  theo_(&nn1,tvs,tthick,tvpvs,tqp,tqs,&rt,&din,&gau,&dt,&t0,&nn,rx);
  model.data.rf.tn.clear();
  model.data.rf.rfn.clear();
  //cout<<"compute_rf! nn: "<<nn<<endl;
  for(i=0;i<nn;i++)
	{
	  model.data.rf.tn.push_back(i*1./model.data.rf.rt);
          if (isnan(rx[i]))
          model.data.rf.rfn.push_back(-999.);
          else
  	  model.data.rf.rfn.push_back(rx[i]);
	}//for i

  delete[] tvs;
  delete[] tvpvs;
  delete[] tqs;
  delete[] tqp;
  delete[] trho;
  delete[] tthick;
  delete[] rx;
  return 1;

}//compute rf
