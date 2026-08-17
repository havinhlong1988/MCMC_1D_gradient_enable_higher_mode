// ====================================================================
//  MODIFIED 2026-08-14: enable higher-mode surface wave (fund + 1st).
//    initmodel: init hp* fields.
//    readdisp1: tflag==5 loads .hph (flag 5 -> keeps 2=group, 3=ellip).
//    readindata: parse OPTIONAL in.data line 10 (hpflag) & 11 (hpw).
//    compute_misfit: fold higher-mode phase into the weighted misfit
//      (tempv1h/phmisfit/hpL/tmfhp, weight hpw, flag ihp).
//    Also: <1.0 s period filter disabled (commented, re-enableable).
//    See modify_to_higher.log (sections C, D).
// ====================================================================
/*#include<fstream>
#include<iostream>
#include<vector>
#include<math.h>
#include<algorithm>
*/
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
//#include"INITstructure.h"
//#include"CALgroup.C"

//#include"/home/jiayi/progs/jy/HEAD/head_c++/string_split.C"
using namespace std;
/*=====CONTENT=======
int initmodel(modeldef &model)
int readdisp(modeldef &model,vector<string> name, int surflag)
int readrf(modeldef &model,char *name)
int readcov (modeldef &model, char *name)
int readmod(modeldef &model,const char *name)
int updatemodel(modeldef &model)
int compute_misfit(modeldef &model,double inp,double nn)
int write_model(modeldef &model,char *name, char *outdir)
int goodmodel(modeldef &model, vector<int> vmono, vector<int> vgrad)
//===================
*/
//class modelcal{
//public:
	int initmodel(modeldef &model)
	{
    // fprintf(stderr,"> CALmodel: initmodel\n");
    // sleep(1);
	  model.laym0.nlayer=0;

	  model.data.rf.nrfo=0;
	  model.data.rf.rt=0;
	  model.data.rf.L=0.;
	  model.data.rf.misfit=0.;
	  model.data.rf.gaus=-1.;

	  model.data.disp.npper=0;
	  model.data.disp.pL=0.;
	  model.data.disp.pmisfit=0.;
	  model.data.disp.fphase=0;

	  model.data.disp.ngper=0;
	  model.data.disp.gL=0.;
	  model.data.disp.gmisfit=0.;
	  model.data.disp.fgroup=0;

    model.data.disp.neper = 0;
    model.data.disp.eL = 0;
    model.data.disp.emisfit = 0;
    model.data.disp.fellip = 0;

	  model.data.disp.nhpper=0;
	  model.data.disp.phmisfit=0.;
	  model.data.disp.hpL=0.;
	  model.data.disp.fhphase=0;

	  model.data.p=0.;
	  model.data.L=0.;
	  model.data.misfit=0.;
	  model.ngroup=0;
	  model.flag=0;
	  model.cc=2;
	  model.tthick = 0.;  // added by weisen 10.11; total thickness; computed when readmodel;
	  model.mohodepth=0.; //added by EMB
      //for loop via EMB
      /*for(int j=0;j<10;j++){
           //cout<<"np?? j="<<j<<" model.groups[i].value[j]="<<model.groups[i].value[j]<<endl;
            model.pdBs0.push_back(0);
            cout<<"pdBs"<<j<<" = "<<model.pdBs0[j]<<endl;
       }*/

      model.Qmodel.Qn = 0;
	  return 1;
	}
//-----------------------------------------------------  read the Qmodel
	int readQmodel(modeldef &model, char *name, vector<int> Ql) {
    // fprintf(stderr,"> CALmodel: readQmodel\n");
    // sleep(1);
	  fstream mff;
          double tmpto,tmprfo,tmpunrfo;
          int i=0,size;
          string line;
          vector<string> v;
          if(model.Qmodel.Qn>0)
                {
                  cout<<"#########data existed in Qmodel already!! exit!\n";
                  exit(0);
                }
          mff.open(name);
          if ( not mff.is_open()){cout<<"#########rf file"<<name<<"does not exist!\n exit!!\n";exit(0);}
          model.Qmodel.Qmod.clear();
          cout<<name<<endl;
	  while(getline(mff,line))
          {
            v.clear();
            Split(line,v," ");
            size=v.size();
            if(size>0)
            {
                //cout<<"read Q here "<<atof(v[0].c_str())<<endl;
                model.Qmodel.Qmod.push_back(atof(v[0].c_str()));
                i = i + 1;
            }
	  }
          model.Qmodel.Qn = i;
          mff.close();
          v.clear();
          return 1;
        }
//-----------------------------------------------------	

	int readdisp1(modeldef &model, int surtype, int surn, vector<int> surflag, vector<string> name) {

    // fprintf(stderr,"> CALmodel: readdisp1\n");
    // sleep(1);

	  ifstream mff;
          double tmpper,tmpvel,tmpunc;
          int i,j,size,type=0,tflag;
          string line,tmpnm;
          char ll[150];
          vector<string> v;
          vector<int> sflag;
          vector<double> cv1,cv2,cv3;

	  if(model.data.disp.npper>0 or model.data.disp.ngper>0 or model.data.disp.neper >0) {
             cout<< "#########disp data existed already!\n exit!"<<endl;
             exit(0);
             }
          if (model.data.disp.fgroup>0 or model.data.disp.fphase>0 or model.data.disp.fellip>0) {
             cout<< "#########disp data existed already!\n exit!"<<endl;
             exit(0);
             }
          model.data.disp.pper.clear(); 
          model.data.disp.pvelo.clear(); 
          model.data.disp.unpvelo.clear();
          //
          model.data.disp.gper.clear(); 
          model.data.disp.gvelo.clear(); 
          model.data.disp.ungvelo.clear();
          //
          model.data.disp.eper.clear();
          model.data.disp.evelo.clear();
          model.data.disp.unevelo.clear();
          //
          model.data.disp.hpper.clear();
          model.data.disp.hpvelo.clear();
          model.data.disp.unhpvelo.clear();

          //
          for (j=0;j<surn;j++) {
             cout<<j<<endl;
             tflag = surflag[j]; 
             // read in file name[j] to pper, pvelo and unpvelo;
             mff.open(name[j].c_str());
             fprintf(stderr,"disp phase file : %s !\n",name[j].c_str());
             if ( not mff.is_open())
             {
              cout<<"###########disp phase file"<<name[j]<<"does not exist!\n exit!!\n";
              fprintf(stderr,"disp phase file : %s does not exist!\n",name[j].c_str());
              exit(0);
              }
             cv1.clear(); cv2.clear(); cv3.clear();
             i = 0;
             cout<<tflag<<" "<<name[j].c_str()<<endl;
             while(getline(mff,line)) {
                v.clear();Split(line,v," ");
                size=v.size();
                if (size <1) continue;
                // cv1.push_back(atoi(v[0].c_str()));
                
                double this_per = atof(v[0].c_str());

                // Force Vph only uses periods < 1.0 sec
                // if (this_per >= 1.0) continue; // remove all period > 1.0
                // ^ period limit disabled; re-enable this line to test the <1.0 s cut

                cv1.push_back(this_per);
                if (size == 1) {cv2.push_back(-1.); cv3.push_back(-1.);}
                if (size > 1 and size <2) {cv2.push_back(atof(v[1].c_str())); cv3.push_back(-1.);}
                if (size > 2) 
                {
                  cv2.push_back(atof(v[1].c_str())); 
                  cv3.push_back(atof(v[2].c_str()));
                  }
                i = i + 1;
                }
	           mff.close();
             cout<<i<<endl;
             // read all file into cv1, cv2, cv3;
             if (tflag == 1) {  // phase velocity
                   model.data.disp.pper = cv1;
                   model.data.disp.pvelo = cv2;
		  // for (auto &x : cv3) {
    		  //   x *= 5.0;} // auto multiple 5 times
                   model.data.disp.unpvelo = cv3;
                   model.data.disp.npper = i;
                   model.data.disp.fphase = 1;
                   }
             else if (tflag ==2) { // group velocity
                   model.data.disp.gper = cv1;
                   model.data.disp.gvelo = cv2;
                   model.data.disp.ungvelo = cv3;
                   model.data.disp.ngper = i;
                   model.data.disp.fgroup = 1;
                   }
             else if (tflag ==3 ) { // ellipticity velocity
                   model.data.disp.eper = cv1;
                   model.data.disp.evelo = cv2;
                   model.data.disp.unevelo = cv3;
                   model.data.disp.neper = i;
                   model.data.disp.fellip = 1;
                   }
             else if (tflag ==5 ) { // higher-mode (1st) phase velocity
                   model.data.disp.hpper = cv1;
                   model.data.disp.hpvelo = cv2;
                   model.data.disp.unhpvelo = cv3;
                   model.data.disp.nhpper = i;
                   model.data.disp.fhphase = 1;
                   }
             else {
                   cout<<"cannot recognize this flag!!!"<<endl; exit(0);
                   }
             }
          cv1 = model.data.disp.pper;
          cv2 = model.data.disp.gper;
          cv3 = model.data.disp.eper;
          sort(cv1.begin(),cv1.end());
          sort(cv2.begin(),cv2.end());
          sort(cv3.begin(),cv3.end());

          vector<double> v1(200);
          vector<double> v2(200);
          vector<double>::iterator it,it1,it2;
          int ictp;
          it1 = set_union (cv1.begin(), cv1.end(), cv2.begin(), cv2.end(), v1.begin());
          it2 = set_union (cv3.begin(), cv3.end(), v1.begin(), it1, v2.begin());
          for (it=v2.begin();it<it2;it=it+1) { if (*it >0.) model.data.disp.period1.push_back(*it-0.);
             //fprintf(stderr,"it : %g %g %g\n",*it,it,it1);
             }
          //for (i=0;i<model.data.disp.period1.size();i++) fprintf(stderr,"T: %d %g\n",i,model.data.disp.period1[i]);
          v1.clear();
          v2.clear();
          cv1.clear();
          cv2.clear();
          cv3.clear();
          fprintf(stderr,"disp check n-phase, n-group, p-ellip: %d %d %d \n",model.data.disp.npper,model.data.disp.ngper,model.data.disp.neper);
          for (int i;i<model.data.disp.npper;++i)
          {
            ictp = i; 
            // fprintf(stderr,"phase list: %d %g %g %g\n",i,model.data.disp.pper[i],model.data.disp.pvelo[i],model.data.disp.unpvelo[i]);

          }
          // fprintf(stderr,"> CALmodel: readdisp1 ok!\n");
          // sleep(5);
	  return 1;
	}
//-------------------------------------------------------
	int readdisp(modeldef &model,vector<string> name, int surflag) //use vector for name, since the name could be 1 name or two names
	{
	  //FILE *ff;

    // fprintf(stderr,"> CALmodel: readdisp\n");
    // sleep(1);
          
	  ifstream mff;
	  double tmpper,tmpvel,tmpunc;
	  int i,size,flag=0;
	  string line,tmpnm;
 	  char ll[150];
          vector<string> v;
	  vector<double> cv1,cv2;
          //cout<<<<endl;
	  if(model.data.disp.npper>0 or model.data.disp.ngper>0)
		{
		  cout<< "#########disp data existed already!\n exit!"<<endl;
		  exit(0);
		}
	  if(model.data.disp.fgroup>0 or model.data.disp.fphase>0)
		{ cout<< "##########disp data existed already!\n exit!"<<endl;
                  exit(0);}

	  model.data.disp.pper.clear(); model.data.disp.pvelo.clear();model.data.disp.unpvelo.clear();
	  if(surflag==1 or surflag==3) //surflag==1: open phase only. surfalg ==3 open phase and group, surflag==2: open group only
	  {
//	    cout<<"read phase disp!\n";
	    i=0;
	    mff.open(name[0].c_str());
            if ( not mff.is_open()){cout<<"###########disp phase file"<<name[0]<<"does not exist!\n exit!!\n";exit(0);}
            while(getline(mff,line))
                {
		  v.clear();
                  Split(line,v," ");
		  size=v.size();
		 //cout<<"now size is "<<size<<endl;
	         if(size>2)
		 {
		 model.data.disp.pper.push_back(atof(v[0].c_str()));
		 model.data.disp.pvelo.push_back(atof(v[1].c_str()));
		 model.data.disp.unpvelo.push_back(atof(v[2].c_str()));
		 }
		 else if(size >1)
		 {
		 model.data.disp.pper.push_back(atof(v[0].c_str()));
     model.data.disp.pvelo.push_back(atof(v[1].c_str()));
		 model.data.disp.unpvelo.push_back(0.);
		 }
		 else if(size>0)
		 {
		 model.data.disp.pper.push_back(atof(v[0].c_str()));
		 model.data.disp.pvelo.push_back(0.);
		 model.data.disp.unpvelo.push_back(0.);
		 }
		 else
		 {
		 model.data.disp.pper.push_back(0.);
		 model.data.disp.pvelo.push_back(0.);
     model.data.disp.unpvelo.push_back(0.);
		 }
		  i=i+1;
		}//while
	    model.data.disp.npper=i;
	    model.data.disp.fphase=1;
	    mff.close();
	  }// 1 and 3
//          cout<<"read phase ok!"<<endl;
	  //============= read group ================
	  i=0;
//          cout<<"surflag "<<surflag<<endl;
	  model.data.disp.gper.clear(); 
    model.data.disp.gvelo.clear();
    model.data.disp.ungvelo.clear();

//    cout<<"surflag "<<surflag<<endl;
//	  cout<<name[0]<<endl;
	  if(surflag==2){tmpnm=name[0];flag=1;}
	  else if(surflag==3){tmpnm=name[1];flag=1;}

//  cout<<flag<<endl;
	  if(flag)
	  {
           //cout<<tmpnm.c_str()<<endl;
 	   mff.open(tmpnm.c_str());
//	   cout<<"read group diisp file! \n";
	   if ( not mff.is_open()){cout<<"####disp group file"<<tmpnm<<"does not exist!\n exit!!\n";exit(0);}
	   while(getline(mff,line))
		{
		 //cout<<line<<endl;
		 v.clear();
		 Split(line,v," ");
		 size=v.size();
	         if(size>2)
		 {
		 model.data.disp.gper.push_back(atof(v[0].c_str()));
		 model.data.disp.gvelo.push_back(atof(v[1].c_str()));
		 model.data.disp.ungvelo.push_back(atof(v[2].c_str()));
		 }
		 else if(size >1)
		 {
		 model.data.disp.gper.push_back(atof(v[0].c_str()));
     model.data.disp.gvelo.push_back(atof(v[1].c_str()));
		 model.data.disp.ungvelo.push_back(0.);
		 }
		 else if(size>0)
		 {
		 model.data.disp.gper.push_back(atof(v[0].c_str()));
		 model.data.disp.gvelo.push_back(0.);
		 model.data.disp.ungvelo.push_back(0.);
		 }
		 else
		 {
		 model.data.disp.gper.push_back(0.);
		 model.data.disp.gvelo.push_back(0.);
     model.data.disp.ungvelo.push_back(0.);
		 }
      i=i+1;
		}//while
	    model.data.disp.ngper=i;
	    model.data.disp.fgroup=1;

	    mff.close();
	  }//if flag
          //for (i=0;i<model.data.disp.npper;i++)
            //  {
              //fprintf(stderr,"check phase here: %g %g %g \n",model.data.disp.pper[i],model.data.disp.pvelo[i],model.data.disp.unpvelo[i]);
            //  }
          //for (i=0;i<model.data.disp.ngper;i++)
            //  {
              //fprintf(stderr,"check group here: %g %g %g \n",model.data.disp.gper[i],model.data.disp.gvelo[i],model.data.disp.ungvelo[i]);
            //  }

	  cv1 = model.data.disp.pper;
    cv2 = model.data.disp.gper;
    sort(cv1.begin(),cv1.end());
    sort(cv2.begin(),cv2.end());

//          for (i=0;i<model.data.disp.npper;i=i+1) cout<<i<<" p check "<<cv1[i]<<endl;
//          for (i=0;i<model.data.disp.ngper;i=i+1) cout<<i<<" g check "<<cv2[i]<<endl;

	  vector<double> v1(200);
    vector<double>::iterator it,it1;

	  it1 = set_union (cv1.begin(), cv1.end(), cv2.begin(), cv2.end(), v1.begin());
	  for (it=v1.begin();it<it1;it=it+1) { if (*it >0.) model.data.disp.period1.push_back(*it-0.); 
             //fprintf(stderr,"it : %g %g %g\n",*it,it,it1);
             }

	  v1.clear();
	  cv1.clear();
	  cv2.clear();
          

          //cout<<"nper   cvcvcvvc "<<model.data.disp.npper<<endl;
//          cout<<model.data.disp.ngper<<endl;
          //cout<<model.data.disp.gper.size()<<endl;
	  return 1;
	}//readdisp
//-----------------------------------------------------	
	int readrf(modeldef &model,char *name, double gauss)
	{

    // fprintf(stderr,"> CALmodel: readrf\n");
    // sleep(1);

	  //FILE *ff;
	  fstream mff;
	  double tmpto,tmprfo,tmpunrfo;
	  int i=0,size;
	  string line;
	  vector<string> v;
	  if(model.data.rf.nrfo>0)
		{
		  cout<<"#########data existed in rf already!! exit!\n";
		  exit(0);
		}
	  mff.open(name);
	  if ( not mff.is_open()){cout<<"#########rf file"<<name<<"does not exist!\n exit!!\n";exit(0);}
	  model.data.rf.to.clear();model.data.rf.rfo.clear(); model.data.rf.unrfo.clear();
	while(getline(mff,line))
	{
	    v.clear();
	    Split(line,v," ");
	    size=v.size();
	    if(size>2)
	    {
		model.data.rf.to.push_back(atof(v[0].c_str()));
		model.data.rf.rfo.push_back(atof(v[1].c_str()));
		model.data.rf.unrfo.push_back(atof(v[2].c_str())*1);
	    }
	    else if(size>1)
	    {
		model.data.rf.rfo.push_back(atof(v[1].c_str()));
		model.data.rf.to.push_back(atof(v[0].c_str()));
		model.data.rf.unrfo.push_back(0.5);
	    }
	    else if(size>0)
	    {
		model.data.rf.rfo.push_back(atof(v[1].c_str()));
    model.data.rf.to.push_back(0.);
    model.data.rf.unrfo.push_back(0.5);
	    }
	   else
	   {
		model.data.rf.rfo.push_back(0.);
    model.data.rf.to.push_back(0.);
    model.data.rf.unrfo.push_back(0.1);
 	   }
	    i=i+1;	    
	  }//while
	  model.data.rf.rt=int(1./(model.data.rf.to[1]-model.data.rf.to[0]));
	  model.data.rf.nrfo=i;
    fprintf(stderr,"RF!!!! %d %d\n", model.data.rf.rt, model.data.rf.nrfo);
          //abort();
          if ( gauss > 0. ) model.data.rf.gaus = gauss;
          else model.data.rf.gaus = -1.;
	  mff.close();
          model.data.rf.nflag = 0; // compute misfit from the uncertainties give 
	  return 1;
	}//readrf
//-----------------------------------------------------	
	int readcov(modeldef &model,const char *name)
	{

    // fprintf(stderr,"> CALmodel: readcov\n");
    // sleep(1);

	  fstream mff;
	  string line;
    char *tname;
	  char *tname1;
 	  int iid,size,i=0,j=0,tnp;
	  vector<string> v;
	  mff.open(name);
	  getline(mff,line);
	  Split(line,v," ");
	  mff.close();  
	  model.data.rf.nflag = atoi(v[0].c_str());
    if (model.data.rf.nflag == 0 ) {    // input nn, rather than covariance matrix
    model.data.rf.nn = atof(v[1].c_str());
    }
	  else {   // read in covariance matrix
    cout<<"now read covariance matrix"<<endl;
    tname = (char*)(v[1].c_str());
    model.data.rf.nc = atof(v[2].c_str());
    tname1 = (char*)(malloc(sizeof(char)*300));
    strcpy(tname1,tname);
    
    mff.open(tname);
    getline(mff,line);
    mff.close();
    v.clear();
    Split(line,v," ");
    size=v.size();
    model.data.rf.covariance = (invcov**)malloc(sizeof(invcov*) * size);
    model.data.rf.nlen = size;
    
    cout<<tname1<<endl;
    mff.open(tname1);
    i = 0;
    while (getline(mff,line)) {
		v.clear();  		
		Split(line,v," ");
		size = v.size();
//		cout<<"now size "<<size<<endl;
		model.data.rf.covariance[i] = (invcov*)malloc(sizeof(invcov) * size);
		for (j=0;j<size;j++) {
      model.data.rf.covariance[i][j].cov0 = atof(v[j].c_str());
      }
		i++;
    }
    mff.close();
    }
//	cout<<model.data.rf.covariance[0][0].cov0;
	  return 1;
	}
//-----------------------------------------------------	
	int readmod(modeldef &model,const char *name)
	{

    // fprintf(stderr,"> CALmodel: readmod\n");
    // sleep(1);

	  string line;
	  int iid,size,i=0,tnp;
	  ifstream mff;
	  vector<string> v;
	  mff.open(name);
	  while(getline(mff,line))
		{
		 groupdef gp;
		 model.groups.push_back(gp);
     model.groups[i].np=0; //HVL: number of period?
		 model.groups[i].flag=-1; 
		 model.groups[i].thick=0.; //HVL: layers thichness?
		 model.groups[i].flagBs=-1; // flag for not use Bspline (-1)| use if = 1 later on updategroup
     model.groups[i].nlay=20; //Number of layers? Why set to 20
		 model.groups[i].vpvs=1.75; //Initial vpvs?
		 i++;
		}
	  model.ngroup=i;
	  mff.close(); //ANY OTHER WAY FOR GOING TO THE BEGINNING OF FILE?
	  mff.open(name);
	  if ( not mff.is_open())
    {
      cout<<"########model file"<<name<<"does not exist!\n exit!!\n";exit(0);
    }
	  while(getline(mff,line))
		{
		  v.clear();
		  Split(line,v," ");
		  size=v.size();
		  iid=atoi(v[0].c_str());
		  model.groups[iid].flag=atoi(v[1].c_str());//type of model: 1--layered 2--Bspline 4--gradient ! Can be capture here!
		  //model.groups[iid].thick=max(atof(v[2].c_str()),0.5);//h of each group
		  model.groups[iid].thick=atof(v[2].c_str()); // group thickness
		  tnp=atoi(v[3].c_str()); // number of layer or parameter
      model.groups[iid].np=tnp;//#of parameters
		  int tnp = model.groups[iid].np;
      model.groups[iid].ddep=atof(v[4].c_str()); // group depth step
		  //=========check the column # of input file====
		  if(model.groups[iid].flag == 4 and tnp!=2) //grad
      {
        cout<<"##########for gradient model, ONLY 2 values!\n";exit(0);
      }
		  else if ( (model.groups[iid].flag==1 and size != 5+2*tnp+11 ) or (model.groups[iid].flag==2 and size!=5+tnp+11)) //Bs or lay
			{
			  cout<<"#########wrong input for model: "<<line<< "\t\t size: "<<size<<" flag: "<<model.groups[iid].flag<<endl;
			  //for(int k=0;k<size;k++)cout<<"//  "<<v[k]<<" "<<k;
			  exit(0);
			}
		  else if(model.groups[iid].flag==5 and tnp !=1)//water layer
			{
        cout<<"########for water model, ONLY 1 values!\n";exit(0);
      }
		  model.groups[iid].value.clear(); // velocity
      model.groups[iid].ratio.clear(); // thickness
		  for(i=0;i< tnp;i++)
			{
			  model.groups[iid].value.push_back(atof(v[5+i].c_str())); // velocity values
			  if( model.groups[iid].flag==1)//layered
        {
          model.groups[iid].ratio.push_back(atof(v[5+i+tnp].c_str())); // push the ratio (total equal 1.0)
        }
      }//for
        model.groups[iid].fdvs1=(atof(v[size-1].c_str()));
        model.groups[iid].fdvs=(atof(v[size-2].c_str()));
        model.groups[iid].dvs1=(atof(v[size-3].c_str()));
        // model.groups[iid].fvs=(atof(v[size-2].c_str()));
        model.groups[iid].dvs=(atof(v[size-4].c_str()));
        //
        model.groups[iid].drho1=(atof(v[size-5].c_str()));
        // model.groups[iid].frho=(atof(v[size-2].c_str()));
        model.groups[iid].drho=(atof(v[size-6].c_str()));
        //
        model.groups[iid].vpvs1=(atof(v[size-7].c_str()));
        // model.groups[iid].fvpvs=(atof(v[size-5].c_str()));
        model.groups[iid].vpvs=(atof(v[size-8].c_str()));
        
        model.groups[iid].p_flag=(atoi(v[size-9].c_str())); // p_flag = vp flag
        model.groups[iid].q_flag=(atoi(v[size-10].c_str()));// p_flag = 0,1,2,3,4 for input, water, sediment, crust, mantle
        model.groups[iid].r_flag=(atoi(v[size-11].c_str()));// r_floag = 
        // fprintf(stderr,"vpvs1? %d \n",(atoi(v[size-7].c_str())));
        // fprintf(stderr,"vpvs? %d \n",(atoi(v[size-8].c_str())));
        // fprintf(stderr,"p_flag? %d \n",(atoi(v[size-9].c_str())));
        // fprintf(stderr,"q_flag? %d \n",(atoi(v[size-10].c_str())));
        // fprintf(stderr,"r_flag? %d \n",(atoi(v[size-11].c_str())));
        // model.groups[iid].vpvs=(atof(v[size-5].c_str()));
        
        model.groups[iid].frho = model.groups[iid].fdvs;
        model.groups[iid].fvpvs = model.groups[iid].fdvs;
        //
        model.groups[iid].frho1 = model.groups[iid].fdvs1;
        model.groups[iid].fvpvs1 = model.groups[iid].fdvs1;
        //
        model.tthick = model.tthick + model.groups[iid].thick;
        // fprintf(stderr,"model thickness update: %f %f \n",model.tthick,model.groups[iid].thick);
        // sleep(1);
        cout<<line<<endl;
		} //while getline
	  mff.close();
	  if(model.ngroup>1)
		{
      cout<<"flag0: "<<model.groups[0].flag<<"   flag1:"<<model.groups[1].flag<<endl;
    }
    fprintf(stderr,"model total thickness: %f \n",model.tthick);
    sleep(1);
	  return 1;
	}//readmod
//----------------------------------------------------- 
  int updatemodel(modeldef &model,int invtype)
  {
    /*
    Push the data of layer into the laym0
    */
    // fprintf(stderr,"> CALmodel: updatemodel\n");
    // sleep(1);
    // fprintf(stderr,"Check inversion style: %d\n",invtype);
    // sleep(1);

    int i,j,tnlay=0;
    double tvalue,tthick,tvpvs,tdep=0.,tvs,tvp,trho,tqs,tqp;
    double tvaluep,tvaluer;
	  
    model.laym0.vs.clear();
    model.laym0.vpvs.clear();
    model.laym0.vp.clear();
    model.laym0.rho.clear();
    model.laym0.qs.clear();
    model.laym0.qp.clear();
    model.laym0.thick.clear();
    model.laym0.cdep.clear();
    model.laym0.q_flag.clear();   // <-- IMPORTANT FIX

    // cout<<"updatemodel.. :( .. model.ngroup:"<<model.ngroup<<endl;
    for(i=0;i<model.ngroup;i++)
    {
      updategroup(model.groups[i]);

      for(j=0;j<model.groups[i].nlay;j++)
      // from each group, reading the value and join them into 1 layerized model
      {
        // tvpvs = model.groups[i].vpvs;
			  tvalue=model.groups[i].value1[j];
        // if (i==0) cout<<tvalue<<endl;
        tvaluep=model.groups[i].value1p[j];
        tvaluer=model.groups[i].value1r[j];
			  tthick=model.groups[i].thick1[j];
        // guard thickness
        if (!(tthick > 0.0) || std::isnan(tthick) || std::isinf(tthick)) 
        {
          tthick = 1e-6;
        }
        // explicit depth bookkeeping
        double ztop = tdep;
        tdep+=tthick;
        double zbot = tdep;
        double zmid = 0.5 * (ztop + zbot);

        // fprintf(stderr,"group depth: %f | group vs: %f | group flag %d! \n",tdep,tvalue,model.groups[i].flag);
        // sleep(1);
        // Q assignment
        double tqs = 500.0, tqp = 1000.0;
			  if(model.groups[i].q_flag == 1)//water layer
        {
          tqs=0.;
          tqp=57822.;
          tvalue = 0.;
        } 
        else if( model.groups[i].q_flag == 2) //layer2 sediment
        {
          tqp=160.;
          tqs=80.;
        }
			  else if (model.groups[i].q_flag == 0)     // input Q model
        {  
          // safer indexing than (int)(tdep)
          int idx = (int)std::floor(zmid); // or zbot, but zmid is smoother
          idx = std::max(0, std::min(idx, (int)model.Qmodel.Qmod.size() - 1));
          tqs = model.Qmodel.Qmod[idx];
          tqp = tqs * 2.0;
        }
        else if (model.groups[i].q_flag == 3 ) // crust Q
        {
          tqs = 599.99;
          tqp = 1368.02;
        }
        else if ( model.groups[i].q_flag == 4 ) // mantle Q
        {
          tqs = 417.59;
          tqp = 1008.71; 
        }
        else if ( model.groups[i].q_flag == 5 ) // mantle Q (scaled)
        {
          int idx = (int)std::floor(zmid);
          idx = std::max(0, std::min(idx, (int)model.Qmodel.Qmod.size() - 1));
          tqs = model.Qmodel.Qmod[idx] / 2.0;
          tqp = tqs * 2.0;
        }
        // vp/vs safe
        tvpvs = (tvalue > 1e-12) ? (tvaluep/tvalue) : 0.0;  
			  // 
			  model.laym0.vpvs.push_back(tvpvs);
			  model.laym0.vs.push_back(tvalue);
			  model.laym0.vp.push_back(tvaluep);
			  model.laym0.rho.push_back(tvaluer);
			  model.laym0.qs.push_back(tqs);
			  model.laym0.qp.push_back(tqp);
			  model.laym0.thick.push_back(tthick);
        model.laym0.q_flag.push_back(model.groups[i].q_flag);

        // CHOOSE what you want cdep to mean:
        model.laym0.cdep.push_back(zbot);   // bottom depth (your old behavior)
        // model.laym0.cdep.push_back(zmid); // midpoint depth (often better for plotting)
        // model.laym0.cdep.push_back(ztop); // top depth
        tnlay++;

			}//for j
		  
    }//for i
    
      model.laym0.nlayer=tnlay;
      // for(i=0;i<model.laym0.nlayer;i++)
      // {
      //   fprintf(stderr,"%d | vs: %f | thick: %f | depth: %f | qflag: %f\n",i,model.laym0.vs[i],model.laym0.thick[i],model.laym0.cdep[i],model.laym0.q_flag[i]);
      // }
      // fprintf(stderr,">> ------------------------ <<\n");
      // sleep(2);
      model.flag=1;//updated, from layered vs,h get other layered para
      //cout<<"got to end of updatemodel, so return1"<<endl;
	  return 1;
  }// updatemodel
//-----------------------------------------------------	
// ---------------------------------------------------
  int readindata(modeldef &model, const char* ifname)
  {
    fstream mff;
    string line;
    int kid;
    vector<string> v;
    // ------------------------------------------------
    mff.open(ifname);
    if ( not mff.is_open()) {cout<<"########indata file "<<ifname<<" does not exist!\n exit!!\n";exit (0);}
    kid=0;
    model.indata.gvflag=0; model.indata.phflag=0; model.indata.hvflag=0; model.indata.rfflag=0;
    model.indata.hpflag=0; model.indata.hpw=0.;  // higher-mode phase (default off; lines 10-11 optional)
    while(getline(mff,line))
		{
      v.clear();  		
		  Split(line,v," ");
      kid++;
      if (kid==1) // 1st line phase file
      {
        model.indata.phflag=(atoi(v[0].c_str()));
      }
      if (kid==2) // 2nd line group file
      {
        model.indata.gvflag=(atoi(v[0].c_str()));
      }
      if (kid==3) // 3rd line ellip file
      {
        model.indata.hvflag=(atoi(v[0].c_str()));
      }
      if (kid==4) // 4th line RF file
      {
        model.indata.rfflag=(atoi(v[0].c_str()));
      }
      if (kid==5) // 5th line weighting flag
      {
        model.indata.isequalweight=(atoi(v[0].c_str()));
      }
      if (kid==6) // 6th line weighting for ph (if not equal)
      {
        model.indata.phw=(atof(v[0].c_str()));
      }
      if (kid==7) // 7th line weighting for GV (if not equal)
      {
        model.indata.gvw=(atof(v[0].c_str()));
      }
      if (kid==8) // 8th line weighting for HV (if not equal)
      {
        model.indata.hvw=(atof(v[0].c_str()));
      }
      if (kid==9) // 9th line weighting for RF (if not equal)
      {
        model.indata.rfw=(atof(v[0].c_str()));
      }
      if (kid==10) // 10th line higher-mode phase flag (optional)
      {
        model.indata.hpflag=(atoi(v[0].c_str()));
      }
      if (kid==11) // 11th line higher-mode phase weight (optional)
      {
        model.indata.hpw=(atof(v[0].c_str()));
      }
    }
    mff.close();
    fprintf(stderr,"data avaiable phase - group - hv - rf %d %d %d %d\n",model.indata.phflag,model.indata.gvflag,model.indata.hvflag,model.indata.rfflag);
    // sleep(10);
    return 1;
  } // readindata
  
///-----------------------------------------------------	
	int compute_misfit(modeldef &model,double inp)
	{
    // fprintf(stderr,"> CALmodel: compute_misfit\n");
    // sleep(1);
    /*
    Calculate the models misfit of synthetic model and observed data.
    Model parameters:
    + tempv1: temporary phase velocity misfit
    + tempv2: temporary group velocity misfit
    + tempv3: temporary ellipticity misfit
    + tmisfit1: misfit of phase velocity and ellipticity
    + tmisfit2:
    + tS:
    + tL1: model likelihood
    + td
    + p:
    + nn:
    + tmfe: mitfit for ellipticity
    + tmfp: mitfit for phase velocity
    + Se: = tempv3
    + Sp: = tempv1
    + S1: = tempv1 + tempv2 + tempv3
    + S2: = tempv2
    */
    int i,j,k;
    int igv,iph,ihv,irf,ihp,nidata=0,nidatadisp=0;
    int isew=1; // Equal weighting flag
    double phw,gvw,hvw,rfw,hpw; // Weighting for each
	  double tempv1=0.,tempv2=0.,tempv3=0.,tmisfit1,tmisfit2,tS,tL1,tp,td,p,nn;
	  double tmfe,tmfp,tmfg,Se,Sp,Sg;//EMB add's for misfit of H/V and phase vel
	  double tmfhp,Sph;//higher-mode (1st) phase velocity misfit / weighted residual
    double deltad[1000],un[1000],S1,S2;
	  S1 = 0.;
	  S2 = 0.;
	  model.data.disp.gL=1.;
	  model.data.disp.pL=1.;
    model.data.disp.eL=1.;
    // fprintf(stderr,"data avaiable in compute_misfit phase - group - hv - rf %d %d %d %d\n",model.indata.phflag,model.indata.gvflag,model.indata.hvflag,model.indata.rfflag);
    // sleep(10);
    // take the file avaiable contritution weighting.
    iph=model.indata.phflag;
    igv=model.indata.gvflag;
    ihv=model.indata.hvflag;
    irf=model.indata.rfflag;
    ihp=model.indata.hpflag;
    nidata=nidata+iph+igv+ihv+irf+ihp; //number of avaiable data
    nidatadisp=nidatadisp+iph+igv+ihv+ihp; // number of dispersion data (without RF)
    // fprintf(stderr,"%d %d\n",model.data.disp.npper,model.data.disp.ngper,model.data.disp.neper);
    isew=model.indata.isequalweight;
    phw=model.indata.phw;
    gvw=model.indata.gvw;
    hvw=model.indata.hvw;
    rfw=model.indata.rfw;
    hpw=model.indata.hpw;
    // fprintf(stderr,"Weighting flag: %d %f %f %f %f\n",isew,phw,gvw,hvw,rfw);
    // sleep(1);
    
    // misfit for Phase velocity dispersion //////////////////////
    tempv1 = 0.;
	  for(i=0;i<model.data.disp.npper;i++)
		{
      // fprintf(stderr,"misfit: %g\n",tempv1);
      tempv1=tempv1+pow((model.data.disp.pvelo[i]-model.data.disp.pvel[i])/model.data.disp.unpvelo[i],2);
      } // for
    if(tempv1>0.)
    {
      model.data.disp.pmisfit=sqrt(tempv1/model.data.disp.npper);
      tS=tempv1;
      // while(tS>50.)tS=sqrt(tS*30.);
      if(tS>30.) tS=sqrt(tS*30.);
      if(tS>30.) tS=sqrt(tS*30.);
      // if(tS>30.) tS=sqrt(tS*30.);
      // if(tS>30.) tS=sqrt(tS*30.);
      // if(tS>30.) tS=sqrt(tS*30.);
      // fprintf(stderr,"tS: %g\n",tS);
      model.data.disp.pL=exp(-0.5*tS);
      // fprintf(stderr,"tS: %g pL: %g\n",tS,model.data.disp.pL);
      }//if
    else
    {
      model.data.disp.pL=1.;
    }

    // misfit for Group velocity dispersion //////////////////////
    tempv2 = 0.;
    for(i=0;i<model.data.disp.ngper;i++)
    {
      tempv2=tempv2+pow((model.data.disp.gvelo[i]-model.data.disp.gvel[i])/model.data.disp.ungvelo[i],2);
      } //for
    if(tempv2>0.)
    {
      model.data.disp.gmisfit=sqrt(tempv2/model.data.disp.ngper);
      tS=tempv2;
      if(tS>30) tS=sqrt(tS*30);
      if(tS>30) tS=sqrt(tS*30);
      model.data.disp.gL=exp(-0.5*tS);
      }
    else
    {
      model.data.disp.gL=1.;
      } //if tempv2
    // misfit for Elipticity ////////////// //////////////////////
    tempv3 = 0.;
    for (i=0;i<model.data.disp.neper;i++) 
    {
      tempv3=tempv3+pow((model.data.disp.evelo[i]-model.data.disp.evel[i])/model.data.disp.unevelo[i],2);
      }// for
    if (tempv3>0.)
    {
      model.data.disp.emisfit=sqrt(tempv3/model.data.disp.neper);
      tS = tempv3;
      if (tS>30.) tS = sqrt(tS*30);
      if (tS>30.) tS = sqrt(tS*30);
      model.data.disp.eL=exp(-0.5*tS);
      }
	  else
    {
      model.data.disp.eL=1.;
      }//if
    // misfit for Higher-mode (1st) phase velocity ////////////////
    double tempv1h=0.;
    for (i=0;i<model.data.disp.nhpper;i++)
    {
      tempv1h=tempv1h+pow((model.data.disp.hpvelo[i]-model.data.disp.hpvel[i])/model.data.disp.unhpvelo[i],2);
    }//for
    if (tempv1h>0.)
    {
      model.data.disp.phmisfit=sqrt(tempv1h/model.data.disp.nhpper);
      tS=tempv1h;
      if (tS>30.) tS=sqrt(tS*30.);
      if (tS>30.) tS=sqrt(tS*30.);
      model.data.disp.hpL=exp(-0.5*tS);
    }
    else
    {
      model.data.disp.hpL=1.;
    }
    //////////////////////////////////////////////////////////////
    // fprintf(stderr,"pL (b0): %g\n",model.data.disp.pL);
    tmisfit1=sqrt(((tempv1*iph)+(tempv2*igv)+(tempv3*ihv))/(model.data.disp.npper+model.data.disp.neper+model.data.disp.ngper)); //original, no gv
    tmfe=sqrt((tempv3*ihv)/model.data.disp.neper);//H/V misfit only
    tmfp=sqrt((tempv1*iph)/model.data.disp.npper);// Vph misfit only
    tmfg=sqrt((tempv2*igv)/model.data.disp.ngper);// Vph misfit only
    tmfhp=sqrt((tempv1h*ihp)/model.data.disp.nhpper);// higher-mode phase misfit only

    if (std::isnan(tmisfit1)) tmisfit1=999.;
    if (std::isnan(tmfp)) tmfp=999.;
    if (std::isnan(tmfg)) tmfg=999.;
    if (std::isnan(tmfe)) tmfe=999.;
    if (std::isnan(tmfhp)) tmfhp=999.;

    
    //tmisfit1=sqrt((tempv1+tempv2+tempv3)/(model.data.disp.npper+model.data.disp.ngper+model.data.disp.neper)); //original, with group velocity
	  tL1=model.data.disp.gL*model.data.disp.pL*model.data.disp.eL*model.data.disp.hpL;//likelihood

	  S1 = tempv1 + tempv2 + tempv3 + tempv1h;//original,no weights
    // S1 = fabs(tempv1) + fabs(tempv2) + fabs(tempv3);// HVL try absolute values only
    Sp=tempv1;//phase velocity
    Sg=tempv2;//group velocity
    Se=tempv3;//ellipticity
    Sph=tempv1h;//higher-mode (1st) phase velocity
    //
    if (std::isnan(S1)) S1=999.;
    if (std::isnan(Sp)) Sp=999.;
    if (std::isnan(Sg)) Sg=999.;
    if (std::isnan(Se)) Se=999.;
    if (std::isnan(Sph)) Sph=999.;
    //
	  model.data.disp.L=tL1;//not terribly important, model.data.L calculated below
	  // model.data.disp.misfit=(0.667)*tmfe+(0.333)*tmfp;//EMB edit for weighted misfit
    // model.data.disp.misfit=(1.0*tmfp);//HVL edit for Vph misfit only
    if (isew==0)
    {
      model.data.disp.misfit=((phw*tmfp))+((gvw*tmfg))+((hvw*tmfe))+((hpw*tmfhp));//HVL edit for misfit based on data avaiable
    }
    else
    {
      model.data.disp.misfit=((phw*tmfp/(phw+gvw+hvw+hpw)))+((gvw*tmfg/(phw+gvw+hvw+hpw)))+((hvw*tmfe/(phw+gvw+hvw+hpw)))+((hpw*tmfhp/(phw+gvw+hvw+hpw)));//HVL edit for misfit based on given weight
    }
    
	  // model.data.disp.misfit=(0.1)*tmfe+(0.9)*tmfp;
	  //model.data.disp.misfit=tmisfit1; //original,non-weighted, misfit
	  tempv2=0.;
	  k=0;
    tp=0.;
	  tmisfit2 = 0.;
    
    // compute misfit and L for RF ///////////////////////////////////////////////////////////////////////////////////
	  if (model.data.rf.nflag < 0.5 )  
    {
      S2 = 0.;
      for(i=0;i<model.data.rf.nrfo;i++)
      {
        for(j=0;j<model.data.rf.tn.size();j++)
        {
          if(model.data.rf.to[i]==model.data.rf.tn[j] and model.data.rf.to[i]>=0)
          {
            tempv2 = tempv2 + pow((model.data.rf.rfo[i]-model.data.rf.rfn[j])/model.data.rf.unrfo[i],2);
            k=k+1;
            }//if
          }//forj
      }//fori
	    S2 = tempv2;
      tmisfit2=sqrt(tempv2/k);
      if (std::isnan(tmisfit2)) tmisfit2 = 999.;
      if (std::isnan(S2)) S2 = 999.;
      //while (tempv2 > 50.) tempv2 = sqrt(tempv2*30);
      if (tempv2 > 30.) tempv2 = sqrt(tempv2*30);
      if (tempv2 > 30.) tempv2 = sqrt(tempv2*30);
      // if (tempv2 > 30.) tempv2 = sqrt(tempv2*30);
      // if (tempv2 > 30.) tempv2 = sqrt(tempv2*30);
      // if (tempv2 > 30.) tempv2 = sqrt(tempv2*30);
      if (tempv2>0.)
      {
      model.data.rf.L=exp(-0.5*tempv2);
      }
      else
      {
        model.data.rf.L=1.;
      }
      model.data.rf.misfit=tmisfit2;
      // cout<<model.data.rf.L<<" "<<tempv2<<" "<<endl;
      }
    else
    {
      // compute L by using the the covariance matrix which is stored in model.data.rf.covariance
      // cout<<"ok here! now compute misfit for RF"<<endl;
      S2 = 0.;
      for(i=0;i<model.data.rf.nrfo;i++)
      {
        for(j=0;j<model.data.rf.tn.size();j++)
        {
          if(fabs(model.data.rf.to[i]-model.data.rf.tn[j])<0.001 and model.data.rf.to[i]<=10 and model.data.rf.to[i]>=0)
          {
            deltad[k] = (model.data.rf.rfo[i]-model.data.rf.rfn[j]);  // delta d
            un[k] = model.data.rf.unrfo[i];
            k=k+1;
            // break;
            }//if
          }//forj
        }//fori 
        // cout<<"compute over! k "<<k<<endl;
	    if (k != model.data.rf.nlen) 
      {
        fprintf(stderr,"length problem!!! exit!!!! %d %d\n",k,model.data.rf.nlen);
        abort();
        }
      tmisfit2 = 0.;
      for (i=0;i<k;i++) 
      {
        tmisfit2 = tmisfit2 + pow(deltad[i]/un[i],2);
        for (j=0;j<k;j++) 
        {
          tempv2 = tempv2 + deltad[i]* model.data.rf.covariance[i][j].cov0 * deltad[j] * model.data.rf.nc;
          }
          // cout<<tempv2<<" check tempv2 "<<i<<endl;
          }   
          // abort();
      if (tempv2 > 30.) tempv2 = sqrt(tempv2*30);
      if (tempv2 > 30.) tempv2 = sqrt(tempv2*30);
      if (tempv2>0.)
      {
        model.data.rf.L=exp(-0.5*tempv2);
      }
      else
      {
        model.data.rf.L=1.;
      }
      
      cout<<"misfitfor RF"<<tempv2<<" "<<model.data.rf.L<<endl;
      
      tmisfit2 = sqrt(tmisfit2/k);
      if (std::isnan(tmisfit2)) tmisfit2 = 999.;
      model.data.rf.misfit= tmisfit2;
      // S2 = fabs(tempv2); // HVL 20231016 trigger abs value
      S2 = tempv2;
      if (S2<0.) 
      {
        cout<<"something is wrong !!! now need to write down the deltad"<<S2<<endl;
        FILE *ffd;
        ffd = fopen("deltad.dat","w");
        for (i=0;i<k;i++) 
        {
          fprintf(ffd,"%g %g %g %g\n", 0.025*i,deltad[i],model.data.rf.rfo[i],model.data.rf.rfn[i]);
          }
          fclose(ffd);
          abort();
        }
      // cout<<"compute for RF over!"<<endl;
    }
	  // combine misfit and liability
    // fprintf(stderr,"%lf %lf\n",S1,S2);
    // fprintf(stderr,"inp check 1: %f \n",inp);
    tS=0;
    if (inp>1) inp=1.; // needs to be changed
    // fprintf(stderr,"inp check 2: %f \n",inp);
    // sleep(5);
    // tS = 0.5*(inp*S1+(1-inp)*S2); // equal weighting; for joint inversion, inp should be 0.5; for RF only, inp == 0.
    //tS = 0.5*((0.67*Se)+(0.33*Sp)); // EMB edited, non-equal weighting (not using RF)
    // tS = 0.5*((1.0*Sp)); 
    if (isew==0)
    {
      tS = 0.5*((hvw*Se)+(phw*Sp)+(gvw*Sg)+(rfw*S2)+(hpw*Sph)); // HVL none-equal weighting and based on data avaiable
    }
    else
    {
      tS = 0.5*((hvw*Se/(phw+gvw+hvw+rfw+hpw))+(phw*Sp/(phw+gvw+hvw+rfw+hpw))+(gvw*Sg/(phw+gvw+hvw+rfw+hpw))+(rfw*S2/(phw+gvw+hvw+rfw+hpw))+(hpw*Sph/(phw+gvw+hvw+rfw+hpw))); // HVL equal weighting and based on data avaiable
    }
    // tS = 0.5*((0.6*Sp)+(0.2*Se)+(0.2*S2)); // HVL edited, non-Vph
    // while (tS>50.)tS=sqrt(tS*30);
    if(tS>30.) tS=sqrt(tS*30);
    if(tS>30.) tS=sqrt(tS*30);
	  // if(tS>30.)tS=sqrt(tS*30);
	  // if(tS>30.)tS=sqrt(tS*30);
	  // if(tS>30.)tS=sqrt(tS*30);
    model.data.L = exp(-tS);
    
    // model.data.misfit=(0.9*tmfp)+(0.1*tmfe);//EMB edited for non-equal weighting
    // model.data.misfit=(1.0*tmfp);
    // model.data.misfit=(0.33*tmfp)+(0.67*tmfe);
    // model.data.misfit=(0.5*tmfe)+(0.5*tmisfit2);//EMB edited for non-equal weighting
    if (isew==0)
    {
      model.data.misfit=(phw*tmfp)+(hvw*tmfe)+(gvw*tmfg)+(rfw*tmisfit2)+(hpw*tmfhp);// HVL none-equal weighting and based on data avaiable
    }
    else
    {
      if (abs(phw+gvw+hvw+rfw+hpw-1) > 1e-3)
      {
        fprintf(stderr,"Total weighting less than 1, stop!");
        exit(EXIT_FAILURE);
      }
      model.data.misfit=(phw*tmfp/(phw+gvw+hvw+rfw+hpw))+(hvw*tmfe/(phw+gvw+hvw+rfw+hpw))+(gvw*tmfg/(phw+gvw+hvw+rfw+hpw))+(rfw*tmisfit2/(phw+gvw+hvw+rfw+hpw))+(hpw*tmfhp/(phw+gvw+hvw+rfw+hpw));// HVL equal weighting and based on data avaiable
    }
    
    // model.data.misfit=(0.34*tmfp)+(0.33*tmfe)+(0.33*tmisfit2);// HVL more weight for Vph and RF
	  // model.data.misfit=inp*tmisfit1+(1-inp)*tmisfit2;//original, equal weighting as applicable
    // fprintf(stderr,"%lf %lf %lf %lf %lf %lf\n",model.data.misfit, model.data.disp.misfit, tmfp, tmfg, tmfe, tmisfit2);
    // sleep(2);
    if (std::isnan(model.data.misfit)) model.data.misfit=9999.;
    if ((model.data.misfit) < 0.)
    {
      model.data.misfit=9999.;
    }  
      
    // cout<<"compute all misfit over!"<<endl;
    return 1;
	}//compute misfit

//-----------------------------------------------------	
	int write_model(modeldef &model,char *name, char *outdir, int flagw)
	{

    // fprintf(stderr,"> CALmodel: write_model\n");
    // sleep(1);

	  char namerf[300],namedispp[300],namedispg[300],namemod[300];
    char namemod1[300],namedispe[300],nameinfo[300],namemod2[300];
    char namedisphp[300];
	  FILE *ff;
	  int i,j,k,nlayer;
	  double dep=0.;
	  sprintf(namerf,"%s/%s.rf",outdir,name);
	  sprintf(namedispp,"%s/%s.p.disp",outdir,name);
	  sprintf(namedispg,"%s/%s.g.disp",outdir,name);
	  sprintf(namedispe,"%s/%s.e.disp",outdir,name);
	  sprintf(namedisphp,"%s/%s.hp.disp",outdir,name); // higher-mode (1st) phase
	  sprintf(namemod,"%s/%s.mod",outdir,name);
	  sprintf(namemod1,"%s/%s.mod.q",outdir,name);
    sprintf(namemod2,"%s/%s.mod.group",outdir,name); // similar input model for mat
	  sprintf(nameinfo,"%s/%s.info",outdir,name);

	  
	  ff=fopen(namemod2,"w");
          //cout<<model.ngroup<<endl;
          //cout<<model.groups.size()<<endl;
	  //for (i=0;i<3;i++) {
	  for (i=0;i<model.groups.size();i++) {
	    fprintf(ff,"%d %d %lf %d ",i,model.groups[i].flag,model.groups[i].thick,model.groups[i].np);
	    for (j=0;j<model.groups[i].np;j++) {
              fprintf(ff,"%lf ",model.groups[i].value[j]);
              }
            if(model.groups[i].flag==1) { // layerized, add ratio
	      fprintf(ff,"%lf ",model.groups[i].ratio[j]);
              }
            fprintf(ff,"%d %d %d %g %g %g %g %g %g %g %g\n",model.groups[i].r_flag,model.groups[i].q_flag,model.groups[i].p_flag,model.groups[i].vpvs,model.groups[i].vpvs1,model.groups[i].drho,model.groups[i].drho1,model.groups[i].dvs,model.groups[i].dvs1,model.groups[i].fdvs,model.groups[i].fdvs1);             
            }
	  fclose(ff);
	  

	  ff=fopen(nameinfo,"w");
	  fprintf(ff,"# ALL %lf %lf\n",model.data.misfit,model.data.L);
	  fprintf(ff,"# DISP %lf %lf\n",model.data.disp.misfit,model.data.disp.L);
	  if (model.data.disp.npper>0) fprintf(ff,"# PHASE %lf %lf\n",model.data.disp.pmisfit,model.data.disp.pL);
	  if (model.data.disp.ngper>0) fprintf(ff,"# GROUP %lf %lf\n",model.data.disp.gmisfit,model.data.disp.gL);
	  if (model.data.disp.neper>0) fprintf(ff,"# ELLIP %lf %lf\n",model.data.disp.emisfit,model.data.disp.eL);
	  if (model.data.disp.nhpper>0) fprintf(ff,"# HMODE %lf %lf\n",model.data.disp.phmisfit,model.data.disp.hpL);
	  if (model.data.rf.nrfo>0) fprintf(ff,"# RF %lf %lf\n",model.data.rf.misfit,model.data.rf.L);
	  fclose(ff);

	  ff=fopen(namerf,"w");
	  
	  for (i=0;i<min(model.data.rf.rfo.size(),model.data.rf.rfn.size());i++)
		{
		  if(model.data.rf.tn.size()>0)
			fprintf(ff,"%g %g %g %g %g\n",model.data.rf.tn[i],model.data.rf.rfn[i],model.data.rf.to[i],model.data.rf.rfo[i],model.data.rf.unrfo[i]);
		}//for i
	 fclose(ff);

	 if(model.data.disp.npper>0)
	 {		
	 ff=fopen(namedispp,"w");
	 for(i=0;i<model.data.disp.npper;i++)	
		{
		  if(model.data.disp.pvel.size()>0)
			fprintf(ff,"%g %g %g %g\n",model.data.disp.pper[i],model.data.disp.pvel[i],model.data.disp.pvelo[i],model.data.disp.unpvelo[i]);
		}//fori  
	  fclose(ff);
	  }//if

	
	 if(model.data.disp.ngper>0)
	 {		
	 ff=fopen(namedispg,"w");
	 for(i=0;i<model.data.disp.ngper;i++)	
		{
		  if(model.data.disp.gvel.size()>0)
			fprintf(ff,"%g %g %g %g\n",model.data.disp.gper[i],model.data.disp.gvel[i],model.data.disp.gvelo[i],model.data.disp.ungvelo[i]);
		}//fori  
	  fclose(ff);
	  }//if
          //cout<<"neper: "<<model.data.disp.neper<<endl;
          if (model.data.disp.neper>0)
            {
            ff=fopen(namedispe,"w");
            for(i=0;i<model.data.disp.neper;i++)                {
               if(model.data.disp.evel.size()>0)
                  fprintf(ff,"%g %g %g %g\n",model.data.disp.eper[i],model.data.disp.evel[i],model.data.disp.evelo[i],model.data.disp.unevelo[i]);
               }//fori
            fclose(ff);
            }
          // higher-mode (1st) phase velocity: per  pred  obs  unc
          if (model.data.disp.nhpper>0)
            {
            ff=fopen(namedisphp,"w");
            for(i=0;i<model.data.disp.nhpper;i++)                {
               if(model.data.disp.hpvel.size()>0)
                  fprintf(ff,"%g %g %g %g\n",model.data.disp.hpper[i],model.data.disp.hpvel[i],model.data.disp.hpvelo[i],model.data.disp.unhpvelo[i]);
               }//fori
            fclose(ff);
            }

          // write continuous model rather than layered model
	  ff=fopen(namemod,"w");
          dep = 0.;
          //fprintf(stderr,"how many groups??? %d groups\n",model.ngroup);
          
          for (i=0;i<model.ngroup;i++) {
            // fprintf(stderr,"model.groups[i].flag: %d, %d \n",i,model.groups[i].flag);
            // sleep(2);
              if (model.groups[i].flag == 5) // water layer
                  {
                  fprintf(ff,"%g 0.\n",dep);
                  dep = dep + model.groups[i].thick;
                  fprintf(ff,"%g 0.\n",dep);
                  }
              if (model.groups[i].flag == 4) // grad
                  {
                  // fprintf(ff,"%g %g\n",dep,model.groups[i].value[0]);
                  fprintf(ff,"%g %g %d\n",dep,model.groups[i].value[0],model.groups[i].q_flag);
                  dep = dep + model.groups[i].thick;
                  // fprintf(ff,"%g %g\n",dep,model.groups[i].value[1]);
                  fprintf(ff,"%g %g %d\n",dep,model.groups[i].value[1],model.groups[i].q_flag);
                  }
              if (model.groups[i].flag == 1) // layered
                  {
                  k = model.groups[i].np;
                  nlayer = k;
                  //fprintf(stderr,"how many layers??? %d layers %d\n",nlayer,k);
                  for (j=0;j<nlayer;j++) {
                    //  fprintf(ff,"%g %g\n",dep,model.groups[i].value[j]);
                    fprintf(ff,"%g %g %d\n",dep,model.groups[i].value[j],model.groups[i].q_flag);
                     dep = dep + model.groups[i].thick*model.groups[i].ratio[j];
                    //  fprintf(ff,"%g %g\n",dep,model.groups[i].value[j]);
                    fprintf(ff,"%g %g %d\n",dep,model.groups[i].value[j],model.groups[i].q_flag);
                     }
                  }
              if (model.groups[i].flag == 2) // B-splines
                  {
                  fprintf(ff,"%g %g %d\n",dep,model.groups[i].value1[0],model.groups[i].q_flag);
                  nlayer = model.groups[i].nlay;
                  //dep = dep + model.groups.thick1[0]/2;
                  for (j=0;j<nlayer-1;j++) {
                     dep = dep + model.groups[i].thick1[j];
                    //  fprintf(ff,"%g %g\n",dep,(model.groups[i].value1[j]+model.groups[i].value1[j+1])/2);
                    fprintf(ff,"%g %g %d\n",dep,(model.groups[i].value1[j]+model.groups[i].value1[j+1])/2,model.groups[i].q_flag);
                     }
                  dep = dep + model.groups[i].thick1[nlayer-1];
                  // fprintf(ff,"%g %g\n",dep,model.groups[i].value1[nlayer-1]);
                  fprintf(ff,"%g %g %d\n",dep,model.groups[i].value1[nlayer-1],model.groups[i].q_flag);
                  }
              }
          /*
	  for(i=0;i<model.laym0.nlayer;i++) //write layered model in stair-step
		{
		  fprintf(ff,"%g %g\n",dep,model.laym0.vs[i]);
		  dep=dep+model.laym0.thick[i];
		  fprintf(ff,"%g %g\n",dep,model.laym0.vs[i]);
		}
          */
	  fclose(ff);

          if(flagw == -1) return 1;           

	  ff=fopen(namemod1,"w");
          dep = 0.;
          for(i=0;i<model.laym0.nlayer;i++) //write layered model in stair-step
                {
                  fprintf(ff,"%g %g %g %g %g %g\n",dep,model.laym0.vs[i],model.laym0.vp[i],model.laym0.rho[i],model.laym0.qs[i], model.laym0.thick[i]);
                  dep=dep+model.laym0.thick[i];
                  fprintf(ff,"%g %g %g %g %g %g\n",dep,model.laym0.vs[i],model.laym0.vp[i],model.laym0.rho[i],model.laym0.qs[i], model.laym0.thick[i]-model.laym0.thick[i]);
                }
          fclose(ff);

	  return 1; 
	}//write model
//-----------------------------------------------------	
	int goodmodel(modeldef &model, vector<int> vmono, vector<int> vgrad, int invtype)
	{	

    // fprintf(stderr,"> CALmodel: goodmodel\n");
    // sleep(1);
// ==========================================================================================================
	  int i,j,ij;
    int ised,jm,js,jc;
    double depthcheck=0;
    double sed_thick=0;
    double crust_thick=0;
	  vector<int>::iterator id;
    constexpr double eps = 1e-9;  // tolerance for float noise
    const char *invtypechar;
// ==========================================================================================================
// Check the forward data first. If any wrong data then discard the model
// ==========================================================================================================
    
    // Phase velocity dispersion check//////////////////////
    if (model.indata.phflag==1)
    {
      for(i=0;i<model.data.disp.npper;i++)
		  {
      // fprintf(stderr,"misfit: %g\n",tempv1);
      double v = model.data.disp.pvelo[i];
      if ( (!std::isfinite(v)) || v < 0.2 || v > 4.9 ) return 0;

      } // for
    }
    // Group velocity dispersion check//////////////////////
    if (model.indata.gvflag==1)
    {
      for(i=0;i<model.data.disp.ngper;i++)
		  {
      // fprintf(stderr,"misfit: %g\n",tempv1);
      double v = model.data.disp.gvelo[i];
      if ( (!std::isfinite(v)) || v < 0.2 || v > 4.9 ) return 0;

      } // for
    }
    // Elipticity check//////////////////////
    if (model.indata.hvflag==1)
    {
      for(i=0;i<model.data.disp.neper;i++)
		  {
      // fprintf(stderr,"misfit: %g\n",tempv1);
      double v = model.data.disp.evelo[i];
      if ( (!std::isfinite(v)) || v < 0.5 || v > 4.0 ) return 0;

      } // for
    }
// ==========================================================================================================
// Check crustal layer is increasing Vs with depth (if not, discard model) - this specific for CHT model only
// ===================== CRUST ONLY: Vs monotonic increase =====================
  for (int ig = 0; ig < model.ngroup; ++ig)
  {
    if (model.groups[ig].q_flag != 3) continue;  // Crust only

    const std::vector<double> &vs = model.groups[ig].value1;
    if (vs.size() < 2) continue;

    for (size_t k = 1; k < vs.size(); ++k)
    {
      if (vs[k] + eps < vs[k - 1])
      {
        fprintf(stderr,"CRUST Vs not monotonic (within group)! group=%d k=%zu  %.6f > %.6f\n",ig, k, vs[k - 1], vs[k]);
        return 0;
      }
    }
  }
  // If crust is split into multiple consecutive groups (all q_flag==3),
  // also enforce no negative jump BETWEEN crust groups:
  for (int ig = 1; ig < model.ngroup; ++ig)
  {
    if (model.groups[ig - 1].q_flag == 3 && model.groups[ig].q_flag == 3)
    {
      if (!model.groups[ig - 1].value1.empty() && !model.groups[ig].value1.empty())
      {
        double v_prev = model.groups[ig - 1].value1.back();
        double v_curr = model.groups[ig].value1.front();
        if (v_curr + eps < v_prev)
        {
          fprintf(stderr,"CRUST Vs negative jump (between crust groups)! g%d->g%d  %.6f > %.6f\n",ig - 1, ig, v_prev, v_curr);
          return 0;
        }
      }
    }
  }
// ==========================================================================================================
// Check mantle layer is increasing Vs with depth (if not, discard model)
// ===================== MANLE ONLY: Vs monotonic increase =====================

  for (int ig = 0; ig < model.ngroup; ++ig)
  {
    if (model.groups[ig].q_flag != 4) continue;  // Mantle only

    const std::vector<double> &vs = model.groups[ig].value1;
    if (vs.size() < 2) continue;

    for (size_t k = 1; k < vs.size(); ++k)
    {
      if (vs[k] + eps < vs[k - 1])
      {
        fprintf(stderr,"MANTLE Vs not monotonic (within group)! group=%d k=%zu  %.6f > %.6f\n",ig, k, vs[k - 1], vs[k]);
        return 0;
      }
    }
  }
  // If mantle is split into multiple consecutive groups (all q_flag==4),
  // also enforce no negative jump BETWEEN mantle groups:
  for (int ig = 1; ig < model.ngroup; ++ig)
  {
    if (model.groups[ig - 1].q_flag == 4 && model.groups[ig].q_flag == 4)
    {
      if (!model.groups[ig - 1].value1.empty() && !model.groups[ig].value1.empty())
      {
        double v_prev = model.groups[ig - 1].value1.back();
        double v_curr = model.groups[ig].value1.front();
        if (v_curr + eps < v_prev)
        {
          fprintf(stderr,"MANTLE Vs negative jump (between crust groups)! g%d->g%d  %.6f > %.6f\n",ig - 1, ig, v_prev, v_curr);
          return 0;
        }
      }
    }
  }
  // ============================================================================ 
    // end of check crust-mantle layer is increasing Vs with depth
// ===========================================================================================================
//  fprintf(stderr, "check inversion style: %d \n", invtype); // check point for the inversion type
// ******************************************** //
      for (i=0;i<model.ngroup;i++) // all group
      {   // abs vel
        if (model.groups[i].q_flag==4) // none mantle (all other group)
        {
          continue;
        }
        else
        {
        for (j=0;j<model.groups[i].nlay;j++)
        {
          // if ((model.groups[i].value1[j] > 4.9)) 
          if ((model.groups[i].value1[j] > 4.9) || (model.groups[i].value1[j] < 0.1)) 
        {
          cout<<"4.9 crit... i = "<<i<<" model.groups[i].value1[j]= "<<model.groups[i].value1[j]<<endl;
          fprintf(stderr,"4.9 wrong!!! %g %d %d\n",model.groups[i].value1[j],i,j);
          return 0;
          }
        } 
        }
      } //end of non-manlte Vs < 4.9km/s
// ******************************************** //
      for (i=0;i<model.ngroup;i++) // all group
      {   // abs vel
        if (model.groups[i].q_flag==2) // sediment
        {
          fprintf(stderr,"Check sediment velocity!!! \n");
          for (j=0;j<model.groups[i].nlay;j++)
          {
            // if ((model.groups[i].value1[j] > 4.9))
            if ((model.groups[i].value1[j] > 4.9) || (model.groups[i].value1[j] < 0.1)) 
            {
              cout<<"4.9 crit... i = "<<i<<" model.groups[i].value1[j]= "<<model.groups[i].value1[j]<<endl;
              fprintf(stderr,"4.9 wrong!!! %g %d %d\n",model.groups[i].value1[j],i,j);
              return 0;
            }

          } 
        }
      }
// ******************************************** //
        //check bottom of overlying layer Vs (e.g., crust) < top of deeper layer (e.g., mantle) Vs
        //cout<<"now do between layers!"<<endl;
        // for(i=1;i<model.ngroup-1;i++)//between crust and mantle
        // for(i=0;i<model.ngroup-1;i++)//between all layers
        for(i=1;i<model.ngroup;i++)//between sed and crust layers only
		    {
          if (model.groups[i].q_flag==3 && model.groups[i-1].q_flag==2) // crust and sed
          {
            if( ((model.groups[i].value1.front() < (model.groups[i-1].value1.back() + eps)) || 
            (model.groups[i].value1p.front() < (model.groups[i-1].value1p.back() + eps)) || 
            (model.groups[i].value1r.front() < (model.groups[i-1].value1r.back() + eps))))
            // if( (model.groups[i+1].value1[0]<model.groups[i].value1.back() || model.groups[i+1].value1p[0]<model.groups[i].value1p.back() || model.groups[i+1].value1r[0]<model.groups[i].value1r.back()) && ((i==0 && model.groups[i].flag !=5 ) || i!=0))
            {
              fprintf(stderr,"bottom sediment greater than top of crust!!! %d  %g %d %g \n",i-1,model.groups[i-1].value1.back(),i,model.groups[i].value1.front());
              return 0;
              }
          }
        }//end of check positive jump from above into below layer


        for(i=1;i<model.ngroup;i++)//between crust and mantle
        // for(i=0;i<model.ngroup-1;i++)//between all layers
        // for(i=0;i<model.ngroup-2;i++)//between sed and crust layers only
		    {
          if (model.groups[i].q_flag==4 && model.groups[i-1].q_flag==3) // crust and mantle
          {
            if( (model.groups[i].value1.front() < model.groups[i-1].value1.back() || model.groups[i].value1p.front() < model.groups[i-1].value1p.back() || model.groups[i].value1r.front() < model.groups[i-1].value1r.back()))
          // if( (model.groups[i+1].value1[0]<model.groups[i].value1.back() || model.groups[i+1].value1p[0]<model.groups[i].value1p.back() || model.groups[i+1].value1r[0]<model.groups[i].value1r.back()) && ((i==0 && model.groups[i].flag !=5 ) || i!=0))
            {
            fprintf(stderr,"Bottom crust velocity greater mantle first layer!!! %d  %g %d %g \n",i-1,model.groups[i-1].value1.back(),i,model.groups[i].value1[0]);
            return 0;
            }

          }
          
        }//end of check positive jump from above into below layer

// ******************************************** //
    //check 1st crust spline < 2nd crust spline, EMB
      // HVL: Check the crustal index by q_flag
      if (invtype == 2) // bspline 
      {
        
        for(ij=0;ij<=model.ngroup;ij++)
        {

          if (model.groups[ij].q_flag==3) // crustal flag
          {
            if(model.groups[ij].value1[1]<model.groups[ij].value1[0])
            {
              cout<<"crust 1 = "<<model.groups[ij].value1[1]<<" < crust 0 = "<<model.groups[ij].value1[0]<<endl;
              fprintf(stderr,"2nd crust spline < 1st crust spline: %g < %g | group q_flag: %d\n",model.groups[ij].value1[1],model.groups[ij].value1[0],model.groups[ij].q_flag);
              return 0;
            } //end of 1st crust spline < 2nd crust spline
          }
        }
  // ******************************************** //
      //check 1st mantle spline < 2nd mantle spline
  // /*        
        // HVL: Check the crustal index by q_flag
        for(ij=0;ij<=model.ngroup;ij++)
        {
          if (model.groups[ij].q_flag==4) // mantle flag
          {
          if(model.groups[ij].value1[1]<model.groups[ij].value1[0])
            {
                cout<<"mantle 1 = "<<model.groups[ij].value1[1]<<" < mantle 0 = "<<model.groups[ij].value1[0]<<endl;
                fprintf(stderr,"2nd mantle spline < 1st mantle spline: %g < %g | group q_flag: %d\n",model.groups[ij].value1[1],model.groups[ij].value1[0], model.groups[ij].q_flag);
                return 0;
            } //end of 1st mantle spline < 2nd mantle spline  
          }
        }
      }

// ******************************************** //
//Check sediment layer is increasing Vs with depth (if not, discard model)
      for(ij=0;ij<=model.ngroup;ij++)
      {
        if (model.groups[ij].q_flag==2) // sediment flag
        {
          fprintf(stderr,"Sediment monotonic checking!\n");   
          fprintf(stderr,"N-lay sed = %d\n",(model.groups[ij].nlay));
          if (model.groups[ij].nlay < 2)
          {
            fprintf(stderr,"Error, sediment have only 1 layer!\n");
            return 0;
          }
          else if (model.groups[ij].nlay == 2)
          {
            if (model.groups[ij].value1[1] < (model.groups[ij].value1[0] + eps))
            {
              fprintf(stderr,"Sediment not monotonic increasing! %d %d %lf > %lf \n",1,ij,model.groups[ij].value1[0],model.groups[ij].value1[1]);
              return 0;
            }
            
          }
          else
          {
            for(int ised=1;ised<model.groups[ij].nlay;ised++)
            {
              const double prev = model.groups[ij].value1[ised - 1];
              const double curr = model.groups[ij].value1[ised];

              //cout<<"Model values in layer "<<j<<" < values in layer "<<i<<"(specifically) "<<(model.groups[j].value1[i+1])<<" <  "<<(model.groups[j].value1[i])<<endl;
              if (curr + eps < prev) 
              {
                fprintf(stderr,"Sediment not monotonic increasing! %d %d %lf > %lf \n",ised,ij,model.groups[ij].value1[ised-1],model.groups[ij].value1[ised]);
              // fprintf(stderr,"monotonic wrong!!! %d %d %lf > %lf \n",i,j,model.groups[j].value1[i],model.groups[j].value1[i+1]); 
              // cout<<"BAD Model values in layer "<<j<<" < values in layer "<<i<<"(specifically) "<<(model.groups[j].value1[i+1])<<" <  "<<(model.groups[j].value1[i])<<endl;
                    return 0;
              }
            }
          }
        }
      }
    //end of check sed layer is increasing Vs with depth
// ******************************************** //
// checking for postitive gradient at top of layers
          //cout<<"now do vgrad!"<<endl;
          //cout<<*(vgrad.end()-1)<<endl;
	  for(id=vgrad.begin();id<vgrad.end();id++)//gradient check for the 1st two velue in group vgrad[?]

	  {
            if (*id > 0)
	    j=*id;
            else
            j=-*id;

            //cout<<j<<" "<<*id<<endl;
            //cout<<model.groups[j].value1[1]<<" "<<model.groups[j].value1[0]<<endl;
	    if( (model.groups[j].value1[1]<model.groups[j].value1[0] && *id>0) || (*id<0 && model.groups[j].value1[1]>model.groups[j].value1[0])) 
		{
                  fprintf(stderr,"vgrad wrong! %d\n",j);
                  return 0;}
	  } 
	  //if(model.groups[3].value1[1]>model.groups[3].value1[0]) {  // negative gradient  forced!!
          //    return 0;
	  //    }
       //end of check positive gradient at top of layers
// ******************************************** //

    //passed all tests, so can now return 1
      return 1;
	  }

//-----------------------------------------------------	
	
//};

        int write_model1(modeldef &model,char *name, char *outdir)
        {

          // fprintf(stderr,"> CALmodel: write_model1\n");
          // sleep(1);

        char namerf[300],namedispp[300],namedispg[300],namemod[300];
          FILE *ff;
          int i;
          double dep=0.;
          sprintf(namerf,"%s/%s.rf",outdir,name);
          sprintf(namedispp,"%s/%s.p.disp",outdir,name);
          sprintf(namedispg,"%s/%s.g.disp",outdir,name);
          sprintf(namemod,"%s/%s.mod",outdir,name);
          ff=fopen(namerf,"w");
          for (i=0;i<model.data.rf.rfn.size();i++)
                {
                  if(model.data.rf.tn.size()>0 and model.data.rf.tn[i]<=10)
                        fprintf(ff,"%g %g %g %g %g\n",model.data.rf.tn[i],model.data.rf.rfn[i],model.data.rf.to[i],model.data.rf.rfo[i],model.data.rf.unrfo[i]);
                }//for i
         fclose(ff);

         if(model.data.disp.npper>0)
         {
         ff=fopen(namedispp,"w");
         for(i=0;i<model.data.disp.npper;i++)
                {
                  if(model.data.disp.pvel.size()>0)
                        fprintf(ff,"%g %g %g %g\n",model.data.disp.pper[i],model.data.disp.pvel[i],model.data.disp.pvelo[i],model.data.disp.unpvelo[i]);
                }//fori  
          fclose(ff);
          }//if


         if(model.data.disp.ngper>0)
         {
         ff=fopen(namedispg,"w");
         for(i=0;i<model.data.disp.ngper;i++)
                {
                  if(model.data.disp.gvel.size()>0)
                        fprintf(ff,"%g %g %g %g\n",model.data.disp.gper[i],model.data.disp.gvel[i],model.data.disp.gvelo[i],model.data.disp.ungvelo[i]);
                }//fori  
          fclose(ff);
          }//if
          /*
          ff=fopen(namemod,"w");
          for(i=0;i<model.laym0.nlayer;i++) //write layered model in stair-step
                {
                  fprintf(ff,"%g %g\n",dep,model.laym0.vs[i]);
                  dep=dep+model.laym0.thick[i];
                  fprintf(ff,"%g %g\n",dep,model.laym0.vs[i]);
                }
          fclose(ff);
          return 1;
          */
          ff=fopen(namemod,"w");
          dep = 0.;
          int k,j,nlayer;
          fprintf(stderr,"how many groups??? %d groups\n",model.ngroup);
          for (i=0;i<model.ngroup;i++) {
              if (model.groups[i].flag == 5) // water layer
                  {
                  fprintf(ff,"%g 0.\n",dep);
                  dep = dep + model.groups[i].thick;
                  fprintf(ff,"%g 0.\n",dep);
                  }
              if (model.groups[i].flag == 4) // grad
                  {
                  fprintf(ff,"%g %g %d\n",dep,model.groups[i].value[0],model.groups[i].q_flag);
                  dep = dep + model.groups[i].thick;
                  fprintf(ff,"%g %g %d\n",dep,model.groups[i].value[1],model.groups[i].q_flag);
                  }
              if (model.groups[i].flag == 1) // layered
                  {
                  k = model.groups[i].np;
                  nlayer = k;
                  fprintf(stderr,"how many layers??? %d layers %d\n",nlayer,k);
                  for (j=0;j<nlayer;j++) {
                     fprintf(ff,"%g %g %d\n",dep,model.groups[i].value[j],model.groups[i].q_flag);
                     dep = dep + model.groups[i].thick*model.groups[i].ratio[j];
                     fprintf(ff,"%g %g %d\n",dep,model.groups[i].value[j],model.groups[i].q_flag);
                     }
                  }
              if (model.groups[i].flag == 2) // B-splines
                  {
                  fprintf(ff,"%g %g %d\n",dep,model.groups[i].value1[0],model.groups[i].q_flag);
                  nlayer = model.groups[i].nlay;
                  //dep = dep + model.groups.thick1[0]/2;
                  for (j=0;j<nlayer-1;j++) {
                     dep = dep + model.groups[i].thick1[j];
                     fprintf(ff,"%g %g %d\n",dep,(model.groups[i].value1[j]+model.groups[i].value1[j+1])/2,model.groups[i].q_flag);
                     }
                  dep = dep + model.groups[i].thick1[nlayer-1];
                  fprintf(ff,"%g %g %d\n",dep,model.groups[i].value1[nlayer-1],model.groups[i].q_flag);
                  }
              }
           return 1;
        }
