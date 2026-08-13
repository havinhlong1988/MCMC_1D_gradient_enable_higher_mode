#include<iostream>
#include<algorithm>
#include<vector>
#include<cmath>
#include<fstream>
#include<string>
#include<omp.h>
//#include"/home/jiayi/progs/jy/HEAD/head_c++/string_split.C"
//#include"/home/jiayi/progs/jy/HEAD/head_c++/generate_Bs.C"
//#include"/home/jiayi/progs/jy/HEAD/head_c++/gen_random.C"
#include"./head_c++/string_split.C"
#include"./head_c++/generate_Bs.C"
#include"./head_c++/gen_random.C"
#include "read_in.C"
#include"INITstructure.h"
#include"CALgroup.C"
#include"CALmodel.C"
#include"CALforward.C"
#include"CALpara.C"
#include"BIN_rw.C"
#include"MC.C"
#include<map>
#include<boost/array.hpp>
using namespace std;
using namespace boost;
void print( vector <string> & v )
{
  for (size_t n = 0; n < v.size(); n++)
    cout << "\"" << v[ n ] << "\"\n";
  cout << endl;
}

int main( int argc, char *argv[])
{
  int surtype,surn,ninput;
  //cout<<"ok here"<<endl;
  if(argc<2)
  {
    cout<<"input [control.file]"<<endl;
    exit(0);
  }//if argc
  //cout<<"ok here"<<endl;

  //cout<<"now begin to define parameters !! "<<argv[1]<<endl;
  modeldef model,tmodel,ttmodel;
  paradef para0,para,para1;
  indatadef indata0;
  string line;
  vector<string> vargv,vdisp;
  vector<int> vmono,vgrad,Ql,surflag;
  vector<paradef> paras;
  vector<modeldef> models;
  vector<boost::array<double, 10> > infos;
  boost::array<double, 10> tinfo;
  char outdir[300],wmodelname1[300],wmodelname2[300],outfnmbin[300],outfnmbin1[300],outfnm[300],strnn[300],t_str[300];
  char scheck[300],inmodelnm[300],inrfnm[300],inparanm[300],Qname[300],outname[300],indatanm[300];
  double pp,oldL,oldmisfit,newL,newmisfit,prob,prandom;
  int monoc,iiter=0,iaccp=0,ibad=0,i,j,k,key,fnn,jumpn,flagw,flag_newp,flag_misfit,n_core;
  time_t start,end;
  double tdif,TTT,gaus;
  double mc_range, dstep;
  int plottype, invtype, selstyle=-1;
  FILE *outf;
  FILE *ff0; //EMB
  int Qflag = 0, maxj;
  double tmp; //EMB
  double mohothick=0; //Added by EMB
  ofstream outbin;
  ofstream outbin1;
  // --------------- Said something for everyone -----------------------
  fprintf(stderr,"running do_MC_Para!\n");
  sleep(2);
  // -------------------------------------------------------------------
  //--model intialization -------------
  initmodel(model);
  initmodel(tmodel);
  initmodel(ttmodel);

  initpara(para0);
  initpara(para);
  initpara(para1);

  //--read in data --------------------

  if (read_in (argv[1],inmodelnm,inparanm, &surtype, &surn, surflag, vdisp, \
  inrfnm, &gaus, Qname, &Qflag, Ql, &pp, &monoc, vmono,\
  vgrad, &TTT, &jumpn, outdir, outname, &flagw, &n_core,\
  indatanm, &mc_range, &dstep, &plottype, &invtype, &selstyle)==0)
  {
    fprintf(stderr,"read_in failed! \n");
    return 0;
    }


  maxj = (int)((TTT-1)/n_core)+1; // TTT: total jump times;
  cout<<"maxj = "<<maxj<<" TTT "<<TTT<<endl;

  readdisp1 (model,surtype, surn, surflag, vdisp);
  readrf(model,inrfnm,gaus);

  readmod(model,inmodelnm);
  cout<<"model.tthick?? "<<model.tthick<<endl;
  cout<<"model.ngroup?? "<<model.ngroup<<endl;
  // fprintf(stderr,"Read model finish, check the moho_thick! \n");
  // get moho thickness,added by EMB
  for(i=0;i<model.ngroup;i++){
      cout<<"mohothick if i: "<<i<<endl;
      
      if (model.groups[i].q_flag != 4)//not in mantle
      { mohothick=mohothick+model.groups[i].thick;
      }
      else {continue;}
    }
  model.mohodepth=mohothick;
  cout<<"mohothick?? "<<model.mohodepth<<endl;
  // fprintf(stderr,"Now read input parameters! \n");
  readpara(para0,inparanm);
  readindata(model,indatanm);

  if ( Qflag > 0) {
    readQmodel(model,Qname,Ql);
    }

  cout<<"read all information ok!"<<endl;
  //--------------------------------------------
  //-- layerize the model and compute rf, disp--
  //-- compute disp for R-P, R-G, R-E ----------
  fprintf(stderr,"layerize the model...\n");
  updatemodel(model,invtype);
  cout<<"update ok!"<<endl;
  fprintf(stderr,"update model ok \n");
//note, still on input model
  fprintf(stderr,"compute start RF\n");
  compute_rf(model,invtype);
  cout<<"compute rf ok!"<<endl;
  fprintf(stderr,"compute RF ok\n");
  
  fprintf(stderr,"compute start DISP\n");
  compute_disp(model,invtype);
  cout<<"compute rf and DISP ok!"<<endl;
  fprintf(stderr,"compute DISP ok\n");
  
  if (pp >=0 && pp<=1) {
    fprintf(stderr,"now compute initial misfit !!\n");
    compute_misfit(model,pp);
    }
  else {
    compute_misfit(model,pp);
    cout<<"do prior sampling, no misfit computed needed, but we still compute initial"<<endl;
    }

  fprintf(stderr,"misfit done...\n");
  fprintf(stderr,"initial misfit values: %f\n",model.data.misfit);
  sleep(1);
  //EMB add, output rf, phase velocity and H/V for the initial model
  sprintf(t_str,"mkdir %s\0",outdir);
  system(t_str);
  // HVL try to print the intial model as well
  sprintf(wmodelname1,"%s/Initial.mod",outdir);
  ff0=fopen(wmodelname1,"w");
  double tdep=0.;
  for(i=0;i<model.laym0.nlayer;i++)
    {
      tdep=tdep+model.laym0.thick[i];
      fprintf(ff0," %g %g %g %g %g %g %g %g\n",tdep,model.laym0.vs[i],model.laym0.vp[i],model.laym0.vpvs[i],model.laym0.rho[i],model.laym0.qs[i],model.laym0.qp[i],model.laym0.q_flag[i]);
      // fprintf(ff0," %g %g \n",tdep,model.laym0.vs[i]);
    }
  fclose(ff0);

  sprintf(wmodelname1,"%s/Initial.hv",outdir);
  ff0=fopen(wmodelname1,"w");
  if (model.data.disp.eper.size() >0) {
      for (i=0;i<model.data.disp.eper.size();i++){
        fprintf(ff0," %g %g \n",model.data.disp.eper[i],model.data.disp.evel[i]);//only need period and ellipticity
      }
      // for (i=0;i<model.data.disp.eper.size();i++) {
        // fprintf(ff0," %g",model.data.disp.evel[i]);//only need period and ellipticity
        //fprintf(ff0,"%g %g %g %g\n",tmodel.data.disp.eper[j],tmodel.data.disp.evel[j],tmodel.data.disp.evelo[j],tmodel.data.disp.unevelo[j]);
        // }
       
        fprintf(ff0,"\n");
    }  
  fclose(ff0);

  sprintf(wmodelname1,"%s/Initial.rf",outdir);
  ff0=fopen(wmodelname1,"w");
  if (model.data.rf.nrfo > 0) {
      for (i=0;i<model.data.rf.nrfo;i++) {
        fprintf(ff0,"%g %g\n",model.data.rf.to[i],model.data.rf.rfn[i]);
        }
      }
  fclose(ff0);

  sprintf(wmodelname1,"%s/Initial.ph",outdir);
  
  ff0=fopen(wmodelname1,"w");
  if (model.data.disp.pper.size() >0) {
      for (i=0;i<model.data.disp.pper.size();i++){
        fprintf(ff0," %g %g \n",model.data.disp.pper[i],model.data.disp.pvel[i]);//only need period and ellipticity
      }
      // for (i=0;i<model.data.disp.pper.size();i++) {
      //   fprintf(ff0," %g",model.data.disp.pvel[i]);//only need period and ellipticity
      //   //fprintf(ff0,"%g %g %g %g\n",tmodel.data.disp.eper[j],tmodel.data.disp.evel[j],tmodel.data.disp.evelo[j],tmodel.data.disp.unevelo[j]);
      //   }
      
        // fprintf(ff0,"\n");
    }  
  fclose(ff0);
//End of EMB adds to output initial model's forward model result
    fprintf(stderr,"Calculate right and left bound\n");
    // HVL 20251112 - calculate the model space from initial model and parameters
    paradef paramid,paraleft,pararight; // min and max value for parameter space - calculate from initial mod.STA and in.para_STA
    modeldef modelleft,modelright;

    // initpara(paramid);
    initpara(paraleft);
    initpara(pararight);

    initmodel(modelleft);
    initmodel(modelright);

    mod2para(model,para0,para,invtype); // now parameters sotred in "paraleft.parameters" and their coresponding space stored in "paraleft.space1"
    char wparaname[300];
    sprintf(wparaname,"Initial.paras");
    write_paras0({para},wparaname,outdir,invtype); // write the initial data
    getminmaxpara(para,paraleft,pararight);
    // fprintf(stderr,"get min max paras done.. \n");
    sprintf(wmodelname1,"%s/Initial.paras_left",outdir);
    ff0=fopen(wmodelname1,"w");
    for(int i=0;i<paraleft.parameter.size();i++)
      {
        fprintf(ff0,"%g ",paraleft.parameter[i]);
      }
    fclose(ff0);
    sprintf(wmodelname1,"%s/Initial.paras_right",outdir);
    ff0=fopen(wmodelname1,"w");
    for(int i=0;i<pararight.parameter.size();i++)
      {
        fprintf(ff0,"%g ",pararight.parameter[i]);
      }
    fclose(ff0);
    fprintf(stderr,"Write out the boundaries\n");

    // ============================ left model ============================
    para2mod(paraleft,model,modelleft,invtype);
    updatemodel(modelleft,invtype);
    sprintf(wmodelname1,"%s/Initial_left.mod",outdir);
    ff0=fopen(wmodelname1,"w");
    tdep=0.;
    for(int i=0;i<modelleft.laym0.nlayer;i++)
      {
        tdep=tdep+modelleft.laym0.thick[i];
        fprintf(ff0,"%g %g %g %g %g %g %g %g\n",tdep,modelleft.laym0.vs[i],modelleft.laym0.vp[i],modelleft.laym0.vpvs[i],modelleft.laym0.rho[i],modelleft.laym0.qs[i],modelleft.laym0.qp[i],modelleft.laym0.q_flag[i]);
      }
    fclose(ff0);
    // ============================ right model ============================
    para2mod(pararight,model,modelright,invtype);
    updatemodel(modelright,invtype);
    sprintf(wmodelname1,"%s/Initial_right.mod",outdir);
    ff0=fopen(wmodelname1,"w");
    tdep=0.;
    for(int i=0;i<modelright.laym0.nlayer;i++)
      {
        tdep=tdep+modelright.laym0.thick[i];
        fprintf(ff0,"%g %g %g %g %g %g %g %g\n",tdep,modelright.laym0.vs[i],modelright.laym0.vp[i],modelright.laym0.vpvs[i],modelright.laym0.rho[i],modelright.laym0.qs[i],modelright.laym0.qp[i],modelright.laym0.q_flag[i]);
      }
    fclose(ff0);
    // fprintf(stderr," Iam going to sleep\n");
    // sleep(10);
    fprintf(stderr,"Write initial bin\n");
  //-- write model here -----------------------
  sprintf(wmodelname1,"MC.%s",outname);
  sprintf(outfnmbin,"%s/BIN.%s",outdir,outname);
  sprintf(outfnmbin1,"%s/BIN.%s.1",outdir,outname);
  sprintf(t_str,"mkdir %s\0",outdir);
  system(t_str);

  outbin.open(outfnmbin,ios::out | ios::binary); 
  write_model(model,wmodelname1,outdir,flagw);
  fprintf(stderr,"writing done, mod2para\n");
  // This overlape the above additional (20251112)
  mod2para(model,para0,para,invtype);//model para0 => para. para.parameter[]/space1[]/flag is changed based on model.groups[].value[]/thick/vpvs
  fprintf(stderr,"mod2para done.. \n");
  tinfo[0] = -2;
  tinfo[1] = 0;
  tinfo[2] = 0;
  tinfo[3] = 0;
  tinfo[4] = model.data.L;
  tinfo[5] = model.data.misfit;
  tinfo[6] = model.data.rf.misfit;
  tinfo[7] = model.data.disp.misfit;
  tinfo[8] = double(para.npara);
  tinfo[9] = double(model.ngroup);
  infos.push_back(tinfo);
  paras.push_back(para);
  write_bin(model,outbin,para,-2,0,0,0);
  outbin.close();
  cout<<"write bin over: "<<outfnmbin<<endl;
  cout<<"everything is fine after reading and writing\n";
  fprintf(stderr,"everything is fine after reading and writing\n");
  
  //-----l and misfit based on the read in model and readin disp/rf-------
  oldL=model.data.L;
  oldmisfit=model.data.misfit;
  cout<<"original l and misfit:"<<oldL<<" "<<oldmisfit<<endl;
  fprintf(stderr,"original l and misfit: %f - %f\n",oldL,oldmisfit);
  const char *invtypechar;
  if (invtype == 1)
  {
      invtypechar = "layercake";
  }
  else if (invtype == 2)
  {
      invtypechar = "b-spline";
  }
  else
  {
      invtypechar = "unknown";
  }
  fprintf(stderr, "check inversion style: %d (%s)\n", invtype, invtypechar);
  // int gm;
  // gm=goodmodel(model,vmono,vgrad,invtype);
  // fprintf(stderr,"the good model is %d \n",gm);
  // abort();
  //----------------------------------
  //mod2para(model,para0,para);//model para0 => para. para.parameter[]/space1[]/flag is changed based on model.groups[].value[]/thick/vpvs

  // now prepare for the loop --------------------------------------------
  sprintf(outfnm,"%s/%s.out",outdir,wmodelname1);
  fprintf(stderr,"open file %s\n",outfnm);
  cout<<outfnm<<endl;
  outf=fopen(outfnm,"w");
  time (&start); 
  cout<<"do_MC_Para pp "<<pp<<endl;
  fprintf(stderr,"do_MC_Para pp: %.2f\n",pp);
  sleep(2);
  if (pp>1.) return 1;

  if (pp<0.) {
   if (do_prior(outf, para, model, monoc, vmono, vgrad, pp, jumpn, iaccp, iiter, wmodelname1, outdir, flagw, paras, invtype) ==0) {
     fprintf(stderr,"cannot do the prior!\n");
     return 0;
     }
   else
     fprintf(stderr,"finish!!! prior\n");
     return 1;
   }  
  cout<<n_core<<endl;
  //abort();
  j=0;
  while (j<maxj) {
  #pragma omp parallel for default(none) shared (paras,infos,i,outbin,outf,para,model,monoc,vmono,vgrad,pp,jumpn,iaccp,iiter,wmodelname1,outdir,flagw,n_core , invtype)
  for (i=0;i<n_core;i++) {
   do_one_jump(outf, para, model, monoc, vmono, vgrad, pp, jumpn, iaccp, iiter, wmodelname1, outdir, flagw, i, paras, infos , invtype);
   }
  cout<<j<<" th done \n"<<wmodelname1<<endl;
  fprintf(stderr,"do the jump %d th done!\n",j);
  j++;
  }
  fclose(outf);
  //outbin.close();
 
  outbin1.open(outfnmbin1,ios::out | ios::binary);
  for (i=0;i<paras.size();i++) {
    write_bin1(outbin1,paras[i],infos[i]);
    }
  outbin1.close();
  // HVL Check the forward paras.size(). If less than 2 then redo the process.
  if (paras.size()<2)
  {
    fprintf(stderr,"Parameter size too small. abort, abort, abort,!!! \n");
    fprintf(stderr,"paras.size(): %d\n",paras.size());
    return 0;
  }
  else
  {
    cout<<"do_MC_Para done!!!"<<endl;
    fprintf(stderr,"do_MC_ Para is done for this station! %s\n",wmodelname1);
    return 1;
  }

}//main

