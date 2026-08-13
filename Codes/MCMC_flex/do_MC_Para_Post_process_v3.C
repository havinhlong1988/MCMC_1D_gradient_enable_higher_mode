#include<iostream>
#include<algorithm>
#include<vector>
#include<cmath>
#include<fstream>
#include<string>
#include<omp.h>
#include<iterator>
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
#include <boost/accumulators/statistics/stats.hpp>
using namespace std;
using namespace boost;

void print( vector <string> & v )
{
  for (size_t n = 0; n < v.size(); n++)
    cout << "\"" << v[ n ] << "\"\n";
  cout << endl;
}

  // ======================================================================
  // Helper function
  // ======================================================================

// ============= layer-cake like value for each depth =============================
double get_vs_at_depth_common_grid(const modeldef& mod, double z)
{
    if (mod.laym0.nlayer <= 0) return -1.0;
    if ((int)mod.laym0.vs.size() != mod.laym0.nlayer) return -1.0;
    if ((int)mod.laym0.thick.size() != mod.laym0.nlayer) return -1.0;

    double ztop = 0.0;
    for (int i = 0; i < mod.laym0.nlayer; ++i)
    {
        double zbot = ztop + mod.laym0.thick[i];

        if (z >= ztop && z <= zbot)
            return mod.laym0.vs[i];

        ztop = zbot;
    }

    return -1.0;
}
// ============= Gradient like value for each depth =============================
double get_vs_at_depth_common_grid_gradient(const modeldef& mod, double z)
{
    if (mod.laym0.nlayer <= 0) return -1.0;
    if (mod.laym0.vs.size() != mod.laym0.thick.size()) return -1.0;

    double ztop = 0.0;

    for (int i = 0; i < mod.laym0.nlayer; ++i)
    {
        double zbot = ztop + mod.laym0.thick[i];

        if (z >= ztop && z <= zbot)
        {
            if (mod.laym0.thick[i] <= 0.0)
                return mod.laym0.vs[i];

            // choose endpoint values
            double v_top, v_bot;

            if (i == 0) {
                v_top = mod.laym0.vs[i];
                v_bot = mod.laym0.vs[i];
            } else {
                v_top = mod.laym0.vs[i - 1];
                v_bot = mod.laym0.vs[i];
            }

            double dtop = z - ztop;
            double dbot = zbot - z;
            double h = zbot - ztop;

            return (dbot / h) * v_top + (dtop / h) * v_bot;
        }

        ztop = zbot;
    }

    return mod.laym0.vs.back();
}
// ============= Stepwise sampling for averaging =============================
double get_vs_at_depth_common_grid_stepwise(const modeldef& mod, double z)
{
    if (mod.laym0.nlayer <= 0) return -1.0;
    if (mod.laym0.vs.size() != mod.laym0.thick.size()) return -1.0;

    double ztop = 0.0;

    for (int i = 0; i < mod.laym0.nlayer; ++i)
    {
        double h = mod.laym0.thick[i];
        if (!(h > 0.0) || !std::isfinite(h)) continue;

        double zbot = ztop + h;

        // left-closed, right-open; last layer includes bottom
        if ((z >= ztop && z < zbot) ||
            (i == mod.laym0.nlayer - 1 && z >= ztop && z <= zbot))
        {
            double v = mod.laym0.vs[i];
            if (!std::isfinite(v)) return -1.0;
            return v;
        }

        ztop = zbot;
    }

    // allow tiny rounding overshoot at bottom
    if (mod.laym0.nlayer > 0 && z >= 0.0)
    {
        double vlast = mod.laym0.vs.back();
        if (std::isfinite(vlast)) return vlast;
    }

    return -1.0;
}

void write_initial_model_bounds(
    modeldef &model,
    paradef &para0,
    paradef &para,
    char *outdir,
    int invtype
)
// Recalculate model-space bounds from the initial model and parameter file.
// This fixes incorrect bounds from do_MC_para without rerunning the full workflow.
{
    fprintf(stderr,"Calculate right and left bound\n");

    paradef paraleft, pararight;
    modeldef modelleft, modelright;

    initpara(paraleft);
    initpara(pararight);
    initmodel(modelleft);
    initmodel(modelright);

    char wname[300];
    FILE *ff0;
    double tdep;

    mod2para(model, para0, para, invtype);

    // sprintf(wname, "Initial.paras");
    // write_paras0(para, wname, outdir, invtype);
    sprintf(wname, "Initial.paras");

    vector<paradef> vpara0;
    vpara0.push_back(para);

    write_paras0(vpara0, wname, outdir, invtype);

    getminmaxpara(para, paraleft, pararight);

    sprintf(wname, "%s/Initial.paras_left", outdir);
    ff0 = fopen(wname, "w");
    for(int i=0; i<paraleft.parameter.size(); i++)
        fprintf(ff0, "%g ", paraleft.parameter[i]);
    fprintf(ff0, "\n");
    fclose(ff0);

    sprintf(wname, "%s/Initial.paras_right", outdir);
    ff0 = fopen(wname, "w");
    for(int i=0; i<pararight.parameter.size(); i++)
        fprintf(ff0, "%g ", pararight.parameter[i]);
    fprintf(ff0, "\n");
    fclose(ff0);

    para2mod(paraleft, model, modelleft, invtype);
    updatemodel(modelleft, invtype);

    sprintf(wname, "%s/Initial_left.mod", outdir);
    ff0 = fopen(wname, "w");

    tdep = 0.0;
    for(int i=0; i<modelleft.laym0.nlayer; i++)
    {
        tdep += modelleft.laym0.thick[i];

        fprintf(ff0,"%g %g %g %g %g %g %g %g\n",
            tdep,
            modelleft.laym0.vs[i],
            modelleft.laym0.vp[i],
            modelleft.laym0.vpvs[i],
            modelleft.laym0.rho[i],
            modelleft.laym0.qs[i],
            modelleft.laym0.qp[i],
            modelleft.laym0.q_flag[i]
        );
    }
    fclose(ff0);

    para2mod(pararight, model, modelright, invtype);
    updatemodel(modelright, invtype);

    sprintf(wname, "%s/Initial_right.mod", outdir);
    ff0 = fopen(wname, "w");

    tdep = 0.0;
    for(int i=0; i<modelright.laym0.nlayer; i++)
    {
        tdep += modelright.laym0.thick[i];

        fprintf(ff0,"%g %g %g %g %g %g %g %g\n",
            tdep,
            modelright.laym0.vs[i],
            modelright.laym0.vp[i],
            modelright.laym0.vpvs[i],
            modelright.laym0.rho[i],
            modelright.laym0.qs[i],
            modelright.laym0.qp[i],
            modelright.laym0.q_flag[i]
        );
    }
    fclose(ff0);

    fprintf(stderr,"Write initial model parameter bounds done.\n");
}

// ==============================================================================
// ======================= Main fucntion start here =============================
// ==============================================================================
int main( int argc, char *argv[])
{
  // --------------- Said something for everyone -----------------------
  fprintf(stderr,"running do_MC_Para_Post_process_v3!\n");
  fprintf(stderr,"Do the post process analysis!\n");
  sleep(2);
// -------------------------------------------------------------------
  int surtype,surn,ninput;
  int option=0;
  //cout<<"ok here"<<endl;
  if(argc<2)
  {
    cout<<"input [control.file]"<<endl;
    exit(0);
  }//if argc
  //cout<<"ok here"<<endl;

  //cout<<"now begin to define parameters !! "<<argv[1]<<endl;
  modeldef model,tmodel,ttmodel;
  modeldef minmodel,maxmodel;
  paradef para0,para,para1,tpara,ave_para;
  paradef minparaemb,maxparaemb;
  paradef minpara,maxpara,inipara;
  indatadef indata0;
  vector<paradef> paras,paras1;
  vector<modeldef> models,models1,models2,models3;
  vector<int> flags, signs1; 
  string line;
  vector<string> vargv,vdisp;
  vector<int> vmono,vgrad,Ql,surflag,accps;
  vector<double> kai_disp,kai_rf,kai_joint;
  vector<double> seddepth,sedperc;//EMB
  vector<boost::array<double, 10> > infos;
  boost::array<double, 10> tinfo;  
  int sign;

  double kai_cri,kai_min;
  char outdir[300],wmodelname1[300],wmodelname2[300],outfnmbin[300],outfnm[300],strnn[300],t_str[300];
  char scheck[300],inmodelnm[300],inrfnm[300],inparanm[300],Qname[300],outname[300],indatanm[300];
  double pp,oldL,oldmisfit,newL,newmisfit,prob,prandom;
  double minmisfit,minmisfit_rf,minmisfit_disp;
  double minmisfit_o,minmisfit_rf_o,minmisfit_disp_o;
  int monoc,iiter=0,iaccp=0,ibad=0,i,j,k,key,fnn,jumpn,flagw,flag_newp,flag_misfit,n_core,nall;
  time_t start,end;
  double tdif,TTT,gaus;
  double mc_range, dstep;
  int plottype,invtype;
  FILE *outf;
  int Qflag = 0, maxj;
  //ofstream outbin;
  ifstream inbin;

  //Added by EMB
  FILE *itmf0emb;
  char itmffEMB[300];
  double mohothick=0; //Added by EMB
  double moddepth=0;
  int ihv,iph,igv,irf;
  // End of adds

  //--model intialization -------------
  initmodel(model);
  initmodel(tmodel);
  initmodel(ttmodel);
  initmodel(minmodel); // add new 20251103
  initmodel(maxmodel); // add new 20251103

  initpara(para0);
  initpara(para);
  initpara(para1);
  initpara(tpara);
  initpara(ave_para);
  initpara(inipara); // add new 20251103
  initpara(minpara); // add new 20251103
  initpara(maxpara); // add new 20251103
  
  //--read in data --------------------
  fprintf(stderr,"Reading the input for post process ...\n");

  if (read_in (argv[1],inmodelnm,inparanm, &surtype, &surn, surflag, vdisp, \
               inrfnm, &gaus, Qname, &Qflag, Ql, &pp, &monoc, vmono,\
               vgrad, &TTT, &jumpn, outdir, outname, &flagw, &n_core,indatanm,\
               &mc_range, &dstep, &plottype, &invtype)==0) 
    {
      fprintf(stderr,"read_in failed! \n");
      return 0;
    }
  fprintf(stderr,"Report new test paras: mcrange - dstep - plottype %f %f %d\n",mc_range, dstep, plottype);
  sleep(5);

  maxj = (int)((TTT-1)/n_core)+1; // TTT: total jump times;

  readdisp1 (model,surtype, surn, surflag, vdisp);
  readrf(model,inrfnm,gaus);
  cout<<"model rf: "<<model.data.rf.nrfo<<endl;
  readmod(model,inmodelnm);
    cout<<"model.tthick?? "<<model.tthick<<endl;
      cout<<"model.ngroup?? "<<model.ngroup<<endl;

    // get moho thickness,added by EMB
  for(i=0;i<model.ngroup;i++){
      cout<<"mohothick if i: "<<i<<endl;
      // get model thickness, added by HVL
      moddepth=moddepth+model.groups[i].thick;
      if (model.groups[i].q_flag != 4)//not in mantle
      { mohothick=mohothick+model.groups[i].thick;
      }
      else {continue;}
    }
  model.mohodepth=mohothick;
  cout<<"mohothick?? "<<model.mohodepth<<endl;


  readpara(para0,inparanm);
  cout<<model.groups[1].fvpvs<<" "<<model.groups[1].vpvs1<<endl;

  // ================== Check the parameter (space1 is not update) and model here =======================
// fprintf(stderr,"Check para!!\n");
// sleep(1);
//   for(i=0;i<para0.npara;i++)
//   {
//     // for (size_t j = 0; j < para0.para0[i].size(); ++j) 
//     // {
//     //   fprintf(stderr, "Check para[%zu][%zu]: %s\n",i, j, para0.para0[i][j].c_str());
//     // }
// }

// fprintf(stderr, "Check model!!\n");
// sleep(1);

// for (size_t i = 0; i < model.groups.size(); ++i) {
//     for (size_t j = 0; j < model.groups[i].value.size(); ++j) {
//         // value[j] is double → use %g (or %f / %.6f as you like)
//         fprintf(stderr, "Check model[%zu][%zu]: %g\n",
//                 i, j, model.groups[i].value[j]);
//     }
// }


  // ==========================================================================
  readindata(model,indatanm);
  iph=model.indata.phflag;
  igv=model.indata.gvflag;
  ihv=model.indata.hvflag;
  irf=model.indata.rfflag;
  // Define the option to switch case
  // if (iph+igv+ihv+irf=4) option=0; // 4 data types
  // if (iph+igv+ihv+irf=1) // only 1 data type
  // {
  //   if (iph==1) option=1;
  //   if (igv==1) option=2;
  //   if (ihv==1) option=3;
  //   if (irf==1) option=4;
  // }
  // if (iph+igv+ihv+irf=2) // 2 data types
  // {
  //   if ((iph==1)&&(igv==1)) option=5;
  //   if ((iph==1)&&(ihv==1)) option=6;
  //   if ((iph==1)&&(irf==1)) option=7;
  //   //
  //   if ((igv==1)&&(ihv==1)) option=8;
  //   if ((igv==1)&&(irf==1)) option=9;
  //   //
  //   if ((ihv==1)&&(irf==1)) option=10;
  // }
  // if (iph+igv+ihv+irf=3) // 2 data types
  // {
  //   if ((iph==1)&&(igv==1)&&(ihv==1)) option=11;
  //   if ((iph==1)&&(igv==1)&&(irf==1)) option=12;
  //   if ((iph==1)&&(irf==1)&&(ihv==1)) option=13;
  //   if ((igv==1)&&(irf==1)&&(ihv==1)) option=14;
  // }
  if ((irf==1)&&(igv+ihv+iph>0)) option=1; // HV and disp
  if ((irf==1)&&(igv+ihv+iph==0)) option=2; // RF only
  if ((irf==0)&&(igv+ihv+iph>0)) option=3; // DISP only



  if ( Qflag > 0) {
    readQmodel(model,Qname,Ql);
    }

  cout<<"read all information ok!"<<endl;
  //--------------------------------------------
  //-- layerize the model and compute rf, disp--
  //-- compute disp for R-P, R-G, R-E ----------

  updatemodel(model,invtype); // initial data transfered to the model values (take a look for laym0 vector)

  // fprintf(stderr, "Check model!!\n");
  // sleep(1);
  // double tdep=0.;
  // fprintf(stderr," %g %g \n",tdep,model.laym0.vs[0]);
  // for (int i=0;i<model.laym0.nlayer;++i)
  // {
  //   tdep=tdep+model.laym0.thick[i];
  //   // fprintf(stderr," %g %g %g %g %g %g %g\n",tdep,model.laym0.vs[i],model.laym0.vp[i],model.laym0.vpvs[i],model.laym0.rho[i],model.laym0.qs[i],model.laym0.qp[i]);
  //   fprintf(stderr," %g %g \n",tdep,model.laym0.vs[i]);
  // }

  mod2para(model,para0,inipara, invtype); // get the initial data and parameter space for intital model (dont know why space1 is not update previously)
  // for (int i=0;i<inipara.parameter.size();i++)
  // {
  //   fprintf(stderr,"%d %d %g\n",para0.nparav,i,inipara.parameter[i]);
  // }
  // sleep(100);
  sprintf(wmodelname2,"Initial.para");
  fprintf(stderr,"Write the initial para to file:\n");
  write_paras0({inipara},wmodelname2,outdir,invtype);

  //   fprintf(stderr,"Check para space again!!\n");
  //   sleep(1);
  //   fprintf(stderr, "npara=%d, space1.size()=%zu\n",
  //             inipara.npara, inipara.space1.size());
  //   for(i=0;i<inipara.npara;i++)
  //   {
  //     for (size_t j = 0; j < inipara.space1[i].size(); ++j) 
  //     {
  //       fprintf(stderr, "Check para[%zu][%zu]: %g \n",i, j, inipara.space1[i][j]);
  //     }
  // }
  // for(i=0;i<inipara.npara;i++)
  // {
  //   fprintf(stderr, "Check para[%zu]: %g+%g-%g\n",i,inipara.parameter[i], inipara.space1[i][0],inipara.space1[i][1]);
  // }
  // ======================== Now we may safety get a min and max model space from the inipara ==============================
  minpara=inipara; // copy
  maxpara=inipara; // copy
  for(i=0;i<inipara.npara;i++)
  {
    minpara.parameter[i] = inipara.space1[i][0]; // min bound;
    maxpara.parameter[i] = inipara.space1[i][1]; // max bound;
  }
  para2mod(minpara,model,minmodel,invtype);
  updatemodel(minmodel, invtype);

  para2mod(maxpara,model,maxmodel, invtype);
  updatemodel(maxmodel, invtype);
// If need to recalculate the left and right bound, then open this option

  // write_initial_model_bounds(model, para0, para, outdir, invtype);


  // double tdep=0.;
  // for (int i=0;i<minmodel.laym0.nlayer;++i)
  // {
  //   tdep=tdep+minmodel.laym0.thick[i];
  //   // fprintf(stderr," %g %g %g %g %g %g %g\n",tdep,model.laym0.vs[i],model.laym0.vp[i],model.laym0.vpvs[i],model.laym0.rho[i],model.laym0.qs[i],model.laym0.qp[i]);
  //   fprintf(stderr," %g %g %g\n",tdep,minmodel.laym0.vs[i],maxmodel.laym0.vs[i]);
  // }
  //-- write model here -----------------------
  sprintf(wmodelname1,"MC.%s",outname);
  sprintf(outfnmbin,"%s/BIN.%s.1",outdir,outname);
  sprintf(t_str,"mkdir %s\0",outdir);
  system(t_str);

  int n_min,n_min_disp,n_min_rf;

  inbin.open(outfnmbin,ios::in | ios::binary);
  fprintf(stderr,"now begin to read file! %s\n",outfnmbin);
  inbin.seekg(0,ios::end);
  int SSIZE = (int)inbin.tellg();
  cout<<SSIZE<<endl;
  inbin.seekg(0,ios::beg);

  tmodel = model;
  mod2para(model,para0,para, invtype);
  tpara=para;
  cout<<" size "<<para.parameter.size()<<endl;
  read_bin1(inbin,tpara,iaccp,tinfo);
  cout<<iaccp<<endl;
  //for(i=0;i<tpara.npara;i++) fprintf(stderr,"%lf \n",tpara.parameter[i]);


  i=0;
  cout<<"read ok here!"<<endl;
  fprintf(stderr,"read ok here! ...\n");
  paras.clear();
  models.clear();
  flags.clear();
  kai_disp.clear();
  kai_rf.clear();
  kai_joint.clear();

  //Added by EMB:
  sprintf(itmffEMB,"%s/%s.itermisfit",outdir,wmodelname1);
  itmf0emb=fopen(itmffEMB,"w");
  //end of adds

  while (inbin.tellg()<SSIZE) {
    initpara(tpara);
    initmodel(tmodel);
    tpara = para;
    read_bin1(inbin,tpara,iaccp,tinfo);

    tmodel.data.L = tinfo[4];
    tmodel.data.misfit = tinfo[5];
    tmodel.data.rf.misfit = tinfo[6];
    tmodel.data.disp.misfit = tinfo[7];
    tmodel.ngroup = tinfo[9];

    sign = (int)tinfo[0];

    paras.push_back(tpara);
    models.push_back(tmodel);
    accps.push_back(iaccp);
    signs1.push_back(sign);
    flags.push_back(1);
    kai_disp.push_back(tmodel.data.disp.misfit);
    kai_rf.push_back(tmodel.data.rf.misfit);
    kai_joint.push_back(tmodel.data.misfit);//EMB edit
    //kai_joint.push_back(999.);
    
    //Added by EMB.. output iteration vs misfits to file
   // cout<<"testing.. ? :( iiter = "<<tinfo[1]<<" tmodel.data.disp.misfit? "<<tmodel.data.disp.misfit<<endl;
    //fprintf(itmf0emb,"%g %g %g %g\n",tinfo[1],tmodel.data.disp.misfit,tmodel.data.rf.misfit,tmodel.data.misfit);
    //end of adds

    i=i+1;
    }
    fclose(itmf0emb);
  nall=i;
  cout<<"read all over!"<<endl;

  minmisfit = 100.;
  cout<<nall<<endl;
  for (i=0;i<nall;i++) {

    if (paras[i].misfit < minmisfit ) {
      minmisfit = paras[i].misfit;
      j=i;
      }
    }

  //added by EMB, get min/max of each depth using all models (see range of model space)
  FILE *ff0emb;
  FILE *ff0;
  char minnameEMB[300];
  char filevsnameEMB[300];
  double temp_demb;
  // int embnlays=151;//150 km (0km to 150km) -- change this if you want to go deeper or shallower
  int embnlays=int(moddepth/dstep)+1;//HVL: pre-condition of number of layers, change later
  double minmodelsemb[embnlays]; double maxmodelsemb[embnlays];
  double printvstmpmodelemb[embnlays]; double printetmpmodelemb[embnlays]; double printphtmpmodelemb[embnlays];
  int jemb;
  // Check point
  fprintf(stderr,"model depth and number of layer: %f %d\n",moddepth,embnlays);
  sleep(2);
  //
  models3.clear();
  cout<<"paras.size?? "<<paras.size()<<endl;
  fprintf(stderr,"paras.size check: %d\n",paras.size());
  sleep(1);

  if (paras.size()<2)
  {
    fprintf(stderr,"Not enough parameter to estimate the finnal model. abort, abort, abort,!!! \n");
    fprintf(stderr,"paras.size(): %d\n",paras.size());
    return 0;
  }
  // fprintf(stderr,"tpara.npara: %d\n",tpara.npara);
  // sleep(1);
  
  // note: Collect the information from whole model-wise
    for (i=0;i<paras.size();i++) 
    {
        tpara=paras[i];
        //EMB addition for min/max of paras
        if (i==0) // first value free difine
        {
            minparaemb=tpara; 
            maxparaemb=tpara;
        }
        else {
          // note: find min and max range of the parameters
            for (jemb=0;jemb<tpara.npara;jemb++) {
                //crust depth instead of thickness
                if (tpara.parameter[jemb]<minparaemb.parameter[jemb]) {
                    minparaemb.parameter[jemb]=tpara.parameter[jemb];
                }
                if (tpara.parameter[jemb]>maxparaemb.parameter[jemb]) {
                    maxparaemb.parameter[jemb]=tpara.parameter[jemb];
                }
            }
        }
        para2mod(tpara,model,tmodel, invtype);
        updatemodel(tmodel, invtype);
        models3.push_back(tmodel);
    }
    cout<<"models.size?? "<<models3.size()<<endl;
    fprintf(stderr,"models3.size(): %d\n",models3.size());
    sleep(1);
  //initialize minmodels and maxmodels
    for (jemb=0;jemb<embnlays;jemb++) 
    {
      temp_demb=(double)jemb;//depth (km)
      minmodelsemb[jemb]=get_v(models3[0],temp_demb);
      maxmodelsemb[jemb]=minmodelsemb[jemb];
      // fprintf(stderr,"see min max(): %f %f\n",minmodelsemb[jemb],maxmodelsemb[jemb]);
      // sleep(1);
    }
 //
 //EMB
 //Loop through sedimentary thickness parameter and input into seddepth/sedperc
 double ddepthemb=dstep; //delta depth for value
 double percemb; // % per model
 int nvalsemb,tmpval;
 // HVL add the flag of initial model have no sediment setting
 int i_sed=0;
 //
  // fprintf(stderr,"nvalsemb,percemb: %d %f %d %d\n",nvalsemb,percemb,minparaemb.npara, maxparaemb.npara);
  // sleep(3);
 // HVL: Check the q_flag to capture the existance of sediment layers
 for (i=0;i<model.ngroup;i++)
 {
  if (model.groups[i].q_flag==2) // sediment flag
  {
    i_sed=1;
    nvalsemb=int((maxparaemb.parameter[maxparaemb.npara-1]-minparaemb.parameter[minparaemb.npara-1])/ddepthemb);
    percemb=(1.0/(paras.size()))*100.0;
  }
  else
  {
    i_sed=0;
    nvalsemb=1; // no sediment
    percemb=0; // no sediment
  }
 }
// 
 for (i=0;i<nvalsemb;i++) 
 {
     seddepth.push_back(minparaemb.parameter[minparaemb.npara-1]+(i*ddepthemb));
     sedperc.push_back(0);
 }
 //cout<<"nvalsemb "<<nvalsemb<<" percemb "<<percemb<<" paras.size()"<<paras.size()<<endl;
 for (i=0;i<paras.size();i++) {
     tpara=paras[i];
     tmpval=trunc((tpara.parameter[tpara.npara-1]-minparaemb.parameter[minparaemb.npara-1])/ddepthemb);
     //cout<<"tmpval?? "<<tmpval<<" sedperc[tmpval]= "<<sedperc[tmpval]+percemb<<endl;
     sedperc[tmpval]=sedperc[tmpval]+percemb;
 }
 //

double tmpvemb;
double tmpq; // temporary q value
 //Loop through each model, updating min/max as necessary
  for (i=0;i<models3.size();i++){
      for (jemb=0;jemb<embnlays;jemb++) {
          
          temp_demb=double(jemb);
          tmpvemb=get_v(models3[i],temp_demb);
          tmpq=get_q_flag(models3[i],temp_demb);
          //cout<<"checking.. "<<"i "<<i<<" j "<<j<<" "<<temp_demb<<" "<<tmpvemb<<endl;
          if (tmpq !=4) // non-mantle value can greater than 5
          {
            if (tmpvemb<minmodelsemb[jemb] && tmpvemb>0.2 && tmpvemb<5.0) {
                minmodelsemb[jemb]=tmpvemb;
            }
            if (tmpvemb>maxmodelsemb[jemb] && tmpvemb>0.2 && tmpvemb<5.0)
            {
                maxmodelsemb[jemb]=tmpvemb;
            }
          }
          else // mantle value can greater than 5
          {
            if (tmpvemb<minmodelsemb[jemb] && tmpvemb>0.2) {
                minmodelsemb[jemb]=tmpvemb;
            }
            if (tmpvemb>maxmodelsemb[jemb] && tmpvemb>0.2)
            {
                maxmodelsemb[jemb]=tmpvemb;
            }
          }
      }
  }

/////////////////////////////////////////////////////////////
// Go through and output information for every model (not just posterior) //
/////////////////////////////////////////////////////////////
// double dep;
 /*sprintf(filevsnameEMB,"%s/%s.ALLmodels_vs",outdir,wmodelname1);
  ff0emb=fopen(filevsnameEMB,"w");
//New way, write similar to write_model in CALmodel.C for namemod1
  // Outputs ALL layers...
    for (i=0;i<models3.size();i++) {
        for (j=0;j<models3[i].laym0.nlayer;j++) { //write layered model in stair-step..
            //cout<<"testing... j "<<j<<" models1[i].laym0.vs[j] "<<models1[i].laym0.vs[j]<<endl;
            fprintf(ff0emb,"%g %g ",models3[i].laym0.vs[j],models3[i].laym0.vs[j]);
        }
        fprintf(ff0emb,"\n");
    }
    fclose(ff0emb);

    sprintf(filevsnameEMB,"%s/%s.ALLmodels_depth",outdir,wmodelname1);
    ff0emb=fopen(filevsnameEMB,"w");
    for (i=0;i<models3.size();i++) {
        dep=0.;
        for (j=0;j<models3[i].laym0.nlayer;j++) { //write layered model in stair-step..
            fprintf(ff0emb,"%g ",dep);
            dep=dep+models3[i].laym0.thick[j];
            fprintf(ff0emb,"%g ",dep);
        }   
        fprintf(ff0emb,"\n");
    }   
    fclose(ff0emb);

sprintf(filevsnameEMB,"%s/%s.ALLmodels_vp",outdir,wmodelname1);
    ff0emb=fopen(filevsnameEMB,"w");
    for (i=0;i<models3.size();i++) {
        dep=0.;
        for (j=0;j<models3[i].laym0.nlayer;j++) { //write layered model in stair-step..
            fprintf(ff0emb,"%g %g ",models3[i].laym0.vp[j],models3[i].laym0.vp[j]);
        }   
        fprintf(ff0emb,"\n");
    }   
    fclose(ff0emb);

 sprintf(filevsnameEMB,"%s/%s.ALLmodels_rho",outdir,wmodelname1);
    ff0emb=fopen(filevsnameEMB,"w");
    for (i=0;i<models3.size();i++) {
        dep=0.;
        for (j=0;j<models3[i].laym0.nlayer;j++) { //write layered model in stair-step..
            fprintf(ff0emb,"%g %g ",models3[i].laym0.rho[j],models3[i].laym0.rho[j]);
        }   
        fprintf(ff0emb,"\n");
    }   
    fclose(ff0emb);
*/

 /*ff0emb=fopen(filevsnameEMB,"w");
  for (i=0;i<models3.size();i++){
      for (jemb=0;jemb<embnlays;jemb++) {
          
          temp_demb=double(jemb);
          tmpvemb=get_v(models3[i],temp_demb);
          printvstmpmodelemb[jemb]=tmpvemb;
          //cout<<"checking.. "<<"i "<<i<<" j "<<j<<" "<<temp_demb<<" "<<tmpvemb<<endl;
          if (tmpvemb<minmodelsemb[jemb] && tmpvemb>0.2 && tmpvemb<5.0) {
              minmodelsemb[jemb]=tmpvemb;
          }
          if (tmpvemb>maxmodelsemb[jemb] && tmpvemb>0.2 && tmpvemb<5.0)
          {
              maxmodelsemb[jemb]=tmpvemb;
      
          }
      }
      //print out tmpvemb
      fprintf(ff0emb,"%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", i, printvstmpmodelemb[0], printvstmpmodelemb[1], printvstmpmodelemb[2], printvstmpmodelemb[3], printvstmpmodelemb[4], printvstmpmodelemb[5], printvstmpmodelemb[6], printvstmpmodelemb[7], printvstmpmodelemb[8], printvstmpmodelemb[9], printvstmpmodelemb[10], printvstmpmodelemb[11], printvstmpmodelemb[12], printvstmpmodelemb[13], printvstmpmodelemb[14], printvstmpmodelemb[15], printvstmpmodelemb[16], printvstmpmodelemb[17], printvstmpmodelemb[18], printvstmpmodelemb[19], printvstmpmodelemb[20], printvstmpmodelemb[21], printvstmpmodelemb[22], printvstmpmodelemb[23], printvstmpmodelemb[24], printvstmpmodelemb[25], printvstmpmodelemb[26], printvstmpmodelemb[27], printvstmpmodelemb[28], printvstmpmodelemb[29], printvstmpmodelemb[30], printvstmpmodelemb[31], printvstmpmodelemb[32], printvstmpmodelemb[33], printvstmpmodelemb[34], printvstmpmodelemb[35], printvstmpmodelemb[36], printvstmpmodelemb[37], printvstmpmodelemb[38], printvstmpmodelemb[39], printvstmpmodelemb[40], printvstmpmodelemb[41], printvstmpmodelemb[42], printvstmpmodelemb[43], printvstmpmodelemb[44], printvstmpmodelemb[45], printvstmpmodelemb[46], printvstmpmodelemb[47], printvstmpmodelemb[48], printvstmpmodelemb[49], printvstmpmodelemb[50]);

  }
  fclose(ff0emb);
*/
 /*
  // Get ellipticity and phase for each model..
  cout<<" domcv3 getting disp etc??? "<<endl;
  models2.clear();
  for (i=0;i<models3.size();i++){
      tmodel=models3[i];
      compute_disp(tmodel);
      compute_rf(tmodel);
      models2.push_back(tmodel);
      }
cout<<" end of domcv3 get all "<<endl;
//print ellipticity to file
    sprintf(filevsnameEMB,"%s/%s.ALLmodels_HV",outdir,wmodelname1);
    ff0emb=fopen(filevsnameEMB,"w");
    for (i=0;i<models2.size();i++) {
    tmodel = models2[i];
    if (tmodel.data.disp.eper.size() > 0) {
      fprintf(ff0emb,"M %d",i);
      for (jemb=0;jemb<tmodel.data.disp.eper.size();jemb++) {
        fprintf(ff0emb," %g",tmodel.data.disp.eper[jemb]);//only need period and ellipticity
      }
      for (jemb=0;jemb<tmodel.data.disp.eper.size();jemb++) {
        fprintf(ff0emb," %g",tmodel.data.disp.evel[jemb]);//only need period and ellipticity
        //fprintf(ff0,"%g %g %g %g\n",tmodel.data.disp.eper[j],tmodel.data.disp.evel[j],tmodel.data.disp.evelo[j],tmodel.data.disp.unevelo[j]);
        }
      } 
        fprintf(ff0emb,"\n");
    }  
  fclose(ff0emb);

//print phase to file
    sprintf(filevsnameEMB,"%s/%s.ALLmodels_phase",outdir,wmodelname1);
    ff0emb=fopen(filevsnameEMB,"w");
    for (i=0;i<models2.size();i++) {
    tmodel = models2[i];
    if (tmodel.data.disp.pper.size() > 0) {
      fprintf(ff0emb,"M %d",i);
      for (jemb=0;jemb<tmodel.data.disp.pper.size();jemb++) {
        fprintf(ff0emb," %g",tmodel.data.disp.pper[jemb]);
        }
      for (jemb=0;jemb<tmodel.data.disp.pper.size();jemb++) {
        fprintf(ff0emb," %g",tmodel.data.disp.pvel[jemb]);
        }
      } 
        fprintf(ff0emb,"\n");
    } 
  fclose(ff0emb);
*/
// End of output every model info //
/////////////////////////////////////////////////////////////


  //write out results of search area
 sprintf(minnameEMB,"%s/%s.minmodelvalue",outdir,wmodelname1);
  ff0emb = fopen(minnameEMB,"w");
  if (plottype==1) // gradient style
  {
    for (jemb=0;jemb<embnlays;jemb++)
    {
      fprintf(ff0emb,"%d %g\n",jemb,minmodelsemb[jemb]);
    }
  }
  else if (plottype==0)
  {
   for (jemb=0;jemb<embnlays-1;jemb++)
   {
    fprintf(ff0emb,"%d %g\n",jemb,minmodelsemb[jemb]);
    fprintf(ff0emb,"%d %g\n",jemb+1,minmodelsemb[jemb]); 
   }
   fprintf(ff0emb,"%d %g\n",embnlays,minmodelsemb[embnlays]);
  }
    fclose(ff0emb);
// -------------------
  sprintf(minnameEMB,"%s/%s.maxmodelvalue",outdir,wmodelname1);
  ff0emb = fopen(minnameEMB,"w");
  if (plottype==1) // gradient style
  {
    for (jemb=0;jemb<embnlays;jemb++)
    {
      fprintf(ff0emb,"%d %g\n",jemb,maxmodelsemb[jemb]);
    }
  }
  else if (plottype==0) // gradient style
  {
    for (jemb=0;jemb<embnlays-1;jemb++)
   {
    fprintf(ff0emb,"%d %g\n",jemb,maxmodelsemb[jemb]);
    fprintf(ff0emb,"%d %g\n",jemb+1,maxmodelsemb[jemb]); 
   }
   fprintf(ff0emb,"%d %g\n",embnlays,maxmodelsemb[embnlays]);
  }
    fclose(ff0emb);
  //
  //  End of EMB addition
 
 //Extra EMB addition, write min/max of all parameters in model
//output average model parameters
// loop through and find min/max of all paras (whatever goes into para2mod..)
// if (i_sed==1)
// {
    sprintf(minnameEMB,"%s/%s.minparas",outdir,wmodelname1);
    ff0emb = fopen(minnameEMB,"w");
      for (jemb=0;jemb<minparaemb.npara;jemb++){
                fprintf(ff0emb,"%g ",minparaemb.parameter[jemb]);
                  }
          fclose(ff0emb);
// }
    sprintf(minnameEMB,"%s/%s.maxparas",outdir,wmodelname1);
    ff0emb = fopen(minnameEMB,"w");
      for (jemb=0;jemb<maxparaemb.npara;jemb++){
                fprintf(ff0emb,"%g ",maxparaemb.parameter[jemb]);
                  }
          fclose(ff0emb);

//output percent of models in 0.1km/s from min to max
    sprintf(minnameEMB,"%s/%s.percseddepth",outdir,wmodelname1);
    ff0emb = fopen(minnameEMB,"w");
      for (jemb=0;jemb<nvalsemb;jemb++)
      {
        fprintf(ff0emb,"%g %g\n",seddepth[jemb],sedperc[jemb]);
      }
    fclose(ff0emb);

// ========================================================================
// Now we search whole the posterior model to get the final subset models 
// ========================================================================

  fprintf(stderr,"Now process the posterior \n");
  minmisfit_rf=99999.;
  minmisfit_disp=99999.;
  k= 0;
  paras1.clear();
  models1.clear();
  // n_min = 100000;
  // n_min_disp = 100000;
  // n_min_rf = 100000;

  cout<<"minmisfit?"<<endl;
  // Condition check
  if (kai_disp.empty() || kai_rf.empty() || kai_joint.empty()) {
    fprintf(stderr,"Empty misfit array!\n");
    return 0;
  }
  minmisfit_disp  = *std::min_element(kai_disp.begin(),kai_disp.end());
  minmisfit_rf = *std::min_element(kai_rf.begin(),kai_rf.end());
  kai_min = *std::min_element(kai_joint.begin(),kai_joint.end());
  if ((kai_min<1000.0) && (kai_min>0.0))
  {
    if ((mc_range > 100.) || (mc_range <= 0.))
    {
      fprintf(stderr,"Wrong in model selection range %f\n",mc_range);
      return(0);
    }
  // Ranging the disp by percentage
  if ((ihv+iph+igv==0) || (kai_disp.size() < 200))  // no dispersion information or number of posterior less than 2000
  {
    // minmisfit_disp_o = minmisfit_disp;
    auto maxmisfit_disp = *std::max_element(kai_disp.begin(),kai_disp.end());
    minmisfit_disp_o = maxmisfit_disp;
  }
  else
  {
    size_t index_range_disp = static_cast<size_t>((mc_range/100.) * (kai_disp.size()));
    if (index_range_disp >= kai_disp.size()) index_range_disp = kai_disp.size() - 1;
    // size_t index_range_disp = static_cast<size_t>(2 * (mc_range/100.) * (kai_disp.size()));
    // Create another vector to copy sorted elements
    vector<float> arr0(kai_disp.size());  // Allocate the same size for arr
    // Copy sorted vec into arr
    copy(kai_disp.begin(), kai_disp.end(), arr0.begin());
    // Sort the vector in ascending order
    sort(arr0.begin(), arr0.end());
    minmisfit_disp_o = arr0[index_range_disp];
    // minmisfit_disp_o = kai_disp[index_range_disp]; // same above line
    if (minmisfit_rf_o == 0)
    {
      minmisfit_disp_o = minmisfit_disp*2.5;
    }
    fprintf(stderr,"Disp misfit out: %f %d %f\n",minmisfit_disp_o,index_range_disp,minmisfit_disp);
  }
  
  // Ranging the receiver function by percentage
  if ((irf==0) || (kai_rf.size() < 200)) 
  {
    // minmisfit_rf_o = minmisfit_rf;
    auto maxmisfit_rf = *std::max_element(kai_rf.begin(),kai_rf.end());
    minmisfit_rf_o = maxmisfit_rf;
  }
  else
  {
    size_t index_range_rf = static_cast<size_t>(int((mc_range/100.) * (kai_rf.size())));
    if (index_range_rf >= kai_rf.size()) index_range_rf = kai_rf.size() - 1;
    // size_t index_range_rf = static_cast<size_t>(int(2 * (mc_range/100.) * (kai_rf.size())));
    // Create another vector to copy sorted elements
    vector<float> arr1(kai_rf.size());  // Allocate the same size for arr
    // Copy sorted vec into arr
    copy(kai_rf.begin(), kai_rf.end(), arr1.begin());
    // Sort the vector in ascending order
    sort(arr1.begin(), arr1.end());
    minmisfit_rf_o = arr1[index_range_rf];
    if (minmisfit_rf_o == 0)
    {
      minmisfit_rf_o = minmisfit_rf*2.5;
    }
    fprintf(stderr,"RF misfit out: %f %d\n",minmisfit_rf_o,index_range_rf);
  }
  if (kai_joint.size() >= 200) // Here, if all posterior model < 200 then take all
  {
    // Ranging the model misfit by percentage
  size_t index_range = static_cast<size_t>(int((mc_range/100.) * (kai_joint.size())));
  if (index_range >= kai_joint.size()) index_range = kai_joint.size() - 1;
  fprintf(stderr,"index_range %zu \n",index_range);
    // Create another vector to copy sorted elements
    vector<float> arr(kai_joint.size());  // Allocate the same size for arr
    // Copy sorted vec into arr
    copy(kai_joint.begin(), kai_joint.end(), arr.begin());
    // Sort the vector in ascending order
    sort(arr.begin(), arr.end());
  // take the misfit value at element xxx
  kai_cri = arr[index_range];
  }
  else // take all 100%
  {
  size_t index_range = static_cast<size_t>((kai_joint.size()))-1;
  fprintf(stderr,"index_range %zu \n",index_range);
    // Create another vector to copy sorted elements
    vector<float> arr(kai_joint.size());  // Allocate the same size for arr
    // Copy sorted vec into arr
    copy(kai_joint.begin(), kai_joint.end(), arr.begin());
    // Sort the vector in ascending order
    sort(arr.begin(), arr.end());
  // take the misfit value at element xxx
  kai_cri = arr.back(); // Take all model
  }

  }
  else
  {
    auto maxmisfit_disp = *std::max_element(kai_disp.begin(),kai_disp.end());
    minmisfit_disp_o = maxmisfit_disp;

    auto maxmisfit_rf = *std::max_element(kai_rf.begin(),kai_rf.end());
    minmisfit_rf_o = maxmisfit_rf;

    // auto maxmisfit_joint = *std::max_element(kai_joint.begin(),kai_joint.end());
    double misfit_joint_criterion = *std::min_element(kai_joint.begin(),kai_joint.end());
    kai_cri = 2*misfit_joint_criterion;
  }
  // if (pp<0.00001) minmisfit_disp_o = 99999999.;
  // if (pp>0.99999) minmisfit_rf_o = 99999999.;
   
    fprintf(stderr,"kai_min %g - kai_cri %g - msifit_cri_disp %g - msifit_cri_rf %g \n",kai_min, kai_cri,minmisfit_disp_o, minmisfit_rf_o);
    sleep(5);

    cout<<"kai_cri "<<kai_cri<<endl;
    // cout<<"kai_cri "<<kai_cri<<" kai_min "<<kai_min<<" misfit_disp_min"<<minmisfit_disp_o<<" misfit_rf_min"<<minmisfit_rf_o<<endl;
  // if (model.indata.gvflag==0 and model.indata.phflag==0 and model.indata.hvflag==0) // no dispersion data
  
  if (std::isnan(minmisfit_disp_o)) minmisfit_disp_o= 99999999.;
  if (std::isnan(minmisfit_rf_o)) minmisfit_rf_o= 99999999.;
  
  // fprintf(stderr,"kai_cri %g kai_min %g misfit_disp_min %g misfit_rf_min %g\n",kai_cri,kai_min,minmisfit_disp_o,minmisfit_rf_o);
  // search a new min misfit for  new model
  // minmisfit_rf=minmisfit_rf_o; // fake values
  // minmisfit_disp=minmisfit_disp_o; // fake values

  // check the option
  // fprintf(stderr,"switch option %d \n",option);
  // sleep(5);
  // registing the temporaly minminfits for a new selection group:
  double minmisfit_rf_new, minmisfit_disp_new, kai_min_new;

  // HVL: 20250104 find the possible for model 
  
  bool found_joint = false;
  bool found_disp  = false;
  bool found_rf    = false;

  n_min_rf = -1;
  n_min_disp = -1;
  n_min = -1;
  // We search for whole population. But only n_min wave use later. 

  // copy copy copy
  minmisfit_rf_new = minmisfit_rf_o;
  minmisfit_disp_new = minmisfit_disp_o;
  kai_min_new = kai_cri;

  for (i=0;i<nall;i++) 
  {
     switch (option) 
     {
      case 1: // Both RF and DISP
      
      if (kai_joint[i] < kai_cri and kai_disp[i] <= minmisfit_disp_o and kai_rf[i] <= minmisfit_rf_o and signs1[i] > 0) 
      {
        fprintf(stderr," >> kai_joint,kai_rf,kai_disp: %f %f %f\n",kai_joint[i],kai_rf[i],kai_disp[i]);
        // find misfit index if model

        // flags[i] = 1; // flags on binary file, no need to change
        if (!found_joint || kai_joint[i] < kai_min_new)
        {
          kai_min_new = kai_joint[i];
          n_min = i;
          found_joint = true;
          // fprintf(stderr,"> > > new n_min_disp: %d\n",n_min);
        }
        // find misfit index of receiver functions
        if (!found_rf || kai_rf[i] < minmisfit_rf_new)
        {
          minmisfit_rf_new = kai_rf[i];
          n_min_rf = i;
          found_rf = true;
          // fprintf(stderr,"> new n_min_rf: %d\n",n_min_rf);
        }
        // find misfit index of dispersion 
        if (!found_disp || kai_disp[i] < minmisfit_disp_new)
        {
          minmisfit_disp_new = kai_disp[i];
          n_min_disp = i;
          found_disp = true;
          // fprintf(stderr,"> > new n_min_disp: %d\n",n_min_disp);
        }
      }   
      else 
      {
        flags[i] = 0;
      }
      fprintf(stderr,"kai_cri,kai_rf_cri,kai_disp_cri: %f %f %f\n",kai_min_new,minmisfit_rf_new,minmisfit_disp_new);
      fprintf(stderr,"nmin,nminrf,nmindisp: %d %d %d\n",n_min,n_min_rf,n_min_disp);
      break; // end case 1

      case 2: // no DISP

      if (kai_joint[i] < kai_cri and  kai_rf[i] <= minmisfit_rf_o and signs1[i] > 0) 
      {
        // find misfit index if model
        // flags[i] = 1; // flags already writen in binary file
        if (!found_joint || kai_joint[i] < kai_min_new)
        {
          kai_min_new = kai_joint[i];
          n_min = i;
          found_joint = true;
          // fprintf(stderr,"> > > new n_min_disp: %d\n",n_min);
        }
        // find misfit index of receiver functions
        if (!found_rf || kai_rf[i] < minmisfit_rf_new)
        {
          minmisfit_rf_new = kai_rf[i];
          n_min_rf = i;
          found_rf = true;
          // fprintf(stderr,"> new n_min_rf: %d\n",n_min_rf);
        }
      }   
      else 
      {
        flags[i] = 0;
      }
      break; // end case 2

      case 3: // No RF
      
      if (kai_joint[i] < kai_cri && kai_disp[i] <= minmisfit_disp_o && signs1[i] > 0) 
      {
        fprintf(stderr," >> kai_joint,kai_rf,kai_disp: %f %f %f\n",kai_joint[i],kai_rf[i],kai_disp[i]);
        // find misfit index if model

        // flags[i] = 1; // flags already writen in binary file
        if (!found_joint || kai_joint[i] < kai_min_new)
        {
          kai_min_new = kai_joint[i];
          n_min = i;
          found_joint = true;
          // fprintf(stderr,"> > > new n_min_disp: %d\n",n_min);
        }
        // find misfit index of dispersion 
        if (!found_disp || kai_disp[i] < minmisfit_disp_new)
        {
          minmisfit_disp_new = kai_disp[i];
          n_min_disp = i;
          found_disp = true;
          // fprintf(stderr,"> > new n_min_disp: %d\n",n_min_disp);
        }
      }   
      else 
      {
        flags[i] = 0;
      }
      break; // end case 3

    }    
  } // end nall loop
  /// report the output
  // fprintf(stderr,"nmin,nminrf,nmindisp: %d %d %d\n",n_min,n_min_rf,n_min_disp);
  // sleep(5);
  if (!found_joint || n_min < 0)
  {
      fprintf(stderr, "No posterior model passed the joint selection criteria. Stop.\n");
      inbin.close();
      return 0;
  }
// ========================================================================
// Finish conduct final subset models 
// ========================================================================
     
  // cout<<k<<" "<<nall<<" "<<n_min_disp<<" "<<n_min_rf<<" "<<n_min<<endl;
  cout<<k<<" "<<nall<<" "<<n_min_disp<<" "<<n_min_rf<<" "<<n_min<<endl;

  fprintf(stderr,"nall: %d %d %d %d\n",nall,n_min_disp,n_min_rf,n_min);
  sleep(2);
  paras1.clear();
  j=0;
  for(i=0;i<nall;i++) {
    if (flags[i]==0) continue;
    j++;  
    paras1.push_back(paras[i]);
    fprintf(stderr,"pushbackparas no %d! id: %d \n",j, i);
    cout<<"posterior:"<<" "<<i<<" "<<j<<endl;
    }
  fprintf(stderr,"getting average model?\n");
  cout<<"getting average model.."<<endl;
  if (paras1.empty()) {
    fprintf(stderr, "No accepted posterior models after selection for average model.\n");
    return 0;
  } 
  get_ave_para(paras1,ave_para);

  sprintf(wmodelname2,"%s.acc.paras",wmodelname1);
  write_paras0(paras1,wmodelname2,outdir,invtype);
  
  //output average model parameters
  sprintf(minnameEMB,"%s/%s.ave.paras",outdir,wmodelname1);
  ff0emb = fopen(minnameEMB,"w");
  for (jemb=0;jemb<ave_para.npara;jemb++)
  {
    fprintf(ff0emb,"%g ",ave_para.parameter[jemb]);
  }
  fclose(ff0emb);

  tpara=ave_para;
  para2mod(tpara,model,tmodel, invtype);
  updatemodel(tmodel, invtype);
  sprintf(wmodelname2,"%s.acc.average",wmodelname1);
  // sprintf(wmodelname2,"%s.mean",wmodelname1);
  get_misfit1(tpara,tmodel,pp,invtype);
  write_model(tmodel,wmodelname2,outdir,flagw);
  cout<<"done writing average model.."<<endl;

  //output average model parameters
  if (j == 0) {
    fprintf(stderr, "No posterior model passed selection criteria.\n");
    return 0;
  }

  tpara=paras[n_min];
  //for (i=0;i<tpara.npara;i++) {fprintf(stderr,"pp: %d %g \n",i,tpara.parameter[i]);}

  para2mod(tpara,model,tmodel, invtype);
  updatemodel(tmodel, invtype);
  sprintf(wmodelname2,"%s.minmisfit",wmodelname1);
  get_misfit1(tpara,tmodel,pp, invtype);
  write_model(tmodel,wmodelname2,outdir,flagw);

  //forward model joint minimum misfit model...
  // tpara=paras[n_min_disp];
  tpara=paras[n_min]; //should be n_min here!
  para2mod(tpara,model,tmodel, invtype);
  updatemodel(tmodel, invtype);
  sprintf(wmodelname2,"%s.minmisfit_disp",wmodelname1);
  get_misfit1(tpara,tmodel,pp, invtype);
  write_model(tmodel,wmodelname2,outdir,flagw);

  // tpara=paras[n_min_rf];
  tpara=paras[n_min]; //should be n_min here!
  para2mod(tpara,model,tmodel, invtype);
  updatemodel(tmodel, invtype);
  sprintf(wmodelname2,"%s.minmisfit_rf",wmodelname1);
  get_misfit1(tpara,tmodel,pp, invtype);
  write_model(tmodel,wmodelname2,outdir,flagw);
  inbin.close();

  //cout<<"now let's do this!"<<endl;
// ======================================================================================
// Calculate the average model for all selected posterior
// ======================================================================================
  models1.clear();
  for (i=0; i<paras1.size(); i++) {
    tpara = paras1[i];
    para2mod(tpara, model, tmodel, invtype);
    updatemodel(tmodel, invtype);
    models1.push_back(tmodel);
  }

  cout << "now let's do this! " << models1.size() << endl;

  // ============================================================
  // New average model on common depth grid
  // Avoid old stitched group-mean method, which can create
  // artificial negative jumps even if each model is monotonic.
  // ============================================================
  fprintf(stderr,"Compute average model on common depth grid...\n");

  double max_depth_all = 0.0;
  for (i = 0; i < (int)models1.size(); i++) {
    double tdep = 0.0;
    for (j = 0; j < models1[i].laym0.nlayer; j++) {
      tdep += models1[i].laym0.thick[j];
    }
    if (tdep > max_depth_all) max_depth_all = tdep;
  }

  // number of depth samples
  int NN = int(max_depth_all / dstep) + 1;

  vector<double> depth_out, mean_out, std_out, min_out, max_out, n1_out, p1_out;
  depth_out.clear();
  mean_out.clear();
  std_out.clear();
  min_out.clear();
  max_out.clear();
  n1_out.clear();
  p1_out.clear();

  double temp_d, temp_v;
  double temp_v_mean, temp_v_std, temp_v_min, temp_v_max, temp_v_n1std, temp_v_p1std;
  vector<double> t_values;

  for (i = 0; i < NN; i++) {
    temp_d = i * dstep;
    t_values.clear();

    // collect Vs from all posterior models at this depth
    for (j = 0; j < (int)models1.size(); j++) {
      temp_v = get_vs_at_depth_common_grid(models1[j], temp_d);
      if (temp_v > 0.0) t_values.push_back(temp_v);
    }

    if (t_values.size() < 2) {
      temp_v_mean  = -1.0;
      temp_v_std   = -1.0;
      temp_v_min   = -1.0;
      temp_v_max   = -1.0;
      temp_v_n1std = -1.0;
      temp_v_p1std = -1.0;
    }
    else {
      temp_v_mean  = cv_mean(t_values);
      temp_v_std   = cv_std(t_values);
      temp_v_min   = cv_min(t_values);
      temp_v_max   = cv_max(t_values);
      temp_v_n1std = temp_v_mean - temp_v_std;
      temp_v_p1std = temp_v_mean + temp_v_std;
    }

    depth_out.push_back(temp_d);
    mean_out.push_back(temp_v_mean);
    std_out.push_back(temp_v_std);
    min_out.push_back(temp_v_min);
    max_out.push_back(temp_v_max);
    n1_out.push_back(temp_v_n1std);
    p1_out.push_back(temp_v_p1std);
  }

  // ============================================================
  // Enforce monotonic increase on final average model
  // This prevents artificial negative movement of the mean curve.
  // ============================================================
  bool FORCE_MONOTONIC_MEAN = true;

  if (FORCE_MONOTONIC_MEAN) {
    for (i = 1; i < (int)mean_out.size(); i++) {
      if (mean_out[i] > 0.0 && mean_out[i-1] > 0.0 && mean_out[i] < mean_out[i-1]) {
        mean_out[i] = mean_out[i-1];
      }
    }

    // optional: keep -1std and +1std also non-decreasing
    for (i = 1; i < (int)n1_out.size(); i++) {
      if (n1_out[i] > 0.0 && n1_out[i-1] > 0.0 && n1_out[i] < n1_out[i-1]) {
        n1_out[i] = n1_out[i-1];
      }
    }

    for (i = 1; i < (int)p1_out.size(); i++) {
      if (p1_out[i] > 0.0 && p1_out[i-1] > 0.0 && p1_out[i] < p1_out[i-1]) {
        p1_out[i] = p1_out[i-1];
      }
    }

    // keep min/max physically consistent with corrected mean if needed
    for (i = 0; i < (int)mean_out.size(); i++) {
      if (mean_out[i] > 0.0 && min_out[i] > 0.0 && min_out[i] > mean_out[i]) {
        min_out[i] = mean_out[i];
      }
      if (mean_out[i] > 0.0 && max_out[i] > 0.0 && max_out[i] < mean_out[i]) {
        max_out[i] = mean_out[i];
      }
    }
  }

  fprintf(stderr,"now write average model in uniform step!\n");
  char minname[300];
  sprintf(minname,"%s/%s.ave.mod",outdir,wmodelname1);
  ff0 = fopen(minname,"w");
  cout << minname << endl;

  for (i = 0; i < (int)depth_out.size(); i++) {
    fprintf(ff0,"%g %g %g %g %g %g %g\n",
            depth_out[i],
            mean_out[i],
            std_out[i],
            min_out[i],
            max_out[i],
            n1_out[i],
            p1_out[i]);
  }
  fclose(ff0);

  fprintf(stderr,"now write avegrage model !\n");

  fprintf(stderr,"now write posterior !\n");
  //double temp1, temp2, temp4, temp8, temp12, temp16, temp20, temp24, temp28, temp32, temp36, temp40, temp44, temp50, cb, mt, mm5, mp5;
    double cb, mt, mm5, mp5;
    double temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25, temp26, temp27, temp28, temp29, temp30, temp31, temp32, temp33, temp34, temp35, temp36, temp37, temp38, temp39, temp40, temp41, temp42, temp43, temp44, temp45, temp46, temp47, temp48, temp49, temp50;
    int nng = 1;
double dep;

/////////////////////////////////////////////////////////////
// Output Posterior model information (paramters, Vp, Vs, depth, Vp, density). //
/////////////////////////////////////////////////////////////

//Output all parameters for posterior models
  sprintf(minname,"%s/%s.all.parameters",outdir,wmodelname1);
  ff0=fopen(minname,"w");
  for (i=0;i<paras1.size();i++) {
        tpara=paras1[i];
        for (j=0;j<tpara.npara;j++) {
            fprintf(ff0,"%g ",tpara.parameter[j]);
           // cout<<"j "<<j<<" para "<<tpara.parameter[j]<<endl;//<<" para "<<tpara.para0[j]<<endl;
        }
        fprintf(ff0,"\n");
  }
  fclose(ff0);


  // Outputs ALL layers Vs from posterior models...
  sprintf(minname,"%s/%s.all.value",outdir,wmodelname1);
  ff0 = fopen(minname,"w");
    for (i=0;i<models1.size();i++) {
        for (j=0;j<models1[i].laym0.nlayer;j++) { //write layered model in stair-step..
            //cout<<"testing... j "<<j<<" models1[i].laym0.vs[j] "<<models1[i].laym0.vs[j]<<endl;
            fprintf(ff0,"%g %g ",models1[i].laym0.vs[j],models1[i].laym0.vs[j]);
        }
        fprintf(ff0,"\n");
    }
    fclose(ff0);

    sprintf(minname,"%s/%s.all.value_depth",outdir,wmodelname1);
    ff0=fopen(minname,"w");
    for (i=0;i<models1.size();i++) {
        dep=0.;
        for (j=0;j<models1[i].laym0.nlayer;j++) { //write layered model in stair-step..
            fprintf(ff0,"%g ",dep);
            dep=dep+models1[i].laym0.thick[j];
            fprintf(ff0,"%g ",dep);
        }   
        fprintf(ff0,"\n");
    }   
    fclose(ff0);

   sprintf(minname,"%s/%s.all.value_vp",outdir,wmodelname1);
    ff0=fopen(minname,"w");
    for (i=0;i<models1.size();i++) {
        dep=0.;
        for (j=0;j<models1[i].laym0.nlayer;j++) { //write layered model in stair-step..
            fprintf(ff0,"%g %g ",models1[i].laym0.vp[j],models1[i].laym0.vp[j]);
        }   
        fprintf(ff0,"\n");
    }   
    fclose(ff0);

 sprintf(minname,"%s/%s.all.value_rho",outdir,wmodelname1);
    ff0=fopen(minname,"w");
    for (i=0;i<models1.size();i++) {
        dep=0.;
        for (j=0;j<models1[i].laym0.nlayer;j++) { //write layered model in stair-step..
            fprintf(ff0,"%g %g ",models1[i].laym0.rho[j],models1[i].laym0.rho[j]);
        }   
        fprintf(ff0,"\n");
    }   
    fclose(ff0);


 
/////////////////////////////////////////////////////////////
// FORWARD MODEL THE POSTERIOR! //
/////////////////////////////////////////////////////////////
    fprintf(stderr,"now do forward computation for all %d models!\n",models1.size());

  models2.clear();
  for (i=0;i<models1.size();i++) {
    fprintf(stderr,"> forward computation model %d \n",i);
    tmodel = models1[i];
    compute_disp(tmodel,invtype);
    compute_rf(tmodel,invtype);
    models2.push_back(tmodel);
    }
  fprintf(stderr,"compute over! now write!\n",models1.size());
  char ename[300];


  sprintf(minname,"%s/%s.all.e",outdir,wmodelname1);
  ff0 = fopen(minname,"w");
  for (i=0;i<models2.size();i++) {
    tmodel = models2[i];
    if (tmodel.data.disp.eper.size() > 0) {
      fprintf(ff0,"M %d",i);
      for (j=0;j<tmodel.data.disp.eper.size();j++) {
        fprintf(ff0," %g",tmodel.data.disp.eper[j]);//only need period and ellipticity
      }
      for (j=0;j<tmodel.data.disp.eper.size();j++) {
        fprintf(ff0," %g",tmodel.data.disp.evel[j]);//only need period and ellipticity
        //fprintf(ff0,"%g %g %g %g\n",tmodel.data.disp.eper[j],tmodel.data.disp.evel[j],tmodel.data.disp.evelo[j],tmodel.data.disp.unevelo[j]);
        }
      } 
        fprintf(ff0,"\n");
    }  
  fclose(ff0);
  //  ===================================================================================
  // Do it again - dont know why if add the loop, can not write data into text file
  //  ===================================================================================
  int npts = -1;
  std::vector<double> mean_eper,mean_evel;   // mean of rfn[j]
  long long nused = 0;
  for (i = 0; i < (int)models2.size(); i++) 
  {
    tmodel = models2[i];
  
    if (tmodel.data.disp.eper.size() > 0) 
    {
      // allocate once from the first valid RF model
      if (npts < 0) 
      {
        npts = tmodel.data.disp.eper.size();
        mean_eper.assign(npts, 0.0);
        mean_evel.assign(npts, 0.0);
      }
      // skip inconsistent-length RFs
      if (tmodel.data.disp.eper.size() != npts) continue;
      // hard safety: skip if any vector is missing/short (prevents segfault)
      if ((int)tmodel.data.disp.eper.size() < npts) continue;

      for (j = 0; j < npts; j++) 
      {
        // time axes: store from first usable model (or overwrite—either is fine)
        mean_eper[j] = tmodel.data.disp.eper[j];

        // accumulate sums (we will divide later)
        mean_evel[j]   += tmodel.data.disp.evel[j];
      }
      ++nused;   // count models (not samples)
    }
  }

  // convert sums -> means
  for (j = 0; j < npts; j++) 
  {
    mean_evel[j]   /= (double)nused;
  }

  sprintf(minname, "%s/%s.mean.e",outdir,wmodelname1);
  fprintf(stderr,"Open the mean file: %s\n",minname);
  fprintf(stderr,">> \n");
  ff0 = fopen(minname, "w");
  for (i = 0; i < npts; i++) 
  {
    fprintf(ff0, "%g %g\n",
            mean_eper[i],mean_evel[i]);
  }
  fclose(ff0);

  sprintf(minname,"%s/%s.all.rf",outdir,wmodelname1);
  ff0 = fopen(minname,"w");
  for (i=0;i<models2.size();i++) {
    tmodel = models2[i];
    if (tmodel.data.rf.nrfo > 0) 
    {
      fprintf(ff0,">> M %d\n",i);
      for (j=0;j<tmodel.data.rf.nrfo;j++) 
      {
        fprintf(ff0,"%g %g %g %g\n",tmodel.data.rf.to[j],tmodel.data.rf.rfn[j],tmodel.data.rf.rfo[j],tmodel.data.rf.unrfo[j]);
      }
    } 
  } 
  fclose(ff0);
//  ===================================================================================
// Do it again - dont know why if add the loop, can not write data into text file
//  ===================================================================================
  npts = -1;
  std::vector<double> mean_rfn,mean_rfo,mean_unrfo;   // mean of rfn[j]
  std::vector<double> mean_to,mean_tn;    // time axis (copied from first valid model)
  nused = 0;
  for (i = 0; i < (int)models2.size(); i++) 
  {
    tmodel = models2[i];
  
    if (tmodel.data.rf.nrfo > 0) 
    {

      // allocate once from the first valid RF model
      if (npts < 0) 
      {
        npts = tmodel.data.rf.nrfo;

        mean_to.assign(npts, 0.0);
        mean_tn.assign(npts, 0.0);
        mean_rfn.assign(npts, 0.0);
        mean_rfo.assign(npts, 0.0);
        mean_unrfo.assign(npts, 0.0);
      }
      // skip inconsistent-length RFs
      if (tmodel.data.rf.nrfo != npts) continue;

      // hard safety: skip if any vector is missing/short (prevents segfault)
      if ((int)tmodel.data.rf.to.size()    < npts) continue;
      if ((int)tmodel.data.rf.rfn.size()   < npts) continue;
      if ((int)tmodel.data.rf.rfo.size()   < npts) continue;
      if ((int)tmodel.data.rf.unrfo.size() < npts) continue;
      for (j = 0; j < npts; j++) 
      {
        // time axes: store from first usable model (or overwrite—either is fine)
        mean_to[j] = tmodel.data.rf.to[j];
        // tn may not exist in some models -> protect it
        if ((int)tmodel.data.rf.tn.size() >= npts) 
        {
          mean_tn[j] = tmodel.data.rf.tn[j];
        }
        else
        {
          mean_tn[j] = tmodel.data.rf.to[j];
        }
        // accumulate sums (we will divide later)
        mean_rfn[j]   += tmodel.data.rf.rfn[j];
        mean_rfo[j]   += tmodel.data.rf.rfo[j];
        mean_unrfo[j] += tmodel.data.rf.unrfo[j];
      }
      ++nused;   // count models (not samples)
    }
  }
  // convert sums -> means
  for (j = 0; j < npts; j++) {
    mean_rfn[j]   /= (double)nused;
    mean_rfo[j]   /= (double)nused;
    mean_unrfo[j] /= (double)nused;
  }

  sprintf(minname, "%s/%s.mean.rf",outdir,wmodelname1);
  fprintf(stderr,"Open the mean file: %s\n",minname);
  fprintf(stderr,">> \n");
  ff0 = fopen(minname, "w");
  for (i = 0; i < npts; i++) 
  {
    // output: tn mean_rfn to mean_rfo mean_unrfo  (same format you used)
    fprintf(ff0, "%g %g %g %g %g\n",
            mean_tn[i], mean_rfn[i], mean_to[i], mean_rfo[i], mean_unrfo[i]);
  }
  fclose(ff0);
  // ==========================================================================
  sprintf(minname,"%s/%s.all.ph",outdir,wmodelname1);
  ff0 = fopen(minname,"w");
  for (i=0;i<models2.size();i++) {
    tmodel = models2[i];
    if (tmodel.data.disp.pper.size() > 0) {
      fprintf(ff0,"M %d",i);
      for (j=0;j<tmodel.data.disp.pper.size();j++) {
        fprintf(ff0," %g",tmodel.data.disp.pper[j]);
        }
      for (j=0;j<tmodel.data.disp.pper.size();j++) {
        fprintf(ff0," %g",tmodel.data.disp.pvel[j]);
        }
      } 
        fprintf(ff0,"\n");
    } 
  fclose(ff0);
  //  ===================================================================================
  // Do it again - dont know why if add the loop, can not write data into text file
  //  ===================================================================================
  npts = -1;
  std::vector<double> mean_pper,mean_pvel;   // mean of rfn[j]
  nused = 0;
  for (i = 0; i < (int)models2.size(); i++) 
  {
    tmodel = models2[i];
  
    if (tmodel.data.disp.pper.size() > 0) 
    {
      // allocate once from the first valid RF model
      if (npts < 0) 
      {
        npts = tmodel.data.disp.pper.size();
        mean_pper.assign(npts, 0.0);
        mean_pvel.assign(npts, 0.0);
      }
      // skip inconsistent-length RFs
      if (tmodel.data.disp.pper.size() != npts) continue;
      // hard safety: skip if any vector is missing/short (prevents segfault)
      if ((int)tmodel.data.disp.pper.size() < npts) continue;

      for (j = 0; j < npts; j++) 
      {
        // time axes: store from first usable model (or overwrite—either is fine)
        mean_pper[j] = tmodel.data.disp.pper[j];

        // accumulate sums (we will divide later)
        mean_pvel[j]   += tmodel.data.disp.pvel[j];
      }
      ++nused;   // count models (not samples)
    }
  }

  // convert sums -> means
  for (j = 0; j < npts; j++) 
  {
    mean_pvel[j]   /= (double)nused;
  }

  sprintf(minname, "%s/%s.mean.ph",outdir,wmodelname1);
  fprintf(stderr,"Open the mean file: %s\n",minname);
  fprintf(stderr,">> \n");
  ff0 = fopen(minname, "w");
  for (i = 0; i < npts; i++) 
  {
    fprintf(ff0, "%g %g\n",
            mean_pper[i],mean_pvel[i]);
  }
  fclose(ff0);
/*
  sprintf(minname,"%s/%s.all.gr",outdir,wmodelname1);
  ff0 = fopen(minname,"w");
  for (i=0;i<models2.size();i++) {
    tmodel = models2[i];
    if (tmodel.data.disp.gper.size() > 0) {
      fprintf(ff0,">> M %d \n",i);
      for (j=0;j<tmodel.data.disp.gper.size();j++) {
        fprintf(ff0,"%g %g %g %g\n",tmodel.data.disp.gper[j],tmodel.data.disp.gvel[j],tmodel.data.disp.gvelo[j],tmodel.data.disp.ungvelo[j]);
        }
      }
    }
  fclose(ff0);
  */
/////////////////////////////////////////////////////////////
// FORWARD MODEL THE VERAGE (previous seem incorrect)! //
/////////////////////////////////////////////////////////////
// models2.clear();

  tpara=ave_para;
  para2mod(tpara,model,tmodel, invtype);
  updatemodel(tmodel, invtype);
  compute_disp(tmodel,invtype);
  compute_rf(tmodel,invtype);
  // models2.push_back(tmodel);
  // ============================================================================================
  sprintf(minname,"%s/%s.final.rf",outdir,wmodelname1);
  ff0 = fopen(minname,"w");

  if (tmodel.data.rf.nrfo > 0) 
  {
    for (int i=0;i<tmodel.data.rf.nrfo;i++) 
    {
      // fprintf(ff0,"%g %g %g %g\n",tmodel.data.rf.to[j],tmodel.data.rf.rfn[j],tmodel.data.rf.rfo[j],tmodel.data.rf.unrfo[j]);
      fprintf(ff0,"%g %g %g %g %g\n",tmodel.data.rf.tn[i],tmodel.data.rf.rfn[i],tmodel.data.rf.to[i],tmodel.data.rf.rfo[i],tmodel.data.rf.unrfo[i]);
    }
  } 

  fclose(ff0);
  // ============================================================================================
  sprintf(minname,"%s/%s.final.e",outdir,wmodelname1);
  ff0 = fopen(minname,"w");
  for(int i=0;i<tmodel.data.disp.neper;i++)
  {
    if (tmodel.data.disp.eper.size() > 0) 
    {
      fprintf(ff0,"%g %g %g %g\n",tmodel.data.disp.eper[i],tmodel.data.disp.evel[i],tmodel.data.disp.evelo[i],tmodel.data.disp.unevelo[i]);                
    } 
  }
  fprintf(ff0,"\n");     
  fclose(ff0);
  // ============================================================================================
  sprintf(minname,"%s/%s.final.ph",outdir,wmodelname1);
  ff0 = fopen(minname,"w");

  for(i=0;i<tmodel.data.disp.npper;i++)	
  {
    if (tmodel.data.disp.pper.size() > 0) 
    {
      fprintf(ff0,"%g %g %g %g\n",tmodel.data.disp.pper[i],tmodel.data.disp.pvel[i],tmodel.data.disp.pvelo[i],tmodel.data.disp.unpvelo[i]);
    } 
  }
  fprintf(ff0,"\n");     
  fclose(ff0);
// ============================================================================================

/////////////////////////////////////////////////////////////
  fprintf(stderr,"write over check the finish condition!\n");
/////////////////////////////////////////////////////////////
// HVL add the signal for the script auto repeat if have no good result
  // if (paras.size()<1)
  // {
  //   cout<<"process finish but seem it not good result! \n"<<endl;
  //   cout<<"change the search range for kai_cri, kai_disp, kai_rf"<<endl;
  // }
  // else
  // {
    cout<<"post process finish!!! ready to plot!\n"<<endl;
    fprintf(stderr,"post process finish!!! ready to plot!\n");
  // }
  return 1;
}//main

