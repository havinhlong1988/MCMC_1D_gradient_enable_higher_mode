// get_new_model
// get_misfit
#include<fstream>
#include<vector>
#include<math.h>
#include<omp.h>
#include<algorithm>
#include<iostream>

//*************** get_new_model ******************************************
//  from a given parameter, you perturb it, and get a new para
//  new para is stored in para1
//  new model is stored in model1
//  if it is successful, return 1; otherwise 0.


int write_paras0(vector<paradef> paras, char *outname, char*outdir, int invtype) 
{

  // fprintf(stderr,"> MC.C: write_paras\n");
  // sleep(1);
  // Change the parameter if it use the layercake style [add the group ratio]
  FILE *ff;
  int i,j;
  char nameparas[300];
  sprintf(nameparas,"%s/%s",outdir,outname);
  fprintf(stderr,"para file name: %s \n",nameparas);
  ff=fopen(nameparas,"w");
  for (i=0;i<paras.size();i++)
  {
    // fprintf(stderr,"paras.npara? %d \n",paras[i].npara);
    // sleep(2);
    for (j=0;j<paras[i].npara;j++)
    {
      fprintf(ff,"%lf ",paras[i].parameter[j]);
    }
    // if (invtype==1) // layer cake or gradient
    // {
    //   for (int k=0;k<paras[i].nparav;k++)
    //   {
    //     fprintf(ff,"%lf ",paras[i].parameter[k+paras[i].npara]);
    //   }
    // }
    fprintf(ff,"\n");
    }
  fclose(ff);
  return 1;
  }

int get_ave_para(vector<paradef> paras, paradef &ave_para)
{

  // fprintf(stderr,"> MC.C: get_ave_para\n");
  // sleep(1);

  fprintf(stderr,"get_ave_para - paras.size()? %d \n",paras.size());
  if (paras.size()<1) return 0;
  if (paras.size()==1)
  {
    ave_para=paras[0]; 
    return 1;
    }
  
  int i,j,k;
  double tvalue;
  paradef tpara;
  ave_para=paras[0];
  int nvalid;
  for (i=0;i<ave_para.npara;i++)
  {
    tvalue=0.;
    nvalid=0;
    for (j=0;j<paras.size();j++)
    {
      // (optional safety) if some para vectors are shorter
      if (i >= (int)paras[j].parameter.size()) continue;

      double v = paras[j].parameter[i];
      if (!std::isfinite(v)) continue;   // skip NaN/Inf

      tvalue += v;
      nvalid++;

    }
    if (nvalid>0) ave_para.parameter[i] = tvalue / nvalid;
    else ave_para.parameter[i] = paras[0].parameter[i]; // keep a safe value
    }
  return 1;
  }

double get_ml(vector<double> values)
{

  // fprintf(stderr,"> MC.C: get_ml\n");
  // sleep(1);

  double tvalue;
  if (values.size()<1) return -999.;
  auto result = std::minmax_element(values.begin(), values.end());
  int NN = 50;
  int TN,fN,ffN;
  double step = (*result.second-*result.first+0.00001)/NN; 
  vector<int> Nums(NN,0);
  int i;
  cout<<step<<" "<<*result.first<<" "<<*result.second<<endl;
  for (i=0;i<values.size();i++)
  {
    tvalue = values[i];
    //cout<<tvalue<<endl;
    TN = int((tvalue-*result.first)/step);
    //cout<<tvalue<<" "<<TN<<endl;
    Nums[TN]++;
    }
  //auto result1 = std::minmax_element(Nums.begin(),Nums.end());
  //fN=*result1.second;
  int n_max=0;
  for (i=0;i<NN;i++)
  {
    if (Nums[i] > n_max)
    {
      ffN=i;
      n_max=Nums[i];
      //break;
      }
      }
  cout<<*result.first + ffN*step<<" "<<ffN<<" "<<step<<endl;
  return *result.first + (ffN)*step + 0.5*step;
  }

int get_ml_para(vector<paradef> paras, paradef &ml_para)
{

  // fprintf(stderr,"> MC.C: get_ml_para\n");
  // sleep(1);

  if (paras.size()<1) return 0;
  if (paras.size()==1)
  {
    ml_para=paras[0];
    return 1;
    }
  int i,j,k;
  double tvalue;
  paradef tpara;
  vector<double> values; 
  ml_para=paras[0];
  for (i=0;i<ml_para.npara;i++)
  {
    tvalue=0.;
    values.clear();
    for (j=0;j<paras.size();j++)
    {
      values.push_back(paras[j].parameter[i]);
      }
    auto result = std::minmax_element(values.begin(), values.end());
    /*
    for (j=0;j<paras.size();j++) {
      tvalue=tvalue+paras[j].parameter[i];
      }
    tvalue = tvalue/paras.size();
    */
    // tvalue is the most likely value.
    tvalue = get_ml(values);
    ml_para.parameter[i] = tvalue;
    }
  return 1;
  }

int get_new_model (paradef &para, paradef &para2, modeldef &model, modeldef &model1, int flag_newp, \
                  vector<int> vmono, vector<int> vgrad, int invtype)
  {

    // fprintf(stderr,"> MC.C: get_new_model\n");
    // sleep(1);

  // fprintf(stderr,"get_new_model outside\n");
  //initpara(para1);
  paradef para1;
  initpara(para1);
  int ibad=0;
  
  // fprintf(stderr,"check inversion style: %d\n", invtype);

  gen_newpara(para,para1,flag_newp,invtype);
  para2mod(para1,model,model1,invtype); // this work for b-spline only
  updatemodel(model1,invtype);
  int ok_m = 0;
  // while ((goodmodel(model1,vmono,vgrad,invtype))==0)
  while (!ok_m)
  {
    // cout<<"get_new_model inside"<<endl;
    // initpara(para1);
    gen_newpara(para,para1,flag_newp,invtype);
    para2mod(para1,model,model1,invtype);
    updatemodel(model1,invtype);
    // cout<<"finish updata model"<<endl;
    ok_m = goodmodel(model1,vmono,vgrad,invtype);

    ibad++;
    if (ibad%2000 == 1999 ) 
    {
      fprintf(stderr,"lots of ibad: %d\n",ibad);
      }
    if (ibad>10000) 
      {
      fprintf(stderr,"too many bad models here, why?\n"); 
      return 0;
      }
    }
  // fprintf(stderr, "goodmodel==1 (accepted) after %d bad tries for get_new_model\n", ibad);
  //cout<<"finish get new para/model "<<ibad<<endl;
  // fprintf(stderr," >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"); 
  // fprintf(stderr," >>>>> good model!\n"); 
  // fprintf(stderr," >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"); 
  para2 = para1;
  return 1;
  }

//************** start with a random model/para *********************
// input is a model, and para
// output is a model that stores disp/rf/misfit 
// flag_newp : 0: brand new para; 1 perturbed based on input para

int start_random(paradef &para, paradef &para1, modeldef &model, modeldef &model1, int flag_newp, \
vector<int> vmono, vector<int> vgrad, double pp, int id_thread, int invtype) 
{
  
  // fprintf(stderr,"> MC.C: start_random\n");
  // sleep(1);

  fprintf(stderr,"now do get_new_model! idthread: %d\n",id_thread);
  if (get_new_model(para,para1,model,model1,flag_newp,vmono,vgrad,invtype) == 0) 
  {
    printf("cannot get random_newmodel to start with!!\n");
    cout<<"cannot get random_newmodel??"<<endl;
    return 0;
    }
  {
    compute_rf(model1,invtype);
  }
  #pragma omp critical(f1)
  {
    compute_disp(model1,invtype);
  }
  compute_misfit(model1,pp);
  
  return 1;
  }

int get_new_para_misfit(paradef &para, paradef &para1, modeldef &model, int monoc, int flag_newp, \
vector<int> vmono, vector<int> vgrad, double pp, int flag_misfit, int invtype)
{

  // fprintf(stderr,"> MC.C: get_new_para_misfit\n");
  // sleep(1);

  modeldef ttmodel;
  paradef tpara;

  gen_newpara(para,tpara,flag_newp, invtype);
  para2mod(tpara,model,ttmodel, invtype);
  updatemodel(ttmodel, invtype);
  int ibad = 0;
  if (monoc==1)
  {
    int ok_p = 0;
    // while ((goodmodel(ttmodel,vmono,vgrad,invtype))==0)
    while (!ok_p)
    {
      gen_newpara(para,tpara,flag_newp, invtype);
      para2mod(tpara,model,ttmodel, invtype);
      updatemodel(ttmodel, invtype);
      ok_p = goodmodel(ttmodel,vmono,vgrad,invtype);
      ibad++;
      if (ibad%2000 == 1999 )
      {
        fprintf(stderr,"lots of ibad: %d\n",ibad);
        }
      if (ibad>10000)
      {
        fprintf(stderr,"too many bad models here, why?\n");
        return 0;
        }
      }
    }
    // fprintf(stderr, "goodmodel==1 (accepted) after %d bad tries\n", ibad);

  if (flag_misfit > 0)
  {
    // fprintf(stderr,"now get misfit !!!\n");
    // #pragma omp critical(forward) 
    {
      get_misfit1(tpara,ttmodel,pp,invtype);
      }
    // fprintf(stderr,"finish get misfit !!!\n");
    }
  model = ttmodel;
  para1 = tpara;
  return 1;
  } 

int write_paras(FILE *outf, int i1, int i2, int i3, paradef &para1, modeldef &ttmodel, double newL, double newmisfit, int id_thread) 
{

  // fprintf(stderr,"> MC.C: write_paras\n");
  // sleep(1);

  int j;
  if (outf == NULL )
  {
    fprintf(stderr,"the file is not open!!!\n");
    return 0;
    }
  fprintf(outf,"%d %d %d %d ",i1,id_thread,i2,i3);
  for(j=0;j<para1.npara;j++)fprintf(outf,"%g",para1.parameter[j]);
  fprintf(outf,"\n > %g %g %g %g %g %g %g %g %g\n",newL,newmisfit,ttmodel.data.rf.L,ttmodel.data.rf.misfit,ttmodel.data.disp.L,ttmodel.data.disp.misfit,ttmodel.data.disp.pmisfit,ttmodel.data.disp.gmisfit,ttmodel.data.disp.emisfit);
  return 1;
  }


int do_one_jump(FILE *outf, paradef &para0, modeldef &model0,  \
int monoc, vector<int> vmono, vector<int> vgrad, double pp, \
int jumpn,int &iaccp, int &iiter, char *wmodelname1, char *outdir, \
int flagw, int id_thread, vector<paradef> &paras, \
std::vector<boost::array<double, 10> > & infos, int invtype)
{
  // fprintf(stderr,"> MC.C: do_one_jump\n");
  // sleep(1);

  cout<<"iaccp: "<<iaccp<<endl;
  // fprintf(stderr,"iaccp: %d \n",iaccp);
  // sleep(3);
  int flag_misfit,flag_newp;
  int i,j,k,iac,id_name;
  double oldL,oldmisfit;
  double newL,newmisfit;
  double prob,prandom,tdif;
  char wmodelname2[300];
  time_t start,end;
  boost::array<double, 10> tinfo;

  paradef para,para1;
  para1 = para0;  
  para = para0;

  modeldef model,model1;
  model = model0;
  model1 = model0;

  /********************************************************************************/
  flag_newp = 0;
  //#pragma omp critical(write_ini) 
  {
    if (start_random(para,para1,model,model1,flag_newp,vmono,vgrad,pp,id_thread, invtype) ==0)
    {
      fprintf(stderr,"initialization failed when genrate random starting point!!!\n");
      return 1;
      }
    }
  //cout<<"get random start ok! "<<id_thread<<endl;
  /*******************************************************************************/
  oldL=model1.data.L;
  oldmisfit=model1.data.misfit;
  if (std::isnan(oldmisfit)) oldmisfit=999; //HVL add condition for -NaN value
  para=para1;
  
  //fprintf(stderr,"now ready to do a jump!!!\n");
  //abort();
  time( &start);

  // fprintf(stderr," -- jumpn : %d\n", jumpn);
  // sleep(3);

  for (i=0;i<jumpn;i++) {
     if (i%500 == 1) {fprintf(stderr,"%d th sampled !\n",i);}
     flag_newp = 1;
     flag_misfit = 1;
     //fprintf(stderr,"now get new para_misfit ! %d %d \n",i,id_thread);
     if (get_new_para_misfit(para, para1, model, monoc, flag_newp, vmono, vgrad, pp, flag_misfit, invtype) == 0) {
        fprintf(stderr,"cannot get new para with misfit! \n");
        continue;
        }
     //cout<<"get_new para ok! "<<i<<" for thread "<<id_thread<<endl;

     #pragma omp critical(write_para_p1)
      {
      iiter = iiter + 1;
      }

     newL=para1.L;
     newmisfit=para1.misfit;
     
     if (std::isnan(newmisfit)) newmisfit=999.; //HVL add condition for -NaN value
     prob=(oldL-newL)/oldL;
     //#pragma omp critical(write_para_pr)
     prandom=gen_random_unif01();
     //cout<<newL<<" "<<oldL<<endl;
     if (newL<oldL and prandom<prob) {
       //cout<<"not accpet : "<<i<<" th model! "<<id_thread<<endl;
      //  fprintf(stderr," %d th model for %d th thread not accepted: %lf %lf\n",i,id_thread,newmisfit,oldmisfit);
       iac=0;
       
       #pragma omp critical(write_para_p2) 
         {
         if (write_paras(outf,iac,iiter,iaccp,para1,model,newL,newmisfit,id_thread) == 0) {
           fprintf(stderr,"cannot write the unacc model to the file , why?\n");
           }
         }    
       //fprintf(stderr," %d th model for %d th thread not accepted: %lf %lf now continue\n",i,id_thread,newmisfit,oldmisfit);
       continue;
       }
     if (newmisfit>900){
         fprintf(stderr,"newmisfit > 900 | {iteration_id, id_name, id_thread: %d | %d |%d }\n",iiter,iaccp,id_thread);
         continue;
     }
     if (newmisfit<0){
         fprintf(stderr,"newmisfit < 0. | {iteration_id, id_name, id_thread: %d | %d |%d }\n",iiter,iaccp,id_thread);
         continue;
     }
     iac=1;
     #pragma omp critical(write_para_p3)
       {
       iaccp = iaccp + 1;
       id_name = iaccp;
       
       time ( &end);
       tdif = difftime(end,start);
       if (fmod(double(iaccp),100.)==1.)  // just print each 10 accepted model due to large amount of this
       {
        //  fprintf(stderr," + model accepted!!{iteration_id - id_name - id_thread, new_misfit, old_misfit, misfit_diff, time_calculation, model_name}:\n");
        //  fprintf(stderr,"   %d - %d - %d - %lf - %lf - %lf - %lf - %s\n",iiter,iaccp,id_thread,newmisfit,oldmisfit,oldmisfit-newmisfit,tdif,wmodelname1);
         fprintf(stderr," + model accepted!!{iteration_id, id_name, id_thread, new_misfit, old_misfit, misfit_diff}: %d | %d | %d | %lf | %lf | %lf | %s\n",iiter,iaccp,id_thread,newmisfit,oldmisfit,oldmisfit-newmisfit,wmodelname1);
       }
       if (write_paras(outf,iac,iiter,iaccp,para1,model,newL,newmisfit,id_thread) == 0 ) {
         fprintf(stderr,"cannot write the acc-ed model to the file , why?\n");
         }
       sprintf(wmodelname2,"%s.%d",wmodelname1,id_name);
       if (flagw != -2) {
         write_model(model,wmodelname2,outdir,flagw);
         }
       //write_bin(model,outbin,para1,iac,iiter,id_name,id_thread);
       paras.push_back(para1);
       tinfo[0] = double(iac);
       tinfo[1] = double(iiter);
       tinfo[2] = double(id_name);
       tinfo[3] = double(id_thread);
       tinfo[4] = model.data.L;
       tinfo[5] = model.data.misfit;
       tinfo[6] = model.data.rf.misfit;
       tinfo[7] = model.data.disp.misfit;
       tinfo[8] = para1.npara;
       tinfo[9] = model.ngroup;
       infos.push_back(tinfo);

       //added by EMB :(
       cout<<"test MC.C i= "<<i<<" iiter "<<iiter<<" misfit "<<tinfo[7]<<" id_thread "<<tinfo[3]<<" iac "<<iac<<endl;
      //  fprintf(stderr,"test MC.C i= %d + iiter: %d + misfit: %lf + id_thread: %lf + iac: %d\n",i,iiter,tinfo[7],tinfo[3],iac);
       }
     
     //cout<<"write acc model ok! "<<id_thread<<endl;
     para=para1;
     oldL=newL;
     oldmisfit=newmisfit;
     }
  fprintf(stderr,"finish the jump! id_thread = %d \n",id_thread);
  return 1;
  }


int do_prior(FILE *outf, paradef &para0, modeldef &model0,  \
            int monoc, vector<int> vmono, vector<int> vgrad, double pp, \
            int jumpn,int iaccp, int iiter, char *wmodelname1, char *outdir, \
            int flagw, vector<paradef> &paras, int invtype) {
  
  // fprintf(stderr,"> MC.C: do_prior\n");
  // sleep(1);

  int flag_misfit,flag_newp;
  int i,j,k,iac;
  int id_thread = 0;
  double oldL,oldmisfit;
  double newL,newmisfit;
  double prob,prandom;
  char wmodelname2[300];
  paradef para,para1;
  para1 = para0;
  para = para0;

  modeldef model,model1;
  model = model0;
  model1 = model0;

  flag_newp = 0;
  if (start_random(para,para1,model,model1,flag_newp,vmono,vgrad,pp,id_thread,invtype) ==0) {
     fprintf(stderr,"initialization failed when genrate random starting point!!!\n");
     return 0;
     }
  oldL=model1.data.L;
  oldmisfit=model1.data.misfit;
  para=para1;

  fprintf(stderr,"now ready to do a jump!!!\n");
  for (i=0;i<jumpn;i++) {
     if (fmod(i,500.)==1.) {fprintf(stderr,"finish %dth sampling!\n",i);}
     flag_newp = 1;
     flag_misfit = 0;
     if (get_new_para_misfit(para, para1, model, monoc, flag_newp, vmono, vgrad, pp, flag_misfit,invtype) == 0) {
        fprintf(stderr,"cannot get new para with misfit! \n");
        return 0;
        }
     #pragma omp critical(updateflag) 
       {
       iiter++;
       iaccp ++;
       }
     
     fprintf(stderr," - model accepted!!{iteration_id, id_name, id_thread, new_misfit, old_misfit, misfit_diff}: %d | %d | %d | %lf | %lf | %lf\n",iiter,iaccp,id_thread,newmisfit,oldmisfit,oldmisfit-newmisfit);
     newmisfit = 0;
     newL = 1;
     iac = 1;
     #pragma omp critical(write_model_1) 
       {
       if (write_paras(outf,iac,iiter,iaccp,para1,model,newL,newmisfit,id_thread) == 0) {
         fprintf(stderr,"cannot write the acc-ed model to the file , why?\n");
         }
       }
     /*
     if(not outbin.is_open()) {
         printf("cannot write the acc-ed model to the BINary file , why?\n");
         return 0;
         }
     */
     sprintf(wmodelname2,"%s.%d",wmodelname1,iaccp);
     // RACE FIX: this push_back must take the SAME critical as the other
     // paras.push_back() above (write_para_p3). OpenMP named criticals are
     // independent mutexes, so guarding the two writes to the shared
     // std::vector `paras` under DIFFERENT names allowed two threads to push
     // concurrently -> vector reallocation race -> heap corruption ->
     // intermittent "Bus error: 10" / "Segmentation fault: 11" inside
     // pthread_mutex_lock. No deadlock: the other site takes only this lock;
     // this site takes updateflag -> write_model_1 -> write_para_p3 (no cycle).
     #pragma omp critical(write_para_p3) 
       {
       //write_bin(model,outbin,para1,iac,iiter,iaccp,0);
       paras.push_back(para1);
       }
     para=para1;
     }
  return 1;
  }

double get_v_tb1(modeldef tmodel, int n_group, double *v1, double *v2) {

  // fprintf(stderr,"> MC.C: get_v_tb1\n");
  // sleep(1);

  *v1 = -999.;
  *v2 = -999.;
  int i,j,k,n;
  double d1,d2,d3,d4,dd,tdep;
  d1 = 0.;
  d2 = 0.;
  for (i=0;i<=n_group;i++) {
    d2 = d1;
    d1 = d1 + tmodel.groups[i].thick;
    }
  tdep = d1 - 5.;
  dd = tmodel.groups[n_group].thick/tmodel.groups[n_group].nlay;
  for (i=0;i<tmodel.groups[n_group].nlay-1;i++) {
    d3 = d2+(i+0.5)*dd;
    d4 = d2+(i+1.5)*dd;
    if (tdep > d3 && tdep < d4+0.00000001) {
       *v1 = tmodel.groups[n_group].value1[i] + (tmodel.groups[n_group].value1[i+1] - tmodel.groups[n_group].value1[i])*(tdep-d3)/(d4-d3);
       }
    }
  
  
  d1 = 0.;
  d2 = 0.;
  for (i=0;i<=n_group+1;i++) {
    d2 = d1;
    d1 = d1 + tmodel.groups[i].thick;
    }
  tdep = d2 + 5.;
  dd = tmodel.groups[n_group+1].thick/tmodel.groups[n_group+1].nlay;
  for (i=0;i<tmodel.groups[n_group+1].nlay-1;i++) {
    d3 = d2+(i+0.5)*dd;
    d4 = d2+(i+1.5)*dd;
    if (tdep > d3 && tdep < d4+0.00000001) {
       *v2 = tmodel.groups[n_group+1].value1[i] + (tmodel.groups[n_group+1].value1[i+1] - tmodel.groups[n_group+1].value1[i])*(tdep-d3)/(d4-d3);
       }
    }
  return 0.;
  }


double get_v_tb(modeldef tmodel, int n_group, double *v1, double *v2) {

  // fprintf(stderr,"> MC.C: get_v_tb\n");
  // sleep(1);

  //if (tdep < 0.) return -1.;
  double d1,d2,tdep;
  double tv0;
  int i,j,k,n;
  d1 = 0.;
  for (i=0;i<tmodel.ngroup;i++)  d1 = d1 + tmodel.groups[i].thick;
  //if (tdep > d1) return -1.;

  //
  d1 = 0.;
  d2 = 0.;
  for (i=0;i<=n_group;i++) {
    d2 = d1;
    d1 = d1 + tmodel.groups[i].thick;
    }
  tv0 = 0.;
  n = 0;
  for (k=0;k<tmodel.groups[n_group].nlay;k++) {
    tdep = d2 + ((float)(k+0.5)/(float)(tmodel.groups[n_group].nlay))*(d1-d2);
    if (tdep > d1-4.) {
       //fprintf(stderr,"%g %g %g %g\n",d1,d2,tdep,tmodel.groups[n_group].value1[k]);
       tv0 = tv0 + tmodel.groups[n_group].value1[k];
       n = n+1;
       }
    }
  *v1 = tv0/n;
  //fprintf(stderr,"%g %d\n",tv0/n,n);
  //
  d1 = 0.;
  d2 = 0.;
  for (i=0;i<=n_group+1;i++) {
    d2 = d1;
    d1 = d1 + tmodel.groups[i].thick;
    }
  tv0 = 0.;
  n = 0;
  for (k=0;k<tmodel.groups[n_group+1].nlay;k++) {
    tdep = d2 + ((float)(k+0.5)/(float)(tmodel.groups[n_group+1].nlay))*(d1-d2);
    if (tdep < d2+4.) {
       tv0 = tv0 + tmodel.groups[n_group+1].value1[k];
       n = n+1;
       }
    }
  *v2 = tv0/n;
  //abort();
  return 0.;
  }


double get_v_1(modeldef tmodel,double tdep, int n_group) {

  // fprintf(stderr,"> MC.C: get_v_1\n");
  // sleep(1);

  if (tdep < 0.) return -1.;
  double d1,d2;
  int i,j,k;
  d1 = 0.;
  for (i=0;i<tmodel.ngroup;i++)  d1 = d1 + tmodel.groups[i].thick;
  if (tdep > d1) return -1.;

  d1 = 0.;
  d2 = 0.;
  for (i=0;i<=n_group;i++) {
    d2 = d1;
    d1 = d1 + tmodel.groups[i].thick;
    }
  //cout<<d2<<" "<<d1<<endl;
  if (tdep >= d2 and tdep < d1){
    k = int((tdep-d2)/(d1-d2) * (tmodel.groups[n_group].nlay));
    //cout<<k<<" "<<(tdep-d2)/(d1-d2)<<" "<<tmodel.groups[n_group].value1.size()<<" "<<tmodel.groups[n_group].value1[k]<<" "<<(tdep-d2)/(d1-d2) * (tmodel.groups[n_group].nlay+1)<<endl;
    return tmodel.groups[n_group].value1[k];
    }
  return -1.;
  }

double get_v(modeldef tmodel,double tdep) 
{

  // fprintf(stderr,"> MC.C: get_v\n");
  // sleep(1);

  if (tdep < 0.) return -1.;
  double d1,d2;
  int i,j,k;
  d1 = 0.;
  for (i=0;i<tmodel.ngroup;i++)  d1 = d1 + tmodel.groups[i].thick;
  if (tdep > d1) return -1.;
  d1 = 0.;
  d2 = 0.;
  for (i=0;i<tmodel.ngroup;i++) 
  {
   d1 = d1 + tmodel.groups[i].thick;
   if (tdep >= d2 and tdep < d1)
   {
     k = int((tdep-d2)/(d1-d2) * (tmodel.groups[i].nlay));
     return tmodel.groups[i].value1[k];
     }
   d2 = d1;
   }
  return -1.;
}

double get_q_flag(modeldef tmodel,double tdep) 
{

  // fprintf(stderr,"> MC.C: get_v\n");
  // sleep(1);

  if (tdep < 0.) return -1.;
  double d1,d2;
  int i,j;
  d1 = 0.;
  for (i=0;i<tmodel.ngroup;i++)  d1 = d1 + tmodel.groups[i].thick;
  if (tdep > d1) return -1.;
  d1 = 0.;
  d2 = 0.;
  for (i=0;i<tmodel.ngroup;i++) 
  {
   d1 = d1 + tmodel.groups[i].thick;
   if (tdep >= d2 and tdep < d1)
   {
     return tmodel.groups[i].q_flag;
     }
   d2 = d1;
   }
  return -1.;
}

//Added getout_vp and getout_rho by EMB to get posterior distribution
double getout_vp(modeldef tmodel,double tdep) {

  // fprintf(stderr,"> MC.C: getout_vp\n");
  // sleep(1);

  if (tdep < 0.) return -1.;
  double d1,d2;
  int i,j,k;
  d1 = 0.;
  for (i=0;i<tmodel.ngroup;i++)  d1 = d1 + tmodel.groups[i].thick;
  if (tdep > d1) return -1.;
  d1 = 0.;
  d2 = 0.;
  for (i=0;i<tmodel.ngroup;i++) {
   d1 = d1 + tmodel.groups[i].thick;
   if (tdep >= d2 and tdep < d1){
     k = int((tdep-d2)/(d1-d2) * (tmodel.groups[i].nlay));
     return tmodel.groups[i].value1p[k];
     }
   d2 = d1;
   }
  return -1.;
  }


double getout_rho(modeldef tmodel,double tdep) {

  // fprintf(stderr,"> MC.C: getout_rho\n");
  // sleep(1);

  if (tdep < 0.) return -1.;
  double d1,d2;
  int i,j,k;
  d1 = 0.;
  for (i=0;i<tmodel.ngroup;i++)  d1 = d1 + tmodel.groups[i].thick;
  if (tdep > d1) return -1.;
  d1 = 0.;
  d2 = 0.;
  for (i=0;i<tmodel.ngroup;i++) {
   d1 = d1 + tmodel.groups[i].thick;
   if (tdep >= d2 and tdep < d1){
     k = int((tdep-d2)/(d1-d2) * (tmodel.groups[i].nlay));
     return tmodel.groups[i].value1r[k];
     }
   d2 = d1;
   }
  return -1.;
  }


double cv_max(vector<double> v1) {

  // fprintf(stderr,"> MC.C: cv_max\n");
  // sleep(1);

  int i;
  double max_v;
  max_v = v1[0];
  for (i=0;i<v1.size();i++) {
    if (max_v < v1[i]) max_v = v1[i];
    }
  return max_v;
  }

double cv_min(vector<double> v1) {
  int i;
  double min_v;
  min_v = v1[0];
  for (i=0;i<v1.size();i++) {
    if (min_v > v1[i]) min_v = v1[i];
    }
  return min_v;
  }


double cv_mean(vector<double> v1) {

  // fprintf(stderr,"> MC.C: cv_mean\n");
  // sleep(1);

  double sum_of_elems=0.;
  for(std::vector<double>::iterator j=v1.begin();j!=v1.end();++j)
     sum_of_elems += *j;
  return sum_of_elems/v1.size();
  }

double cv_std(vector<double> v1) {

  // fprintf(stderr,"> MC.C: cv_std\n");
  // sleep(1);

  double v_mean = cv_mean(v1);
  double var=0.;;
  for(std::vector<double>::iterator j=v1.begin();j!=v1.end();++j)
     var += (*j-v_mean)*(*j-v_mean);
  return sqrt(var/v1.size());
  }

