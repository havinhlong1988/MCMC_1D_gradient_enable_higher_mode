/*#include<iostream>
#include<vector>
#include<fstream>
#include<algorithm>
#include"INITstructure.h"
#include"/home/jiayi/progs/jy/HEAD/head_c++/gen_random.C"
*/
#include<iostream>
using namespace std;
/*========CONTENT============
int initpara(paradef &para)
int readpara(paradef &para, const char* fname)
int gen_newpara ( paradef &inpara, paradef &outpara,int pflag)
int check_para(paradef &inpara, int invtype) 
int mod2para(modeldef &model, paradef &inpara, paradef &outpara)
int getminmaxpara(paradef &inpara, paradef &paraleft, paradef &pararight)
int para2mod(paradef &para, modeldef &inmodel, modeldef &outmodel)
int get_misfit(paradef &para, modeldef &inmodel,modeldef &outmodel,double p, int nn)
//============================
*/
//class paracal{
//public:
	int initpara(paradef &para)
	{
	
	// fprintf(stderr,"> CALpara: initpara\n");
    // sleep(1);;

	  para.npara=0;
	  para.nparav=0;
	  para.npara_tt=0; // total thickness pertubation (sediment and crust)
	  para.npara_tr=0; // sublayer thickness ratio pertubation
	  para.L=0.;
	  para.misfit=0.;
	  para.flag=0;
	  return 1;
	}
//-----------------------------------------------------	 
	int readpara(paradef &para, const char* fname)
	{  
	/*
	Read the in.paras_STA file
	*/
	// fprintf(stderr,"> CALpara: readpara\n");
    // sleep(1);;

	  fstream mff;
	  int k,i;
	  string line;
	  vector<string> v;
	  mff.open(fname);
	  if ( not mff.is_open()){cout<<"########para file "<<fname<<" does not exist!\n exit!!\n";exit (0);}
	  k=0;
	  int nparav=0;
	  int npara_tt=0;
	  int npara_tr=0;
	  para.para0.clear();
	  while(getline(mff,line))
	  {
		v.clear();
		Split(line,v," ");
		// for (size_t ii = 0; ii < v.size(); ++ii) 
		// {
		// 	fprintf(stderr, " [%zu]=%s", ii, v[ii].c_str());
  		// }
  		// fprintf(stderr, "\n");
		if (std::atoi(v[0].c_str()) == 0) 
		{
			nparav++;
		}
		if (std::atoi(v[0].c_str()) == 1) // thickness pertubation 
		{
			if (v.size() == 5) // line total thickness pertubation have 5 paras
			{
				npara_tt++;
			}
			else if (v.size() > 5) // line of ratio pertubation have 6 paras
			{
				npara_tr++;
			}
		}
		k=k+1;
		para.para0.push_back(v);
	  }
	  mff.close();
	  para.npara=k;
	  para.nparav=nparav; // since layercake model give number depth-ratio = number of thickness
	  para.npara_tt=npara_tt; // use later to arrange the parameter
	  para.npara_tr=npara_tr; // use later to perturbe the sublayer thickness
	  fprintf(stderr, " readpara > npara: %d | npara velocity: %d | npara total thicness: %d | npara thickness ratio: %d\n", k,nparav,npara_tt,npara_tr);
	  cout<<"finish reading para! npara="<<para.npara<<endl;
	  return 1;
	}	

//-----------------------------------------------------	 
	int perturb_thickness_ratios_weighted(const std::vector<double>& rin,
										std::vector<double>& rout,
										const std::vector<double>& w_in,
										double eps = 1e-10)
	{
		const int n = (int)rin.size();
		if (n <= 0) return 0;
		if ((int)w_in.size() != n) return 0;

		// Draw random simplex u ~ Dirichlet(1,...,1)
		std::vector<double> u(n, 0.0);
		double su = 0.0;
		for (int i = 0; i < n; ++i) {
			double r = gen_random_unif01();
			if (r < 1e-15) r = 1e-15;
			u[i] = -std::log(r);   // Exp(1)
			su += u[i];
		}
		if (su <= 0.0) return 0;
		for (int i = 0; i < n; ++i) u[i] /= su;

		// Blend per-component with weight w[i] in [0,1]
		rout.assign(n, 0.0);
		double s = 0.0;
		for (int i = 0; i < n; ++i) {
			double a = rin[i];
			if (a < 0.0) a = 0.0;

			double wi = w_in[i];
			if (wi < 0.0) wi = 0.0;
			if (wi > 1.0) wi = 1.0;

			rout[i] = (1.0 - wi) * a + wi * u[i];
			if (rout[i] < eps) rout[i] = eps;
			s += rout[i];
		}

		if (s <= 0.0) return 0;
		for (int i = 0; i < n; ++i) rout[i] /= s;
		return 1;
	}

	int perturb_thickness_ratios_inplace_weighted(std::vector<double>& params,
												int offset, int nr,
												const std::vector<double>& w)
	{
		if (nr <= 0) return 0;
		std::vector<double> rin(nr), rout;
		for (int i = 0; i < nr; ++i) rin[i] = params[offset + i];

		if (!perturb_thickness_ratios_weighted(rin, rout, w)) return 0;

		for (int i = 0; i < nr; ++i) params[offset + i] = rout[i];
		return 1;
	}



// randomly generating new parameters, uniform or normal distribution
	int gen_newpara ( paradef &inpara, paradef &outpara,int pflag, int invtype)
	{

	// fprintf(stderr,"> CALpara: gen_newpara\n");
    // sleep(1);;

	  int i,flag,ibad;
	  double newv,sigma,mean,cvr;
	  outpara=inpara; 
	  int offset = (inpara.nparav+inpara.npara_tt);  // ratios start here in your code
	  int nr = (inpara.npara-offset);                // number of ratios in this group
	  if (offset < 0 || nr < 0 || offset + nr > inpara.npara) 
	  {
        fprintf(stderr, "bad ratio indexing: npara=%d offset=%d nr=%d\n",
                inpara.npara, offset, nr);
        return 0;
      }  
      const int n_main = offset; // velocity + total thickness part
	  
	//   fprintf(stderr,"parameter check!! %d %d %d %d\n",inpara.npara,offset,nr,n_main);
	  
	  if(pflag==0) // uniform genration when old model remain forward
	  {
		for(i=0;i<n_main;i++) // for velocity only
			{
			  cvr = gen_random_unif01(); // Crusial point to make a new model based on a random value
			  newv=cvr*(inpara.space1[i][1]-inpara.space1[i][0])+inpara.space1[i][0];//gen_random_unif01 => [0,1)
			//   fprintf(stderr,"check genpara!! cv: %g | val: %g | min: %g | max: %g\n",cvr,newv,inpara.space1[i][0],inpara.space1[i][1]);
			  if (newv > inpara.space1[i][1] || newv < inpara.space1[i][0])
			  	{fprintf(stderr,"something is wrong here!! %g %g %g %g\n",cvr,newv,inpara.space1[i][0],inpara.space1[i][1]);
				return 0;
				};
				outpara.parameter[i]=newv;  
			}
		// -------- Ratios: simplex, per-component random cvr --------
        if (nr > 1) 
		{
            std::vector<double> w(nr);
            for (int k = 0; k < nr; ++k) w[k] = gen_random_unif01(); // per ratio
            if (!perturb_thickness_ratios_inplace_weighted(outpara.parameter, offset, nr, w))
                return 0;
        } 
		else if (nr == 1) 
		{
        	outpara.parameter[offset] = 1.0; // only one ratio => must be 1
        }
	  }//if		  	  
	  else //normal distribution genration when new model need to forward
	  {
		for(i=0;i<n_main;i++)
		{
                //cout<<"which edit? "<<inpara.parameter[i]<<endl; //EMB
		  flag=0;
		  mean=inpara.parameter[i];
		  sigma=inpara.space1[i][2];
		  while(flag<1)
		  {
			newv=gen_random_normal(mean,sigma); // normal distribution
			ibad++;
			if(newv>inpara.space1[i][1] || newv<inpara.space1[i][0]) {
			    if(ibad%5000==4999) 
				{
					fprintf(stderr,"cannot find good random ! %g %g %g %g %g %g %g %d %d\n",inpara.parameter[i],mean,sigma,newv,inpara.space1[i][0],inpara.space1[i][1],inpara.parameter[i],i,ibad); 
					// abort();
					exit(1);
					fprintf(stderr,"??? Program not terminated ???");
					sleep(2);
					} 
			    if(ibad > 9997 && ibad <10002) { 
					fprintf(stderr,"wrong here ,needs check!!! %g %g %g %g\n",sigma,mean,inpara.space1[i][0],inpara.space1[i][1]); sigma = sigma * 2; if (mean > inpara.space1[i][1]) mean = inpara.space1[i][1]; if (mean < inpara.space1[i][0]) mean = inpara.space1[i][0];
					abort();
					}
			    continue;}
			else {ibad = 0;flag=2;}
		  }//while flag
		  outpara.parameter[i]=newv;
		}//fori
		// -------- Ratios: still simplex (NOT normal) --------
        if (nr > 1) 
		{
			std::vector<double> w(nr);
            for (int k = 0; k < nr; ++k) w[k] = gen_random_unif01(); // per ratio

            if (!perturb_thickness_ratios_inplace_weighted(outpara.parameter, offset, nr, w))
                return 0;
        } 
		else if (nr == 1) 
		{
            outpara.parameter[offset] = 1.0;
        }
	  }
	// Print to check
	double total_tr = 0.0;

	// // print all params
	// for (size_t i = 0; i < outpara.parameter.size(); ++i) {
	// 	fprintf(stderr, "outpara.parameter[%zu] = %.8f\n", i, outpara.parameter[i]);
	// }

	// // sum only ratios: [offset, offset+nr)
	// for (int k = 0; k < nr; ++k) {
	// 	total_tr += outpara.parameter[offset + k];
	// }

	// fprintf(stderr, "ratio block sum (offset=%d, nr=%d) = %.12f\n", offset, nr, total_tr);

	  return 1;
	}//gen new para
	
// ----------------------------------------------------
	int check_para(paradef &inpara, int invtype) 
	{

	// fprintf(stderr,"> CALpara: check_para\n");
    // sleep(1);;

	int i=0;
	for (i=0;i<inpara.npara;i++) 
	  {
	  if (inpara.parameter[i] > inpara.space1[i][1] || inpara.parameter[i] < inpara.space1[i][0])
	    {
	    fprintf(stderr,"this para is wrong ! %d %g %g %g\n",i,inpara.parameter[i],inpara.space1[i][0],inpara.space1[i][1]);
	    return 0;
	    }
	  }
	return 1;
	}
//-----------------------------------------------------	 
	int mod2para(modeldef &model, paradef &inpara, paradef &outpara, int invtype)
	{

	// fprintf(stderr,"> CALpara: mod2para\n");
    // sleep(1);;
	// ========================================================
	// Convert the model from mod.{STA} to the parameter by combining with in.para_{STA}
	// ========================================================
	  int i,ng,nv;
	  double tv,tmin,tmax,sigma;
	  int p0,p1,p3,p4,p5;
	  double p2;
	  vector<double> vtr(3,0.);
	  std::vector<double> tv_ratio;
	  outpara=inpara;
	  outpara.parameter.clear();
	  outpara.L=model.data.L;
	  outpara.misfit=model.data.misfit;
	  //cout<<"mod2para .. inpara.npara? "<<inpara.npara<<endl;
      for(i=0;i<inpara.npara;i++)
	  {
		// for (size_t j = 0; j < inpara.para0[i].size(); ++j) {
		// 	fprintf(stderr, "%s%s", inpara.para0[i][j].c_str(),(j + 1 == inpara.para0[i].size()) ? "" : " ");
		// }
		// fprintf(stderr, "\n");
	  	p0=atoi(inpara.para0[i][0].c_str()); // id_flag
		p1=atoi(inpara.para0[i][1].c_str()); // perturb_flag
		p2=atof(inpara.para0[i][2].c_str()); // perturb range
		ng=atoi(inpara.para0[i][4].c_str()); // layer id

		// fprintf(stderr,"layer id check point : %d !\n",ng);
		if (p0==0) //value. vs or Bspline coeff
		{
			p5=atoi(inpara.para0[i][5].c_str()); // sub-layer id
			nv = p5;
			tv=model.groups[ng].value[nv];
			// if (invtype == 1) 
			// {
			// 	tv_ratio.push_back(model.groups[ng].ratio[nv]);
			// }
		}
		else if(p0==1)// thickness of group
		{
			if (inpara.para0[i].size() == 5) // total thickness pertubation have 5 paras
			{
				tv=model.groups[ng].thick;
				// fprintf(stderr,"model.groups[ng].thick check point 1: %d %f !\n",ng,tv);
				// sleep(2);
			}
			else if (inpara.para0[i].size() > 5) // now thickness ratio perturbation
			{
				nv=atoi(inpara.para0[i][5].c_str()); // sub-layer id
				tv=model.groups[ng].ratio[nv];
			}
		}
		else if(p0==-1)// vpvs
		{
			tv=model.groups[ng].vpvs1;
		}
        else if (p0==-2)
		{
            tv=model.groups[ng].vpvs;
		}
        else if (p0==-3)
		{
            tv=model.groups[ng].drho1;
		}
        else if (p0 == -4)
		{
            tv=model.groups[ng].drho;
		}
        else if (p0 == -5)
		{
            tv=model.groups[ng].dvs1;
		}
    	else if (p0 == -6)
		{
			tv=model.groups[ng].dvs;
		}
    	else if (p0 == -7) 
		{
            tv=model.groups[ng].fdvs;                        
        }
        else if (p0 == -8) 
		{
			tv=model.groups[ng].fdvs1;
        }
		else
			{cout<<"#####Wrong!!! mode2para\n";
			exit (0);}

		outpara.parameter.push_back(tv); // pushback parameter for other model
		
		// fprintf(stderr,"paracheck point 2: %d %d %f %d %d!\n",p0,p1,p2,p3,p4,p5);
		// sleep(1);
//		cout<<outpara.parameter[i]<<endl;
		if(inpara.flag<1) //initial flag;
		{
				// fprintf(stderr,"inpara.flag check point 2: %d !\n",inpara.flag);
				// sleep(1);
			sigma=atof(inpara.para0[i][3].c_str());
		  	if(p1==1) // asolute pertubation
			{
			  tmin=tv-p2;
			  tmax=tv+p2;
			}
		  else
		  // percentage pertubation
			{
			  tmin=tv-tv*p2/100.;
			  tmax=tv+tv*p2/100.;
			}
            if (p0 != -4 && p0 != -3  && p0 != -5 && p0 != -6 &&(p0 != -2 && p0 != -1)) 
			{ // for velocity, vpvs, and thickness
		     tmin=max(0.,tmin);
		     tmax=max(tmin+0.0001,tmax);
            }
            if (p0 == -7 or p0 == -8) 
			{
				tmax = min(1.,tmax);
				tmin = min(tmax-0.0001,tmin);

                tmin=max(0.,tmin);
                tmax=max(tmin+0.0001,tmax);
            }
			if(p0 == 0 and p0+i+p5==0)
			{
				tmin=max(0.2,tmin); tmax=max(0.5,tmax);
			} 
		  	vtr[0]=tmin; vtr[1]=tmax; vtr[2]=sigma;
		  	// fprintf(stderr,"space now get: %d - %g - %g - %g - %g\n",i,tmin,(tmin+tmax)/2,tmax,sigma);
			// sleep(1);
		  	outpara.space1.push_back(vtr);		
		  	outpara.flag=1;
			}//if flag
	  }//for
	//   for (size_t ir=0; ir<tv_ratio.size(); ir++) outpara.parameter.push_back(tv_ratio[ir]);
//         cout<<outpara.parameter[0]<<endl;
	  // Print and check all parameter values
	//   for (size_t i = 0; i < outpara.parameter.size(); ++i) 
	//   {
	// 	fprintf(stderr, "outpara.parameter[%zu] = %.8f\n", i, outpara.parameter[i]);
	//   }
	 return 1;
	}//mod2para

	int getminmaxpara(paradef &inpara, paradef &paraleft, paradef &pararight)
	{	
		//------------------------------------------------------------------
		// get the minimum and maximum of the parameter based on the space
		// HVL added 20251112
		//------------------------------------------------------------------
		paraleft=inpara;
		pararight=inpara;
	  	/*
		For a new setting of parameter the parameters was arranged as:
		 1. velocity pertubation (para.nparav) 
		 2. total thickness variation (para.npara_tt)
		 3. thickness ratio variation (para.npara_tr)
		 So we can find the min max of each based on order of 1 - 2 - 3
		*/
		// fprintf(stderr, " getminmaxpara > npara: %d | npara velocity: %d | npara total thicness: %d | npara thickness ratio: %d\n",inpara.npara,inpara.nparav,inpara.npara_tt,inpara.npara_tr);
		for(int i=0;i<inpara.npara;i++)
	  	{
			// fprintf(stderr," parameter: %zu %g %g %g\n",i,inpara.parameter[i],inpara.space1[i][0],inpara.space1[i][1]);
			 if ( (i >= inpara.nparav) && (i < inpara.nparav + inpara.npara_tt) )
			{
				// Total thickness variation controls boundary depth.
				// For model-space bounds:
				// left model  = maximum depth / deepest boundary
				// right model = minimum depth / shallowest boundary
				paraleft.parameter[i] = inpara.space1[i][1];
				pararight.parameter[i] = inpara.space1[i][0];

				paraleft.space1[i][0] = inpara.space1[i][1];
				paraleft.space1[i][1] = inpara.space1[i][1];

				pararight.space1[i][0] = inpara.space1[i][0];
				pararight.space1[i][1] = inpara.space1[i][0];
				// fprintf(stderr," >> parameter: %zu %g %g %g\n",i,paraleft.parameter[i],paraleft.space1[i][0],paraleft.space1[i][1]);
				// fprintf(stderr," >> parameter: %zu %g %g %g\n",i,pararight.parameter[i],pararight.space1[i][0],pararight.space1[i][1]);
			}
			else if (i >= inpara.nparav + inpara.npara_tt) // thickness ratio variation
			{
				paraleft.parameter[i] = inpara.parameter[i];
				pararight.parameter[i] = inpara.parameter[i];

				paraleft.space1[i][0] = inpara.parameter[i];
				paraleft.space1[i][1] = inpara.parameter[i];
				
				pararight.space1[i][0] = inpara.parameter[i];
				pararight.space1[i][1] = inpara.parameter[i];

			}
			else // All other values (velocity)
			{
				paraleft.parameter[i] = inpara.space1[i][0];
				pararight.parameter[i] = inpara.space1[i][1];

				paraleft.space1[i][0] = inpara.space1[i][0];
				paraleft.space1[i][1] = inpara.space1[i][0];
				
				pararight.space1[i][0] = inpara.space1[i][1];
				pararight.space1[i][1] = inpara.space1[i][1];
			}


	  	} // end for

	  return 1;
	}
	
//-----------------------------------------------------	 
	int para2mod(paradef &para, modeldef &inmodel, modeldef &outmodel, int invtype)
	{

	// fprintf(stderr,"> CALpara: para2mod\n");
    // sleep(1);
	
	  double newv,tthick0=0.;
	  int ng,nv,p0;
	  int kflag = 0; // flag for group thickness pertubation
	  outmodel=inmodel;
	  for(int i=0;i<para.npara;i++)
	  {
	  	p0=atoi(para.para0[i][0].c_str()); // Line 1
		newv=para.parameter[i]; // center value coressponding with para0
	 	ng=atoi(para.para0[i][4].c_str()); // Group ID
        //cout<<":( p0: "<<p0<<" newv: "<<newv<<" ng: "<<ng<<" i: "<<i<<endl;
        if ( p0==0) // velocity
		{
			// fprintf(stderr,"para2mod > vel paras? %d | %g ",i,para.parameter[i]);
			// for (size_t j = 0; j < para.para0[i].size(); ++j) {
			// fprintf(stderr, " %s", para.para0[i][j].c_str());
			// }
			// fprintf(stderr, "\n");

		  	nv=atoi(para.para0[i][5].c_str()); // sublayer id

          	//cout<<":( check nv: "<<nv<<endl;
		  	outmodel.groups[ng].value[nv]=newv;

		  	// fprintf(stderr,"\npara2mod > add value v: %d %d %d | %g \n",i,ng,nv,outmodel.groups[ng].value[nv]);
			
          //EMB, average in between crustal params
          // ng==1,crust ; nv=0,2,4,6,8 crustparams
         // if ( (ng==1) && (nv>=2) && (nv<8) ) {
           //   outmodel.groups[ng].value[nv-1]=(0.5*((outmodel.groups[ng].value[nv-2])+(outmodel.groups[ng].value[nv]))); //average values 1, 3, 5, 7
          //}
		}//if
		else if (p0==1) { // thickness pertubation
			
			// nv is OPTIONAL now
			if ((int)para.para0[i].size() > 5) {
				// fprintf(stderr,"para2mod > thickness paras? %d %g ",i,para.parameter[i]);
				// for (size_t j = 0; j < para.para0[i].size(); ++j) {
				// fprintf(stderr, "| %s", para.para0[i][j].c_str());
				// }
				// fprintf(stderr,"\n");
				// sub-layer thickness (alway less-equal 1)
				nv=atoi(para.para0[i][5].c_str());
				outmodel.groups[ng].ratio[nv]=newv;

				// fprintf(stderr,"\npara2mod > add value: %d %d %d %g ",i,ng,nv,outmodel.groups[ng].ratio[nv]);
			}
			else if ((int)para.para0[i].size() == 5)
			{
				// fprintf(stderr,"para2mod > Total thickness paras? %d %g ",i,para.parameter[i]);
				// for (size_t j = 0; j < para.para0[i].size(); ++j) {
				// fprintf(stderr, "| %s", para.para0[i][j].c_str());
				// }
				// fprintf(stderr,"\n");
				//group thickness
				outmodel.groups[ng].thick=newv;
				// fprintf(stderr,"\npara2mod > add value: %d %d %g ",i,ng,outmodel.groups[ng].thick);
				kflag = 1;
			}

			// fprintf(stderr, "\n");
			}

		else if(p0==-1)// vpvs
			outmodel.groups[ng].vpvs1=newv;
                else if(p0==-2) 
                        outmodel.groups[ng].vpvs=newv;
                else if (p0==-3)
                        outmodel.groups[ng].drho1=newv;
                else if (p0 == -4)
                        outmodel.groups[ng].drho = newv;
                else if (p0 == -5)
                        outmodel.groups[ng].dvs1 = newv;
                else if (p0 == -6)
                        outmodel.groups[ng].dvs = newv;
                else if (p0 == -7) {
                        outmodel.groups[ng].frho = newv;
                        outmodel.groups[ng].fvpvs = newv;
                        outmodel.groups[ng].fdvs = newv;
                        }
                else if (p0 == -8) {
                        outmodel.groups[ng].frho1 = newv;
                        outmodel.groups[ng].fvpvs1 = newv;
                        outmodel.groups[ng].fdvs1 = newv;
                        }
		else
			{cout<<"########wrong!!!!! in para2mode"<<endl;
			 exit (0);}
	  }//for i
          if (kflag) {
             for (int i=0;i<outmodel.ngroup-1;i++) {
		tthick0 = tthick0 + outmodel.groups[i].thick;
		}
		
        outmodel.groups[outmodel.ngroup-1].thick = outmodel.tthick - tthick0;
		// fprintf(stderr,"thickness check outmodel.tthick: %f, tthick0: %f \n",outmodel.tthick,tthick0);
		// // sleep(1);;    
		// for (int i=0;i<outmodel.ngroup;++i)
		// {
		// 	for (int j=0;j<para.nparav;++j)
		// 	{
		// 		fprintf(stderr,"Layer value: %d %d | %g\n",i,j,outmodel.groups[i].value[j]);
		// 		fprintf(stderr,"Layer thick: %d %d | %g\n",i,j,outmodel.groups[i].ratio[j]);
		// 	}
		// 	fprintf(stderr,"Group thick: %d %g\n",i, outmodel.groups[i].thick);
		// }  
		}

	  return 1;
	}//para2mod

//-----------------------------------------------------	 
	int get_misfit(paradef &para, modeldef &inmodel,modeldef &outmodel,double p, int invtype)
	{

	// fprintf(stderr,"> CALpara: get_misfit\n");
    // sleep(1);;

	  paradef temppara;
	  temppara=para;
	  para2mod(para,inmodel,outmodel,invtype);
//	  cout<<outmodel.ngroup<<endl;
	  updatemodel(outmodel,invtype);
	  compute_rf(outmodel,invtype);
	  compute_disp(outmodel,invtype);
	  compute_misfit(outmodel,p);
	  mod2para(outmodel,temppara,para,invtype);
	  return 1;
	}//get misfit
//};
	int get_misfit1(paradef &para, modeldef &outmodel, double p, int invtype) {

		// fprintf(stderr,"> CALpara: get_misfit1\n");
    	// sleep(1);;

        paradef temppara;
        temppara = para;
        //#pragma omp critical(f2)
        //if (p!=1)
        compute_rf(outmodel,invtype);
	//if (p!=0)
        #pragma omp critical(f2)//taken out by EMB
        {
        compute_disp(outmodel,invtype);
        }
        //cout<<"now compute misfit"<<endl;
        compute_misfit(outmodel,p);
        //cout<<"now convert back"<<endl;
        mod2para(outmodel,temppara,para, invtype);
        return 1;
        }

