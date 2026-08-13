// #include<iostream>
// #include <stdio.h>
using namespace std;

/*========CONTENT===========
int updategroup1(groupdef &group)//layer the model from layered model
int updategroupBs(groupdef &group) //calculate the Bspline function
int updategroup2(groupdef &group)//layer the model from Bspline model
int updategroup3(groupdef &group)//layer the model from gradient model
int updategroup4(groupdef &group)//..water layer
int updategroup(groupdef &group)
//==========================
*/
//----------------------------------------------------- 
        double get_vp(double ts, int flag, double vpvs) 
          {

            // fprintf(stderr,"> CALgroup: get_vp\n");
            // sleep(1);;
          double tp;
          // flag: 1(Brocher, 2005); 
          // flag: 2 AK135, 120km, vpvs=1.789; 
          // flag: 3 vpvs = vpvs
          // flag: 4 Huang 2014 relationship
          if (flag==1) {

              tp = 0.9409 + 2.0947*ts - 0.8206*ts*ts + \
                  0.2683*ts*ts*ts -0.0251*ts*ts*ts*ts;
              }
          else if (flag==2) {
              tp = ts*1.789;
              //fprintf(stderr,"flag vp is 2 \n");
          }
          else if (flag==3) tp = ts*vpvs;
          else if (flag==4) tp = ts*(1.373+0.699/ts);
          else tp = ts*1.75;
            //fprintf(stderr,"vp flag: %d \n",flag);
          return tp;
          }
//-----------------------------------------------------
        double get_rho(double tp,double ts,int flag, double drho) 
          {
            // fprintf(stderr,"> CALgroup: get_rho\n");
            // sleep(1);;
          double trho,trho1;
            //fprintf(stderr,"rho flag: %d \n",flag);
          //flag=1; //added by EMB
          // flag: 1(Brocher, 2005); flag: 2(Kaban, M. K et al. (2003), Density of the continental roots)
          // flag: 1(Brocher, 2005, with relationship between rho and vs), regression made by Weisen Shen 2014 Jan.;
          if (flag == 1) {
            //fprintf(stderr,"rho flag 1\n");
            trho = 1.22679 + 1.53201*ts -0.83668*ts*ts + 0.20673*ts*ts*ts -0.01656*ts*ts*ts*ts;//Weisen's
            trho1=trho;
            if (trho1 < 0.) trho1 = 0.;
            return trho1;
            }
          else if (flag == 2) {
            //fprintf(stderr,'rho flag==2\n');
            //fprintf(stderr,"rho flag 2\n");
            trho = 3.42+0.01*100*(ts-4.5)/4.5; //10 kg/m3 per 1%
            trho1 = trho + drho;
            if (trho1 < 0.) trho1 = 0.;
            return trho1;
            }
          /*if (flag == 3) {
            trho = 1.05; 
            trho1 = trho + drho;
            if (trho1 < 0.) trho1 = 0.;
            return trho1;
            //return trho+drho;
            }
          if (flag == 4) {
            trho=2.7;
            trho1=trho+drho;
            if (trho1 < 0.) trho1 = 0.;
            return trho1;
            }
          if (flag == 5) {
            trho=3.38;
            trho1=trho+drho;
            if (trho1 < 0.) trho1 = 0.;
            return trho1;
            }
          if (flag == 6) {   // http://en.wikipedia.org/wiki/Gardner's_relation 
            trho = 0.31 * pow (tp*1000., 0.25);
            trho1 = trho + drho;
            if (trho1 < 0.) trho1 = 0.;
            return trho1;
            }*/
          else {
            fprintf(stderr,"trho is not specified! using Brocher's formula\n");
            //trho = 1.6612*tp - 0.4721*tp*tp + 0.0671 *(tp*tp*tp) \
                         -0.0043 *(tp*tp*tp*tp) + 0.000106 *(tp*tp*tp*tp*tp);
            trho = 1.22679 + 1.53201*ts -0.83668*ts*ts + 0.20673*ts*ts*ts -0.01656*ts*ts*ts*ts;
            if (trho1 < 0.) trho1 = 0.;
            return trho1;
            }
          }
// =================================================================================================
int updategroup1(groupdef &group)  // layer the model from layered model
{
    int tn0;
    double dd;
    double tp, ts, trho, tvpvs, tdrho;
    double cst; // current sub-block thickness
    int r_flag;
    tdrho = 0.;

    // Only flag==1 is valid for layered model update
    if (group.flag != 1) {
        cout << "Wrong flag: updating layered model requires flag=1, but group.flag="
             << group.flag << endl;
        return 0;
    }

    r_flag = group.r_flag;
    dd = group.ddep;      // target sub-layer thickness

    const double eps = 1e-8;  // tolerance for floating-point remainder

    // ----------------------------------------------------------------------
    // 1) Precompute Vs, Vp/Vs control and density increment per coarse block
    //    (block-level properties, before vertical subdivision)
    // ----------------------------------------------------------------------
    std::vector<double> vs_block(group.np);
    std::vector<double> tvpvs_block(group.np);
    std::vector<double> tdrho_block(group.np);

    for (int i = 0; i < group.np; ++i) {
        ts = group.value[i];

        if (group.fvpvs1 <= 0.) {
            // Vs refinement
            if (int(group.fvpvs) <= i + 1)
                ts = ts + group.dvs1;
            else
                ts = ts + group.dvs;

            // Vp/Vs ratio
            tvpvs = group.vpvs;
            if (int(group.fvpvs) <= i + 1 && group.fvpvs != 0)
                tvpvs = group.vpvs1;

            // Density increment
            tdrho = group.drho;
            if (int(group.frho) <= i + 1 && group.frho != 0)
                tdrho = group.drho1;
        }
        else {
            // Vs refinement for fvpvs1 case
            if (int(group.fvpvs1) <= i + 1 && int(group.fvpvs) > i + 1)
                ts = ts + group.dvs1;
            else
                ts = ts + group.dvs;

            // Vp/Vs ratio
            tvpvs = group.vpvs;
            if (int(group.fvpvs1) <= i + 1 && int(group.fvpvs) > i + 1 && group.fvpvs != 0)
                tvpvs = group.vpvs1;

            // Density increment
            tdrho = group.drho;
            if (int(group.fvpvs1) <= i + 1 && int(group.frho) > i + 1 && group.frho != 0)
                tdrho = group.drho1;
        }

        vs_block[i]     = ts;      // block Vs after refinement
        tvpvs_block[i]  = tvpvs;   // block Vp/Vs
        tdrho_block[i]  = tdrho;   // block density increment param
    }

    // ----------------------------------------------------------------------
    // 2) Build the fine-layered model with Vs scaled between block tops/bottoms
    // ----------------------------------------------------------------------

    // Clear previous layered arrays
    group.thick1.clear();
    group.value1.clear();
    group.value1p.clear();
    group.value1r.clear();

    tn0 = 0;
    int laycount = 0;
    double l_thick_check = 0.0;

    for (int i = 0; i < group.np; i++)
    {
        // Total thickness of coarse block i
        cst = group.thick * group.ratio[i];

        // Decide number of sub-layers based on cst and dd
        int tn = 1;
        if (cst > dd) {
            // at least 2 sub-layers if thicker than dd
            tn = std::max(2, int(std::round(cst / dd)));
        }

        double sub_layer_thick = cst / double(tn);          // equal thickness per sub-layer
        double rem = cst - sub_layer_thick * tn;            // numerical remainder

        // Vs at top and bottom of this block (after refinement)
        double vs_top    = vs_block[i];
        double vs_bottom = vs_top;
        if (i < group.np - 1) {
            vs_bottom = vs_block[i + 1];
        }

        // Use block i's Vp/Vs and density increment for all sub-layers
        tvpvs = tvpvs_block[i];
        tdrho = tdrho_block[i];

        // Push tn sub-layers with Vs scaled from top to bottom
        for (int k = 0; k < tn; k++) {
            tn0++;
            laycount++;

            double alpha = (tn > 1) ? (double)k / (tn - 1) : 0.0;   // 0 at top, 1 at bottom
            double ts_k = (vs_top * (1.0 - alpha)) + (vs_bottom * alpha); // this give the gradient value
            // double ts_k = vs_top; // this give the uniform value
            // double ts_k = 0.5*(vs_top+vs_bottom); // this give the mean uniform value a

            double tp_k   = get_vp(ts_k, group.p_flag, tvpvs);
            double trho_k = get_rho(tp_k, ts_k, r_flag, tdrho);

            group.thick1.push_back(sub_layer_thick);
            group.value1.push_back(ts_k);
            group.value1p.push_back(tp_k);
            group.value1r.push_back(trho_k);

            l_thick_check += sub_layer_thick;
        }

        // Optionally add the small remainder rem into the last sub-layer thickness
        if (std::fabs(rem) > eps * dd && !group.thick1.empty()) {
            group.thick1.back() += rem;
            l_thick_check += rem;
        }
    }

    // Update total number of layers
    group.nlay = tn0;

    // fprintf(stderr, "TOTAL tn0 = %d, laycount = %d\n", tn0, laycount);
    // // fprintf(stderr,"Total thickness check = %f (should match group.thick)\n",
    // //         l_thick_check);
    // for (int i=0;i<tn0;i++)
    // {
    //   fprintf(stderr, "Layer = %d | thickness = %f | vs = %f\n", i, group.thick1[i],group.value1[i]);
    // }
    return 1;
}
//updategroup1
//----------------------------------------------------- 
        int updategroupBs(groupdef &group) //calculate the Bspline function     
        {
          // fprintf(stderr,"> CALgroup: updategroupBs\n");
          // sleep(1);;

          int nBs,order;
          double factor;
          double dd;
          /*
          // HVL We can layerization the model here!
          // This is default setting from EMD
          if(group.thick>=150)group.nlay=60;
	        else if(group.thick<5)group.nlay=5;
	        else if(group.thick<10)group.nlay=10;
	        else if(group.thick<20 && group.thick>10) group.nlay=20;
          else group.nlay=60;
          
          */
          // if(group.thick>=150)group.nlay=30;
          // else if(group.thick<5) group.nlay=int(group.thick/1.0);
          // else group.nlay=int(group.thick/2.5); // HVL: try layerized the model each 1.5 km
          // //Number of layers
          // if(group.thick>=150)group.nlay=60;
	        // else if(group.thick<5)group.nlay=5;
	        // else if(group.thick<10)group.nlay=10;
	        // else if(group.thick<20 && group.thick>10) group.nlay=20;
          // else group.nlay=150;
          dd = group.ddep;
          if (group.thick > dd) group.nlay=int(group.thick/dd);
          if (group.thick < dd) group.nlay=2;



          nBs=group.np;
          if(nBs<4) order=3;
          else order=4; //cubic Bspline
          factor = 1.35; // factor = 1
	        gen_B_spline(nBs,order,0.,group.thick,factor,group.nlay,group.Bsplines);
          group.flagBs=1;
	  return 1;
      }//updategroupBs
//----------------------------------------------------- 
        int updategroup2(groupdef &group)//layer the model from Bspline model   
        {
          
          // fprintf(stderr,"> CALgroup: updategroup2\n");
          // sleep(1);;

	        int i,j,nnlay,nBs,r_flag;
	        double tmp,tvpvs,tdrho;
          double tp,ts,trho;
          tdrho = 0.;
          if(group.flag!=2){cout<<"wrong flag, here updating Bspline model, flag=2, group.flag="<<group.flag<<endl;return 0;}
	  //if (group.thick1.size()+group.value1.size()!=0){cout<<"problem! the thick1 or value1 is not empty!"<<endl;return 0;}
	        group.thick1.clear();
          group.value1.clear();
      /*cout<<"D: gv0 "<<group.value[0]<<endl;
      cout<<"D: gv1 "<<group.value[1]<<endl;
      cout<<"D: gv2 "<<group.value[2]<<endl;
      cout<<"D: gv3 "<<group.value[3]<<endl;
      cout<<"D: gv4 "<<group.value[4]<<endl;*/
          if(group.flagBs!=1)updategroupBs(group);
	        nnlay=group.nlay;
          // fprintf(stderr,"nnlay: %d & group thick: %f \n",nnlay,group.thick);
          // sleep(3);
          r_flag = group.r_flag;
          group.value1p.clear();
          group.value1r.clear();
	  for (i=0;i<nnlay;i++)
		{
		  nBs=group.np;
		  tmp=0.;
		  for (j=0;j<nBs;j++){tmp=tmp+group.Bsplines[j*nnlay+i]*group.value[j];} //???
		  //group.value1.push_back(tmp);
      //ts = tmp + group.dvs;
      //

      //cout<<"fvpvs1"<<group.fvpvs1<<endl;
      if (group.fvpvs1 <= 0.) {
        if (group.fvpvs != 0 and i < group.fvpvs * nnlay) ts = tmp + group.dvs1;
        else ts = tmp + group.dvs;

        tvpvs = group.vpvs;
        if (group.fvpvs != 0 and i < group.fvpvs * nnlay) tvpvs = group.vpvs1;
        tp = get_vp(ts,group.p_flag,tvpvs);

        tdrho = group.drho;
        if (group.frho != 0 and i < (int)(group.frho * nnlay)) tdrho = group.drho1;
        trho = get_rho(tp,ts,r_flag,tdrho);
        }
      else {
        if (group.fvpvs != 0 and i >= group.fvpvs * nnlay and i < (int)((group.fvpvs1+group.fvpvs) * nnlay)) ts = tmp + group.dvs1;
        else ts = tmp + group.dvs;

        tvpvs = group.vpvs;
        if (group.fvpvs != 0 and i >= group.fvpvs * nnlay and i < (int)((group.fvpvs1+group.fvpvs) * nnlay)) tvpvs = group.vpvs1;
        tp = get_vp(ts,group.p_flag,tvpvs);

        tdrho = group.drho;
        if (group.frho != 0 and i >= (int)(group.frho * nnlay) and (int)((group.fvpvs1+group.fvpvs) * nnlay)) tdrho = group.drho1;
        trho = get_rho(tp,ts,r_flag,tdrho);
        }
      
      group.value1.push_back(ts);
      group.value1p.push_back(tp);
      group.value1r.push_back(trho);
      group.thick1.push_back(group.thick/nnlay);
		}
	  return 1;		
        }//updategroup2

//----------------------------------------------------- 
        int updategroup3(groupdef &group) //layer the model from gradient model
        {

          // fprintf(stderr,"> CALgroup: updategroup3\n");
          // sleep(1);;

	        int tn,r_flag;
	        double dh,dv;
          double dd;
          double tp,ts,trho,tvpvs,tdrho;
          tdrho = 0.;
	        if(group.flag!=4){cout<<"wrong flag, here updating gradient model, flag=4, group.flag="<<group.flag<<endl;return 0;}
	//  if (group.thick1.size()+group.value1.size()!=0){cout<<"problem! the thick1 or value1 is not empty!"<<endl;return 0;}
	        group.thick1.clear();group.value1.clear();
          group.value1p.clear(); group.value1r.clear();

	        // if(group.thick>2.5)tn=min(int(group.thick/0.5),10);
	        // else if(group.thick<0.5 && group.thick>0.1)tn=5;
          // else if (group.thick<0.1000001) tn=2;
	        // else tn=5;
          dd=group.ddep;
          if (group.thick > dd) tn=int(group.thick/dd);
          if (group.thick < dd) tn=2;

	        dh=group.thick/double(tn);
	        dv=(group.value[1]-group.value[0])/(tn-1.);//by default, only two values,[0] is smaller than [1]
	        r_flag = group.r_flag;
	        for(int j=0;j<tn;j++)
          {
		      //group.value1.push_back(group.value[0]+j*dv);
		      group.thick1.push_back(dh);
                  //ts = group.value[0]+j*dv;
                  if (group.fvpvs != 0 and j< tn*group.fvpvs) ts = group.value[0]+j*dv+group.dvs1;
                  else  ts = group.value[0]+j*dv+group.dvs;

                  //cout<<ts<<endl;

                  tvpvs = group.vpvs;
                  if (group.fvpvs != 0 and j< tn*group.fvpvs) tvpvs = group.vpvs1;
                  //tp = tvpvs*ts;
                  tp = get_vp(ts,group.p_flag,tvpvs);

                  tdrho = group.drho;

                  group.value1.push_back(ts);
                  trho = get_rho(tp,ts,r_flag,tdrho);
                  group.value1p.push_back(tp);
                  group.value1r.push_back(trho);
		}
	  group.nlay=tn;
	  return 1;
        }//updategroup3
//----------------------------------------------------- 
	int updategroup4(groupdef &group)//layer the water layer
	{

      // fprintf(stderr,"> CALgroup: updategroup4\n");
      // sleep(1);;

	  if(group.flag!=5){cout<<"wrong flag, here updating water model, flag=5, group.flag="<<group.flag<<endl;return 0;}
	  group.thick1.clear();group.value1.clear();
    group.value1p.clear();group.value1r.clear();
	  //group.value1.push_back(group.value[0]-0.);
	  group.value1.push_back(0.);   // inpute value of 0km/sec
	  group.thick1.push_back(group.thick-0.);
    group.value1p.push_back(group.value[0]-0.);
    group.value1r.push_back(1.02);
	  group.nlay=1;
	  return 1;
	}//updategroup4
//----------------------------------------------------- 
        int updategroup(groupdef &group)
        {
          // fprintf(stderr,"> CALgroup: updategroup\n");
          // sleep(1);
	      if(group.flag==1)updategroup1(group);//lay
	      else if(group.flag==2)updategroup2(group);//Bs
	      else if(group.flag==4)updategroup3(group);//grad
	      else if (group.flag==5)updategroup4(group);//water layer
	      else cout<<"wrong value in group.flag : "<<group.flag<<endl;
	      return 1;
        }//updategroup
