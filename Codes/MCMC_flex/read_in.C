/////////////////////////////////////////////////////////////////////////
/////////////////////////// read_in /////////////////////////////////////
// Read the control file to grasp the informations
//
//
///////////////////////////////////////////////////////////////////////// 
#include<map>

int read_in( char *fn, char *inmodelnm, char *inparanm , int *surtype, int *surn, \
vector <int> &surflag, vector <string> &vdisp, \
char *inrfnm, double *gaus, char *Qname, int *Qflag, \
vector <int> &Ql, double *pp, int *monoc, vector <int> &vmono, \
vector <int> &vgrad, double *TTT, int *jumpn, \
char *outdir, char *outname, int *flagw, int *n_core, char *indatanm,\
double *mc_range, double *dstep, int *plottype, int *invtype) 
{
  
  fprintf(stderr,"> read_in: read in\n");
  sleep(1);

  ifstream mf;
  string line;
  int i,j,k;
  char scheck[300];

  mf.open(fn);
  if(not mf.is_open())
    {
      cout<<"#####control file does not exist!\n"<<fn<<endl;
      fprintf(stderr,"#####control file %s does not exist!\n",fn);
      exit(0);
      }

  // initialize numbers for switch // 
  std::map<std::string,int> numbers;

  numbers["model"] = 1;
  numbers["para"] = 2;
  numbers["disp"] = 3;
  numbers["rf"] = 4;
  numbers["Qmodel"] = 5;
  numbers["p"] = 6;
  numbers["mono"] = 7;
  numbers["monol"] = 8;
  numbers["gradl"] = 9;
  numbers["tt"] = 10;
  numbers["jt"] = 11;
  numbers["outdir"] = 12;
  numbers["n_thread"] = 13;
  numbers["indata"] = 14;
  numbers["mc_range"] = 15;
  numbers["dstep"] = 16;
  numbers["plottype"] = 17;
  numbers["invtype"] = 18;
  numbers["end"] = 0;

  *flagw = -1;
  
  vector <string> vargv;
  
  while (getline(mf,line))
  {
    cout<<line<<endl;
    vargv.clear();
    Split(line,vargv," ");
    sprintf(scheck,"%s",vargv[0].c_str());
  
    switch (numbers[scheck]) 
    {
      case 1:
        if (vargv.size() != 2) 
        {
          cout<<"wrong input mmodel!"<<endl; 
          return 0;
          }
          sprintf(inmodelnm,"%s",vargv[1].c_str());
          break;

      case 2:
        if (vargv.size() != 2)
        {
          cout<<"wrong input mpara!"<<endl; return 0;
          }
          sprintf(inparanm,"%s",vargv[1].c_str());
          break;

      case 3:
        if (vargv.size() < 2) 
        {
          cout<<"wrong input surf!"<<endl; return 0;
          }
          *surtype = atoi(vargv[1].c_str());
          *surn = atoi(vargv[2].c_str());
          if (vargv.size() != 3+2*(*surn)) 
          {
            cout<<"wrong input surf!"<<endl; 
            return 0;
            }
          for (i=0;i<*surn;i++) 
          {
            surflag.push_back(atoi(vargv[3+i*2].c_str()));
            vdisp.push_back(vargv[4+i*2].c_str());
            cout<<vdisp[i]<<endl;
            }
          cout<<*surn<<endl;
          break;

         case 4:
           if (vargv.size() < 2) {cout<<"wrong input rf!"<<endl; return 0;}
           sprintf(inrfnm,"%s",vargv[1].c_str());
           if (vargv.size() > 2) *gaus = atof(vargv[2].c_str());
           else *gaus = -1.;
         break;

         case 5:
           if (vargv.size() < 2) {cout<<"wrong input Qmodel!"<<endl; return 0;}
           sprintf(Qname,"%s",vargv[1].c_str());
           *Qflag = 1;
           //j = vargv.size()-2;
           //for (i=0;i<j;i++) Ql.push_back(atoi(vargv[2+i].c_str()));
         break;

         case 6: *pp = atof(vargv[1].c_str());
         break;

         case 7: *monoc = atoi(vargv[1].c_str());
         break;

         case 8:
           j = vargv.size()-1;
           for (i=0;i<j;i++) vmono.push_back(atoi(vargv[1+i].c_str()));
         break;

         case 9:
           j = vargv.size()-1;
           for (i=0;i<j;i++) vgrad.push_back(atoi(vargv[1+i].c_str()));
         break;

         case 10:
           *TTT = atof(vargv[1].c_str());
         break;

         case 11:
           *jumpn = atoi(vargv[1].c_str());
         break;

         case 12:
           sprintf(outdir,"%s\0",vargv[1].c_str());
           sprintf(outname,"%s\0",vargv[2].c_str());
           if (vargv.size() > 3) {*flagw = atoi(vargv[3].c_str());}
         break;

         case 13:
           *n_core = atoi(vargv[1].c_str());
        
         case 14:
          if (vargv.size() != 2)
          {
            cout<<"wrong input indata!"<<endl; return 0;
          }
            sprintf(indatanm,"%s",vargv[1].c_str());
          break;
        
        case 15:
        if (vargv.size() != 2)
          {
            cout<<"wrong mc_range input!"<<endl; return 0;
          }
           *mc_range = atof(vargv[1].c_str());
        break;

        case 16:
        if (vargv.size() != 2)
          {
            cout<<"wrong dstep input!"<<endl; return 0;
          }
           *dstep = atof(vargv[1].c_str());
        break;

        case 17:
        if (vargv.size() != 2)
          {
            cout<<"wrong plot type flag!"<<endl; return 0;
          }
           *plottype = atoi(vargv[1].c_str());
        break;

        case 18:
        if (vargv.size() != 2)
          {
            cout<<"wrong invert type flag!"<<endl; return 0;
          }
           *invtype = atoi(vargv[1].c_str());
        break;

         case 0:
           cout<<"input is over!"<<endl;
         break;

         default:
           cout<<"I don't know waht's this, return"<<endl;
           return 0;
         break;
       }
     }

  mf.close();
  cout<<"########## read over ! ###################"<<endl;
  cout<<"#####      read control.dat ok ! #########\n";

  fprintf(stderr,"> read %s ok !\n",fn);
  // sleep(10);

  return 1;
  }

