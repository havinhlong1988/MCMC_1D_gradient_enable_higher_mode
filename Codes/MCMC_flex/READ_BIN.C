#include<iostream>
#include<algorithm>
#include<vector>
#include<cmath>
#include<fstream>

#include"/home/jiayi/progs/jy/HEAD/head_c++/string_split.C"
#include"/home/jiayi/progs/jy/HEAD/head_c++/generate_Bs.C"
#include"/home/jiayi/progs/jy/HEAD/head_c++/gen_random.C"
#include"INITstructure.h"
#include"CALgroup.C"
#include"CALmodel.C"
//#include"CALforward.C"
//#include"CALpara.C"

#include"BIN_rw.C"
using namespace std;
int main( int argc, char *argv[])
{
  fprintf(stderr,"> READ_BIN: main\n");
  sleep(1);

if (argc != 2 ) {
  cout<<"input [in.BIN]"<<endl;
  return 0;
  }
ifstream inbin;
inbin.open(argv[1],ios::in| ios::binary);
cout<<argv[1]<<endl;
//modeldef mod0;
//paradef para0;

read_bin(inbin);
}
