#include <iostream>
#include <map>
#include <stdlib.h>
#include <fstream>
using namespace std;

bool fncomp (char lhs, char rhs) {return lhs<rhs;}

struct classcomp {
  bool operator() (const char& lhs, const char& rhs) const
  {return lhs<rhs;}
};

int main ()
{
  std::map<std::string,int> first;

  first["ab"]=10;
  first["cd"]=30;
  first["ef"]=50;
  first["gh"]=70;

  cout<<first["ab"]<<endl;
  return 0;
}
