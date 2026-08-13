#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>

using namespace std;

typedef struct data{
  int myint;
  char* mystring;
};

int whatsMyInt(data** arrayPtr, int x, int y){
  return arrayPtr[y][x].myint;
}

int main () {
//initialize
int x,y,w,h;
data** array;

w = 10; //width of array
h = 20; //height of array

//malloc the 'y' dimension
array = (data**) malloc(sizeof(data*) * h);

//iterate over 'y' dimension
for(y=0;y<h;y++){
  //malloc the 'x' dimension
  array[y] = (data*)malloc(sizeof(data) * w);

  //iterate over the 'x' dimension
  for(x=0;x<w;x++){
    //malloc the string in the data structure
    array[y][x].mystring = (char*)malloc(sizeof(char) * 50); //50 chars

    //initialize
    array[y][x].myint = 6;
    array[y][x].mystring = "w00t";
  }
}
printf("My int is %d.\n", whatsMyInt(array, 2, 4));
return 1;
}
