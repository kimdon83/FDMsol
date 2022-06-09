#ifndef _DEFINE_H_
#define _DEFINE_H_

#include <iomanip>
#include <iostream>
#include <string>

#include <omp.h>

#include <time.h>
#include <sys/timeb.h>

using namespace std;

typedef double real;

#define NUM_GRP     2 //CX's energy group number

#define CARDLEN     64  //Maximum character length of cards used in input file
#define ENDSTR      ");" //Ending string for a input section
#define STR_EMPTY   "."

//Indices of Direction
#define XDIR    0
#define YDIR    1
#define ZDIR    2

//Indices of Neighbors
#define NEIGH_NUM 6

#define XL 0
#define XR 1
#define YL 2
#define YR 3
#define ZL 4
#define ZR 5
#define OWN 6

//Returning reaction types of selecting collision type
#define REA_ABS     0
#define REA_SCA     1
#define REA_FIS     2

#endif