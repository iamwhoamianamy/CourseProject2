#include <iostream>
#include "BoundaryValueProblem.h"

using namespace std;

double phi1(double t)
{
   return 2 * (t - 0.5) * (t - 1);
}

double phi2(double t)
{
   return -4 * t * (t - 1);
}

double phi3(double t)
{
   return 2 * t * (t - 0.5);
}

int main()
{
   BoundaryValueProblem bvp;

   bvp.ReadFormGrid("data/grid.txt");
   bvp.ReadBoundaries("data/boundaries.txt");
   bvp.ReadFormTimeGrid("data/time.txt");
   bvp.ReadMatrices();
   bvp.InitializeMemory();

   bvp.FormPortrait();

   ofstream fout("result.txt");

   bvp.TimeIterations(fout);

   fout.close();
}