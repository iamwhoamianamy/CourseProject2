#include <iostream>
#include "BoundaryValueProblem.h"

using namespace std;

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