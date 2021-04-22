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
   bvp.test = Test(1);

   bvp.FormPortrait();
   //bvp.BuildMatrices(0);
   //bvp.AssembleGlobalMatrix();
   //bvp.AccountFirstBound(0);

   //bvp.Solve();

   ofstream fout("result.txt");
   //bvp.PrintSolution(fout, 0.0);
   bvp.TimeIterations(fout);
   fout.close();

   int asdasd;
}