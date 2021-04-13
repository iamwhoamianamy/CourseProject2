#include <iostream>
#include "BoundaryValueProblem.h"

using namespace std;

int main()
{
   BoundaryValueProblem bvp;

   bvp.ReadFormGrid("data/grid.txt");
   bvp.ReadMatrices();
   bvp.InitializeMemory();
   bvp.test = Test(3);

   bvp.FormPortrait();
   bvp.BuildMatrices();
   bvp.AssembleGlobalMatrix();
   bvp.AccountFirstBound();

   bvp.Solve();

   ofstream fout("result.txt");
   bvp.PrintSolution(fout);
   fout.close();

   int asdasd;
}