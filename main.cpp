#include <iostream>
#include "BoundaryValueProblem.h"

using namespace std;

int main()
{
   BoundaryValueProblem bvp;

   bvp.ReadFormGrid("grid.txt");
   bvp.InitializeMemory();

   vector<int> vec(9);
   bvp.FormPortrait();

   int asdasd;
}