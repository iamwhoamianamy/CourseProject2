#pragma once
#pragma once
using namespace std;

class Test
{
public:

   int N;

   Test(const int& t_N) : N(t_N) {};

   Test() : N(0){};

   double f(const double& x, const double& y, const double& t)
   {
      return -divgrad(x, y, t) * lambda() + sigma() * u(x, y, t);
   }

   /*double lambda(const double& x, const double& y)
   {
      return 1;
   }*/

   double lambda()
   {
      return 1;
   }

   double sigma()
   {
      return 1;
   }

   // ������ �������
   double u(const double& x, const double& y, const double& t)
   {
      switch(N)
      {
         case(0): return 2.0;
         case(1): return x;
         case(2): return x * x;
         case(3): return x * x * x;
         case(4): return x * x * x * x;
      };
   }

   double divgrad(const double& x, const double& y, const double& t)
   {
      switch(N)
      {
         case(0): return 0;
         case(1): return 0;
         case(2): return 2;
         case(3): return 6 * x;
         case(4): return 12 * x * x;
      };
   }

   // ����������� ������� ������� �� t
   double dudt(const double& x, const double& t)
   {
      switch(N)
      {
         case(0): return 0;
         case(1): return 3 * t * t;
         case(2): return 0;
         case(3): return 0;
         case(4): return 0;
      };
   }
};