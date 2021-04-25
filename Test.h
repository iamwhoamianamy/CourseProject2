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
      return -1 * divlambdagrad(x, y, t) +
         sigma() * dudt(x, y, t) + chi() * d2udt2(x, y, t);
   }

   double lambda(const double& x, const double& y)
   {
      return 1;
   }

   // Точное решение
   double u(const double& x, const double& y, const double& t)
   {
      switch(N)
      {
         case(0): return 2.0;
         case(1): return x + y;
         case(2): return x * x + y * y;
         case(3): return x * x * x + y * y * y;
         case(4): return x * x * x * x + y * y * y * y;
      };
   }

   double divlambdagrad(const double& x, const double& y, const double& t)
   {
      switch(N)
      {
         case(0): return 0;
         case(1): return 0;
         case(2): return 4;
         case(3): return 4;
         case(4): return 12 * x * x + 12 * y * y;
      };
   }

   double sigma()
   {
      return 1;
   }

   double chi()
   {
      return 1;
   }

   // Производная точного решения по t
   double dudt(const double& x, const double& y, const double& t)
   {
      switch(N)
      {
         case(0): return 0;
         case(1): return 0;
         case(2): return 0;
         case(3): return 0;
         case(4): return 0;
      };
   }

   // Вторая производная точного решения по t
   double d2udt2(const double& x, const double& y, const double& t)
   {
      switch(N)
      {
         case(0): return 0;
         case(1): return 0;
         case(2): return 0;
         case(3): return 0;
         case(4): return 0;
      };
   }
};