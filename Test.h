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
      return x * x;
   }

   // ������ �������
   double u(const double& x, const double& y, const double& t)
   {
      switch(N)
      {
         case(0): return 2.0;
         case(1): return x + t * t * t;
         case(2): return x * x + y * y;
         case(3): return x * x * x + y * y * y;
         case(4): return x * x * x * x;
      };
   }

   double divlambdagrad(const double& x, const double& y, const double& t)
   {
      switch(N)
      {
         case(0): return 0;
         case(1): return 2 * x;
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

   // ����������� ������� ������� �� t
   double dudt(const double& x, const double& y, const double& t)
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

   // ������ ����������� ������� ������� �� t
   double d2udt2(const double& x, const double& y, const double& t)
   {
      switch(N)
      {
         case(0): return 0;
         case(1): return 6 * t;
         case(2): return 0;
         case(3): return 0;
         case(4): return 0;
      };
   }
};