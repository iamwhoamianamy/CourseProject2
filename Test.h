#pragma once
#pragma once
using namespace std;

class Test
{
public:
   int test_n = 0;            // Номер теста
   int lambda_n = 0;          // Номер формулы для коэффициента lambda
   int sigma_n = 0;           // Номер формулы для коэффициента sigma
   int chi_n = 0;             // Номер формулы для коэффициента chi

   Test(const int& t_N) : test_n(t_N) {};
   Test() { };

   double f(const double& x, const double& y, const double& t)
   {
      return -1 * divlambdagrad(x, y, t) +
         sigma() * dudt(x, y, t) + chi() * d2udt2(x, y, t);
   }

   double lambda(const double& x, const double& y)
   {
      switch(lambda_n)
      {
         case 0: return 1;
         case 1: return 2;
      }
   }

   double sigma()
   {
      switch(sigma_n)
      {
         case 0: return 1;
         case 1: return 4;
      }
   }

   double chi()
   {
      switch(chi_n)
      {
         case 0: return 1;
         case 1: return 0.3;
      }
   }

   // Точное решение
   double u(const double& x, const double& y, const double& t)
   {
      switch(test_n)
      {
         case 0: return 2.0;
         case 1: return x + y;
         case 2: return x * x + y * y;
         case 3: return x * x * x + y * y * y;
         case 4: return x * x * x * x + y * y * y * y;

         case 5: return t * t * t * t;
         case 6: return t * t * t;
         case 7: return x * x * t * t * t * t;
         
      };
   }

   double divlambdagrad(const double& x, const double& y, const double& t)
   {

      switch(test_n)
      {
         case 0: return 0;
         case 1: return 0;
         case 2: return 4;
         case 3: return 6 * x + 6 * y;
         case 4: return 12 * x * x + 12 * y * y;

         case 5: return 0;
         case 6: return 0;
         case 7: return 2 * t * t * t * t;
      };

   }


   // Производная точного решения по t
   double dudt(const double& x, const double& y, const double& t)
   {
      switch(test_n)
      {
         case 0: return 0;
         case 1: return 0;
         case 2: return 0;
         case 3: return 0;
         case 4: return 0;
              
         case 5: return 4 * t * t * t;
         case 6: return 3 * t * t;
         case 7: return 4 * t * t * t * x * x;
      };
   }

   // Вторая производная точного решения по t
   double d2udt2(const double& x, const double& y, const double& t)
   {
      switch(test_n)
      {
         case 0: return 0;
         case 1: return 0;
         case 2: return 0;
         case 3: return 0;
         case 4: return 0;
              
         case 5: return 12 * t * t;
         case 6: return 6 * t;
         case 7: return 12 * t * t * x * x;
      };
   }
};