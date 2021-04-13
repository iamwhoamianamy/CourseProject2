#pragma once
#include <vector>
using namespace std;
typedef double real;

class Test
{
public:
   vector<real> f(real x, real y)
   {
      return  { 2 * (sin(x)),
                2 * (sin(x)),
                2 * (sin(x)) };
   }

   vector<real> lambda()
   {
      return { 1, 1, 1 };
   }

   vector<real> gamma()
   {
      return { 1, 1, 1 };
   }

   vector<real> ug(real x, real y)
   {
      return { sin(x),
               sin(x),
               sin(x) };
   }

   vector<real> u(real x, real y)
   {
      return { sin(x),
               sin(x),
               sin(x) };
   }
};

//class Test
//{
//public:
//   vector<real> f(real x, real y)
//   {
//      return  { -6 * (x + y) + (pow(x, 3) + pow(y, 3)),
//                -6 * (x + y) + (pow(x, 3) + pow(y, 3)),
//                -6 * (x + y) + (pow(x, 3) + pow(y, 3)) };
//   }
//
//   vector<real> lambda()
//   {
//      return { 1,
//               1,
//               1 };
//   }
//
//   vector<real> gamma()
//   {
//      return { 1,
//               1,
//               1 };
//   }
//
//   vector<real> ug(real x, real y)
//   {
//      return { pow(x, 3) + pow(y, 3),
//               pow(x, 3) + pow(y, 3),
//               pow(x, 3) + pow(y, 3) };
//   }
//
//   vector<real> u(real x, real y)
//   {
//      return { pow(x, 3) + pow(y, 3),
//               pow(x, 3) + pow(y, 3),
//               pow(x, 3) + pow(y, 3) };
//   }
//};

//class Test
//{
//public:
//   vector<real> f(real x, real y)
//   {
//      return  { 1 * (x * x + y * y),
//                1 * (x * x + y * y),
//                1 * (x * x + y * y) };
//   }
//
//   vector<real> lambda()
//   {
//      return { 0,
//               0,
//               0 };
//   }
//
//   vector<real> gamma()
//   {
//      return { 1,
//               1,
//               1 };
//   }
//
//   vector<real> ug(real x, real y)
//   {
//      return { x * x + y * y,
//               x * x + y * y,
//               x * x + y * y };
//   }
//
//   vector<real> u(real x, real y)
//   {
//      return { x * x + y * y,
//               x * x + y * y,
//               x * x + y * y };
//   }
//};

//class Test
//{
//public:
//   vector<real> f(real x, real y)
//   {
//      return  { -4 + 1 * (x * x + y * y),
//                -4 + 1 * (x * x + y * y),
//                -4 + 1 * (x * x + y * y) };
//   }
//
//   vector<real> lambda()
//   {
//      return { 1,
//               1,
//               1 };
//   }
//
//   vector<real> gamma()
//   {
//      return { 1,
//               1,
//               1 };
//   }
//
//   vector<real> ug(real x, real y)
//   {
//      return { x * x + y * y,
//               x * x + y * y,
//               x * x + y * y };
//   }
//
//   vector<real> u(real x, real y)
//   {
//      return { x * x + y * y,
//               x * x + y * y,
//               x * x + y * y };
//   }
//};

//class Test
//{
//public:
//   vector<real> f(real x, real y)
//   {
//      return  { 0,
//                0,
//                0 };
//   }
//
//   vector<real> lambda()
//   {
//      return { 1,
//               1,
//               1 };
//   }
//
//   vector<real> gamma()
//   {
//      return { 0,
//               0,
//               0 };
//   }
//
//   vector<real> ug(real x, real y)
//   {
//      return { 5,
//               6,
//               7 };
//   }
//
//   vector<real> u(real x, real y)
//   {
//      return { 5,
//               6,
//               7 };
//   }
//};

//class Test
//{
//public:
//   vector<real> f(real x, real y)
//   {
//      return  { x + y,
//                x + y,
//                x + y };
//   }
//
//   vector<real> lambda()
//   {
//      return { 1, 1, 1 };
//   }
//
//   vector<real> gamma()
//   {
//      return { 1, 1, 1 };
//   }
//
//   vector<real> ug(real x, real y)
//   {
//      return { x + y,
//               x + y,
//               x + y };
//   }
//
//   vector<real> u(real x, real y)
//   {
//      return { x + y,
//               x + y,
//               x + y };
//   }
//};
