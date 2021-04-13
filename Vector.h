#pragma once
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;

// Умножение вектора на число
vector<double> operator * (const double& val, vector<double> vec)
{
   for (size_t i = 0; i < vec.size(); i++)
      vec[i] *= val;

   return vec;
}

// Сложение векторов
vector<double> operator + (vector<double> vec1, const vector<double>& vec2)
{
   for (size_t i = 0; i < vec1.size(); i++)
      vec1[i] += vec2[i];
   return vec1;
}

// Вычитание векторов
vector<double> operator - (vector<double> vec1, const vector<double>& vec2)
{
   for(size_t i = 0; i < vec1.size(); i++)
      vec1[i] -= vec2[i];
   return vec1;
}

double ScalarMult(const vector<double>& vec1, const vector<double>& vec2)
{
   double res = 0;
   for(int i = 0; i < vec1.size(); i++)
      res += vec1[i] * vec2[i];
   return res;
}

double operator * (const vector<double>& vec1, const vector<double>& vec2)
{
   return ScalarMult(vec1, vec2);
}

double Norm(const vector<double>& vec)
{
   return sqrt(ScalarMult(vec, vec));
}

