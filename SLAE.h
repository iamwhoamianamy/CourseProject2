#pragma once
#include "Vector.h"
#include "Matrix.h"
using namespace std;

class SLAE
{
public:
   int maxiter;           // Максимальное количество итераций
   double eps;            // Велечина требуемой относительной невязки

   vector<double> b;      // Вектор правой части
   vector<double> t;      // Вспомогательный вектор для МСГ
   vector<double> tt;     // Вспомогательный вектор для МСГ
   vector<double> rk1;    // Вектор невязки на перд. итерации МСГ
   vector<double> zk1;    // Вектор спуска на пред. итерации МСГ
   vector<double> AtAzk1; // Вспомогательный вектор для МСГ

   SLAE(int size, int _maxiter, double _eps)
   {
      maxiter = _maxiter;
      eps = _eps;

      b.resize(size);

      t.resize(size);
      tt.resize(size);
      rk1.resize(size);
      zk1.resize(size);
      AtAzk1.resize(size);
   }

   SLAE()
   {

   }

   // Метод сопряженных градиентов, возвращает количество итераций
   int ConjGradMethod(vector<double>& xk1, vector<double>& res, Matrix& mat)
   {
      for(int i = 0; i < mat.size; i++)
         res[i] = xk1[i] = 0;

      mat.MatVecMult(xk1, t, mat.bot_tr, mat.top_tr);       // t = A * x0
      mat.MatVecMult(b - t, rk1, mat.top_tr, mat.bot_tr);  // r0 = AT(f - A * x0)
      zk1 = rk1;

      int k = 1;
      while(k < maxiter)
      {
         mat.MatVecMult(zk1, t, mat.bot_tr, mat.top_tr);    // t = A * zk-1

         mat.MatVecMult(t, AtAzk1, mat.top_tr, mat.bot_tr); // AtAzk1 = At * A * zk-1
         double ak = (rk1 * rk1) / (AtAzk1 * zk1);
         xk1 = xk1 + ak * zk1;
         double bk = rk1 * rk1;
         rk1 = rk1 - ak * AtAzk1;
         bk = (rk1 * rk1) / bk;
         zk1 = rk1 + bk * zk1;

         double disc = Norm(rk1) / Norm(b); // Относительная невязка

         if(disc < eps)
            break;
         else
            k++;
      }

      res = xk1;
      return k;
   }

   // Метод сопряженных градиентов с предобусловденной неполной
   // факторизацией матрицей, возвращает количество итераций
   int ConjGradPredMethod(vector<double>& xk1, vector<double>& res, Matrix& mat, SLAE& fac_slae, Matrix& fac_mat)
   {
      for(int i = 0; i < mat.size; i++)
         res[i] = xk1[i] = 0;

      mat.MatVecMult(xk1, t, mat.bot_tr, mat.top_tr);       // t = A * x0
      mat.MatVecMult(b - t, rk1, mat.top_tr, mat.bot_tr);  // r0 = AT(f - A * x0)

      // Решаем z0 = M-1 * r0
      fac_slae.b = rk1;
      fac_slae.ConjGradMethod(t, zk1, fac_mat);

      int k = 1;
      while(k < maxiter)
      {
         // Решаем tt = M-1 * rk-1
         fac_slae.b = rk1;
         fac_slae.ConjGradMethod(t, tt, fac_mat);

         mat.MatVecMult(zk1, t, mat.bot_tr, mat.top_tr);    // t = A * zk-1
         mat.MatVecMult(t, AtAzk1, mat.top_tr, mat.bot_tr); // AtAzk1 = At * A * zk-1

         double ak = (tt * rk1) / (AtAzk1 * zk1);
         xk1 = xk1 + ak * zk1;
         double bk = tt * rk1;
         rk1 = rk1 - ak * AtAzk1;

         // Решаем tt = M-1 * rk
         fac_slae.b = rk1;
         fac_slae.ConjGradMethod(t, tt, fac_mat);

         bk = (tt * rk1) / bk;
         zk1 = tt + bk * zk1;

         double disc = Norm(rk1) / Norm(b);

         if(disc < eps)
            break;
         else
            k++;
      }

      res = xk1;
      return k;
   }
};
