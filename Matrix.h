#pragma once
#include <vector>
#include <fstream>
using namespace std;

class Matrix
{
public:
   int size;                   // Размер матрицы
   int tr_size;                // Количество элементов в треугольнике
   
   vector<int> ind;            // Указатели начала строк
   vector<int> columns_ind;    // Номера столбцов внедиагональных элементов

   vector<double> top_tr;      // Верхний треугольник
   vector<double> bot_tr;      // Нижний треугольник
                               
   vector<double> diag;        // Диагональ

   Matrix(const int& t_size, const int& t_tr_size) : size(t_size), tr_size(t_tr_size) 
   {
      top_tr = vector<double>(tr_size);
      bot_tr = vector<double>(tr_size);
      columns_ind = vector<int>(tr_size);
      diag = vector<double>(size);
      ind = vector<int>(size + 1);
   }

   Matrix(const int& t_size) : size(t_size)
   {
      tr_size = size * (size - 1) / 2;

      top_tr.resize(tr_size);
      bot_tr.resize(tr_size);
      columns_ind.resize(tr_size);
      diag.resize(size);
      ind.resize(size + 1);

      ind[0] = ind[1] = 0;

      diag[0] = 1.0;

      for (int i = 1; i < size; i++)
      {
         int i0 = ind[i + 0];
         int i1 = ind[i + 1] = ind[i + 0] + i;

         for (int j = 0, k = i0; j < i; j++, k++)
            columns_ind[k] = j;
      }
   }

   Matrix(const Matrix& mat)
   {
      size = mat.size;
      tr_size = mat.tr_size;

      top_tr = mat.top_tr;
      bot_tr = mat.bot_tr;
      diag = mat.diag;
      ind = mat.ind;
      columns_ind = mat.columns_ind;
   }

   Matrix()
   {

   }

   // Получение диагональной факторизации матрицы
   void DiagFact(Matrix& fact)
   {
      fact.diag = diag;

      for (int i = 0; i < size + 1; i++)
         fact.ind[i] = 0;
   }

   // Функция умножения матрицы на вектор vec, результат в res
   void MatVecMult(const vector<double>& vec, vector<double>& res,
      const vector<double>& bot_tr, const vector<double>& top_tr)
   {
      for (int i = 0; i < size; i++)
         res[i] = 0;

      for (int i = 0; i < size; i++)
      {
         res[i] += vec[i] * diag[i];

         int prof_len = ind[i + 1] - ind[i];
         for (int k = 0; k < prof_len; k++)
         {
            int i_in_gg = ind[i] + k;
            int j = columns_ind[i_in_gg];
            res[i] += vec[j] * bot_tr[i_in_gg];
            res[j] += vec[i] * top_tr[i_in_gg];
         }
      }
   }

   // Функция считывания диагонали и треугольников из файлов
   void ReadDiTr(const string& file_name)
   {
      ifstream fin;
      fin.open(file_name);

      for(int i = 0; i < size; i++)
         fin >> diag[i];

      for(int i = 0; i < tr_size; i++)
      {
         double t;
         fin >> t;
         top_tr[i] = bot_tr[i] = t;
      }

      fin.close();
   }

   // Зануление всех элементов матрицы
   void ResetValues()
   {
      for(int i = 0; i < size; i++)
         diag[i] = 0;

      for(int i = 0; i < tr_size; i++)
      {
         bot_tr[i] = 0;
         top_tr[i] = 0;
      }
   }

   // Умножение матрицы на число
   void Mult(const double& val)
   {
      for(int i = 0; i < size; i++)
         diag[i] *= val;

      for(int i = 0; i < tr_size; i++)
      {
         bot_tr[i] *= val;
         top_tr[i] *= val;
      }
   }
};

//// Умножение матрицы на число
//Matrix operator * (double val, const Matrix& mat)
//{
//   Matrix res = Matrix(mat);
//
//   res.di = val * res.di;
//   res.ggl = val * res.ggl;
//   res.ggu = val * res.ggu;
//
//   return res;
//}
//
//// Сложение матриц
//Matrix operator + (const Matrix& mat1, const Matrix& mat2)
//{
//   Matrix res = Matrix(mat1);
//
//   res.di = res.di + mat2.di;
//   res.ggl = res.ggl + mat2.ggl;
//   res.ggu = res.ggu + mat2.ggu;
//
//   return res;
//}