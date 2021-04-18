#pragma once
#include "Matrix.h"
#include "SLAE.h"
#include "Test.h"
#include "Region.h"

using namespace std;

class BoundaryValueProblem
{
public:
   Matrix global;                 // Глобальная матрица
   Matrix fac_global;             // Неполная факторизация глобальной матрицы
   Matrix stiff_mat;              // Матрица жесткости
   Matrix weight_mat;             // Матрица масс
   Matrix G1, G2, M;              // Вспомогательные матрицы для вычисления элементов
                                  // локальных матриц и векторов правых частей

   SLAE slae;                     // Решатель системы без предобуславливания
   SLAE fac_slae;                 // Решатель системы c предобуславливанием

   vector<double> b;              // Глобальный вектор правой части
   vector<double> loc_b;          // Локальный вектор правой части
   vector<double> solution;       // Решение
   vector<double> true_solution;  // Точное решение для сравнения

   Test test;                     // Информация о значениях функции,
                                  // парматров лямбда и гамма

   int elems_count;               // Количество конечных элементов
   int nodes_count;               // Количество узлов

   vector<double> x_nodes;        // Координаты сетки по x
   vector<double> y_nodes;        // Координаты сетки по y
                                  
   int x_nodes_count;             // Количество узлов по x
   int y_nodes_count;             // Количество узлов по y

   vector<vector<int>> regions;   // Информация о подобластях (регионах)
   int regions_count;             // Количество регионов

   vector<int> x_cords_i;         // Индексы исходных координатных линий в векторе сетки по x
   vector<int> y_cords_i;         // Индексы исходных координатных линий в векторе сетки по y

   vector<vector<int>> boundaries; // Информация о краевых условиях

    // Считывание и формирование сетки из файла file_name
   void ReadFormGrid(const string& file_name)
   {
      ifstream fin(file_name);

      // Считывание координатных линий по x
      int x_coord_count;
      fin >> x_coord_count;

      vector<double> x_coords(x_coord_count);

      for(int i = 0; i < x_coord_count; i++)
         fin >> x_coords[i];

      // Формирование сетки по x
      x_nodes_count = 1;
      x_nodes = vector<double>(1);
      x_nodes[0] = x_coords[0];
      vector<int> n_x(x_coord_count - 1);

      for(int i = 0; i < x_coord_count - 1; i++)
      {
         fin >> n_x[i];
         x_nodes.resize(x_nodes_count + n_x[i]);

         double h = (x_coords[i + 1] - x_coords[i]) / n_x[i];

         for(int j = 0; j < n_x[i]; j++)
            x_nodes[j + x_nodes_count] = x_nodes[x_nodes_count - 1] + h * (j + 1);

         x_nodes_count += n_x[i];
      }

      // Считывание координатных линий по y
      int y_coord_count;
      fin >> y_coord_count;

      vector<double> y_coords(y_coord_count);

      for(int i = 0; i < y_coord_count; i++)
         fin >> y_coords[i];

      // Формирование сетки по y
      y_nodes_count = 1;
      y_nodes = vector<double>(1);
      y_nodes[0] = y_coords[0];
      vector<int> n_y(y_coord_count - 1);

      for(int i = 0; i < y_coord_count - 1; i++)
      {
         fin >> n_y[i];
         y_nodes.resize(y_nodes_count + n_y[i]);

         double h = (y_coords[i + 1] - y_coords[i]) / n_y[i];

         for(int j = 0; j < n_y[i]; j++)
            y_nodes[j + y_nodes_count] = y_nodes[y_nodes_count - 1] + h * (j + 1);

         y_nodes_count += n_y[i];

      }
      
      // Считывание информации о подобластях
      fin >> regions_count;
      regions = vector<vector<int>>(regions_count, vector<int>(4));

      for(int i = 0; i < regions_count; i++)
         fin >> regions[i][0] >> regions[i][1] >> regions[i][2] >> regions[i][3];

      fin.close();

      // Пересчет индексов границ подобластей
      for(int i = 1; i < n_x.size(); i++)
         n_x[i] += n_x[i - 1];

      for(int i = 1; i < n_y.size(); i++)
         n_y[i] += n_y[i - 1];

      x_cords_i = vector<int>(n_x.size() + 1);

      for(int i = 0; i < n_x.size(); i++)
         x_cords_i[i + 1] = n_x[i];

      y_cords_i = vector<int>(n_y.size() + 1);

      for(int i = 0; i < n_y.size(); i++)
         y_cords_i[i + 1] = n_y[i];

      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         regions[reg_i][0] = x_cords_i[regions[reg_i][0]];
         regions[reg_i][1] = x_cords_i[regions[reg_i][1]];
         regions[reg_i][2] = y_cords_i[regions[reg_i][2]];
         regions[reg_i][3] = y_cords_i[regions[reg_i][3]];
      }

      nodes_count = x_nodes_count * y_nodes_count;
      elems_count += n_x[n_x.size() - 1] / 2 * n_y[n_y.size() - 1] / 2;
   }

   // Считывание информации о краевых условиях из файла file_name
   void ReadBoundaries(const string& file_name)
   {
      ifstream fin(file_name);



      fin.close();
   }

   //// Считывание и формирование сетки из файла file_name
   //void ReadFormGrid(const string& file_name)
   //{
   //   ifstream fin(file_name);

   //   // Координаты границы сетки
   //   double left, right, bot, top;

   //   // Считывание границы границы
   //   fin >> left;
   //   fin >> right;
   //   fin >> bot;
   //   fin >> top;

   //   // Генерация координат узлов по x
   //   int n;
   //   double h, q;

   //   fin >> q >> n;

   //   x_nodes_count = n + 1;
   //   x_nodes.resize(x_nodes_count);

   //   h = right - left;

   //   if(q != 1)
   //      h *= (1 - q) / (1 - pow(q, n));
   //   else
   //      h /= n;

   //   x_nodes[0] = left;

   //   for(int i = 0; i < n; i++)
   //      x_nodes[i + 1] = x_nodes[i] + h * pow(q, i);

   //   elems_count = n / 2;

   //   // Генерация координат узлов по y
   //   fin >> q >> n;

   //   y_nodes_count = n + 1;
   //   y_nodes.resize(y_nodes_count);

   //   h = top - bot;

   //   if(q != 1)
   //      h *= (1 - q) / (1 - pow(q, n));
   //   else
   //      h /= n;

   //   y_nodes[0] = bot;

   //   for(int i = 0; i < n; i++)
   //      y_nodes[i + 1] = y_nodes[i] + h * pow(q, i);

   //   fin.close();

   //   nodes_count = x_nodes_count * y_nodes_count;
   //   elems_count *= n / 2;
   //}

   // Считывание вспомогательных матриц для формирования
   // матриц жесткости и массы
   void ReadMatrices()
   {
      // Чтение вспомогательных матриц для построения матриц жесткости и массы
      G1 = Matrix(9);
      G1.ReadDiTr("data/G1.txt");

      G2 = Matrix(9);
      G2.ReadDiTr("data/G2.txt");

      M = Matrix(9);
      M.ReadDiTr("data/M1.txt");
   }

   // Выделение памяти под массивы
   void InitializeMemory()
   {
      slae = SLAE(nodes_count, 10000, 1e-20);
      fac_slae = SLAE(nodes_count, 10000, 1e-20);

      fac_global = Matrix(nodes_count, 0);

      global.ind.resize(nodes_count + 1);
      b.resize(nodes_count);
      solution.resize(nodes_count);
      true_solution.resize(nodes_count);

      stiff_mat = Matrix(nodes_count);
      weight_mat = Matrix(nodes_count);
   }

   // Заполнение массива global_indices индексами, соответствующими глобальной номерации
   // узлов конечного элемента с номером elem_index(индексация с нуля)
   void CalcGlobalIndices(int elem_index, vector<int>& global_indices)
   {
      int n_coords = x_nodes_count / 2 + 1;
      int k = 2 * floor((elem_index) / (n_coords - 1)) * (2 * n_coords - 1) + 2 * ((elem_index) % (n_coords - 1)) + 1;
      k--;

      global_indices[0] = k + 0;
      global_indices[1] = k + 1;
      global_indices[2] = k + 2;

      global_indices[3] = k + 2 * n_coords - 1;
      global_indices[4] = k + 2 * n_coords - 0;
      global_indices[5] = k + 2 * n_coords + 1;

      global_indices[6] = k + 2 * (2 * n_coords - 1);
      global_indices[7] = k + 2 * (2 * n_coords - 1) + 1;
      global_indices[8] = k + 2 * (2 * n_coords - 1) + 2;
   }

   // Вспомогательная функция для формирования портрета
   void IncertToRow(const int& r, const int& c)
   {
      int i_in_jg = global.ind[r];
      int prof_len = global.ind[r + 1] - global.ind[r];

      bool found = false;

      for(int k = i_in_jg; k < i_in_jg + prof_len; k++)
         if(global.columns_ind[k] == c)
         {
            found = true;
            break;
         }

      if(!found)
      {
         for(int l = r + 1; l < global.ind.size(); l++)
            global.ind[l]++;

         int k = i_in_jg;

         while((k < i_in_jg + prof_len) && global.columns_ind[k] < c)
            k++;

         global.columns_ind.insert(global.columns_ind.begin() + k, c);
      }
   }

   // Формирование портрета глобальной матрицы
   void FormPortrait()
   {
      global.ind[0] = global.ind[1] = 0;

      for(int i = 0; i < elems_count; i++)
      {
         vector<int> global_indices(9);
         CalcGlobalIndices(i, global_indices);
         vector<vector<vector<int>>> help(9);

         for(int i = 0; i < 9; i++)
            help[i].resize(9);

         for(int i = 0; i < 9; i++)
            for(int j = 0; j < 9; j++)
               help[i][j] = { global_indices[i], global_indices[j] };


         for(int i = 0; i < 9; i++)
            for(int j = 0; j < i; j++)
               IncertToRow(help[i][j][0], help[i][j][1]);
      }

      global.size = global.ind.size() - 1;
      global.diag.resize(global.size);
      
      global.tr_size = global.columns_ind.size();
      global.bot_tr.resize(global.tr_size);
      global.top_tr.resize(global.tr_size);

      stiff_mat = global;
      weight_mat = global;
   }

   // Добавление элемента в матрицу
   void AddToMat(Matrix& mat, const int& row_i, const int& col_i, const double& val_l, const double& val_u)
   {
      int beg_prof = mat.ind[row_i];
      int end_prof = mat.ind[row_i + 1];

      for(int i_in_prof = beg_prof; i_in_prof < end_prof; i_in_prof++)
      {
         if(mat.columns_ind[i_in_prof] == col_i)
         {
            mat.bot_tr[i_in_prof] += val_l;
            mat.top_tr[i_in_prof] += val_u;
            break;
         }
      }
   }

   // Находит координаты узлов конечного элемента
   // с номером elem_index(индексация с нуля)
   void CalcElemNodes(int elem_index, vector<double>& x_nodes_elem, vector<double>& y_nodes_elem)
   {
      int n_coords = x_nodes_count / 2 + 1;
      int x0 = (elem_index) % (n_coords - 1) * 2;
      int y0 = floor((elem_index) / (n_coords - 1)) * 2;

      x_nodes_elem[0] = x_nodes[x0];
      x_nodes_elem[1] = x_nodes[x0 + 1];
      x_nodes_elem[2] = x_nodes[x0 + 2];

      y_nodes_elem[0] = y_nodes[y0];
      y_nodes_elem[1] = y_nodes[y0 + 1];
      y_nodes_elem[2] = y_nodes[y0 + 2];
   }

   // Сборка матриц жесткости и массы
   void BuildMatrices()
   {
      vector<int> global_indices(9);

      for(int elem_i = 0; elem_i < elems_count; elem_i++)
      {
         CalcGlobalIndices(elem_i, global_indices);

         vector<double> x_nodes_elem(3);       // Координаты конечного элемента по x
         vector<double> y_nodes_elem(3);       // Координаты конечного элемента по y

         CalcElemNodes(elem_i, x_nodes_elem, y_nodes_elem);

         double hx = x_nodes_elem[2] - x_nodes_elem[0];
         double hy = y_nodes_elem[2] - y_nodes_elem[0];

         vector<double> local_f(9);
         
         for(int i = 0; i < 9; i++)
         {
            stiff_mat.diag[global_indices[i]] += (test.lambda() / 90.0) * (hy / hx * G1.diag[i] + hx / hy * G2.diag[i]);
            weight_mat.diag[global_indices[i]] += (test.sigma() * hx * hy / 900.0) * M.diag[i];

            local_f[i] = test.f(x_nodes_elem[i % 3], y_nodes_elem[floor(i / 3)], 0);
            true_solution[global_indices[i]] = test.u(x_nodes_elem[i % 3], y_nodes_elem[floor(i / 3)], 0);

            int beg_prof = M.ind[i];
            int end_prof = M.ind[i + 1];

            for(int i_in_prof = beg_prof; i_in_prof < end_prof; i_in_prof++)
            {
               int j = M.columns_ind[i_in_prof];

               double val_l = (test.lambda() / 90.0) * (hy / hx * G1.bot_tr[i_in_prof] + hx / hy * G2.bot_tr[i_in_prof]);
               double val_u = (test.lambda() / 90.0) * (hy / hx * G1.top_tr[i_in_prof] + hx / hy * G2.top_tr[i_in_prof]);

               AddToMat(stiff_mat, global_indices[i], global_indices[j], val_l, val_u);

               val_l = (test.sigma() * hx * hy / 900.0) * M.bot_tr[i_in_prof];
               val_u = (test.sigma() * hx * hy / 900.0) * M.top_tr[i_in_prof];

               AddToMat(weight_mat, global_indices[i], global_indices[j], val_l, val_u);
            }
         }

         vector<double> local_b(9);
         M.MatVecMult(local_f, local_b, M.bot_tr, M.top_tr);
         for(int i = 0; i < 9; i++)
            b[global_indices[i]] += hx * hy / 900.0 * local_b[i];
      }
   }

   // Сборка глобальной матрицы
   void AssembleGlobalMatrix()
   {
      for(int i = 0; i < nodes_count; i++)
         global.diag[i] = stiff_mat.diag[i] + weight_mat.diag[i];

      for(int i = 0; i < global.tr_size; i++)
      {
         global.bot_tr[i] = stiff_mat.bot_tr[i] + weight_mat.bot_tr[i];
         global.top_tr[i] = stiff_mat.top_tr[i] + weight_mat.top_tr[i];
      }
   }

   void FirstBoundOnLine(const int& line_i, const double& x, const double& y, const double& t)
   {
      global.diag[line_i] = 1;

      b[line_i] = test.u(x, y, t);

      for(int prof_i = global.ind[line_i]; prof_i < global.ind[line_i + 1]; prof_i++)
      {
         global.bot_tr[prof_i] = 0;
      }

      for(int i = 0; i < nodes_count; i++)
      {
         for(int prof_i = global.ind[i]; prof_i < global.ind[i + 1]; prof_i++)
         {
            if(global.columns_ind[prof_i] == line_i)
               global.top_tr[prof_i] = 0;
         }
      }
   }

   // Учет первых краевых условий
   void AccountFirstBound()
   {
      for(int x_i = 0; x_i < x_nodes_count; x_i++)
      {
         for(int y_i = 0; y_i < y_nodes_count; y_i++)
         {
            if(x_i == 0 || x_i == x_nodes_count - 1 ||  y_i == 0 || y_i == y_nodes_count - 1)
               FirstBoundOnLine(x_i + y_i * x_nodes_count, x_nodes[x_i], y_nodes[y_i], 0);
         }
      }
   }

   // Нахождение решения
   void Solve()
   {
      slae.b = b;
      global.DiagFact(fac_global);

      vector<double> x0(nodes_count);
      cout << slae.ConjGradPredMethod(x0, solution, global, fac_slae, fac_global) << endl;
      //cout << slae.ConjGradMethod(x0, solution, global) << endl;
   }

   void PrintSolution(ofstream& fout)
   {
      fout << setw(14) << "x" << setw(14) << "y";
      fout << setw(14) << "prec" << setw(14) << "calc" << setw(14) << "diff" << setw(5) << "n" << " loc" << endl;

      for(int y_i = 0; y_i < y_nodes_count; y_i++)
      {
         for(int x_i = 0; x_i < x_nodes_count; x_i++)
         {
            int i = x_i + y_i * x_nodes_count;

            fout << scientific;
            fout << setw(14) << x_nodes[x_i];
            fout << setw(14) << y_nodes[y_i];
            fout << setw(14) << true_solution[i];
            fout << setw(14) << solution[i];
            fout << setw(14) << abs(true_solution[i] - solution[i]);
            fout << fixed << setw(5) << i;

            if(x_i == 0 || x_i == x_nodes_count - 1 || y_i == 0 || y_i == y_nodes_count - 1)
               fout << " border";
            else
               fout << " inner";

            fout << endl;
         }
      }
   }
};