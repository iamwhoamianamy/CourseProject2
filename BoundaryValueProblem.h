#pragma once
#include <iomanip>
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
   Matrix sigma_mass_mat;         // Матрица массы для параметра Сигма
   Matrix chi_mass_mat;           // Матрица массы для параметра Хи
   Matrix  M;                     // Вспомогательные матрицы для вычисления элементов
   vector<Matrix> Gl, Gr;         // локальных матриц и векторов правых частей
   
   SLAE slae;                     // Решатель системы без предобуславливания
   SLAE fac_slae;                 // Решатель системы c предобуславливанием

   vector<double> b;              // Глобальный вектор правой части
   vector<double> loc_b;          // Локальный вектор правой части
   vector<double> solution;       // Решение
   vector<double> prev_solution1; // Решение на предыдущей временной итерации
   vector<double> prev_solution2; // Решение на предпредыдущей временной итерации
   vector<double> prev_solution3; // Решение на предпредпредыдущей временной итерации
   vector<double> true_solution;  // Точное решение для сравнения
   vector<int> location;          // Положение на сетке для каждого узла

   int elems_count;               // Количество конечных элементов
   int nodes_count;               // Количество узлов

   vector<double> x_nodes;        // Координаты сетки по x
   vector<double> y_nodes;        // Координаты сетки по y
                                  
   int x_nodes_count;             // Количество узлов по x
   int y_nodes_count;             // Количество узлов по y

   vector<Region> regions;        // Информация о подобластях (регионах)
   int regions_count;             // Количество регионов

   vector<int> x_cords_i;         // Индексы исходных координатных линий в векторе сетки по x
   vector<int> y_cords_i;         // Индексы исходных координатных линий в векторе сетки по y

   vector<vector<int>> boundaries; // Информация о краевых условиях
   int bound_count;                // Количество краевых условий

   vector<double> time_grid;       // Сетка по времени
   int timesteps_count;            // Количество временных слоев

   vector<double> Mnqj_1;          // Вспомогательный вектор для итерации по времени
   vector<double> Mnqj_2;          // Вспомогательный вектор для итерации по времени
   vector<double> Mnqj_3;          // Вспомогательный вектор для итерации по времени

   vector<double> Moqj_1;          // Вспомогательный вектор для итерации по времени
   vector<double> Moqj_2;          // Вспомогательный вектор для итерации по времени
   vector<double> Moqj_3;          // Вспомогательный вектор для итерации по времени

   vector<double> xk;              // Вспомогательный вектор для МСГ

   // Вспомогательный массив для вычисления индексов узлов в глобальной нумерации
   vector<vector<vector<int>>> help;

   // Считывание и формирование сетки по времени из файла file_name
   void ReadFormTimeGrid(const string& file_name)
   {
      ifstream fin(file_name);

      double t_first, t_last;
      fin >> t_first >> t_last >> timesteps_count;
      timesteps_count++;

      time_grid = vector<double>(timesteps_count);
      double h = (t_last - t_first) / (timesteps_count - 1);

      time_grid[0] = t_first;
      for(int i = 1; i < timesteps_count; i++)
         time_grid[i] = h * i;

      fin.close();
   }

   // Считывание и формирование сетки из файла file_name
   void ReadFormGrid(const string& file_name)
   {
      ifstream fin(file_name);
      
      // Переменная для считывания разделителей
      string fake;
      fin >> fake;

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

      fin >> fake;

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

      fin >> fake;
      
      // Считывание информации о подобластях
      fin >> regions_count;
      regions = vector<Region>(regions_count);

      for(int i = 0; i < regions_count; i++)
      {
         fin >> regions[i].left >> regions[i].right >> regions[i].bot >> regions[i].top;

         int test_n;
         fin >> test_n;

         regions[i].test = Test(test_n);
         fin >> regions[i].test.lambda_n >> regions[i].test.sigma_n >> regions[i].test.chi_n;
      }

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
         regions[reg_i].left = x_cords_i[regions[reg_i].left];
         regions[reg_i].right = x_cords_i[regions[reg_i].right];
         regions[reg_i].bot = x_cords_i[regions[reg_i].bot];
         regions[reg_i].top = x_cords_i[regions[reg_i].top];
      }

      nodes_count = x_nodes_count * y_nodes_count;
      elems_count += n_x[n_x.size() - 1] / 2 * n_y[n_y.size() - 1] / 2;
   }

   // Считывание информации о краевых условиях из файла file_name
   void ReadBoundaries(const string& file_name)
   {
      ifstream fin(file_name);

      fin >> bound_count;

      boundaries = vector<vector<int>>(bound_count, vector<int>(5));

      for(int bound_i = 0; bound_i < bound_count; bound_i++)
      {
         for(int i = 0; i < 5; i++)
            fin >> boundaries[bound_i][i];

         boundaries[bound_i][1] = x_cords_i[boundaries[bound_i][1]];
         boundaries[bound_i][2] = x_cords_i[boundaries[bound_i][2]];
         boundaries[bound_i][3] = y_cords_i[boundaries[bound_i][3]];
         boundaries[bound_i][4] = y_cords_i[boundaries[bound_i][4]];
      }

      fin.close();
   }

   // Считывание вспомогательных матриц для формирования
   // матриц жесткости и массы
   void ReadMatrices()
   {
      M = Matrix(9);
      M.ReadDiTr("data/M1.txt");

      Gl = vector<Matrix>(4, Matrix(9));
      Gl[0].ReadDiTr("data/Gl0.txt");
      Gl[1].ReadDiTr("data/Gl1.txt");
      Gl[2].ReadDiTr("data/Gl2.txt");
      Gl[3].ReadDiTr("data/Gl3.txt");

      Gr = vector<Matrix>(4, Matrix(9));
      Gr[0].ReadDiTr("data/Gr0.txt");
      Gr[1].ReadDiTr("data/Gr1.txt");
      Gr[2].ReadDiTr("data/Gr2.txt");
      Gr[3].ReadDiTr("data/Gr3.txt");

   }

   // Выделение памяти под массивы
   void InitializeMemory()
   {
      slae = SLAE(nodes_count, 10000, 1e-20);
      fac_slae = SLAE(nodes_count, 10000, 1e-20);

      global.ind = vector<int>(nodes_count + 1);
      b = vector<double>(nodes_count);

      solution = vector<double>(nodes_count);
      prev_solution1 = vector<double>(nodes_count);
      prev_solution2 = vector<double>(nodes_count);
      prev_solution3 = vector<double>(nodes_count);
      true_solution = vector<double>(nodes_count);

      location = vector<int>(nodes_count);

      Mnqj_1 = vector<double>(nodes_count);
      Mnqj_2 = vector<double>(nodes_count);
      Mnqj_3 = vector<double>(nodes_count);

      Moqj_1 = vector<double>(nodes_count);
      Moqj_2 = vector<double>(nodes_count);
      Moqj_3 = vector<double>(nodes_count);

      stiff_mat = Matrix(nodes_count);
      sigma_mass_mat = Matrix(nodes_count);
      chi_mass_mat = Matrix(nodes_count);
      fac_global = Matrix(nodes_count, 0);

      xk = vector<double>(nodes_count);

      help = vector<vector<vector<int>>>(9);

      for(int i = 0; i < 9; i++)
         help[i].resize(9);
   }

   // Заполнение массива global_indices индексами, соответствующими глобальной номерации
   // узлов конечного элемента с номером elem_i(индексация с нуля)
   void CalcGlobalIndices(int elem_i, vector<int>& global_indices)
   {
      int n_coords = x_nodes_count / 2 + 1;
      int k = 2 * floor((elem_i) / (n_coords - 1)) * (2 * n_coords - 1) + 2 * ((elem_i) % (n_coords - 1));

      global_indices[0] = k + 0;
      global_indices[1] = k + 1;
      global_indices[2] = k + 2;

      global_indices[3] = k + 2 * n_coords - 1;
      global_indices[4] = k + 2 * n_coords - 0;
      global_indices[5] = k + 2 * n_coords + 1;

      global_indices[6] = k + 2 * (2 * n_coords - 1) + 0;
      global_indices[7] = k + 2 * (2 * n_coords - 1) + 1;
      global_indices[8] = k + 2 * (2 * n_coords - 1) + 2;
   }

   // Поиск региона по номеру конечного элемента
   int CalcRegionIndex(const int& elem_i)
   {
      int n_coords = x_nodes_count / 2 + 1;
      int x_i = (elem_i) % (n_coords - 1) * 2 + 1;
      int y_i = floor((elem_i) / (n_coords - 1)) * 2 + 1;

      return CalcRegionIndex(x_i, y_i);
   }

   // Поиск региона по индексам в сетке
   int CalcRegionIndex(const int& x_i, const int& y_i)
   {
      int found_reg_i = -1;

      for(int reg_i = 0; reg_i < regions_count; reg_i++)
      {
         if(regions[reg_i].left <= x_i && x_i <= regions[reg_i].right &&
            regions[reg_i].bot <= y_i && y_i <= regions[reg_i].top)
         {
            found_reg_i = reg_i;
            break;
         }
      }
      return found_reg_i;
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

      for(int elem_i = 0; elem_i < elems_count; elem_i++)
      {
         int reg_i = CalcRegionIndex(elem_i);

         if(reg_i != -1)
         {
            vector<int> global_indices(9);
            CalcGlobalIndices(elem_i, global_indices);

            for(int i = 0; i < 9; i++)
               for(int j = 0; j < 9; j++)
                  help[i][j] = { global_indices[i], global_indices[j] };


            for(int i = 0; i < 9; i++)
               for(int j = 0; j < i; j++)
                  IncertToRow(help[i][j][0], help[i][j][1]);
         }
      }

      global.size = global.ind.size() - 1;
      global.diag.resize(global.size);
      
      global.tr_size = global.columns_ind.size();
      global.bot_tr.resize(global.tr_size);
      global.top_tr.resize(global.tr_size);

      stiff_mat = global;
      sigma_mass_mat = global;
      chi_mass_mat = global;
   }

   // Добавление элемента в матрицу в строку row_i столбец col_i
   // val_l добавляется в нижний треугольник, val_u добавляется в верхний треугольник
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
   // с номером elem_i(индексация с нуля)
   void CalcElemNodes(const int& elem_i, vector<double>& x_nodes_elem, vector<double>& y_nodes_elem)
   {
      int n_coords = x_nodes_count / 2 + 1;
      int x0 = (elem_i) % (n_coords - 1) * 2;
      int y0 = floor((elem_i) / (n_coords - 1)) * 2;

      x_nodes_elem[0] = x_nodes[x0];
      x_nodes_elem[1] = x_nodes[x0 + 1];
      x_nodes_elem[2] = x_nodes[x0 + 2];

      y_nodes_elem[0] = y_nodes[y0];
      y_nodes_elem[1] = y_nodes[y0 + 1];
      y_nodes_elem[2] = y_nodes[y0 + 2];
   }

   // Сборка матриц жесткости и массы
   void BuildMatrices(const double& t)
   {
      vector<int> global_indices(9);

      for(int elem_i = 0; elem_i < elems_count; elem_i++)
      {
         int reg_i = CalcRegionIndex(elem_i);
         CalcGlobalIndices(elem_i, global_indices);

         if(reg_i == -1)
         {
            for(int i = 0; i < 9; i++)
            {
               stiff_mat.diag[global_indices[i]] = 1;
               sigma_mass_mat.diag[global_indices[i]] = 0;
               chi_mass_mat.diag[global_indices[i]] = 0;
               b[global_indices[i]] = 0;
               location[global_indices[i]] = 2;
            }
         }
         else
         {
            Region r = regions[reg_i];
            vector<double> x_nodes_elem(3);       // Координаты конечного элемента по x
            vector<double> y_nodes_elem(3);       // Координаты конечного элемента по y
            CalcElemNodes(elem_i, x_nodes_elem, y_nodes_elem);

            double hx = x_nodes_elem[2] - x_nodes_elem[0];
            double hy = y_nodes_elem[2] - y_nodes_elem[0];

            vector<double> local_f(9);

            vector<double> lambda {
                  r.test.lambda(x_nodes_elem[0], y_nodes_elem[0]),
                  r.test.lambda(x_nodes_elem[0], y_nodes_elem[2]),
                  r.test.lambda(x_nodes_elem[2], y_nodes_elem[0]),
                  r.test.lambda(x_nodes_elem[2], y_nodes_elem[2]) };

            for(int i = 0; i < 9; i++)
            {
               double x = x_nodes_elem[i % 3];
               double y = y_nodes_elem[floor(i / 3)];

               stiff_mat.diag[global_indices[i]] += (1.0 / 90.0) * 
                  (hy / hx * (Gl[0].diag[i] * lambda[0] + Gl[1].diag[i] * lambda[1] + Gl[2].diag[i] * lambda[2] + Gl[3].diag[i] * lambda[3]) +
                   hx / hy * (Gr[0].diag[i] * lambda[0] + Gr[1].diag[i] * lambda[1] + Gr[2].diag[i] * lambda[2] + Gr[3].diag[i] * lambda[3]));

               sigma_mass_mat.diag[global_indices[i]] += (r.test.sigma() * hx * hy / 900.0) * M.diag[i];
               chi_mass_mat.diag[global_indices[i]] += (r.test.chi() * hx * hy / 900.0) * M.diag[i];

               local_f[i] = r.test.f(x, y, t);
               true_solution[global_indices[i]] = r.test.u(x, y, t);

               int beg_prof = M.ind[i];
               int end_prof = M.ind[i + 1];

               for(int i_in_prof = beg_prof; i_in_prof < end_prof; i_in_prof++)
               {
                  int j = M.columns_ind[i_in_prof];

                  double val_l = (1.0 / 90.0) * 
                     (hy / hx * (Gl[0].bot_tr[i_in_prof] * lambda[0] + Gl[1].bot_tr[i_in_prof] * lambda[1] + Gl[2].bot_tr[i_in_prof] * lambda[2] + Gl[3].bot_tr[i_in_prof] * lambda[3]) +
                      hx / hy * (Gr[0].bot_tr[i_in_prof] * lambda[0] + Gr[1].bot_tr[i_in_prof] * lambda[1] + Gr[2].bot_tr[i_in_prof] * lambda[2] + Gr[3].bot_tr[i_in_prof] * lambda[3]));
                  
                  double val_u = (1.0 / 90.0) * 
                     (hy / hx * (Gl[0].top_tr[i_in_prof] * lambda[0] + Gl[1].top_tr[i_in_prof] * lambda[1] + Gl[2].top_tr[i_in_prof] * lambda[2] + Gl[3].top_tr[i_in_prof] * lambda[3]) +
                      hx / hy * (Gr[0].top_tr[i_in_prof] * lambda[0] + Gr[1].top_tr[i_in_prof] * lambda[1] + Gr[2].top_tr[i_in_prof] * lambda[2] + Gr[3].top_tr[i_in_prof] * lambda[3]));

                  AddToMat(stiff_mat, global_indices[i], global_indices[j], val_l, val_u);

                  val_l = (r.test.sigma() * hx * hy / 900.0) * M.bot_tr[i_in_prof];
                  val_u = (r.test.sigma() * hx * hy / 900.0) * M.top_tr[i_in_prof];

                  AddToMat(sigma_mass_mat, global_indices[i], global_indices[j], val_l, val_u);

                  val_l = (r.test.chi() * hx * hy / 900.0) * M.bot_tr[i_in_prof];
                  val_u = (r.test.chi() * hx * hy / 900.0) * M.top_tr[i_in_prof];

                  AddToMat(chi_mass_mat, global_indices[i], global_indices[j], val_l, val_u);
               }
            }

            vector<double> local_b(9);
            M.MatVecMult(local_f, local_b, M.bot_tr, M.top_tr);
            for(int i = 0; i < 9; i++)
               b[global_indices[i]] += hx * hy / 900.0 * local_b[i];
         }
      }
   }

   // Сборка глобальной матрицы
   void AssembleGlobalMatrix()
   {
      for(int i = 0; i < nodes_count; i++)
         global.diag[i] = stiff_mat.diag[i] + sigma_mass_mat.diag[i] + chi_mass_mat.diag[i];

      for(int i = 0; i < global.tr_size; i++)
      {
         global.bot_tr[i] = stiff_mat.bot_tr[i] + sigma_mass_mat.bot_tr[i] + chi_mass_mat.bot_tr[i];
         global.top_tr[i] = stiff_mat.top_tr[i] + sigma_mass_mat.top_tr[i] + chi_mass_mat.top_tr[i];
      }
   }

   // Учет первых краевых условий u в точке (x,y) на временном слое t на строке с номером line_i
   void FirstBoundOnLine(const int& line_i, const double& x, const double& y, const double& t, const double& u)
   {
      global.diag[line_i] = 1;
      true_solution[line_i] = u;
      b[line_i] = u;

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

   // Учет первых краевых условий в СЛАУ на временном слое t
   void AccountFirstBound(const double& t)
   {
      for(int x_i = 0; x_i < x_nodes_count; x_i++)
      {
         for(int y_i = 0; y_i < y_nodes_count; y_i++)
         {
            int reg_i = CalcRegionIndex(x_i, y_i);

            if(reg_i != -1)
            {
               for(int bound_i = 0; bound_i < bound_count; bound_i++)
               {
                  if(boundaries[bound_i][1] <= x_i && x_i <= boundaries[bound_i][2] &&
                     boundaries[bound_i][3] <= y_i && y_i <= boundaries[bound_i][4])
                  {
                     int i = x_i + y_i * x_nodes_count;
                     FirstBoundOnLine(i, x_nodes[x_i], y_nodes[y_i], t, regions[reg_i].test.u(x_nodes[x_i], y_nodes[y_i], t));
                     location[i] = 1;
                     break;
                  }
               }
            }
         }
      }
   }

   // Нахождение решения СЛАУ
   void Solve()
   {
      slae.b = b;
      global.DiagFact(fac_global);

      cout << slae.ConjGradPredMethod(xk, solution, global, fac_slae, fac_global) << endl;
      //slae.ConjGradMethod(xk, solution, global);
   }

   // Итерации по времени, вывод в поток fout
   void TimeIterations(ofstream& fout)
   {
      for(int y_i = 0; y_i < y_nodes_count; y_i++)
      {
         for(int x_i = 0; x_i < x_nodes_count; x_i++)
         {
            int reg_i = CalcRegionIndex(x_i, y_i);
            if(reg_i != -1)
            {
              prev_solution3[x_i + y_i * x_nodes_count] =
                 regions[reg_i].test.u(x_nodes[x_i], y_nodes[y_i], time_grid[0]);

               prev_solution2[x_i + y_i * x_nodes_count] =
                 regions[reg_i].test.u(x_nodes[x_i], y_nodes[y_i], time_grid[1]);

               prev_solution1[x_i + y_i * x_nodes_count] =
                 regions[reg_i].test.u(x_nodes[x_i], y_nodes[y_i], time_grid[2]);
            }
         }
      }
      for(int time_i = 3; time_i < timesteps_count; time_i++)
      {
         // Временные слои
         double t3 = time_grid[time_i - 3];
         double t2 = time_grid[time_i - 2];
         double t1 = time_grid[time_i - 1];
         double t0 = time_grid[time_i - 0];

         // Числители дробей
         double den3 = (t3 - t2) * (t3 - t1) * (t3 - t0);
         double den2 = (t2 - t3) * (t2 - t1) * (t2 - t0);
         double den1 = (t1 - t3) * (t1 - t2) * (t1 - t0);
         double den0 = (t0 - t3) * (t0 - t2) * (t0 - t1);

         // Первые произодные
         double n3 = (3 * t0 * t0 - 2 * t0 * t2 - 2 * t0 * t1 - 2 * t0 * t0 +
                      t2 * t1 + t2 * t0 + t1 * t0) / den3;

         double n2 = (3 * t0 * t0 - 2 * t0 * t3 - 2 * t0 * t1 - 2 * t0 * t0 +
                      t3 * t1 + t3 * t0 + t1 * t0) / den2;

         double n1 = (3 * t0 * t0 - 2 * t0 * t3 - 2 * t0 * t2 - 2 * t0 * t0 +
                      t3 * t2 + t3 * t0 + t2 * t0) / den1;

         double n0 = (3 * t0 * t0 - 2 * t0 * t3 - 2 * t0 * t2 - 2 * t0 * t1 +
                      t3 * t2 + t3 * t1 + t2 * t1) / den0;

         // Вторые производные
         double o3 = (6 * t0 - 2 * t2 - 2 * t1 - 2 * t0) / den3;
         double o2 = (6 * t0 - 2 * t3 - 2 * t1 - 2 * t0) / den2;
         double o1 = (6 * t0 - 2 * t3 - 2 * t2 - 2 * t0) / den1;
         double o0 = (6 * t0 - 2 * t3 - 2 * t2 - 2 * t1) / den0;

         sigma_mass_mat.ResetValues();
         chi_mass_mat.ResetValues();
         stiff_mat.ResetValues();
         global.ResetValues();

         for(int i = 0; i < nodes_count; i++)
            b[i] = 0;

         BuildMatrices(t0);

         sigma_mass_mat.MatVecMult(prev_solution1, Mnqj_1, sigma_mass_mat.bot_tr, sigma_mass_mat.top_tr);
         sigma_mass_mat.MatVecMult(prev_solution2, Mnqj_2, sigma_mass_mat.bot_tr, sigma_mass_mat.top_tr);
         sigma_mass_mat.MatVecMult(prev_solution3, Mnqj_3, sigma_mass_mat.bot_tr, sigma_mass_mat.top_tr);

         chi_mass_mat.MatVecMult(prev_solution1, Moqj_1, chi_mass_mat.bot_tr, chi_mass_mat.top_tr);
         chi_mass_mat.MatVecMult(prev_solution2, Moqj_2, chi_mass_mat.bot_tr, chi_mass_mat.top_tr);
         chi_mass_mat.MatVecMult(prev_solution3, Moqj_3, chi_mass_mat.bot_tr, chi_mass_mat.top_tr);

         for(int i = 0; i < nodes_count; i++)
         {
            Mnqj_1[i] *= n1;
            Mnqj_2[i] *= n2;
            Mnqj_3[i] *= n3;

            Moqj_1[i] *= o1;
            Moqj_2[i] *= o2;
            Moqj_3[i] *= o3;
         }

         sigma_mass_mat.Mult(n0);
         chi_mass_mat.Mult(o0);

         AssembleGlobalMatrix();

         for(int i = 0; i < nodes_count; i++)
            b[i] -= Mnqj_1[i] + Mnqj_2[i] + Mnqj_3[i] + Moqj_1[i] + Moqj_2[i] + Moqj_3[i];

         AccountFirstBound(t0);
         Solve();

         PrintSolution(fout, t0);
         fout << endl;

         prev_solution3 = prev_solution2;
         prev_solution2 = prev_solution1;
         prev_solution1 = solution;
      }
   }

   // Базисные функции для расчета решения не в узле сетки
   double phi1(double t)
   {
      return 2 * (t - 0.5) * (t - 1);
   }

   double phi2(double t)
   {
      return -4 * t * (t - 1);
   }

   double phi3(double t)
   {
      return 2 * t * (t - 0.5);
   }

   // Вывод решения на временном слое t в поток fout 
   void PrintSolution(ofstream& fout, const double& t)
   {
      fout << "t = " << fixed << t << endl;
      fout << setw(14) << "x" << setw(14) << "y";
      fout << setw(14) << "prec" << setw(14) << "calc" << setw(14) << "diff" << setw(5) << "n" << " loc" << endl;

      double norm = 0, norm_u = 0;

      for(int y_i = 0; y_i < y_nodes_count; y_i++)
      {
         for(int x_i = 0; x_i < x_nodes_count; x_i++)
         {
            int i = x_i + y_i * x_nodes_count;

            double prec = true_solution[i];
            double calc = solution[i];

            //if(x_i % 32 == 0 && y_i % 32 == 0)
            {
               fout << scientific;
               fout << setw(14) << x_nodes[x_i];
               fout << setw(14) << y_nodes[y_i];
               fout << setw(14) << prec;
               fout << setw(14) << calc;
               fout << setw(14) << abs(true_solution[i] - solution[i]);
               fout << fixed << setw(5) << i;

               if(location[i] == 2)
                  fout << " outer";
               else if(location[i] == 1)
                  fout << " border";
               else
                  fout << " inner";

               fout << endl;

            }
            norm_u += prec * prec;
            norm += abs(prec - calc) * abs(prec - calc);
         }
      }

     // Блок вывода для вывода решения в произвольной точке расчетной области
     /* double x = 0.25;
      double y = 0.25;

      vector<double> s = solution;

      double calc = 0;
      calc += phi1(x) * phi1(y) * s[0];
      calc += phi2(x) * phi1(y) * s[1];
      calc += phi3(x) * phi1(y) * s[2];
      calc += phi1(x) * phi2(y) * s[3];
      calc += phi2(x) * phi2(y) * s[4];
      calc += phi3(x) * phi2(y) * s[5];
      calc += phi1(x) * phi3(y) * s[6];
      calc += phi2(x) * phi3(y) * s[7];
      calc += phi3(x) * phi3(y) * s[8];

      double prec = regions[0].test.u(2, 2, 3);

      fout << scientific;
      fout << setw(14) << 2.0;
      fout << setw(14) << 2.0;
      fout << setw(14) << prec;
      fout << setw(14) << calc;
      fout << setw(14) << abs(prec - calc);
      fout << fixed << setw(5) << "-";

      fout << " point" << endl;*/

      fout << "||u-u*||/||u*|| = " << scientific << sqrt(norm) / sqrt(norm_u) << endl;
      fout << "||u-u*||" << scientific << sqrt(norm) << endl;
   }
};