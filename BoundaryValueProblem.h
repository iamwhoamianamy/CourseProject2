#pragma once
#include "Matrix.h"
#include "SLAE.h"
#include "Test.h"
using namespace std;

class BoundaryValueProblem
{
public:
   Matrix global;               // ���������� �������
   Matrix fac_global;           // �������� ������������ ���������� �������
   Matrix stiff_mat;            // ������� ���������
   Matrix weight_mat;           // ������� ����
   Matrix G1, G2, M;            // ��������������� ������� ��� ���������� ���������
                                // ��������� ������ � �������� ������ ������

   SLAE slae;                   // �������� ������� ��� ������������������
   SLAE fac_slae;               // �������� ������� c �������������������

   vector<real> b;              // ���������� ������ ������ �����
   vector<real> loc_b;          // ��������� ������ ������ �����
   vector<real> solution;       // �������

   Test test;                   // ���������� � ��������� �������,
                                // ��������� ������ � �����

   int elems_count;              // ���������� �������� ���������
   int nodes_count;              // ���������� �����

   vector<double> x_nodes;      // ���������� ����� �� x
   vector<double> y_nodes;      // ���������� ����� �� y

   int x_nodes_count;           // ���������� ����� �� x
   int y_nodes_count;           // ���������� ����� �� y

   // ������� ���������� �������� �� ����� file_name
   // � ������������ �����
   void ReadFormGrid(const string& file_name)
   {
      ifstream fin(file_name);

      // ���������� ������� �����
      double left, right, bot, top;

      // ���������� ������� �������
      fin >> left;
      fin >> right;
      fin >> bot;
      fin >> top;

      // ��������� ��������� ����� �� x
      int n;
      double h, q;

      fin >> q >> n;

      x_nodes_count = n + 1;
      x_nodes.resize(x_nodes_count);

      h = right - left;

      if(q != 1)
         h *= (1 - q) / (1 - pow(q, n));
      else
         h /= n;

      x_nodes[0] = left;

      for(int i = 0; i < n; i++)
         x_nodes[i + 1] = x_nodes[i] + h * pow(q, i);

      // ��������� ��������� ����� �� y
      fin >> q >> n;

      y_nodes_count = n + 1;
      y_nodes.resize(y_nodes_count);

      h = top - bot;

      if(q != 1)
         h *= (1 - q) / (1 - pow(q, n));
      else
         h /= n;

      y_nodes[0] = bot;

      for(int i = 0; i < n; i++)
         y_nodes[i + 1] = y_nodes[i] + h * pow(q, i);

      fin.close();

      nodes_count = x_nodes_count * y_nodes_count;
      elems_count = ceil(nodes_count / 9.0);
   }

   // ��������� ������ ��� �������
   void InitializeMemory()
   {
      slae = SLAE(nodes_count, 10000, 1e-20);
      fac_slae = SLAE(nodes_count, 10000, 1e-20);

      fac_global = Matrix(nodes_count, 0);

      global.ind.resize(nodes_count + 1);
      b.resize(nodes_count);
      solution.resize(nodes_count);
   }

   // ���������� ������� global_indices ���������, ���������������� ���������� ���������
   // ����� ��������� �������� � ������� elem_index(���������� � ����)
   void CalcGlobalIndices(int elem_index, vector<int>& global_indices)
   {
      elem_index++;
      int n_coords = x_nodes_count / 2 + 1;
      int k = 2 * floor((elem_index - 1) / (n_coords - 1)) * (2 * n_coords - 1) + 2 * ((elem_index - 1) % (n_coords - 1)) + 1;
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

   // ��������������� ������� ��� ������������ ��������
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

   // ������������ �������� ��������� �������
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


         for(int i = 1; i < 9; i++)
            for(int j = 0; j < i; j++)
               IncertToRow(help[i][j][0], help[i][j][1]);
      }

      global.size = global.ind.size() - 1;
      global.diag.resize(global.size);
      
      global.tr_size = global.columns_ind.size();
      global.bot_tr.resize(global.tr_size);
      global.top_tr.resize(global.tr_size);
   }
};