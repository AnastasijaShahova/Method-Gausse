#include <iostream>
#include<vector>
#include<math.h>
#include<omp.h>
#include<chrono>
 
using namespace std;
const double eps = 0.000001;
 
void getmatr(vector<vector<double>>& M, int n) {
 
    M.resize(n);
    for (int i = 0; i < n; ++i) {
        M[i].resize(n);
        for (int j = 0; j < n; ++j)
            M[i][j] = rand() % 100;
    }
 
}
void getb(vector<double>& b, int n) {
    b.resize(n);
    for (int i = 0; i < n; ++i)
 
        b[i] = rand() % 100;
}
int find_max(vector<vector<double>>& M, int t, int n) {
    double max = fabs(M[t][t]);
    int index = t;
#pragma omp parallel for  
    for (int i = t; i < n; ++i) {
        if (fabs(M[i][t]) > max) {
            max = fabs(M[i][t]);
            index = i;
        }
    }
    if (max < eps)
        return -1;
    else  return index;
}
void show_matr(vector<vector<double>>& M, int n) {
    cout << "Значения матрицы:" << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            cout << M[i][j] << " ";
        cout << endl;
    }
}
void showx(vector<double>x, int n)
{
    for (int i = 0; i < n; ++i)
        cout << "X[" << i << "]=" << x[i] << endl;
}
 
vector<double> gauss(vector<vector<double>>& M, vector<double>& b, int n) {
    vector<double> x(n);
    int k;
    for (int i = 0; i < n; ++i) {
       int  m = find_max(M, i, n);
        if (m == -1) {
            cout << "Reshenie polychit neBozmozno";
            exit(1);
        }
        if (m != i) {
            swap(b[m], b[i]);
            swap(M[m], M[i]);
        }
        for (int j = i; j < n; ++j) {
            double temp = M[j][i];
            if (fabs(temp) < eps) continue;
#pragma omp parallel for private(k)
            for ( k = 0; k < n; ++k)
                M[j][k] = M[j][k] / temp;
            b[j] = b[j] / temp;
            if (j == i) continue;
            
#pragma omp parallel for private(k)
            for (k = 0; k < n; ++k)
                M[j][k] = M[j][k] - M[i][k];
            b[j] = b[j] - b[i];
 
        }
    }
   
    for (int k = n - 1; k >= 0; --k) {
        x[k] = b[k];
        int i;
#pragma omp parallel for  private(i)
        for ( i = 0; i < k; ++i)
           b[i] = b[i] - M[i][k] * x[k];
    }
    return x;
 
}
 
int main()
{
       setlocale(LC_ALL, "Russian");
for (int thread = 1; thread <= 16; thread * 2) {
        cout << thread << endl;
        omp_set_num_threads(thread);
        for (int n = 1000; n <= 10000; n += 1000) {
            omp_set_num_threads(threads); 
            vector<vector<double>>M(n);
            for (int i = 0; i < n; ++i) {
                M[i].resize(n);
                for (int j = 0; j < n; ++j)
                    M[i][j] = rand()%100;
            }
            vector<double>b(n);
            for (int i = 0; i < n; ++i)
                b[i] = rand()%100;
 
            auto begin = chrono::steady_clock::now();
            vector<double> x = gauss(M, b, n);
            auto end = chrono::steady_clock::now();
            auto time = chrono::duration_cast<chrono::seconds>(end - begin);
            cout << " Размер= " << n << "  Время:" << time.count() << endl;
        }      
    }
   return 0;
    
}

