#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;
const int N = 1e5 + 10;
const double C = 0.1;
int n, m;
vector<vector<double> > x;
vector<int> y;
vector<double> alp;
double b;

inline double dot(vector<double> &a, vector<double> &b) {
    static double res = 0;
    for (int i = 0; i < m; i++)
        res += a[i] * b[i];
    return res;
}

inline double Kernel(vector<double> &x1, vector<double> &x2) {
    double tmp = 0;
    const double sigma = 2;
    for (int i = 0; i < m; i++)
        tmp += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    return exp((tmp * tmp) / 2 * sigma * sigma);
}

inline double f(vector<double> &_x) {
    double res = 0;
    static vector<double> w;
    w.reserve(m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            w[j] += alp[i] * y[i] * x[i][j];
    return dot(w, _x) - b;
}

inline double eta(int i, int j) {
    return Kernel(x[i], x[i]) + Kernel(x[j], x[j]) - 2 * Kernel(x[i], x[j]);
}

inline double E(int i) {
    return f(x[i]) - y[i];
}

inline double alp2new(int i, int j) {
    static double res = alp[j] + y[j] * (E(i) - E(j)) / eta(i, j);
    static double L, H, alp3_n = 0;

    for (int k = 0; k < n; k++)
        if (k != i && k != j)
            alp3_n += alp[k] * y[k];

    if (y[i] == y[j])
        L = max(0.0,alp3_n-C), H = min(C,alp3_n);
    else
        L = max(0.0,-alp3_n), H = min(C,C-alp3_n);

    if(res>H)
        return H;
    if(res<L)
        return L;
    return res;
}

inline double alp1new(int i,int j,double a2old,double a2new) {
    return alp[i] + y[i]*y[j]*(a2old-a2new);
}

inline double bnew(int i,int j,double a1old,double a1new,double a2old,double a2new) {
    static double b1,b2;
    b1 = b - E(i) - y[i]*(a1new-a1old)*Kernel(x[i],x[i]) - y[j]*(a2new-a2old)*Kernel(x[i],x[j]);
    b2 = b - E(j) - y[i]*(a1new-a1old)*Kernel(x[i],x[j]) - y[j]*(a2new-a2old)*Kernel(x[j],x[j]);
    if(0<a1new && a1new<C && 0<a2new && a2new<C)
        return (b1+b2)/2;
    if(0<a1new && a1new<C)
        return b1;
    if(0<a2new && a2new<C)
        return b2;
}

void init() {
    ifstream ifs("in.txt");
    ifs >> n >> m;
    y.reserve(n);
    x.reserve(n);
    alp.reserve(n);
    for (int i = 0; i < n; i++) {
        x[i].reserve(m);
        for (int j = 0; j < m; j++)
            ifs >> x[i][j];
    }
    for (int i = 0; i < n; i++)
        ifs >> y[i];
    ifs.close();
}


int main() {
    return 0;
}
