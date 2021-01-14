#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <ctime>
/*
 * 0<alpha_i<C
 */
using namespace std;
const int N = 1e5 + 10;
const double C = 100, eps = 1e-3, tolerance = 1e-3;
int n, m;
double x[1000][1000];
int y[1000];
double alp[1000],err[1000];
//vector<vector<double> > x;
//vector<int> y;
//vector<double> alp, err;
double b;

double dot(vector<double> &a, vector<double> &b) {
    static double res = 0;
    for (int i = 0; i < m; i++)
        res += a[i] * b[i];
    return res;
}
double dot(double *a, double *b) {
    double res = 0;
    for (int i = 0; i < m; i++)
        res += a[i] * b[i];
    return res;
}
double Kernel(double*x1, double*x2) {
    return dot(x1,x2);
}
double Kernel(vector<double> &x1, vector<double> &x2) {
//    double tmp = 0;
//    const double sigma = 50;
//    for (int i = 0; i < m; i++)
//        tmp += (x1[i] - x2[i]) * (x1[i] - x2[i]);
//    return exp((tmp * tmp) / 2 * sigma * sigma);
    return dot(x1,x2);
}
double f(double*_x) {//correct
    double res = 0;
    for (int i = 0; i < n; i++)
        res += alp[i] * y[i] * Kernel(x[i], _x);
    res -= b;
    return res;
}
//double f(vector<double> &_x) {//correct
//    static double res = 0;
//    for (int i = 0; i < n; i++)
//        if(abs(alp[i])<eps)
//            continue;
//        else
//            res += alp[i] * y[i] * Kernel(x[i], _x);
//    res -= b;
//    return res;
//}

inline bool update_alp1alp2(int &i, int &j) {
    if (i == j)
        return false;
    int &y1 = y[i], &y2 = y[j], s = y[i] * y[j];
    double &alp1 = alp[i], &alp2 = alp[j], a1new, a2new, bnew;
    double E1, E2, L, H, L_obj, H_obj;
    double K11 = Kernel(x[i], x[i]),
            K22 = Kernel(x[j], x[j]),
            K12 = Kernel(x[i], x[j]);
    double eta = K11 + K22 - 2 * K12;

    E1 = (0 < alp1 && alp1 < C) ? err[i] : f(x[i]) - y1;
    E2 = (0 < alp2 && alp2 < C) ? err[j] : f(x[j]) - y2;
    /*
     * updating for alpha2_new
     */
    if (y1 == y2) {
        if(alp1+alp2>C)
            L = alp1 + alp2 - C,
            H = C;
        else
            L = 0,
            H = alp1 + alp2;
    } else {
        if(alp1-alp2>0)
            L = 0,
            H = C - (alp1-alp2);
        else
            L = -(alp1-alp2),
            H = C;
    }

    if (abs(L - H) < eps)
        return false;

    if (eta > 0) {
        a2new = alp2 + y2 * (E1 - E2) / eta;
        if (a2new < L)
            a2new = L;
        else if (a2new > H)
            a2new = H;
    } else {
        double t1 = -eta / 2, t2 = y2 * (E1 - E2) + eta * alp2;
        L_obj = t1 * L * L + t2 * L;
        H_obj = t1 * H * H + t2 * H;
        if (L_obj - H_obj > eps)
            a2new = L;
        else if (L_obj - H_obj < -eps)
            a2new = H;
        else
            a2new = alp2;
    }

    if (abs(a2new - alp2) < eps)
        return false;
    /*
     * updating for alpha1_new
     */
    a1new = alp1 - s * (a2new - alp2);
    if (a1new < 0)
        a2new += s * a1new, a1new = 0;
    else if (a1new > C)
        a2new += s * (a1new - C), a1new = C;
    /*
     * updating for constant b
     */
    if (0 < a1new && a1new < C)
        bnew = b + E1 + y1 * (a1new - alp1) * K11 + y2 * (a2new - alp2) * K12;
    else {
        double b1, b2;
        b1 = b + E1 + y1 * (a1new - alp1) * K11 + y2 * (a2new - alp2) * K12;
        b2 = b + E2 + y1 * (a1new - alp1) * K12 + y2 * (a2new - alp2) * K22;
        bnew = (b1 + b2) / 2;
    }
    /*
     * updating for err
     */
    double t1 = y1 * (a1new - alp1), t2 = y2 * (a2new - alp2);
    for (int k = 0; k < n; k++)
        if (0 < alp[k] && alp[k] < C)
            err[k] += t1 * Kernel(x[k], x[i]) + t2 * Kernel(x[k], x[j]) - (bnew - b);

    err[i] = err[j] = 0;
    alp1 = a1new, alp2 = a2new, b = bnew;
    return true;
}

bool examine_first(int &i, double &E1) {
    static int j = -1;
    static double tmax = 0, E2;
    for (int k = 0; k < n; k++) {
        if (0 < alp[k] && alp[k] < C) { // 在间隔边界即支持向量上寻找满足条件的第二个点
            E2 = err[k];
            if (abs(E1 - E2) > tmax)
                tmax = abs(E1 - E2), j = k;
        }
    }
    if (j != -1)
        return update_alp1alp2(i, j);
    return false;
}

bool examine_all(int &i) {
    /*
     * 可以增添随机遍历
     * 先遍历支持向量，再遍历其他
     */
    for (int j = 0; j < n; j++) {
        if (0 < alp[j] && alp[j] < C)
            if (update_alp1alp2(i, j))
                return true;
    }
    for (int j = 0; j < n; j++) {
//        if (alp[j] <=  0 || alp[j] >= C) 可能有问题
            if (update_alp1alp2(i, j))
                return true;
    }
    return false;
}

bool examine_example(int &i) {
    int &y1 = y[i];
    double &alp1 = alp[i], E1, r1;
    if (0 > alp1 && alp1 < C)
        E1 = err[i];
    else
        E1 = f(x[i]) - y1;
    r1 = y1 * E1;

    if ((r1 > tolerance && alp1 > 0) || (r1 < -tolerance && alp1 < C)) { // 不满足kkt
        if (examine_first(i, E1))
            return true;
        if (examine_all(i)) //遍历
            return true;
    }
    return false;
}

void work() {
    int num_changed = 0;
    bool examine_all = true;
    while (num_changed > 0 || examine_all) {
        num_changed = 0;
        if (examine_all)
            for (int i = 0; i < n; i++)
                if(examine_example(i))
                    ++num_changed;
        else
            for (int i = 0; i < n; i++)
                if (alp[i] != 0 && alp[i] != C)
                    if(examine_example(i))
                        ++num_changed;

        if (examine_all)
            examine_all = false;
        else if (num_changed == 0)
            examine_all = true;
    }
}

void load_and_init() {
    ifstream ifs("data.txt");
    ifs >> n >> m;
    b = 0;
    for(int i=0;i<n;i++)
        alp[i] = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            ifs >> x[i][j];
    }
    for (int i = 0; i < n; i++)
        ifs >> y[i];
    ifs.close();
    return;
}

void save() {
    ofstream ofs("model.txt");
    ofs << n << " " << m << " " << tolerance << endl;
    for (int i = 0; i < n; i++)
        ofs << alp[i] << " ";
    ofs<<endl;
    double w[5]={0,0,0,0,0};
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            w[j]+=y[i]*alp[i]*x[i][j];
    ofs<<"w:"<<endl;
    for(int i=0;i<m;i++)
        ofs<<w[i]<<" ";
    ofs<<endl<<"b:"<<b<<endl;
    ofs.close();
    return;
}

void rand_data() {
    srand(time(0));
    int _n = 100, _m = m = 5;
    double _b = double(rand()%1000)/(rand()%1000+1);
    vector<double> _x(_m),_w(_m);
    vector<int> _y(_n);

    for(int i=0;i<_m;i++)
        _w[i] = double(rand()%100)/(rand()%100+1);


    ofstream ofs("data.txt");

    ofs<<_n<<" "<<_m<<endl;
    for(int i=0;i<_n;i++){
        for(int j=0;j<_m;j++){
            _x[j] = rand()/double(rand()+1);
            if(rand()&1)
                _x[j] = -_x[j];
            ofs<<_x[j]<<" ";
        }
        ofs<<endl;
        double tmp = dot(_w,_x)-_b;
        if(tmp==0){
            i--;
            continue;
        }
        if(tmp>0)_y[i]=1;
        else _y[i]=-1;
    }
    for(int i=0;i<_n;i++)
        ofs<<_y[i]<<endl;


    ofs<<"w:";
    for(int i=0;i<_m;i++)
        ofs<<_w[i]<<" ";
    ofs<<endl<<"b:"<<_b<<endl;

    ofs.close();
    return;
}
void test(){
    int tot=0;
    double w[1000];
    memset(w,0,sizeof(w));
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++)
            w[j] += alp[i]*y[i]*x[i][j];
    }
    for(int i=0;i<n;i++){
        cout<<i<<":wx-b="<<(dot(w,x[i])-b)<<", y="<<y[i];
        if((dot(w,x[i])-b)*y[i]>0){
            cout<<"yes";tot++;
        }
        else cout<<"no";
        cout<<endl;
    }
    cout<<"rate:"<<double(tot)/n<<endl;
}
int main() {
    rand_data();
double startime = clock();
    load_and_init();
//    for(int k=0;k<n;k++){
//        for(int i=0;i<m;i++)
//            cout<<x[k][i]<<" ";
//        cout<<endl;
//    }
    work();
    save();
    test();
    cout<<"cost:"<<(clock()-startime)/CLOCKS_PER_SEC<<"s"<<endl;
    return 0;
}
