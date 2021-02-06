#include<iostream>
using namespace std;
const int N = 1e2;
int DEP,WID[N];
const double step = 0.2;
double W[N][N][N], b[N][N], X[N][N], Z[N][N];
double dC_dW[N][N][N], dC_db[N][N], dC_dX[N][N];

double act(const double &x){return x<0?0:x;}
double act_(const double &x){return x<0?0:1;}
void fp_cal_layer(const int &l){
    for(int i=0;i<WID[l+1];i++){
        double res = 0;
        for(int j=0;j<WID[l];j++){
            res += W[l][i][j] * X[l][j];
        }
        Z[l][i] = res;
    }
}
void fp_act(const int &l){
    for(int i=0;i<WID[l];i++)
        X[l][i] = act(Z[l][i]);
}
void front_propagation(){
    for(int i=1;i<DEP;i++){
        // Y_L+1 = W_L * X_L + b
        fp_cal_layer(i);
        // X_L = act(Y_L)
        fp_act(i);
    }
}
void bp_cal_dC_dW(const int &l){
    for(int i=0;i<WID[l];i++)
        for(int j=0;j<WID[l-1];j++)
            dC_dW[l-1][i][j] = dC_dX[l][i]*act_(Z[l][i])*X[l-1][j];
}
void bp_cal_dC_db(const int &l){
    for(int j=0;j<WID[l-1];j++){
        dC_db[l-1][j] = 0;
        for(int i=0;i<WID[l];i++)
            dC_db[l-1][j] += dC_dX[l][i]*act_(Z[l][i]);
    }
}
void bp_cal_dC_dX(const int &l){
    for(int j=0;j<WID[l-1];j++){
        dC_dX[l-1][j] = 0;
        for(int i=0;i<WID[l];i++)
            dC_dX[l-1][j] += dC_dX[l][i]*act_(Z[l][i])*act_(Z[l-1][j]);
    }
}
void back_propagation(){
    for(int i=DEP;i>0;i--){
        bp_cal_dC_dW(i);
        bp_cal_dC_db(i);
        bp_cal_dC_dX(i);
    }
}
