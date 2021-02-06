#include<iostream>
using namespace std;
const int N = 1e2;
double input[N];
int DEP,WID[N];
const double alpha = 0.2;
double W[N][N][N], b[N][N], X[N][N], Z[N][N];
double dC_dW[N][N][N], dC_db[N][N], dC_dX[N][N];

double rand_num(double l,double u){
    double x=rand(),y=rand();
    if(x<y)swap(x,y);
    return y/x * (u-l) + l;
}
double loss(){

}

double act(const double &x){return x<0?0:x;}
double act_(const double &x){return x<0?0:1;}
void fp_cal_layer(const int &l){
    for(int i=0;i<WID[l];i++){
        double res = 0;
        for(int j=0;j<WID[l-1];j++){
            res += W[l-1][i][j] * X[l-1][j] + b[l-1][j];
        }
        Z[l][i] = res;
    }
}
void fp_act_layer(const int &l){
    for(int i=0;i<WID[l];i++)
        X[l][i] = act(Z[l][i]);
}
void front_propagation(){
    for(int i=1;i<DEP;i++){
        // Z_L = W_L-1 * X_L-1 + b_L-1
        fp_cal_layer(i);
        // X_L = act(Z_L)
        fp_act_layer(i);
    }
}
void bp_cal_dC_dW(const int &l){
    for(int i=0;i<WID[l+1];i++)
        for(int j=0;j<WID[l];j++)
            dC_dW[l][i][j] = dC_dX[l+1][i]*act_(Z[l+1][i])*X[l][j];
}
void bp_cal_dC_db(const int &l){
    for(int j=0;j<WID[l];j++){
        dC_db[l][j] = 0;
        for(int i=0;i<WID[l+1];i++)
            dC_db[l][j] += dC_dX[l+1][i]*act_(Z[l+1][i]);
    }
}
void bp_cal_dC_dX(const int &l){
    for(int j=0;j<WID[l];j++){
        dC_dX[l][j] = 0;
        for(int i=0;i<WID[l+1];i++)
            dC_dX[l][j] += dC_dX[l+1][i]*act_(Z[l+1][i])*act_(Z[l][j]);
    }
}
void update(const int &l){
    for(int i=0;i<WID[l+1];i++)
        for(int j=0;j<WID[l];j++)
            W[l][i][j] -= alpha * dC_dW[l][i][j];
    for(int i=0;i<WID[l];i++)
        b[l][i] -= alpha * dC_db[l][i];
}
void back_propagation(){
    for(int i=DEP-1;i>0;i--){
        bp_cal_dC_dW(i);
        bp_cal_dC_db(i);
        bp_cal_dC_dX(i);
        update(i);
    }
}
void run_net(int times){
    for(int i=1;i<=times;i++){
        front_propagation();
        cout<<"generation "<<i<<":"<<endl;
        cout<<"the loss = "<<loss()<<endl;
        back_propagation();
    }
}
void init_rand(){
    for(int l=0;l<=DEP;l++)
        for(int i=0;i<WID[l+1];i++)
            for(int j=0;j<WID[l];j++)
                W[l][i][j] = rand_num(-100,100);
    for(int l=0;l<DEP;l++)
        for(int i=0;i<WID[l];i++)
            b[l][i] = rand_num(-100,100);
}
