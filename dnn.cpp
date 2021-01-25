#include <bits/stdc++.h>

using namespace std;
const int N = 1e2;
const double step = 0.2;
double Y[N];
double W[N][N][N], b[N][N], X[N][N][2];
/*
 * x[][][0]未经过激活函数的值(公式中x)，x[][][1]经过激活函数的值(公式中a)
 * W[i][][]的规格为(width[i+1],width[i])
 * 每一层公式: wx+b=z ,g(z) = a
 * 大写为矩阵、向量，小写标量
 */

double W_[N][N],b_[N],X_[N];
/*
 * W_,b_,X_为反向传播时用到的导数变量
 */
int depth, width[N]; // 0为输入层，1~depth为隐藏层，depth+1为输出层

double act(const double &x){return x<0?0:x;}
double act_(const double &x){return x<0?0:1;}
double cost(){
    double res = 0;
    for(int i=0;i<width[depth+1];i++)
        res += (X[depth+1][i][1] - Y[i])*(X[depth+1][i][1] - Y[i]);
    return res/width[depth+1];
}
double dcost_dx1(){
    for(int i=0;i<width[depth+1];i++)
        X_[i] = (X[depth+1][i][1]-Y[i])/width[depth+1];
}
void update_W(int dep){
    for(int i=0;i<width[dep+1];i++)
        for(int j=0;j<width[dep];j++)
            W[dep][i][j] -= W_[i][j]*step;
}
void update_b(int dep){
    for(int i=0;i<width[dep];i++)
        b[dep][i] -= b_[i]*step;
}
void back_propagation_layer(int dep){
    // X_ 存储了dJ/da(当前层a)的值
    if(dep>depth)
        throw std::invalid_argument("depth must lower than the last layer.");
    for(int i=0;i<width[dep];i++)
        X_[i]*=act_(X[dep][i][1]); // dJ/da * da/dZ = dJ/dZ


    for(int i=0;i<width[dep+1];i++)
        for(int j=0;j<width[dep];j++){
            W_[i][j] = X[dep][j][0]; // dZ/dW
            W_[i][j] *= X_[j];  // dJ/dZ * dZ/dW = dJ/dW
        }
    update_W(dep);

    for(int i=0;i<width[dep];i++)
        b_[i] = X_[i]; // dJ/dZ * dZ/db = dJ/db
    update_b(dep);

    for(int j=0;j<width[dep];j++){
        double dz_dx = 0;
        for(int i=0;i<width[dep+1];i++) // zi
            dz_dx += W_[i][j];
        X_[j] *= dz_dx;
    }
}
void cal_layer(int dep) {
    if (dep < 1)
        throw std::invalid_argument("depth must greater than 0.");
    for (int i = 0; i < width[dep]; i++) {
        for (int j = 0; j < width[dep - 1]; j++)
            X[dep][i][0] += W[dep - 1][i][j] * X[dep - 1][j][1] + b[dep - 1][i];
        X[dep][i][1] = act(X[dep][i][0]);
    }
}

void input(double data[], int size) {
    for (int i = 0; i < size; i++)
        X[0][i] = data[i];
}

double *output(int dep) {
    return X[dep];
}

int main() {

    return 0;
}