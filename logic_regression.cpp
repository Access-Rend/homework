#include <bits/stdc++.h>
using namespace std;
//
//
//
inline double sigmoid(const double &x){
    return 1/(1+exp(-x));
}
const int N = 1e5+10;
double x[N],y[N],w[N],n;
double y_[N];
double J(){
    double res = 0;
    for(int i=0;i<n;i++)
        res -= (y[i]*log(y_[i]) + (1-y[i])*log(1-y_[i]) );
    return res / n;
}

inline double J_a(const double &x,const double &c){
    return - (c/x) + ((1-c)/(1-x));
}
inline double a_y(const double &x){
    return sigmoid(x)*(1-sigmoid(x));
}
inline double y_w(const double &x){
    return x;
}

void GD(int times){
    double alpha = 1;
    while(times--){
        for(int i=0;i<n;i++)
            y_[i] = sigmoid(x[i]*w[i]);

        for(int i=0;i<n;i++){
            w[i] -= J_a(y_[i],y[i]) * a_y(x[i]*w[i]) * y_w(x[i]) * alpha;
        }
        cout<<J()<<endl;
    }
}
void init_w(double d){
    for(int i=0;i<n;i++)
        w[i] = d;
}
void rand_data(){
    ofstream ofs("data");
    ofs<<"1000\n";
    for(int i=0;i<1000;i++)
        ofs<<rand() / double(RAND_MAX)<<' ';

    ofs<<endl;

    for(int i=0;i<1000;i++)
        ofs<<rand()%2<<' ';

    ofs.close();
}

int main() {
    srand(time(0));

    ifstream ifs("data");
    ifs>>n;
    for(int i=0;i<n;i++)ifs>>x[i];
    for(int i=0;i<n;i++)ifs>>y[i];
    init_w(0.5);

    GD(100000);

    ofstream ofs("result");
    ofs<<n<<endl;
    for(int i=0;i<n;i++)
        ofs<<w[i]<<' ';

    return 0;
}
