#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
struct mat{
	int row,col;
	vector<vector<double> > element;
	mat() = default;
    mat(int r,int c,double d = 0):row(r),col(c){
        element.reserve(r);
        for(int i=0;i<r;i++)
            element[i].reserve(c),
            element[i].insert(element[i].begin(),c,d);
    }
	mat(const vector<vector<double> > &e) : element(e){row = e.size(),col = e[0].size();}
	vector<double>& operator [](const unsigned int x){return element[x];}
	mat operator +(mat &m){
	    if(this->row!=m.row || this->col!=m.col)
	        throw std::string("型不一样加你妈的矩阵呢");
	    mat res(this->row,this->col,0);
	    for(int i=0;i<this->row;i++)
	        for(int j=0;j<this->col;j++)
	            res[i][j] = this->element[i][j] + m[i][j];
        return res;
	}
	mat operator *(mat &m){
        if(this->col!=m.row)
            throw std::string("型对不上乘你妈的矩阵呢");
        mat res(this->row,m.col,0);
        for(int i=0;i<this->row;i++)
            for(int j=0;j<m.col;j++)
                for(int k=0;k<m.row;k++)
                    res[i][j] += this->element[i][k]*m[k][j];
        return res;
	}
	mat operator = (const mat m){
        this->row = m.row;
        this->col = m.col;
        this->element = m.element;
    }

//	ostream&operator <<(ostream& out,mat &m){
//        out<<'[';
//        for(auto x:m.element){
//            out<<'[';
//            for(auto y:x)
//                out<<y<<' ';
//            out<<']\n';
//        }
//        out<<']';
//    }
};
void print(mat &m){
    for(int i=0;i<m.row;i++){
        for(int j=0;j<m.col;j++)
            cerr<<m.element[i][j]<<' ';
        cerr<<endl;
    }
}

mat readmat(ifstream &ifs,int n,int m){
    mat res(n,m);
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            ifs>>res[i][j];
    return res;
}
int main(){
    ifstream ifs("data");

    mat x;

    mat X[20];
    for(int i=0;i<20;i++)
        X[i] = readmat(ifs,10,20);


	return 0;
}