#include<vector>
#include<string>
class Mat_r{
public:
	std::vector<std::vector<double> > element;
	unsigned int col,row;
	Mat_r () = default;
    Mat_r (unsigned int c,unsigned int r):col(c),row(r){
        element.reserve(c);
        for(auto &x:element)
            x.reserve(r);
    }
	auto operator [] (unsigned int x){
	    return &this->element[x];
	}
	Mat_r operator + (const Mat_r &m){
		if(this->col != m.col || this->row != m.row)
			throw std::string("矩阵形状不一样加你妈呢.");
		Mat_r res(this->row,this->col);
		for(int i=0;i<this->row;i++)
			for(int j=0;j<this->col;j++)
				res.element[i][j] = this->element[i][j] + m.element[i][j];
		return res;
	}
	Mat_r operator * (const Mat_r &m){
		if(this->col != m.row)
			throw std::string("矩阵形状都不对乘你妈呢.");
		Mat_r res(this->row,m.col);
		for(int i=0;i<this->row;i++)
			for(int j=0;j<m.col;j++)
			    for(int k=0;k<this->col;k++)
				    res.element[i][j] += this->element[i][k] * m.element[k][j];
		return res;
	}
	Mat_r T(){
        Mat_r res(this->col,this->row);
        for(int i=0;i<this->col;i++)
            for(int j=0;j<this->row;j++)
                res.element[i][j] = this->element[j][j];
        return res;
    }
    ~Mat_r() = default;
};