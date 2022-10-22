#include <cstdio>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include  <vector>
#include  <fstream>
#include <set>
#include  <array>
#include <tuple>
using namespace std;
vector<vector<double> > Matrix;
double Z;
set<int> P;
size_t cn, bn;
bool verbose = false;
vector<vector<double>> hcat(vector<vector<double>> m, vector<double> n) {
    for (int i = 0; i < m.size(); i++) {
        
            m[i].push_back(n[i]);
   
    }

    return m;

}
vector<double>vec_mul(vector<double> a, vector<vector<double>> b) {

    // 创建一个一维数组
    vector<double> mat;
    // A矩阵的行数，A矩阵的列数，B矩阵的行数，B矩阵的列数
    int col1 = a.size(), row2 = b.size(), col2 = b[0].size();
    // A矩阵的列数和B矩阵的行数不相等时，矩阵无法相乘
    if (col1 != row2)
    {
        cout << "A矩阵的列数和B矩阵的行数不相等！" << endl;
        return mat;
    }

    // 构建出 A矩阵的行数 * B矩阵的列数 的新矩阵
    // 矩阵A的行数
    for (int i = 0; i < col2; i++)
    {
        double ans = 0;
       
        for (int j = 0; j < row2; j++)
        {
            ans = ans + a[j] * b[j][i];
            // 把值存到暂存数组
         
        }
        // 把暂存数组放到矩阵C
        mat.push_back(ans);
    }

    return mat;

}
vector<double>vec_mul2(vector<vector<double>> a, vector<double> b) {

    // 创建一个一维数组
    vector<double> mat;
    // A矩阵的行数，A矩阵的列数，B矩阵的行数，B矩阵的列数
    int row = a.size(), col = a[0].size(), row2 = b.size();
    // A矩阵的列数和B矩阵的行数不相等时，矩阵无法相乘
    if (col != row2)
    {
        cout << "A矩阵的列数和B矩阵的行数不相等！" << endl;
        return mat;
    }

    // 构建出 A矩阵的行数 * B矩阵的列数 的新矩阵
    // 矩阵A的行数
    for (int i = 0; i < row; i++)
    {
        double ans = 0;
        // 创建一个暂存数组
        //vector<double> v;
        // 矩阵B的列数
        for (int j = 0; j < col; j++)
        {
            // 暂存新矩阵 C[i][j] 的值
            
            // 矩阵A的列数

            ans = ans + a[i][j] * b[j];

            // 把值存到暂存数组
           
        }
        // 把暂存数组放到矩阵C
        //mat.push_back(v);
        mat.push_back(ans);
    }

    return mat;

}
vector<vector<double> > PivotalTransform(vector<vector<double> > a, vector<int> bases, int p, int q)//枢轴变换
{

    vector<vector<double> > b(a);
    // printV(b);
    // cout << p << q<<endl;
    for (int i = 0; i < a.size(); i++) {
        if (i != p)
            b[i][q] = a[i][q] / a[p][q];
    }
    for (int m = 0; m < a.size(); m++) {
        if (m == p) {
            for (int j = 0; j < a[0].size(); j++) {

                b[m][j] = a[m][j] / a[p][q];

            }
        }
        if (m != p) {
            for (int j = 0; j < a[0].size(); j++) {
                if (j != q) {
                    b[m][j] = a[m][j] - b[m][q] * a[p][j];
                }


            }
        }



    }
    for (int t1 = 0; t1 < b.size(); t1++) {
        if (t1 != p) {
            for (int j1 = 0; j1 < b[0].size(); j1++) {
                b[t1][q] = 0;



            }
        }

    }
    return b;
}
int SelectP(vector<double> xo, vector<double> xq)//出基
{
    int m = xo.size();
    if (m!=xq.size()) {
        cout << "vector x_o and vector x_q does not has zhe same length!";
    
    }
    int p = -1;
    double pv = INT_MAX;
    for (int i = 0; i < m; i++) {
        if (xq[i]>0) {
            double temp = xo[i] / xq[i];
            if (temp<pv) {
                p = i;
                pv = temp;

            }
        
        }
    
    }
    return p;

}

int SelectQ(vector<double> r)
{
    int n = r.size();
    int q = -1;
    int qv = 0;
    for (int i = 0; i < n;i++) {
        if (r[i]>qv) {
            q = i;
            qv = r[i];
        }
    
    }
    return q;

}
void printV(vector<vector<double >> a ) {
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < a[i].size(); j++) {
            cout << a[i][j] << "  ";

        }
        cout << endl;

    }



}
tuple<vector<vector<double>>, vector<double>, vector<int>, vector<int>>Simplex(
    vector<vector<double>> A_a,
    vector<vector<double>> B_b,
    vector<double> X_x,
    vector<double> C,
    vector<int> basesB,
    vector<int> basesD)
{
    vector<double> X(X_x);
    vector<vector<double> > A(A_a);
    int n = A.size();
    int m = A.size();
    vector<vector<double> > B(B_b);

    for (int v = 0; v < 128; v++) {
       
        vector<double> cb;
        cout << "基为";
        for (int i = 0; i < basesB.size(); i++) {
            cb.push_back(C[basesB[i]]);
           cout << basesB[i] + 1 << "   ";
        }
        cout << endl;
        vector<double> cd;
        for (int i = 0; i < basesD.size(); i++) {
            cd.push_back(C[basesD[i]]);

        }
        vector<vector<double>> D;
        for (int i = 0; i < A.size(); i++) {
            vector<double> d_row;
            for (int j = 0; j < A[i].size(); j++) {
                for (int m = 0; m < basesD.size(); m++) {
                    if (j == basesD[m]) {
                        d_row.push_back(A[i][j]);

                    }
                }

            }
            D.push_back(d_row);
        }
        vector<double> r;
        vector<double> res;
        res = vec_mul(vec_mul(cb, B), A);
        for (int i = 0; i < C.size(); i++) {
            r.push_back(C[i] - res[i]);


        }
        int q = SelectQ(r);
        if (q == -1) {
            double ans = 0;
            for (int k = 0; k < cb.size(); k++)
            {
                // 矩阵A x轴的数 乘 矩阵B y轴的数，累加得到，新矩阵 C[i][j] 的值
                ans = ans + cb[k] * X[k];
            }
            cout << "minimumcost is" << ans;
            tuple<vector<vector<double>>, vector<double>, vector<int>, vector<int>> result =
                make_tuple(B, X, basesB, basesD);
            return result;
        }

        vector<double> d_l;
        for (int l = 0; l < A.size(); l++) {
            d_l.push_back(A[l][q]);
        }
        //vector<double> xq = vec_mul2(B, d_l);
        int p = SelectP(X, d_l);
        if (p == -1) {
            cout << "not bound";

            tuple<vector<vector<double>>, vector<double>, vector<int>, vector<int>> result =
                make_tuple(B, X, basesB, basesD);
            return result;
        }
        vector<vector<double>> A_pivot;
        A_pivot = PivotalTransform(hcat(A, X), basesB,p, q);
        printV(A_pivot);
        for (int i = 0; i < A_pivot.size(); i++) {
           X[i] = A_pivot[i][A_pivot[0].size() - 1];

        }
        
        for (int i = 0; i < A_pivot.size(); i++) {
            //vector<double> b_l1;
            for (int j = 0; j < m ; j++) {
                A[i][j] = A_pivot[i][j];

            }
            //b_l.push_back(b_l1);
        }
        
        
        int temp = basesB[p];
        basesB[p] = q;
       // basesD[q] = temp;
        //vector<vector<double>> b_l;
        for (int i = 0; i < A_pivot.size(); i++) {
            //vector<double> b_l1;
            for (int j = 0; j < basesB.size(); j++) {
                B[i][j] = A_pivot[i][basesB[j]];

            }
            //b_l.push_back(b_l1);
        }
    
        cout << "第" << v + 1 << "次迭代" << endl;
        printV(A);
        // cout << "BJuzhen" << endl;
        printV(B);
        for (int i = 0; i < X.size();i++) {
            cout << X[i] << ",";
        }
        cout << endl;
    }
    cout << "warning!";
    tuple<vector<vector<double>>, vector<double>, vector<int>, vector<int>> result =
        make_tuple(B, X, basesB, basesD);
    return result;

}

int main(int argc, char* argv[])
{
    ifstream fin;
    vector<int> basesB;  
    vector<int> basesD;
    vector<double> x;           //b vector in Ax=b
    
    vector<vector<double> > A;
    vector<vector<double> > b;  //基矩阵
    fin.open("test.txt");       //输入矩阵
    fin >> cn >> bn;            //矩阵大小
    for (int g = 0; g < cn;g++) {//初始基可行解
        int base = 0;
        fin >> base;
        if (base != 0) { basesB.push_back(base-1); }
        
    }
    for (int s = 0; s < cn - 1; s++) {
        bool flag = true;
        for (int w = 0; w < basesB.size(); w++) {
            if (s==basesB[w]) {
                flag = false;
            }
        
        }
        if (flag) {
            basesD.push_back(s);
        }
    }
    for (size_t i = 0; i < bn-2; i++)
    {
        vector<double> vectmp;
        vector<double> b_tmp;
        double b_x = 0;
        for (size_t j = 0; j < cn-1; j++)
        {
            double tmp = 0;
            double b_i = 0;
            fin >> tmp;
            for (int k = 0; k < basesB.size();k++) {
                if (basesB[k] == j){ 
                    b_tmp.push_back(tmp);
                    
                }
                   
            }
            vectmp.push_back(tmp);
           
        }
        b.push_back(b_tmp);
        fin >> b_x;
        x.push_back(b_x);
        A.push_back(vectmp);
    }
    vector<double> c_vector;    //c vector in object function,c'*x
    for (size_t i = 0; i < cn-1; i++) {
    
        double tmp = 0;
        fin >> tmp;
        c_vector.push_back(tmp);
    
    }
    if (verbose) {
        cout << "Original Basics is" ;
        for (set<int>::iterator it = P.begin(); it != P.end(); it++) {
            cout << *it << " ";
        }
    
    }
    Simplex(A,b,x, c_vector,basesB,basesD);
}

/////////////////////////////////////
//glpk input:
///* Variables */
//var x1 >= 0;
//var x2 >= 0;
//var x3 >= 0;
///* Object function */
//maximize z: x1 + 14*x2 + 6*x3;
///* Constrains */
//s.t. con1: x1 + x2 + x3 <= 4;
//s.t. con2: x1  <= 2;
//s.t. con3: x3  <= 3;
//s.t. con4: 3*x2 + x3  <= 6;
//end;
/////////////////////////////////////
//myinput:
//
//7 4 下列矩阵大小
//5 6 0 0 0 0 0 初始基向量列索引
//1 - 2 3 - 1 1 0 15
//2 1 - 1 2 0 1 10//////Ax=b
//- 3 - 2 - 1 1 0 0 0  目标函数系数
/////////////////////////////////////