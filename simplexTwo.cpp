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
vector<vector<double>> hcatI(vector<vector<double>> m, int n) {
    for (int i = 0; i < m.size(); i++) {

        for (int j = 0; j < n;j++) {
            if (i == j) { 
                m[i].push_back(1); 
            }else
            m[i].push_back(0);
        }

    }

    return m;

}
vector<vector<double>> buildI(int n) {

    vector<vector<double>> i_b;
    for (int i = 0; i < n; i++) {
        vector<double> i_row;
        for (int j = 0; j < n; j++) {
            if (i == j) { i_row.push_back(1); }
            else
                i_row.push_back(0);
        }
        i_b.push_back(i_row);
    }

    return i_b;

}
vector<double>vec_mul(vector<double> a, vector<vector<double>> b) {

    // ����һ��һά����
    vector<double> mat;
    // A�����������A�����������B�����������B���������
    int col1 = a.size(), row2 = b.size(), col2 = b[0].size();
    // A�����������B��������������ʱ�������޷����
    if (col1 != row2)
    {
        cout << "A�����������B�������������ȣ�" << endl;
        return mat;
    }

    // ������ A��������� * B��������� ���¾���
    // ����A������
    for (int i = 0; i < col2; i++)
    {
        double ans = 0;

        for (int j = 0; j < row2; j++)
        {
            ans = ans + a[j] * b[j][i];
            // ��ֵ�浽�ݴ�����

        }
        // ���ݴ�����ŵ�����C
        mat.push_back(ans);
    }

    return mat;

}
vector<double>vec_mul2(vector<vector<double>> a, vector<double> b) {

    // ����һ��һά����
    vector<double> mat;
    // A�����������A�����������B�����������B���������
    int row = a.size(), col = a[0].size(), row2 = b.size();
    // A�����������B��������������ʱ�������޷����
    if (col != row2)
    {
        cout << "A�����������B�������������ȣ�" << endl;
        return mat;
    }

    // ������ A��������� * B��������� ���¾���
    // ����A������
    for (int i = 0; i < row; i++)
    {
        int ans = 0;
        // ����һ���ݴ�����
        //vector<double> v;
        // ����B������
        for (int j = 0; j < col; j++)
        {
            // �ݴ��¾��� C[i][j] ��ֵ

            // ����A������

            ans = ans + a[i][j] * b[j];

            // ��ֵ�浽�ݴ�����

        }
        // ���ݴ�����ŵ�����C
        //mat.push_back(v);
        mat.push_back(ans);
    }

    return mat;

}
vector<vector<double> > PivotalTransform(vector<vector<double> > a, vector<int> bases, int p, int q)//����任
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
                if(j!=q){
                    b[m][j] = a[m][j] - b[m][q] * a[p][j];
                }
              

            }
        }

       

    }
    for (int t1 = 0; t1 < b.size();t1++) {
        if (t1!=p) {
            for (int j1 = 0; j1 < b[0].size(); j1++) {
                b[t1][q] = 0;



            }
        }
        
    }
    return b;
}
int SelectP(vector<double> xo, vector<double> xq)//����
{
    int m = xo.size();
    if (m != xq.size()) {
        cout << "vector x_o and vector x_q does not has zhe same length!";

    }
    int p = -1;
    double pv = 10000;
    for (int i = 0; i < m; i++) {
        if (xq[i] > 0) {
            double temp = xo[i] / xq[i];
            if (temp < pv) {
                p = i;
                pv = temp;

            }

        }

    }
    return p;

}

int SelectQ(vector<double> r, vector<vector<double>> a,int minmax)
{
    int n = r.size();
    int q = -1;
    if (minmax ==1) {
        
        int qv = 0;
        int count = 0;
        bool flag = true;
        for (int i = 0; i < n; i++) {
            if (r[i] < qv) {
                q = i;
                count++;
                qv = r[i];
            }

        }
        if (count == 1) {
            for (int k = 0; k < a.size(); k++) {
                if (a[k][q] > 0) {
                    flag = false;
                }

            }
            if (flag) { return -2; }

        }
    }
    else {
       
        int qv = 0;
        int count = 0;
        bool flag = true;
        for (int i = 0; i < n; i++) {
            if (r[i] > qv) {
                q = i;
                count++;
                qv = r[i];
            }

        }
        if (count == 1) {
            for (int k = 0; k < a.size(); k++) {
                if (a[k][q] > 0) {
                    flag = false;
                }

            }
            if (flag) { return -2; }

        }
    }
    
    
    return q;

}
void printV(vector<vector<double >> a) {
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < a[i].size(); j++) {
            cout << a[i][j] << "  ";

        }
        cout << endl;

    }



}
tuple<vector<vector<double>>, vector<vector<double>>, vector<double>, vector<int>, vector<int>>Simplex(
    vector<vector<double>> A_a,
    vector<vector<double>> B_b,
    vector<double> X_x,
    vector<double> C,
    vector<int> basesB,
    vector<int> basesD,
    int minmax)
{
    vector<double> X(X_x);
    vector<vector<double> > A(A_a);
    int n = A.size();
    int m = A.size();
    vector<vector<double> > B(B_b);

    for (int v = 0; v < 128; v++) {
        cout << "��" << v + 1 << "�ε���" << endl;
        printV(A);
        // cout << "BJuzhen" << endl;
        printV(B);
        vector<double> cb;
        cout << "��Ϊ";
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
        int q = SelectQ(r,A, minmax);
        if (q == -2) {
            cout << "not bound";

            tuple<vector<vector<double>>, vector<vector<double>>, vector<double>, vector<int>, vector<int>> result =
                make_tuple(A, B, X, basesB, basesD);
            return result;
        }
        if (q == -1) {
            double ans = 0;
            for (int k = 0; k < cb.size(); k++)
            {
                // ����A x����� �� ����B y��������ۼӵõ����¾��� C[i][j] ��ֵ
                ans = ans + cb[k] * X[k];
            }
            cout << "minimumcost is" << ans<<endl;
            tuple<vector<vector<double>>,vector<vector<double>>, vector<double>, vector<int>, vector<int>> result =
                make_tuple(A,B, X, basesB, basesD);
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

            tuple<vector<vector<double>>, vector<vector<double>>, vector<double>, vector<int>, vector<int>> result =
                make_tuple(A, B, X, basesB, basesD);
            return result;
        }
        vector<vector<double>> A_pivot;
        A_pivot = PivotalTransform(hcat(A, X), basesB, p, q);
        printV(A_pivot);
        for (int i = 0; i < A_pivot.size(); i++) {
            X[i] = A_pivot[i][A_pivot[0].size() - 1];

        }

        for (int i = 0; i < A_pivot.size(); i++) {
            //vector<double> b_l1;
            for (int j = 0; j < A_pivot[0].size() - 1; j++) {
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

       
        for (int i = 0; i < X.size(); i++) {
            cout << X[i] << ",";
        }
        cout << endl;
    }
    cout << "warning!";
    tuple<vector<vector<double>>, vector<vector<double>>, vector<double>, vector<int>, vector<int>> result =
        make_tuple(A, B, X, basesB, basesD);
    return result;

}

int main(int argc, char* argv[])
{
    ifstream fin;
    vector<int> basesB;
    vector<int> basesD;
    vector<double> x;           //b vector in Ax=b

    vector<vector<double> > A;
    vector<vector<double> > b;  //������
    fin.open("test.txt");       //�������
    fin >> cn >> bn;            //�����С

    for (size_t i = 0; i < bn - 1; i++)
    {
        vector<double> vectmp;
        double b_x = 0;
        for (size_t j = 0; j < cn - 1; j++)
        {
            double tmp = 0;
           // double b_i = 0;
            fin >> tmp;
            
            vectmp.push_back(tmp);

        }
        //b.push_back(b_tmp);
        fin >> b_x;
        x.push_back(b_x);
        A.push_back(vectmp);
    }
    vector<double> c_vector;    //c vector in object function,c'*x
    for (size_t i = 0; i < cn - 1; i++) {

        double tmp = 0;
        fin >> tmp;
        c_vector.push_back(tmp);

    }
    if (verbose) {
        cout << "Original Basics is";
        for (set<int>::iterator it = P.begin(); it != P.end(); it++) {
            cout << *it << " ";
        }

    }
    int a_s = A.size();
    vector<double> c1_vector;
    for (int m = 0; m < c_vector.size();m++) {
        c1_vector.push_back(0);
    }
    for (int n = 0; n < a_s;n++) {
    
        c1_vector.push_back(1);
    }
    vector<int> baseb1;
    vector<int> based1;
    for (int n = 0; n < a_s; n++) {

        baseb1.push_back(A[0].size()+n);
    }
    for (int a1 = 0; a1 < A[0].size(); a1++) {

        based1.push_back(a1);
    }
    //paseI
    cout << "paseI"<<endl;
    tuple<vector<vector<double>>, vector<vector<double>>, vector<double>, vector<int>, vector<int>> par = Simplex(hcatI(A,a_s), buildI(a_s), x, c1_vector, baseb1, based1,1);
    //paseII
    cout << "paseII" << endl;
    vector<vector<double>>A_II ;
    vector<vector<double>>B_II;
    vector<double> x_ii;
    vector<int> b_ii;
    vector<int> d_ii;
    for (int a_ii = 0; a_ii < get<0>(par).size(); a_ii++) {
        vector<double> a_temp;
        for (int a_ij = 0; a_ij < get<0>(par)[0].size()-a_s; a_ij++) {
            a_temp.push_back(get<0>(par)[a_ii][a_ij]);
        
        }
        A_II.push_back(a_temp);
    }
    B_II = get<1>(par);
    x_ii = get<2>(par);
    for (int i = 0; i < get<4>(par).size();i++) {
        if (get<4>(par)[i]< get<0>(par)[0].size() - a_s) {
            d_ii.push_back(get<4>(par)[i]);
        }
    }
    b_ii = get<3>(par);

    Simplex(A_II, B_II, x_ii, c_vector, b_ii, d_ii,0);
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
//8 5
//7 4 ���о����С
//5 6 0 0 0 0 0 ��ʼ������������
//1 - 2 3 - 1 1 0 15
//2 1 - 1 2 0 1 10//////Ax=b
//- 3 - 2 - 1 1 0 0 0  Ŀ�꺯��ϵ��
/////////////////////////////////////