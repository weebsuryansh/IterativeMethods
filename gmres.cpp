#include <utility>
#include "Eigen/Dense"
#include <cmath>
#include <iostream>


using namespace std;
using namespace Eigen;

class GMRES{
    MatrixXd _matrix_A;
    MatrixXd _matrix_H; // Hessenberg Matrix
    MatrixXd _matrix_V;
    VectorXd _vector_X;
    VectorXd _vector_b;
    pair<int,int> _dimensions;
    int _max_iter;
    double _epsilon;


public:
    GMRES(const MatrixXd &matrix, const VectorXd &vector_b, const int dimension_X, const int dimension_Y, const double epsilon, const int max_iter){
        _matrix_A = matrix;
        _vector_b = vector_b;
        _epsilon = epsilon;
        _dimensions = make_pair(dimension_X, dimension_Y);
        _max_iter = max_iter;
        _matrix_H = MatrixXd::Zero(_max_iter+1,_max_iter);
        _vector_X =VectorXd::Zero(_dimensions.second);
        _matrix_V = MatrixXd::Zero(_dimensions.second, _max_iter+1);
    }
    MatrixXd getMatrix() {
        return _matrix_A;
    }
    double twoNorm(VectorXd _vector) {
        double sum=0;
        for (int i = 0; i < _vector.rows(); i++) {
            sum+=_vector(i)*_vector(i);
        }
        return sqrt(sum);
    }

    pair<MatrixXd,MatrixXd> GivensRotationQR(MatrixXd A) {
        long const rows = A.rows();
        long const cols = A.cols();
        MatrixXd Q_transpose = MatrixXd::Identity(rows,rows);
        MatrixXd R = A.eval();


        for (long i = 0; i < cols; i++) {
            long const j = i + 1;
            if(R(j,i)==0) {
                continue;
            }
            const double hypot = sqrt((R(i,i)*R(i,i)) + (R(j,i)*R(j,i)));
            const double cos = R(i,i)/hypot;
            const double sin = R(j,i)/hypot;

            for(long k = 0; k < rows; k++) {
                double a_ik;
                double a_jk;
                if (k<cols) {
                    a_ik = R(i,k);
                    a_jk = R(j,k);
                    R(i,k ) = (cos*a_ik - sin*a_jk);
                    R(j,k) = (sin*a_ik + cos*a_jk);
                }
                a_ik = Q_transpose(i,k);
                a_jk = Q_transpose(j,k);
                Q_transpose(i,k ) = (cos*a_ik - sin*a_jk);
                Q_transpose(j,k) = (sin*a_ik + cos*a_jk);
            }
        }

        return make_pair(R,Q_transpose);
    }

    VectorXd solve() {
        VectorXd _vector_r = _vector_b - _matrix_A*_vector_X;
        double r_mod = _vector_r.norm();
        VectorXd _vector_V_current = _vector_r / r_mod;
        _matrix_V.col(0) = _vector_V_current;

        int iter = 0;

        do {
            iter++;
            _vector_V_current = _vector_r / r_mod;
            _matrix_V = MatrixXd::Zero(_dimensions.second, _max_iter+1);
            _matrix_V.col(0) = _vector_V_current;
            _matrix_H = MatrixXd::Zero(_max_iter+1,_max_iter);

            for (int j=0; j<_max_iter; j++) {
                VectorXd _vector_AV = _matrix_A*_matrix_V.col(j);
                VectorXd sum_vector=VectorXd::Zero(_dimensions.second);

                for (int i = 0; i < j+1; i++) {
                    _matrix_H(i,j) = (_vector_AV.transpose()*_matrix_V.col(i));
                    sum_vector+=_matrix_H(i,j)*_matrix_V.col(i);
                }

                VectorXd _vector_V_next = _vector_AV - sum_vector;
                _matrix_H(j+1,j) = twoNorm(_vector_V_next);
                if (_matrix_H(j+1,j) > _epsilon) {
                    _vector_V_next = _vector_V_next / _matrix_H(j+1,j);
                }
                _matrix_V.col(j+1)=_vector_V_next;
            }

            VectorXd _vector_e1 = VectorXd::Zero(_max_iter+1);
            _vector_e1(0)=1;

            pair<MatrixXd,MatrixXd> R_Qt = GivensRotationQR(_matrix_H);
            MatrixXd _matrix_R = R_Qt.first;
            MatrixXd _matrix_Qt = R_Qt.second;

            VectorXd _vector_G = r_mod*_matrix_Qt*_vector_e1;
            //back-substitution
            VectorXd _vector_y = _matrix_R.topRows(_max_iter).triangularView<Upper>().solve(_vector_G.head(_max_iter));

            //multiplying V_m*y_m
            _vector_X = _vector_X + _matrix_V.leftCols(_max_iter)*_vector_y;
            _vector_r = _vector_b - _matrix_A*_vector_X;
            r_mod = _vector_r.norm();
            if (_max_iter<_dimensions.second) {
                _max_iter++;
            }
            cout<<iter<<" "<<r_mod<<endl;
        }while (r_mod > _epsilon);

        cout <<"The solver took "<<iter <<" iterations to complete."<< endl;
        return _vector_X;
    }
};
