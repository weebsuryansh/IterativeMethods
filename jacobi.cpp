#include <utility>
#include "Eigen/Dense"
#include <cmath>
#include <iostream>


using namespace std;
using namespace Eigen;

class Jacobi {
    MatrixXd _matrix_A;
    VectorXd _vector_X;
    VectorXd _vector_b;
    pair<int,int> _dimensions;
    double _epsilon;


public:
    Jacobi(const MatrixXd &matrix, const VectorXd &vector_b, const int dimension_X, const int dimension_Y, const double epsilon){
        _matrix_A = matrix;
        _vector_b = vector_b;
        _epsilon = epsilon;
        _dimensions = make_pair(dimension_X, dimension_Y);
        _vector_X =VectorXd::Zero(_dimensions.second);
    }
    MatrixXd getMatrix() {
        return _matrix_A;
    }
    double twoNorm(VectorXd _vector) {
        double sum=0;
        for (int i = 0; i <_vector.rows(); i++) {
            sum+=_vector(i)*_vector(i);
        }
        return sqrt(sum);
    }

    VectorXd solve() {
        if (_matrix_A.diagonal().array().abs().minCoeff() < std::numeric_limits<double>::epsilon()) {
            throw std::runtime_error("Matrix A has zero diagonal elements; Jacobi may fail.");
        }
        int iteration=0;
        VectorXd _vector_X_next=VectorXd::Zero(_dimensions.second);
        VectorXd _vector_r;
        do {
            iteration++;
            for (int i = 0; i < _dimensions.second; i++) {
                double sum=0;
                for (int j = 0; j < _dimensions.second; j++) {
                    if (j != i) {
                        sum+=(_matrix_A(i,j)*_vector_X(j));
                    }
                }
                _vector_X_next(i) = (_vector_b(i)-sum)/_matrix_A(i,i);
            }
            _vector_r = _vector_b - (_matrix_A*_vector_X_next);
            _vector_X = _vector_X_next;
            cout<<iteration<<" "<<_vector_r.norm()<<endl;
        } while(twoNorm(_vector_r)>_epsilon);

        cout <<"The solver took "<<iteration <<" iterations to complete."<< endl;
        return _vector_X;
    }

};