#include <utility>
#include "Eigen/Dense"
#include <cmath>
#include <iostream>


using namespace std;
using namespace Eigen;

class GMRES{
    MatrixXd _matrix_A;
    MatrixXd _vector_X_current;
    MatrixXd _vector_b;
    pair<int,int> _dimensions;
    double _epsilon;


public:
    GMRES(const MatrixXd &matrix, const MatrixXd &vector_b, const int dimension_X, const int dimension_Y, const double epsilon){
        _matrix_A = matrix;
        _vector_b = vector_b;
        _epsilon = epsilon;
        _dimensions = make_pair(dimension_X, dimension_Y);
        _vector_X_current = MatrixXd::Zero(_dimensions.second, 1);
    }
    MatrixXd getMatrix() {
        return _matrix_A;
    }
    double twoNorm(MatrixXd _vector) {
        double sum=0;
        for (int i = 0; i < _vector.rows(); i++) {
            sum+=_vector(i,0)*_vector(i,0);
        }
        return sqrt(sum);
    }

    MatrixXd solve() {

    }
};
