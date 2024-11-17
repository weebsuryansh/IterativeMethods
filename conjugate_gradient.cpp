#include <utility>
#include "Eigen/Dense"
#include <cmath>
#include <iostream>


using namespace std;
using namespace Eigen;

class ConjugateGradient{
    MatrixXd _matrix_A;
    MatrixXd _vector_X_current;
    MatrixXd _vector_b;
    pair<int,int> _dimensions;
    double _epsilon;


public:
    ConjugateGradient(const MatrixXd &matrix, const MatrixXd &vector_b, const int dimension_X, const int dimension_Y, const double epsilon){
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
        MatrixXd _vector_r_current = _vector_b - (_matrix_A*_vector_X_current);
        MatrixXd _vector_P_current = _vector_r_current.eval();;
        int iteration=0;
        do{
            iteration++;
            const double r_current_mod = (_vector_r_current.transpose()*_vector_r_current).value();
            const double P_current_mod_A = (_vector_P_current.transpose()*_matrix_A*_vector_P_current).value();
            double alpha_current = r_current_mod / P_current_mod_A;

            MatrixXd _vector_X_next = _vector_X_current + alpha_current * _vector_P_current;

            MatrixXd _vector_r_next = _vector_r_current - alpha_current * _matrix_A * _vector_P_current;

            const double r_next_mod = (_vector_r_next.transpose()*_vector_r_next).value();

            double beta_current = r_next_mod / r_current_mod;

            MatrixXd _vector_P_next = _vector_r_next + beta_current * _vector_P_current;

            //updating current vectors
            _vector_r_current = _vector_r_next;
            _vector_P_current = _vector_P_next;
            _vector_X_current = _vector_X_next;
        }while (twoNorm(_vector_r_current)>_epsilon);

        cout <<"The solver took "<<iteration <<" iterations to complete."<< endl;
        return _vector_X_current;
    }
};

