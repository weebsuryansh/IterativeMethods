#include <iostream>
#include "Eigen/Dense"
#include "jacobi.cpp"
#include "gauss_seidel.cpp"

using namespace std;
using namespace Eigen;

int main() {
    string input;
    Matrix4d matrix_A;
    matrix_A << 10,2,1,2,
                2,10,2,1,
                1,2,10,2,
                2,1,2,10;
    Matrix<double,4,1> vector_b;
    MatrixXd vector_x;
    vector_b << 15000,15000,15000,15000;

    // 10^-8
    constexpr double epsilon = 0.00000001;
    int choice;

    cout<<"Choose methos:\n1. Jacobi\n2. Gauss-Seidel\nEnter choice:"<<endl;

    cin>>choice;
    getline(cin,input);

    if (choice == 1) {
        Jacobi jacobi_solver (matrix_A, vector_b, 4, 4, epsilon);
        vector_x = jacobi_solver.solve();
    }
    else if (choice == 2) {
        GaussSeidel gauss_seidel_solver (matrix_A, vector_b, 4, 4, epsilon);
        vector_x = gauss_seidel_solver.solve();
    }
    else {
        cout<<"Wrong choice!"<<endl;
    }

    cout<<"Print the result?"<<endl;
    getline(cin, input);
    cout << vector_x << endl;

    return 0;
}
