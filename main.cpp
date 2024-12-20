#include <iostream>
#include <fstream>
#include <chrono>
#include "Eigen/Dense"
#include "jacobi.cpp"
#include "gauss_seidel.cpp"
#include "conjugate_gradient.cpp"
#include "gmres.cpp"

using namespace std;
using namespace Eigen;

int main() {
    string input;
    MatrixXd matrix_A;
    VectorXd vector_b;
    VectorXd vector_x_true;
    VectorXd vector_x;

    int choice;
    int num_rows;
    int num_cols;
    int num_lines;
    string filename;
    cout<<"Enter the matrix you want to run on: \n1. 3x3\n2. 4x4\n3. 5x5\n4. 6x6"<<endl;
    cin >> choice;
    if (choice==1) {
        filename = "../data/matrix3.mtx";
    }
    else if (choice==2) {
        filename = "../data/matrix4.mtx";
    }
    else if (choice==3) {
        filename = "../data/matrix5.mtx";
    }
    else if (choice==4) {
        filename = "../data/matrix6.mtx";
    }
    else {
        cout<<"Wrong choice!"<<endl;
        exit(0);
    }

    ifstream input_file(filename);
    if (!input_file) {
        cerr << "Error: Could not open the file!" << endl;
        return 1; // Exit with error code
    }

    input_file >> num_rows >> num_cols >> num_lines;

    matrix_A = MatrixXd::Zero(num_rows, num_cols);

    for (int l = 0; l < num_lines; l++)
    {
        getline(input_file, input);
        double data;
        int row, col;
        input_file >> row >> col >> data;
        matrix_A(--row,--col) = data;
        cout<<row<<" "<<col<<" "<<matrix_A(row,col)<<endl;
    }

    vector_x_true = VectorXd::Zero(num_cols);
    vector_x_true.setRandom();
    vector_x_true *= 10;
    vector_x_true+=50*VectorXd::Ones(num_cols);
    vector_b = matrix_A * vector_x_true;

    constexpr double epsilon = 1e-8;
    if ((matrix_A-matrix_A.transpose()).norm()>epsilon) {
        cout<<"This matrix ain't symmetric"<<endl;
    }

    cout<<"True x:"<<endl;
    cout<<vector_x_true<<endl;

    cout<<"Choose method:\n1. Jacobi\n2. Gauss-Seidel\n3. Conjugate Gradient\n4. GMRES\nEnter choice:"<<endl;

    cin>>choice;
    getline(cin,input);

    if (choice == 1) {
        Jacobi jacobi_solver (matrix_A, vector_b, num_rows, num_cols, epsilon);
        auto start = chrono::high_resolution_clock::now();
        vector_x = jacobi_solver.solve();
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end-start);
        cout<<"The iterator took "<<duration.count()<<" microseconds to execute."<<endl;
    }
    else if (choice == 2) {
        GaussSeidel gauss_seidel_solver (matrix_A, vector_b, num_rows, num_cols, epsilon);
        auto start = chrono::high_resolution_clock::now();
        vector_x = gauss_seidel_solver.solve();
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end-start);
        cout<<"The iterator took "<<duration.count()<<" microseconds to execute."<<endl;
    }
    else if (choice == 3) {
        ConjugateGradient conjugate_gradient_solver (matrix_A, vector_b, num_rows, num_cols, epsilon);
        auto start = chrono::high_resolution_clock::now();
        vector_x = conjugate_gradient_solver.solve();
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end-start);
        cout<<"The iterator took "<<duration.count()<<" microseconds to execute."<<endl;
    }
    else if (choice == 4) {
        GMRES gmres_solver(matrix_A, vector_b, num_rows, num_cols, epsilon,num_cols);
        auto start = chrono::high_resolution_clock::now();
        vector_x =gmres_solver.solve();
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end-start);
        cout<<"The iterator took "<<duration.count()<<" microseconds to execute."<<endl;
    }
    else {
        cout<<"Wrong choice!"<<endl;
        exit(0);
    }


    cout<<"Print the result?"<<endl;
    getline(cin, input);
    cout << vector_x<< endl;

    return 0;
}
