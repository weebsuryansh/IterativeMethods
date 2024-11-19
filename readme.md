# Comparison of Stationary and Non-Stationary Iterative Methods

This project implements and compares stationary and non-stationary iterative methods for solving large sparse linear systems. The implemented methods are:

- **Jacobi Method** (Stationary)
- **Gauss-Seidel Method** (Stationary)
- **Conjugate Gradient (CG) Method** (Non-Stationary)
- **Generalized Minimal Residual (GMRES) Method** (Non-Stationary)

The implementation is written in **C++**, utilizing the [Eigen](https://eigen.tuxfamily.org/) library for basic matrix operations and functionalities.

---

## Features

- Implementation of **stationary** and **non-stationary** iterative solvers.
- Comparative analysis of performance metrics:
    - Convergence speed.
    - Accuracy of solutions.
    - Computational efficiency.
- Designed to handle **large sparse matrices**.

---

## Prerequisites

To run this project, ensure the following are installed on your system:

1. **C++ Compiler** (e.g., g++, clang++)
2. **CMake** (For building the project)
3. [**Eigen**](https://eigen.tuxfamily.org/) (Matrix library)

---

## Building and Running

### Clone the Repository

```bash
git clone https://github.com/weebsuryansh/IterativeMethods
cd IterativeMethods
```

### Build the project
```bash
mkdir build
cd build
cmake ..
make
```
### Running the executable
```bash
.\IterativeMethods
```

---

## Input format
The code uses 3x3, 4x4, 5x5 and 6x6 matrices that can be edited by going into the `\data` directory.
The matrices used are in `.mtx` format.

---

## Result and Analysis
Result and analysis can be found in `Iterative_Methods.pdf`.

---

## Project Structure

```
project_root/ 
├── CMakeLists.txt          # Build configuration
├── Eigen/                  # Library used for basic matrix operations 
├── data/                   # Example input files (matrices and vectors) 
├── build/                  # Directory for compiled files (generated) 
├── Iterative_Methods.pdf   # The research paper entailing the result and analysis 
└── README.md               # Project documentation
```