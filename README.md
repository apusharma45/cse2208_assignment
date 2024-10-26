# cse2208_assignment
Numerical assignment of group 26 (2107005, 2107030, 2107099)


## Run Command

g++ main.cpp linear_methods.cpp non_linear_methods.cpp inverseMatrix.cpp runge_kutta.cpp -o main && ./main


## Table of Contents

- [Jacobi Iterative Method](#jacobi-iterative-method)
- [Gauss-Seidel Method](#gauss-seidel-method)
- [Gauss Elimination](#gauss-elimination)
- [Gauss-Jordan Elimination](#gauss-jordan-elimination)
- [LU Factorization](#lu-factorization)
- [Bisection Method](#bisection-method)
- [False Position Method](#false-position-method)
- [Secant Method](#secant-method)
- [Newton-Raphson Method](#newton-raphson-method)
- [Runge-Kutta Method](#runge-kutta-method)
- [Matrix Inversion](#matrix-inversion)

---
### Algorithms for Linear Equations

### Jacobi Iterative Method

The **Jacobi Iterative Method** is an algorithm for solving a system of linear equations using an initial guess and iterative updates. It converges if the matrix is diagonally dominant. Each variable is solved in terms of other variables using the previous iteration's values.

### Gauss-Seidel Method

The **Gauss-Seidel Method** is similar to the Jacobi method but improves convergence by using the updated values as soon as they are computed. This method often converges faster than Jacobi, especially for strictly diagonally dominant matrices.

### Gauss Elimination

**Gauss Elimination** is a method to solve linear systems by reducing the matrix to an upper triangular form. Once the matrix is in this form, back substitution is used to solve for the unknowns.

### Gauss-Jordan Elimination

The **Gauss-Jordan Elimination** method extends Gauss Elimination by reducing the matrix to a reduced row echelon form (diagonalized), allowing for a direct solution without back substitution.

### LU Factorization

**LU Factorization** decomposes a matrix into a lower triangular matrix \(L\) and an upper triangular matrix \(U\) such that \(A = LU\). This decomposition simplifies the solution of linear equations and is particularly useful for solving multiple systems with the same coefficient matrix.

### Algorithms for Non-Linear Equations


### Bisection Method

The **Bisection Method** is a root-finding method that iteratively divides an interval in half and selects the subinterval in which the function changes sign. This method guarantees convergence to a root if the function is continuous and changes signs over the interval.

### False Position Method

The **False Position Method** (or Regula Falsi) is a root-finding algorithm that uses linear interpolation between endpoints to find roots. This method is similar to bisection but often converges more quickly for functions that are not strictly linear.

### Secant Method

The **Secant Method** is an iterative root-finding algorithm that uses a sequence of secants to approximate the root. It is faster than the bisection and false position methods but does not guarantee convergence unless the initial guesses are close to the root.

### Newton-Raphson Method

The **Newton-Raphson Method** is a fast and widely used root-finding algorithm that requires the derivative of the function. Starting from an initial guess, it uses tangent lines to converge to a root. This method is effective for well-behaved functions but may fail if the derivative is zero or the function is not smooth.


### Algorithm for Differential Equations 

### Runge-Kutta Method

The **Runge-Kutta Method** is a numerical integration method for solving ordinary differential equations (ODEs). This repository implements the popular fourth-order Runge-Kutta method, which is widely used for its accuracy and stability.


### Matrix Inversion

**Matrix Inversion** provides the inverse of a square matrix, if it exists. Matrix inversion is a foundational operation in linear algebra and has applications in various fields, including solving systems of linear equations.

---




