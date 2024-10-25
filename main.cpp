//main program to navigate here
#include "linear_methods.h"
#include "non_linear_methods.h"
#include"inverseMatrix.h"
#include <bits/stdc++.h>
using namespace std;

int main() {
    int choice;

    while (true) {
        cout << "\nChoose an option:\n";
        cout << "1. Jacobi Iteration           2. Gauss Seidel Iteration       3. Gauss Elimination\n";
        cout << "4. Gauss-Jordan Elimination   5. LU Factorization              6. Bisection Method\n";
        cout << "7. False Position Method       8. Newton's Method              9. Secant Method\n";
        cout << "10. Other Method 1            11. Other Method 2              12. Exit\n";
        cout << "Enter your choice: ";

        cin >> choice;

        if (cin.fail()) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "Invalid input! Please enter a number between 1 and 12.\n";
            continue;
        }

        switch (choice) {
            case 1:
                jacobi();
                break;
            case 2:
                gauss_seidel();
                break;
            case 3:
                gaussElimination();
                break;
            case 4:
                gaussJordan();
                break;
            case 5:
                lu_factor();
                break;
            case 6:
                biSection();
                break;
            case 7:
                falsePosition();
                break;
            case 8:
                newton();
                break;
            case 9:
                secant();
                break;
            case 10:
                // Call runge kutta method
                break;
            case 11:
                // Call matrix inversion method
                break;
            case 12:
                cout << "Exiting the program.\n";
                return 0;
            default:
                cout << "Invalid choice! Please try again.\n";
                break;
        }
    }

    return 0;
}
