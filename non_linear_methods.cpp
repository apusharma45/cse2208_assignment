//non_linear all methods here
#include "non_linear_methods.h"
#include <bits/stdc++.h>
using namespace std;

void find_ab(double A, double B, double C, double &a, double &b) {
    b = sqrt(pow(B/A, 2) - 2*C/A);
    a = -b;
}

double f(double x, double A, double B, double C) {
    return A * pow(x, 2) + B * x + C;
}

void biSection() {
    double A, B, C;
    cout << "Enter coefficients A, B, C for the equation Ax^2 + Bx + C = 0: ";
    cin >> A >> B >> C;
    
    double dif, x = 0, pre_x, a = 0, b = 0, res;
    find_ab(A, B, C, a, b);
    
    int count = 0;
    bool resultFound = false;

    while (count < 500) {
        ++count;
        pre_x = x;
        double f_a = f(a, A, B, C);
        double f_b = f(b, A, B, C);
        x = (a + b) / 2;
        double f_x = f(x, A, B, C);
        
        if (abs(f_x) < 0.0001) {
            res = x;
            resultFound = true;
            break;
        }
        
        if (f_x * f_a > 0) {
            a = x;
        } else {
            b = x;
        }
    }

    if (resultFound) {
        cout << "Result: " << res << endl;
        cout << "Iterations: " << count << endl;
    } else {
        cout << "Result not found in 500 iterations." << endl;
    }
}

void falsePosition() {
    double A, B, C;
    cout << "Enter coefficients A, B, C for the equation Ax^2 + Bx + C = 0: ";
    cin >> A >> B >> C;

    double dif, x = 0, pre_x, a = 0, b = 0, res;
    find_ab(A, B, C, a, b);

    int count = 0;
    bool resultFound = false;

    while (count < 500) {
        ++count;
        pre_x = x;
        double f_a = f(a, A, B, C);
        double f_b = f(b, A, B, C);
        x = (a * f_b - b * f_a) / (f_b - f_a);
        double f_x = f(x, A, B, C);

        if (abs(f_x) < 0.0001) {
            res = x;
            resultFound = true;
            break;
        }

        if (f_x * f_a > 0) {
            a = x;
        } else {
            b = x;
        }
    }

    if (resultFound) {
        cout << "Result: " << res << endl;
        cout << "Iterations: " << count << endl;
    } else {
        cout << "Result not found in 500 iterations." << endl;
    }
}

double df(double x, double A, double B, double C) {
    return 2 * A * x + B;
}

void newton() {
    double A, B, C;
    cout << "Enter coefficients A, B, C for the equation Ax^2 + Bx + C = 0: ";
    cin >> A >> B >> C;

    double x = 0, pre_x = 4, f_pre_x, df_pre_x;
    int count = 0;
    double res;
    int iteration;
    bool gotResult = false;

    while (count < 500) {
        ++count;
        f_pre_x = f(pre_x, A, B, C);
        if (abs(f_pre_x) < 0.0001 && !gotResult) {
            res = pre_x;
            iteration = count;
            gotResult = true;
            break;
        }

        df_pre_x = df(pre_x, A, B, C);
        if (df_pre_x == 0) {
            cout << "Can't be calculated. Derivative is zero!" << endl;
            return;
        }

        x = pre_x - f_pre_x / df_pre_x;
        pre_x = x;
    }

    if (gotResult) {
        cout << "Result: " << res << "  Iterations: " << iteration << endl;
    } else {
        cout << "Result not found in 500 iterations." << endl;
    }
}


void secant() {
    double A, B, C;
    cout << "Enter coefficients A, B, C for the equation Ax^2 + Bx + C = 0: ";
    cin >> A >> B >> C;

    double x = 0, a = 5, b = 10, f_x, f_a, f_b;
    int count = 0;
    double res;
    int iteration;
    bool gotResult = false;

    while (count < 500) {
        ++count; 

        f_a = f(a, A, B, C);
        f_b = f(b, A, B, C);
        if (abs(f_a - f_b) < 1e-10) {
            cout << "Function values are too close. No solution found." << endl;
            return;
        }

        x = b - (f_b * (b - a)) / (f_b - f_a);
        f_x = f(x, A, B, C);

        if (abs(f_x) < 0.0001 && !gotResult) {
            res = x;
            iteration = count;
            gotResult = true;
            break;
        }
        
        a = b;
        b = x;
    }

    if (gotResult) {
        cout << "Result: " << res << "  Iterations: " << iteration << endl;
    } else {
        cout << "Result not found in 500 iterations." << endl;
    }
}
