//non_linear all methods here
#include <bits/stdc++.h>
#include "non_linear_methods.h"
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
