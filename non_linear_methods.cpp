//non_linear all methods here
#include "non_linear_methods.h"
#include <bits/stdc++.h>
using namespace std;
double maxAmong(const vector<double>& a) {
    double max = -DBL_MAX;
    for (int i = 1; i < a.size(); ++i) {
        if (abs(a[i]) > max) {
            max = abs(a[i]);
        }
    }
    return max;
}

void findAB(const vector<double>& a, double& A, double& B) {
    A = 1 + (1 / abs(a[0])) * maxAmong(a);
    B = -A;
}

double f(const vector<double>& a, double x) {
    double value = 0;
    int n = a.size();
    for (int i = 0; i < n; ++i) {
        value += a[i] * pow(x, n - i - 1);
    }
    return value;
}

double df(const vector<double>& a, double x) {
    double value = 0;
    int n = a.size();
    for (int i = 0; i < n - 1; ++i) { 
        value += a[i] * (n - i - 1) * pow(x, n - i - 2);
    }
    return value;
}

void biSection() {
    int maxPow;
    cout << "Enter maximum power of the polynomial: ";
    cin >> maxPow;

    cout << "Enter " << maxPow + 1 << " coefficients: ";
    vector<double> a(maxPow + 1);
    for (int i = 0; i <= maxPow; ++i) {
        cin >> a[i];
    }

    double A, B;
    double res;
    findAB(a, A, B);

    double f_a = f(a, A);
    double f_b = f(a, B);

    while (f_a * f_b > 0) {
        B += 1.0;
        f_b = f(a, B);
        if (static_cast<int>(A) == static_cast<int>(B)) {
            cout << "The solution may not be found!" << endl;
            return;
        }
    }

    int count = 0;
    bool resultFound = false;
    double x = 0, f_x;

    while (count < 500) {
        ++count;
        x = (A + B) / 2;
        f_x = f(a, x);

        if (abs(f_x) < 0.0001) {
            res = x;
            resultFound = true;
            break;
        }

        if (f_x * f_a > 0) {
            A = x;
            f_a = f_x;
        } else {
            B = x;
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
    int maxPow;
    cout << "Enter maximum power of the polynomial: ";
    cin >> maxPow;

    cout << "Enter " << maxPow + 1 << " coefficients: ";
    vector<double> a(maxPow + 1);
    for (int i = 0; i <= maxPow; ++i) {
        cin >> a[i];
    }

    double A, B;
    findAB(a, A, B);

    double f_a = f(a, A);
    double f_b = f(a, B);

    while (f_a * f_b > 0) {
        B += 1.0;
        f_b = f(a, B);
        if (static_cast<int>(A) == static_cast<int>(B)) {
            cout << "The solution may not be found!" << endl;
            return;
        }
    }

    int count = 0;
    bool resultFound = false;
    double x = 0, f_x, res;

    while (count < 500) {
        ++count;
        f_a = f(a, A);
        f_b = f(a, B);
        
        x = (A * f_b - B * f_a) / (f_b - f_a);
        f_x = f(a, x);

        if (abs(f_x) < 0.0001) {
            res = x;
            resultFound = true;
            break;
        }

        if (f_x * f_a > 0) {
            A = x;  
            f_a = f_x;
        } else {
            B = x;  
        }
    }

    if (resultFound){
        cout << "Result: " << (res == 0 ? 0 : res) << endl;
        cout << "Iterations: " << count << endl;
    }
    else {
        cout << "Result not found in 500 iterations." << endl;
    }
}

void newton() {
    int maxPow;
    cout << "Enter maximum power of the polynomial: ";
    cin >> maxPow;

    cout << "Enter " << maxPow + 1 << " coefficients: ";
    vector<double> a(maxPow + 1);
    for (int i = 0; i <= maxPow; ++i) {
        cin >> a[i];
    }

    double x = 0, pre_x = 4;
    int count = 0;
    double res;
    int iteration;
    bool gotResult = false;

    while (count < 500) {
        ++count;
        double f_pre_x = f(a, pre_x);
        if (abs(f_pre_x) < 0.0001 && !gotResult) {
            res = pre_x;
            iteration = count;
            gotResult = true;
            break;
        }

        double df_pre_x = df(a, pre_x);
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
    int maxPow;
    cout << "Enter maximum power of the polynomial: ";
    cin >> maxPow;

    cout << "Enter " << maxPow + 1 << " coefficients: ";
    vector<double> a(maxPow + 1);
    for (int i = 0; i <= maxPow; ++i) {
        cin >> a[i];
    }

    double a_val = 1; 
    double b_val = 10; 
    double x = 0, f_x, f_a, f_b;
    int count = 0;
    double res;
    int iteration;
    bool gotResult = false;

    while (count < 500) {
        ++count; 

        f_a = f(a, a_val);
        f_b = f(a, b_val);
        if (abs(f_a - f_b) < 1e-10) {
            cout << "Function values are too close. No solution found." << endl;
            return;
        }

        x = b_val - (f_b * (b_val - a_val)) / (f_b - f_a);
        f_x = f(a, x);

        if (abs(f_x) < 0.0001 && !gotResult) {
            res = x;
            iteration = count;
            gotResult = true;
            break;
        }
        
        a_val = b_val;
        b_val = x;
    }

    if (gotResult) {
        cout << "Result: " << res << "  Iterations: " << iteration << endl;
    } else {
        cout << "Result not found in 500 iterations." << endl;
    }
}
