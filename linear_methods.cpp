#include "linear_methods.h"
#include <bits/stdc++.h>
using namespace std;

//Jacobi iterative method
void jacobi() {
    vector<int> a(5), b(5), c(5), d(5), e(5), f(5);

    cout << "Enter the coefficients for 5 equations (a1, b1, c1, d1, e1, f1 for equation 1, ...):\n";
    for (int i = 0; i < 5; ++i) {
        cin >> a[i] >> b[i] >> c[i] >> d[i] >> e[i] >> f[i];
    }

    double x = 0, y = 0, z = 0, w = 0, v = 0;
    double x1, y1, z1, w1, v1;
    int count = 0;

    while (true) {
        count++;
        x1 = (f[0] - b[0] * y - c[0] * z - d[0] * w - e[0] * v) / a[0];
        y1 = (f[1] - a[1] * x - c[1] * z - d[1] * w - e[1] * v) / b[1];
        z1 = (f[2] - a[2] * x - b[2] * y - d[2] * w - e[2] * v) / c[2];
        w1 = (f[3] - a[3] * x - b[3] * y - c[3] * z - e[3] * v) / d[3];
        v1 = (f[4] - a[4] * x - b[4] * y - c[4] * z - d[4] * w) / e[4];

        if (abs(x - x1) <= 0.0001 && abs(y - y1) <= 0.0001 && 
            abs(z - z1) <= 0.0001 && abs(w - w1) <= 0.0001 && 
            abs(v - v1) <= 0.0001) {
            break;
        }

        x = x1; 
        y = y1; 
        z = z1; 
        w = w1; 
        v = v1;
    }

    cout << "x: " << x1 << "  y: " << y1 << "  z: " << z1 << "  w: " << w1 << "  v: " << v1 << endl;
    cout << "Result found in " << count << " iterations" << endl;
}

void gauss_seidel() {
    vector<int> a(5), b(5), c(5), d(5), e(5), f(5);

    cout << "Enter the coefficients for 5 equations (a1, b1, c1, d1, e1, f1 for equation 1, ...):\n";
    for (int i = 0; i < 5; ++i) {
        cin >> a[i] >> b[i] >> c[i] >> d[i] >> e[i] >> f[i];
    }

    double x = 0, y = 0, z = 0, w = 0, v = 0;
    double x1, y1, z1, w1, v1;
    int count = 0;

    while (true) {
        count++;
        x1 = (f[0] - b[0] * y - c[0] * z - d[0] * w - e[0] * v) / a[0];
        y1 = (f[1] - a[1] * x1 - c[1] * z - d[1] * w - e[1] * v) / b[1];
        z1 = (f[2] - a[2] * x1 - b[2] * y1 - d[2] * w - e[2] * v) / c[2];
        w1 = (f[3] - a[3] * x1 - b[3] * y1 - c[3] * z - e[3] * v) / d[3];
        v1 = (f[4] - a[4] * x1 - b[4] * y1 - c[4] * z - d[4] * w) / e[4];

        if (abs(x - x1) <= 0.0001 && abs(y - y1) <= 0.0001 && 
            abs(z - z1) <= 0.0001 && abs(w - w1) <= 0.0001 && 
            abs(v - v1) <= 0.0001) {
            break;
        }

        x = x1; 
        y = y1; 
        z = z1; 
        w = w1; 
        v = v1;
    }

    cout << "x: " << x1 << "  y: " << y1 << "  z: " << z1 << "  w: " << w1 << "  v: " << v1 << endl;
    cout << "Result found in " << count << " iterations" << endl;
}

void gaussElimination() {
    int n;
    cout << "Enter the number of variables: ";
    cin >> n;

    vector<vector<double>> a(n, vector<double>(n));
    vector<double> b(n);

    cout << "Enter the coefficients of the equations (row by row):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> a[i][j];
        }
        cin >> b[i];
    }

    for (int k = 0; k < n - 1; k++) {
        int maxIndex = k;
        for (int i = k + 1; i < n; i++) {
            if (abs(a[i][k]) > abs(a[maxIndex][k])) {
                maxIndex = i;
            }
        }

        if (maxIndex != k) {
            swap(a[k], a[maxIndex]);
            swap(b[k], b[maxIndex]);
        }

        if (a[k][k] != 0.0) {
            for (int i = k + 1; i < n; i++) {
                double lam = a[i][k] / a[k][k];
                for (int j = k + 1; j < n; j++) {
                    a[i][j] -= lam * a[k][j];
                }
                b[i] -= lam * b[k];
            }
        } else {
            cout << "Matrix is singular or nearly singular!" << endl;
            return;
        }
    }

    for (int k = n - 1; k >= 0; k--) {
        for (int j = k + 1; j < n; j++) {
            b[k] -= a[k][j] * b[j];
        }
        b[k] /= a[k][k];
    }

    cout << endl;
    cout << fixed << setprecision(6) << "Solution: ";
    for (double x : b) {
        cout << x << " ";
    }
    cout << endl;
}

void gaussJordan() {
    int n;
    cout << "Enter the number of variables: ";
    cin >> n;

    vector<vector<double>> a(n, vector<double>(n));
    vector<double> b(n);

    cout << "Enter the coefficients of the equations (row by row):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> a[i][j];
        }
        cin >> b[i];
    }

    for (int k = 0; k < n; k++) {
        int maxIndex = k;
        for (int i = k + 1; i < n; i++) {
            if (abs(a[i][k]) > abs(a[maxIndex][k])) {
                maxIndex = i;
            }
        }

        if (maxIndex != k) {
            swap(a[k], a[maxIndex]);
            swap(b[k], b[maxIndex]);
        }

        double pivot = a[k][k];
        if (pivot != 0) {
            for (int j = 0; j < n; j++) {
                a[k][j] /= pivot;
            }
            b[k] /= pivot;
        } else {
            cout << "Matrix is singular or nearly singular!" << endl;
            return;
        }

        for (int i = 0; i < n; i++) {
            if (i != k) {
                double factor = a[i][k];
                for (int j = 0; j < n; j++) {
                    a[i][j] -= factor * a[k][j];
                }
                b[i] -= factor * b[k];
            }
        }
    }

    cout << endl;
    cout << fixed << setprecision(4) << "Solution: ";
    for (double x : b) {
        cout << x << " ";
    }
    cout << endl << endl;
}

void lu_factor() {
    
}