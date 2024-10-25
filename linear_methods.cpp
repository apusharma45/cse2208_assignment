#include "linear_methods.h"
#include <bits/stdc++.h>
using namespace std;


bool isDiagonalDominant(const vector<vector<int>>& a, const vector<int>& b) {
    int n = a.size();
    for (int i = 0; i < n; ++i) {
        int diagonal = abs(a[i][i]);
        int sum = 0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                sum += abs(a[i][j]);
            }
        }
        if (diagonal < sum) {
            return false;
        }
    }
    return true;
}

// Calculate the sum of elements of row i, except the kth element
int rowSumExceptK(const vector<vector<int>>& a, int i, int k) {
    int sum = 0;
    int n = a.size();
    for (int j = 0; j < n; ++j) {
        if (j != k) {
            sum += abs(a[i][j]);
        }
    }
    return sum;
}

bool makeDiagonalDominant(vector<vector<int>>& a, vector<int>& b) {
    int n = b.size();
    for (int k = 0; k < n; ++k) {
        int foundRow = -1;
        if (abs(a[k][k]) < rowSumExceptK(a, k, k)) {
            for (int i = k + 1; i < n; ++i) {
                if (abs(a[i][k]) >= rowSumExceptK(a, i, k)) {
                    foundRow = i;
                    break;
                }
            }
            if (foundRow != -1) {
                swap(a[k], a[foundRow]);
                swap(b[k], b[foundRow]);
            }
        }
    }
    return isDiagonalDominant(a, b);
}

void jacobi() {
    int n;
    cout << "Enter the number of variables: ";
    cin >> n;

    vector<vector<int>> a(n, vector<int>(n));
    vector<int> b(n);
    cout << "Enter the coefficients for the equations (a1, b1, ..., an, f for each equation):\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> a[i][j];
        }
        cin >> b[i];
    }

    if (!isDiagonalDominant(a, b)) {
        if (!makeDiagonalDominant(a, b)) {
            cout << "Can't make diagonally dominant" << endl;
            return;
        }
    }

    vector<double> x(n, 0), x1(n, 0);
    int count = 0;

    while (true) {
        count++;
        for (int i = 0; i < n; ++i) {
            double sum = b[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum -= a[i][j] * x[j];
                }
            }
            x1[i] = sum / a[i][i];
        }
        bool convergence = true;
        for (int i = 0; i < n; ++i) {
            if (abs(x[i] - x1[i]) > 0.0001) {
                convergence = false;
                break;
            }
        }
        if (convergence) break;
        x = x1;
    }

    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << ": " << x1[i] << " ";
    }
    cout << "\nResult found in " << count << " iterations\n";
}

void gauss_seidel() {
    int n;
    cout << "Enter the number of variables: ";
    cin >> n;

    vector<vector<int>> a(n, vector<int>(n));
    vector<int> b(n);
    cout << "Enter the coefficients for the equations (a1, b1, ..., an, f for each equation):\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> a[i][j];
        }
        cin >> b[i];
    }

    if (!isDiagonalDominant(a, b)) {
        if (!makeDiagonalDominant(a, b)) {
            cout << "Can't make diagonally dominant" << endl;
            return;
        }
    }

    vector<double> x(n, 0);
    int count = 0;

    while (true) {
        count++;
        vector<double> x_old = x;
        for (int i = 0; i < n; ++i) {
            double sum = b[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum -= a[i][j] * x[j];
                }
            }
            x[i] = sum / a[i][i];
        }
        bool convergence = true;
        for (int i = 0; i < n; ++i) {
            if (abs(x[i] - x_old[i]) > 0.0001) {
                convergence = false;
                break;
            }
        }
        if (convergence) break;
    }

    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << ": " << x[i] << " ";
    }
    cout << "\nResult found in " << count << " iterations\n";
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