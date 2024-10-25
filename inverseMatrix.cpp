
#include <bits/stdc++.h>
#include"inverseMatrix.h"
using namespace std;

void print_matrix(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}

// LU Decomposition of matrix A
bool lu_decomposition(const vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U) {
    int n = A.size();
    L.assign(n, vector<double>(n, 0));
    U = A;

    for (int i = 0; i < n; i++) {
        L[i][i] = 1;  // Set the diagonal of L to 1
        for (int j = i + 1; j < n; j++) {
            if (U[i][i] == 0) {
                cerr << "Matrix is singular, cannot compute LU decomposition." << endl;
                return false;
            }
            double factor = U[j][i] / U[i][i];
            L[j][i] = factor;
            for (int k = i; k < n; k++) {
                U[j][k] -= factor * U[i][k];
            }
        }
    }
    return true;
}

// Forward substitution to solve LY = B
vector<double> forward_substitution(const vector<vector<double>>& L, const vector<double>& B) {
    int n = L.size();
    vector<double> Y(n, 0);
    for (int i = 0; i < n; i++) {
        Y[i] = B[i];
        for (int j = 0; j < i; j++) {
            Y[i] -= L[i][j] * Y[j];
        }
        Y[i] /= L[i][i];
    }
    return Y;
}

// Back substitution to solve UX = Y
vector<double> back_substitution(const vector<vector<double>>& U, const vector<double>& Y) {
    int n = U.size();
    vector<double> X(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        X[i] = Y[i];
        for (int j = i + 1; j < n; j++) {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }
    return X;
}

// calculating inverse 
void inverse() {
     int n;
    cout << "Enter the dimension : ";
    cin >> n;
    vector<vector<double>> A(n, vector<double>(n));
    cout << "Enter the elements :" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }
   
    vector<vector<double>> L, U, inverse(n, vector<double>(n));
    if (!lu_decomposition(A, L, U)) {
        cerr << "Matrix inversion failed cause of singularity." << endl;
        return ;
    }

    
    for (int i = 0; i < n; i++) {
        vector<double> e(n, 0);
        e[i] = 1; 

        // Solve LY = e 
        vector<double> Y = forward_substitution(L, e);

        // Solve UX = Y 
        vector<double> X = back_substitution(U, Y);

    
        for (int j = 0; j < n; j++) {
            inverse[j][i] = X[j];
        }
    }
    if (!inverse.empty()) {
        cout << "Inverse of the matrix:" << endl;
        print_matrix(inverse);
    } else {
        cout << "Inverse is not possible " << endl;
    }
}
