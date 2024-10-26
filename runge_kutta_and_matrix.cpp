//Differential Equations Solution Using 4th Order Runge-Kutta Method
#include <bits/stdc++.h>
using namespace std;

typedef vector<double> Vec;
typedef function<Vec(double, const Vec&)> DerivativeFunction;

Vec addVectors(const Vec& v1, const Vec& v2) {
    Vec result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i)
        result[i] = v1[i] + v2[i];
    return result;
}

Vec scaleVector(const Vec& vec, double factor) {
    Vec result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        result[i] = vec[i] * factor;
    return result;
}

Vec rungeKutta4Step(DerivativeFunction func, const Vec& state, double time, double stepSize) {
    Vec k1 = func(time, state);
    Vec k2 = func(time + 0.5 * stepSize, addVectors(state, scaleVector(k1, 0.5 * stepSize)));
    Vec k3 = func(time + 0.5 * stepSize, addVectors(state, scaleVector(k2, 0.5 * stepSize)));
    Vec k4 = func(time + stepSize, addVectors(state, scaleVector(k3, stepSize)));

    Vec nextState = state;
    for (size_t i = 0; i < state.size(); ++i) {
        nextState[i] += stepSize / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }

    return nextState;
}

void rungeKutta4() {
    DerivativeFunction diffEq = [](double t, const Vec& y) -> Vec {
        return {-y[0]};
    };

    double initialCondition, startTime, endTime, stepSize;

    cout << "Enter initial condition y(0): ";
    cin >> initialCondition;

    cout << "Enter start time t0: ";
    cin >> startTime;

    cout << "Enter end time tEnd: ";
    cin >> endTime;

    cout << "Enter step size h: ";
    cin >> stepSize;

    Vec currentState = {initialCondition};
    double currentTime = startTime;
    while (currentTime < endTime) {
        if (currentTime + stepSize > endTime) stepSize = endTime - currentTime;
        currentState = rungeKutta4Step(diffEq, currentState, currentTime, stepSize);
        currentTime += stepSize;
    }

    cout << "y(" << endTime << ") = " << currentState[0] << endl;
}

