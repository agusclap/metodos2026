#include <iostream>
#include <cmath>
using namespace std;

// ============================================================================
// Cuadratura de Gauss-Legendre. Solo sirve con F(x) definida (no con tabla de
// datos): necesita evaluar F en puntos muy puntuales (x0,x1,...) que no
// coinciden con ningun punto de una tabla en general.
//
// Aproxima I=Integral[-1,1] F(x)dx como C0*F(x0)+C1*F(x1)+...+Cn*F(xn). Los
// Cx y los coeficientes salen de pedirle al metodo que ajuste EXACTO
// polinomios hasta grado 2*puntos-1 (tabla 22.1 de Chapra, ya calculada).
//
// Como el intervalo de la formula es fijo [-1,1] y el problema puede pedir
// cualquier [a,b], primero hay que cambiar de variable:
//   x = ((b-a)*x' + (b+a)) / 2      dx = (b-a)/2 dx'
// ============================================================================

int main() {
    auto f = [](double x) { return sin(2*x) * exp(-x); }; // <-- ACA se elige el problema

    double a, b;
    cout << "Ingrese los limites de integracion a y b: ";
    cin >> a >> b;

    int puntos;
    cout << "Ingrese el numero de puntos a usar (entre 2 y 6): ";
    cin >> puntos;

    double I;

    switch (puntos) {
        case 2: {
            double C0 = 1.0,       C1 = 1.0;
            double x0 = -0.5773503, x1 = -x0;
            I = (b - a) / 2.0 * (
                C0 * f(((b - a) * x0 + (b + a)) / 2.0) +
                C1 * f(((b - a) * x1 + (b + a)) / 2.0)
            );
            break;
        }
        case 3: {
            double C0 = 0.5555556, C1 = 0.8888889, C2 = 0.5555556;
            double x0 = -0.7745967, x1 = 0.0,       x2 = -x0;
            I = (b - a) / 2.0 * (
                C0 * f(((b - a) * x0 + (b + a)) / 2.0) +
                C1 * f(((b - a) * x1 + (b + a)) / 2.0) +
                C2 * f(((b - a) * x2 + (b + a)) / 2.0)
            );
            break;
        }
        case 4: {
            double C0 = 0.3478548, C1 = 0.6521452, C2 = 0.6521452, C3 = 0.3478548;
            double x0 = -0.8611363, x1 = -0.3399810, x2 = -x1, x3 = -x0;
            I = (b - a) / 2.0 * (
                C0 * f(((b - a) * x0 + (b + a)) / 2.0) +
                C1 * f(((b - a) * x1 + (b + a)) / 2.0) +
                C2 * f(((b - a) * x2 + (b + a)) / 2.0) +
                C3 * f(((b - a) * x3 + (b + a)) / 2.0)
            );
            break;
        }
        case 5: {
            double C0 = 0.2369269, C1 = 0.4786287, C2 = 0.5688889, C3 = 0.4786287, C4 = 0.2369269;
            double x0 = -0.9061798, x1 = -0.5384693, x2 = 0.0, x3 = -x1, x4 = -x0;
            I = (b - a) / 2.0 * (
                C0 * f(((b - a) * x0 + (b + a)) / 2.0) +
                C1 * f(((b - a) * x1 + (b + a)) / 2.0) +
                C2 * f(((b - a) * x2 + (b + a)) / 2.0) +
                C3 * f(((b - a) * x3 + (b + a)) / 2.0) +
                C4 * f(((b - a) * x4 + (b + a)) / 2.0)
            );
            break;
        }
        case 6: {
            double C0 = 0.1713245, C1 = 0.3607616, C2 = 0.4679139,
                   C3 = 0.4679139, C4 = 0.3607616, C5 = 0.1713245;
            double x0 = -0.9324695, x1 = -0.6612094, x2 = -0.2386192,
                   x3 = -x2, x4 = -x1, x5 = -x0;
            I = (b - a) / 2.0 * (
                C0 * f(((b - a) * x0 + (b + a)) / 2.0) +
                C1 * f(((b - a) * x1 + (b + a)) / 2.0) +
                C2 * f(((b - a) * x2 + (b + a)) / 2.0) +
                C3 * f(((b - a) * x3 + (b + a)) / 2.0) +
                C4 * f(((b - a) * x4 + (b + a)) / 2.0) +
                C5 * f(((b - a) * x5 + (b + a)) / 2.0)
            );
            break;
        }
        default:
            cerr << "Cantidad de puntos invalida (tiene que ser entre 2 y 6)." << endl;
            return 1;
    }

    cout << "La integral aproximada es: I = " << I << endl;
    return 0;
}
