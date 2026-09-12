#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

// ============================================================================
// Igual que lagrangeFuncion.cpp (se conoce f(x) real, se puede calcular el
// error exacto), pero para cuando el enunciado da nodos PUNTUALES y no
// equiespaciados (una tabla), en vez de pedir generarlos automaticamente.
//
// Cubre los 3 items de este tipo de problema en una sola corrida:
//   a) arma el polinomio de Lagrange con los nodos de la tabla (implicito:
//      se arma adentro de la funcion lagrange(), no hace falta expandirlo)
//   b) lo evalua en un punto puntual y da el error exacto ahi
//   c) lo evalua en una grilla fina (solo dentro de [x0,xN]) y exporta un
//      CSV con x, f_exacta, P_n, error para graficar e(x)
// ============================================================================

// ---- ACA se elige el problema: la funcion conocida ----
double f(double x) {
    return log(x * x + 1) - sin(x);   // misma f() del Problema 1
}

double lagrange(const vector<double> &x, const vector<double> &y, double xhat) {
    int n = (int)x.size() - 1;
    double suma = 0.0;
    for (int k = 0; k <= n; k++) {
        double prod = 1.0;
        for (int i = 0; i <= n; i++) {
            if (i != k) prod *= (xhat - x[i]) / (x[k] - x[i]);
        }
        suma += y[k] * prod;
    }
    return suma;
}

int main() {
    // ---- Leer los nodos x (uno por linea, no equiespaciados) ----
    string nombreArchivo;
    cout << "Ingrese el archivo con los nodos x (uno por linea): ";
    cin >> nombreArchivo;

    ifstream archivo(nombreArchivo);
    if (!archivo.is_open()) {
        cerr << "No se pudo abrir el archivo." << endl;
        return 1;
    }

    vector<double> xNodos;
    double xi;
    while (archivo >> xi) xNodos.push_back(xi);
    archivo.close();

    int n = (int)xNodos.size() - 1;
    if (n < 1) {
        cerr << "Hacen falta al menos 2 nodos." << endl;
        return 1;
    }

    // ---- a) los y_i salen de evaluar f() en cada nodo dado ----
    vector<double> yNodos(n + 1);
    for (int i = 0; i <= n; i++) yNodos[i] = f(xNodos[i]);

    cout << "Nodos leidos (n=" << n << ", polinomio de grado " << n << "):" << endl;
    for (int i = 0; i <= n; i++) {
        cout << "  x" << i << " = " << xNodos[i] << "   y" << i << " = " << yNodos[i] << endl;
    }

    // ---- b) evaluar en un punto puntual, con error exacto ----
    double xhat;
    cout << "\nIngrese el punto a evaluar (item b): ";
    cin >> xhat;
    double Pxhat = lagrange(xNodos, yNodos, xhat);
    double fxhat = f(xhat);
    cout << "P(" << xhat << ") = " << Pxhat << endl;
    cout << "f(" << xhat << ") = " << fxhat << endl;
    cout << "Error absoluto exacto = |f(x)-P(x)| = " << fabs(fxhat - Pxhat) << endl;

    // ---- c) grilla fina dentro de [x0,xN] para graficar e(x) ----
    string salida;
    cout << "\nNombre del archivo CSV para el grafico de e(x) (item c): ";
    cin >> salida;

    double xIni = xNodos[0], xFin = xNodos[n];
    const int npuntosGrilla = 400;
    ofstream out(salida);
    if (!out.is_open()) {
        cerr << "No se pudo crear " << salida << endl;
        return 1;
    }
    out << setprecision(12);
    out << "x,f_exacta,P_n,error\n";
    for (int i = 0; i < npuntosGrilla; i++) {
        double x = xIni + i * (xFin - xIni) / (npuntosGrilla - 1);
        double fx = f(x);
        double px = lagrange(xNodos, yNodos, x);
        out << x << "," << fx << "," << px << "," << fabs(fx - px) << "\n";
    }
    out.close();
    cout << "Listo: " << salida << " generado (para graficar e(x) y ver que pasa cerca de los nodos)." << endl;

    return 0;
}
