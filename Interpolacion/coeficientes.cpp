#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

// Interpolacion por coeficientes indeterminados: arma la matriz de
// Vandermonde (A[i][j] = x[i]^j) y resuelve con Gauss + pivoteo parcial
// (mismo algoritmo que en Matrices/eliminacionGausseana.cpp) para obtener
// los coeficientes a_i del polinomio Pn(x) = a0 + a1*x + ... + an*x^n.
int main() {
    string nombreArchivo;
    cout << "Ingrese el nombre del archivo con los puntos (x y por linea): ";
    cin >> nombreArchivo;

    ifstream archivo(nombreArchivo);
    if (!archivo.is_open()) {
        cerr << "No se pudo abrir el archivo." << endl;
        return 1;
    }

    vector<double> x, y;
    double xi, yi;
    while (archivo >> xi >> yi) {
        x.push_back(xi);
        y.push_back(yi);
    }
    archivo.close();

    int n = (int)x.size() - 1;
    if (n < 0) {
        cerr << "El archivo no tiene puntos validos." << endl;
        return 1;
    }
    int m = n + 1; // tamano del sistema (n+1 incognitas: a0..an)
    cout << "Se leyeron " << m << " puntos (polinomio de grado " << n << ")." << endl;

    // ---- Armar A y b (matriz de Vandermonde) ----
    vector<vector<double>> A(m, vector<double>(m));
    vector<double> b(m);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            A[i][j] = pow(x[i], j);
        }
        b[i] = y[i];
    }

    // ---- Resolver el sistema: Gauss con pivoteo parcial ----
    for (int i = 0; i < m - 1; i++) {
        int p = i;
        if (fabs(A[i][i]) < 1e-4) {
            double mmax = fabs(A[i][i]);
            for (int l = i + 1; l < m; l++) {
                if (fabs(A[l][i]) > mmax) { mmax = fabs(A[l][i]); p = l; }
            }
            for (int l = 0; l < m; l++) { double aux = A[p][l]; A[p][l] = A[i][l]; A[i][l] = aux; }
            double aux = b[p]; b[p] = b[i]; b[i] = aux;
        }
        for (int j = i + 1; j < m; j++) {
            double factor = -A[j][i] / A[i][i];
            for (int k = i + 1; k < m; k++) A[j][k] += A[i][k] * factor;
            b[j] += b[i] * factor;
            A[j][i] = 0.0;
        }
    }

    vector<double> a(m);
    a[m - 1] = b[m - 1] / A[m - 1][m - 1];
    for (int i = m - 2; i >= 0; i--) {
        double suma = b[i];
        for (int j = i + 1; j < m; j++) suma -= A[i][j] * a[j];
        a[i] = suma / A[i][i];
    }

    cout << "Coeficientes del polinomio (Pn(x) = a0 + a1*x + ... + an*x^n):" << endl;
    for (int i = 0; i < m; i++) cout << "a" << i << " = " << a[i] << endl;

    // ---- Evaluar en un punto ----
    double xhat;
    cout << "Ingrese el valor x a evaluar: ";
    cin >> xhat;

    double suma = 0;
    for (int i = 0; i < m; i++) suma += a[i] * pow(xhat, i);
    cout << "P(" << xhat << ") = " << suma << endl;

    return 0;
}
