#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
using namespace std;

// Regresion por cuadrados minimos con funciones base arbitrarias:
// F(x) = a0*phi0(x) + a1*phi1(x) + ... + ap*phip(x)
//
// Derivando Sr respecto a cada a_k SIEMPRE da el mismo patron, sea cual sea
// la funcion base (es la generalizacion de A_lm=sum(x^(l+m)) de polinomios):
//   A[k][j] = sum( phi_k(x_i) * phi_j(x_i) )
//   b[k]    = sum( y_i * phi_k(x_i) )
//
// Lo UNICO que cambia de problema a problema son las funciones base de abajo
// (linea marcada "<-- ACA se elige el problema"). El resto es generico.

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

    int npuntos = (int)x.size();
    cout << "Se leyeron " << npuntos << " puntos." << endl;

    // ---- ACA se elige el problema: las funciones base ----
    // Ejemplo: F(x) = a*sen(x) + b*cos(x)  (P4 - Parcial 1, 2023)
    vector<function<double(double)>> phi = {
        [](double t) { return sin(t); },
        [](double t) { return cos(t); }
    };
    int m = (int)phi.size();

    if (npuntos < m) {
        cerr << "No hay suficientes puntos para " << m << " funciones base." << endl;
        return 1;
    }

    // ---- Armar A y b (genera para CUALQUIER phi, sin tocar nada de aca para abajo) ----
    vector<vector<double>> A(m, vector<double>(m));
    vector<double> b(m);
    for (int k = 0; k < m; k++) {
        double sumaB = 0;
        for (int i = 0; i < npuntos; i++) sumaB += y[i] * phi[k](x[i]);
        b[k] = sumaB;

        for (int j = 0; j < m; j++) {
            double sumaA = 0;
            for (int i = 0; i < npuntos; i++) sumaA += phi[k](x[i]) * phi[j](x[i]);
            A[k][j] = sumaA;
        }
    }

    // Se muestran A y b ANTES de resolver -- son las que hay que comparar
    // contra la matriz derivada a mano en la hoja.
    cout << "\nMatriz A (comparar contra la derivada a mano):" << endl;
    for (int k = 0; k < m; k++) {
        for (int j = 0; j < m; j++) cout << A[k][j] << "\t";
        cout << "| " << b[k] << endl;
    }

    // ---- Resolver el sistema: Gauss con pivoteo parcial (mismo de siempre) ----
    for (int i = 0; i < m - 1; i++) {
        int piv = i;
        if (fabs(A[i][i]) < 1e-9) {
            double mmax = fabs(A[i][i]);
            for (int l = i + 1; l < m; l++) {
                if (fabs(A[l][i]) > mmax) { mmax = fabs(A[l][i]); piv = l; }
            }
            for (int l = 0; l < m; l++) { double aux = A[piv][l]; A[piv][l] = A[i][l]; A[i][l] = aux; }
            double aux = b[piv]; b[piv] = b[i]; b[i] = aux;
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

    cout << "\nCoeficientes:" << endl;
    for (int i = 0; i < m; i++) cout << "a" << i << " = " << a[i] << endl;

    // ---- Coeficiente de correlacion (F(xi) = suma_k a[k]*phi[k](xi)) ----
    double sumaY = 0;
    for (int i = 0; i < npuntos; i++) sumaY += y[i];
    double yprom = sumaY / npuntos;

    double sumaT = 0, sumaR = 0;
    for (int i = 0; i < npuntos; i++) {
        sumaT += pow(y[i] - yprom, 2);
        double fi = 0;
        for (int k = 0; k < m; k++) fi += a[k] * phi[k](x[i]);
        sumaR += pow(y[i] - fi, 2);
    }
    double r = sqrt((sumaT - sumaR) / sumaT);

    cout << "\nSt = " << sumaT << ", Sr = " << sumaR << endl;
    cout << "Coeficiente de correlacion r = " << r << endl;

    return 0;
}
