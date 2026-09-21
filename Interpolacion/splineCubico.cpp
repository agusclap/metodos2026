#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

// ============================================================================
// Spline cubico NATURAL, resuelto armando el sistema completo de 4n incognitas
// y resolviendolo con Eliminacion Gaussiana -- pseudocodigo tal cual lo dio
// la profesora en clase (15/9). No es el metodo tridiagonal mas eficiente:
// se plantean las 4n incognitas (a_k,b_k,c_k,d_k de cada uno de los n tramos)
// juntas, en una matriz A[4n][4n], y se resuelven todas de una.
//
// F_k(x) = a_k*x^3 + b_k*x^2 + c_k*x + d_k      (tramo k, k=0,...,n-1)
// z[4k]=a_k, z[4k+1]=b_k, z[4k+2]=c_k, z[4k+3]=d_k
//
// Condiciones armadas (en ese orden, filas de A):
//   Bloque 1 (filas 0..2n-1):   F_k(x_k)=y_k, F_k(x_k+1)=y_k+1   (pasa por los puntos)
//   Bloque 2 (filas 2n..3n-2):  F'_k(x_k+1) = F'_k+1(x_k+1)      (derivada 1ra continua)
//   Bloque 3 (filas 3n-1..4n-3): F''_k(x_k+1) = F''_k+1(x_k+1)   (derivada 2da continua)
//   Bloque 4 (filas 4n-2,4n-1): F''_0(x_0)=0, F''_n-1(x_n)=0      (spline NATURAL)
// ============================================================================

bool leerPuntos(vector<double> &x, vector<double> &y) {
    string nombreArchivo;
    cout << "Ingrese el nombre del archivo con los puntos (x y por linea): ";
    cin >> nombreArchivo;

    ifstream archivo(nombreArchivo);
    if (!archivo.is_open()) {
        cerr << "No se pudo abrir el archivo." << endl;
        return false;
    }
    double xi, yi;
    while (archivo >> xi >> yi) {
        x.push_back(xi);
        y.push_back(yi);
    }
    archivo.close();
    return true;
}

// ---- Eliminacion Gaussiana con pivoteo parcial (mismo de siempre) ----
vector<double> eliminacionGaussiana(vector<vector<double>> A, vector<double> b, int m) {
    for (int i = 0; i < m - 1; i++) {
        int piv = i;
        if (fabs(A[i][i]) < 1e-9) {
            cout << "Pivoteo necesario: |A[" << i+1 << "][" << i+1 << "]| = "
                 << fabs(A[i][i]) << " es muy chico (< 1e-9)." << endl;
            double mmax = fabs(A[i][i]);
            for (int l = i + 1; l < m; l++) {
                if (fabs(A[l][i]) > mmax) { mmax = fabs(A[l][i]); piv = l; }
            }
            if (piv != i) {
                cout << "  -> se intercambia la fila " << i+1 << " con la fila " << piv+1 << endl;
            } else {
                cout << "  -> no se encontro una fila mejor, se sigue con el mismo pivote" << endl;
            }
            for (int l = 0; l < m; l++) { double aux = A[piv][l]; A[piv][l] = A[i][l]; A[i][l] = aux; }
            double aux = b[piv]; b[piv] = b[i]; b[i] = aux;
        }
        for (int j = i + 1; j < m; j++) {
            double factor = -A[j][i] / A[i][i];
            for (int k = i; k < m; k++) A[j][k] += A[i][k] * factor;
            b[j] += b[i] * factor;
        }
    }

    vector<double> z(m);
    z[m - 1] = b[m - 1] / A[m - 1][m - 1];
    for (int i = m - 2; i >= 0; i--) {
        double suma = b[i];
        for (int j = i + 1; j < m; j++) suma -= A[i][j] * z[j];
        z[i] = suma / A[i][i];
    }
    return z;
}

int main() {
    vector<double> x, y;
    if (!leerPuntos(x, y)) return 1;

    int npuntos = (int)x.size();
    int n = npuntos - 1; // (n+1) nodos -> n tramos
    if (n < 2) {
        cerr << "Hacen falta al menos 3 puntos (2 tramos) para que el spline tenga sentido." << endl;
        return 1;
    }
    int m = 4 * n; // tamano del sistema: 4n incognitas (a_k,b_k,c_k,d_k por tramo)
    cout << "Se leyeron " << npuntos << " puntos (" << n << " tramos, sistema " << m << "x" << m << ")." << endl;

    vector<vector<double>> A(m, vector<double>(m, 0.0));
    vector<double> b(m, 0.0);

    // ---- Bloque 1: F_k(x_k)=y_k, F_k(x_k+1)=y_k+1 -- filas 0 a 2n-1 ----
    for (int k = 0; k < n; k++) {
        for (int j = 0; j <= 3; j++) {
            A[2*k][4*k + j]     = pow(x[k],   3 - j);
            A[2*k + 1][4*k + j] = pow(x[k+1], 3 - j);
        }
        b[2*k]     = y[k];
        b[2*k + 1] = y[k+1];
    }

    // ---- Bloque 2: F'_k(x_k+1) = F'_k+1(x_k+1) -- filas 2n a 3n-2 ----
    for (int k = 0; k <= n - 2; k++) {
        int i = 2*n + k;
        for (int j = 0; j <= 2; j++) {
            A[i][4*k + j]     = (3 - j) * pow(x[k+1], 2 - j);
            A[i][4*(k+1) + j] = -(3 - j) * pow(x[k+1], 2 - j);
        }
        b[i] = 0.0;
    }

    // ---- Bloque 3: F''_k(x_k+1) = F''_k+1(x_k+1) -- filas 3n-1 a 4n-3 ----
    // (ecuacion ya dividida por 2, tal como la dejo la profesora: coeficientes
    // 3*x y 1 en vez de 6*x y 2 -- son la misma ecuacion, mas simple de escribir)
    for (int k = 0; k <= n - 2; k++) {
        int i = 3*n - 1 + k;
        A[i][4*k]         = 3 * x[k+1];
        A[i][4*k + 1]     = 1.0;
        A[i][4*(k+1)]     = -3 * x[k+1];
        A[i][4*(k+1) + 1] = -1.0;
        b[i] = 0.0;
    }

    // ---- Bloque 4: spline NATURAL, F''=0 en los dos extremos -- filas 4n-2, 4n-1 ----
    A[4*n - 2][0] = 3 * x[0];
    A[4*n - 2][1] = 1.0;
    b[4*n - 2] = 0.0;

    A[4*n - 1][4*n - 4] = 3 * x[n];
    A[4*n - 1][4*n - 3] = 1.0;
    b[4*n - 1] = 0.0;

    // ---- Resolver el sistema llamando a Eliminacion Gaussiana -> z (4n componentes) ----
    vector<double> z = eliminacionGaussiana(A, b, m);

    cout << "\nCoeficientes por tramo (F_k(x) = a_k*x^3 + b_k*x^2 + c_k*x + d_k):" << endl;
    for (int k = 0; k < n; k++) {
        cout << "  Tramo " << k << "  [" << x[k] << ", " << x[k+1] << "]:  "
             << "a=" << z[4*k] << "  b=" << z[4*k+1]
             << "  c=" << z[4*k+2] << "  d=" << z[4*k+3] << endl;
    }

    // ---- Ingresar x_hat y localizar a que subintervalo pertenece (determina k) ----
    double xhat;
    cout << "\nIngrese el valor x a interpolar: ";
    cin >> xhat;

    // OJO: igual que en el pseudocodigo, la comparacion es estricta (> y <),
    // asi que si x_hat cae EXACTO sobre un nodo de la tabla, este loop no lo
    // encuentra y no se imprime nada -- pedir un valor estrictamente adentro
    // de algun tramo, no un nodo exacto.
    for (int i = 0; i < n; i++) {
        if (xhat > x[i] && xhat < x[i+1]) {
            double yhat = z[4*i]*pow(xhat,3) + z[4*i+1]*pow(xhat,2) + z[4*i+2]*xhat + z[4*i+3];
            cout << "El valor interpolado es: " << yhat << endl;
            return 0;
        }
    }
    cerr << "x_hat no cayo estrictamente adentro de ningun tramo (fuera de rango, o cayo justo sobre un nodo)." << endl;
    return 1;
}
