#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

// ============================================================================
// Igual que gaussSeidel.cpp, pero con el criterio de corte del Problema 3:
// el residuo ||A*x - b|| en vez de la diferencia sucesiva ||Xnuevo-Xviejo||.
// ============================================================================

double calcularResiduo(const vector<vector<double>> &A, const vector<double> &b,
                        const vector<double> &X, int n) {
    double suma = 0.0;
    for (int i = 0; i < n; i++) {
        double Ax_i = 0.0;
        for (int j = 0; j < n; j++) {
            Ax_i += A[i][j] * X[j];
        }
        double r_i = Ax_i - b[i];
        suma += r_i * r_i;
    }
    return sqrt(suma);
}

int main() {

    // ---- Definicion de A y b (desde archivo de texto) ----
    string nombreArchivo;
    cout << "Ingrese el nombre del archivo con la matriz aumentada [A|b]: ";
    cin >> nombreArchivo;

    ifstream archivo(nombreArchivo);
    if (!archivo.is_open()) {
        cerr << "No se pudo abrir el archivo." << endl;
        return 1;
    }

    vector<double> datos;
    double valor;
    while (archivo >> valor) {
        datos.push_back(valor);
    }
    archivo.close();

    int total = datos.size();
    int n = static_cast<int>((-1 + sqrt(1.0 + 4.0 * total)) / 2);
    if (n <= 0 || n * (n + 1) != total) {
        cerr << "El archivo no tiene el formato esperado (n filas x n+1 columnas)." << endl;
        return 1;
    }

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);
    int idx = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = datos[idx++];
        }
        b[i] = datos[idx++];
    }

    // ---- Verificar si A es diagonalmente dominante ----
    bool esDD = true;
    for (int i = 0; i < n; i++) {
        double suma = 0;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                suma += fabs(A[i][j]);
            }
        }
        if (fabs(A[i][i]) < suma) {
            cout << "La matriz no es diagonalmente dominante (fila " << i + 1 << ")" << endl;
            esDD = false;
        }
        if (A[i][i] == 0) {
            cerr << "a" << i + 1 << i + 1 << " = 0, no se puede dividir. Abortando." << endl;
            return 1;
        }
    }
    if (!esDD) {
        cout << "Aviso: al no ser diagonalmente dominante, no esta garantizada la convergencia." << endl;
    }

    vector<double> X(n, 0.0);

    double tolerancia;
    cout << "Ingrese la tolerancia: ";
    cin >> tolerancia;

    double w;
    cout << "Ingrese el coeficiente de relajacion w (1 = Gauss-Seidel sin relajacion): ";
    cin >> w;

    double residuo_viejo = 1000;
    double residuo;
    int iter = 0;

    do {
        iter++;

        for (int i = 0; i < n; i++) {
            double suma = b[i];
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    suma = suma - A[i][j] * X[j];
                }
            }
            double xGS = suma / A[i][i];
            double xViejo_i = X[i];
            X[i] = w * xGS + (1 - w) * xViejo_i;
        }

        // ---- Criterio de corte del Problema 3: residuo ||A*X - b|| ----
        residuo = calcularResiduo(A, b, X, n);

        if (residuo > residuo_viejo) {
            cout << "El metodo no converge." << endl;
            return 1;
        }
        residuo_viejo = residuo;

        cout << "Iteracion " << iter << ": residuo = " << residuo << endl;

    } while (residuo > tolerancia);

    cout << "Solucion:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << X[i] << endl;
    }
    cout << "Residuo final: " << residuo << ", iteraciones: " << iter << endl;

    return 0;
}
