#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
using namespace std;
using namespace std::chrono;

// ---- Ingresar A y b (desde archivo de texto) ----
// Mismo formato que los demas: matriz aumentada [A|b], n filas, n+1 columnas.
bool leerSistema(vector<vector<double>> &A, vector<double> &b, int &n) {
    string nombreArchivo;
    cout << "Ingrese el nombre del archivo con la matriz aumentada [A|b]: ";
    cin >> nombreArchivo;

    ifstream archivo(nombreArchivo);
    if (!archivo.is_open()) {
        cerr << "No se pudo abrir el archivo." << endl;
        return false;
    }

    vector<double> datos;
    double valor;
    while (archivo >> valor) {
        datos.push_back(valor);
    }
    archivo.close();

    int total = datos.size();
    n = static_cast<int>((-1 + sqrt(1.0 + 4.0 * total)) / 2);
    if (n <= 0 || n * (n + 1) != total) {
        cerr << "El archivo no tiene el formato esperado (n filas x n+1 columnas)." << endl;
        return false;
    }

    A.assign(n, vector<double>(n));
    b.assign(n, 0.0);
    int idx = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = datos[idx++];
        }
        b[i] = datos[idx++];
    }
    return true;
}

// ---- Ancho de banda de A ----
// p = max{ |i-j| } para todo (i,j) con A[i][j] != 0.
// Funciona para cualquier matriz, no asume nada sobre su estructura.
int anchoBanda(const vector<vector<double>> &A, int n) {
    int p = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (A[i][j] != 0.0) {
                int dist = abs(i - j);
                if (dist > p) p = dist;
            }
        }
    }
    return p;
}

// ---- Gauss-Seidel "denso": recorre TODAS las columnas de cada fila ----
// (lo que ya veniamos haciendo hasta ahora, sin aprovechar la banda)
vector<double> gaussSeidelDenso(const vector<vector<double>> &A, const vector<double> &b,
                                 int n, double tol, int &iterOut) {
    vector<double> X(n, 0.0);
    double error_viejo = 1000, error;
    int iter = 0;
    do {
        iter++;
        vector<double> Xviejo = X;
        for (int i = 0; i < n; i++) {
            double suma = b[i];
            for (int j = 0; j < n; j++) {
                if (j != i) suma -= A[i][j] * X[j];
            }
            X[i] = suma / A[i][i];
        }
        error = 0;
        for (int i = 0; i < n; i++) error += (X[i] - Xviejo[i]) * (X[i] - Xviejo[i]);
        error = sqrt(error);
        if (error > error_viejo) { iterOut = -1; return X; }
        error_viejo = error;
    } while (error > tol);
    iterOut = iter;
    return X;
}

// ---- Gauss-Seidel "de banda": solo recorre columnas dentro de [i-p, i+p] ----
// Fuera de ese rango, A[i][j] es cero por construccion -> ni vale la pena
// mirarlo. Por fila, hace ~(2p+1) cuentas en vez de n.
vector<double> gaussSeidelBanda(const vector<vector<double>> &A, const vector<double> &b,
                                 int n, int p, double tol, int &iterOut) {
    vector<double> X(n, 0.0);
    double error_viejo = 1000, error;
    int iter = 0;
    do {
        iter++;
        vector<double> Xviejo = X;
        for (int i = 0; i < n; i++) {
            double suma = b[i];
            int jIni = max(0, i - p);
            int jFin = min(n - 1, i + p);
            for (int j = jIni; j <= jFin; j++) {
                if (j != i) suma -= A[i][j] * X[j];
            }
            X[i] = suma / A[i][i];
        }
        error = 0;
        for (int i = 0; i < n; i++) error += (X[i] - Xviejo[i]) * (X[i] - Xviejo[i]);
        error = sqrt(error);
        if (error > error_viejo) { iterOut = -1; return X; }
        error_viejo = error;
    } while (error > tol);
    iterOut = iter;
    return X;
}

int main() {
    vector<vector<double>> A;
    vector<double> b;
    int n;
    if (!leerSistema(A, b, n)) return 1;

    int p = anchoBanda(A, n);
    cout << "n = " << n << ", ancho de banda detectado = " << p << endl;

    double tol;
    cout << "Ingrese la tolerancia: ";
    cin >> tol;

    // Se repite muchas veces para que el tiempo sea medible (con n chico,
    // una sola resolucion tarda microsegundos, poco confiable para comparar).
    const int REPS = 500;
    int iterDenso = 0, iterBanda = 0;
    vector<double> xDenso, xBanda;

    auto t1 = high_resolution_clock::now();
    for (int r = 0; r < REPS; r++) xDenso = gaussSeidelDenso(A, b, n, tol, iterDenso);
    auto t2 = high_resolution_clock::now();
    for (int r = 0; r < REPS; r++) xBanda = gaussSeidelBanda(A, b, n, p, tol, iterBanda);
    auto t3 = high_resolution_clock::now();

    double usDenso = duration_cast<microseconds>(t2 - t1).count() / double(REPS);
    double usBanda = duration_cast<microseconds>(t3 - t2).count() / double(REPS);

    cout << "\n--- Solucion (denso), " << iterDenso << " iteraciones ---" << endl;
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << xDenso[i] << endl;

    cout << "\n--- Solucion (banda), " << iterBanda << " iteraciones ---" << endl;
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << xBanda[i] << endl;

    cout << "\n--- Comparacion de tiempo (promedio de " << REPS << " corridas) ---" << endl;
    cout << "Denso: " << usDenso << " us por resolucion" << endl;
    cout << "Banda: " << usBanda << " us por resolucion" << endl;
    cout << "Speedup: " << (usDenso / usBanda) << "x" << endl;

    return 0;
}
