#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

// ============================================================================
// Variante "mejorada" de gaussSeidel.cpp, para casos puntuales donde el
// chequeo de divergencia del pseudocodigo de la profesora (error > error
// de la iteracion anterior) da un falso "no converge".
//
// Eso pasa cuando el radio espectral de la matriz de iteracion esta MUY
// cerca de 1: la convergencia sigue siendo real, pero tan lenta y no
// estrictamente monotona que un repunte chico y puntual del error dispara
// el chequeo original sin que el metodo haya divergido de verdad.
//
// OJO: esta version NO es la que hay que entregar como respuesta de la
// materia -- el algoritmo dado en clase es el de gaussSeidel.cpp, tal cual.
// Esta es solo para explorar/confirmar el comportamiento real de un sistema
// puntual (por ejemplo, cuando gaussSeidel.cpp da "no converge" y da dudas
// si es un verdadero no-converge o una convergencia demasiado lenta).
//
// Cambios respecto a gaussSeidel.cpp:
//   1) Solo declara divergencia si el error explota de verdad (>1e10 o
//      no finito), en vez de comparar contra la iteracion anterior.
//   2) max_iter mucho mas alto (5000), para darle margen a una
//      convergencia lenta.
//   3) Si se llega a max_iter sin alcanzar la tolerancia, avisa que el
//      resultado es parcial (en vez de mentir con un "convergio bien").
// ============================================================================

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

    // ---- Verificar si A es diagonalmente dominante (informativo) ----
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
    vector<double> Xviejo(n, 0.0);

    double tolerancia;
    cout << "Ingrese la tolerancia: ";
    cin >> tolerancia;

    double w;
    cout << "Ingrese el coeficiente de relajacion w (1 = Gauss-Seidel sin relajacion): ";
    cin >> w;

    double error;
    int iter = 0, max_iter = 5000;

    do {
        iter++;
        Xviejo = X;

        for (int i = 0; i < n; i++) {
            double suma = b[i];
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    suma = suma - A[i][j] * X[j];
                }
            }
            double xGS = suma / A[i][i];
            X[i] = w * xGS + (1 - w) * Xviejo[i];
        }

        error = 0;
        for (int i = 0; i < n; i++) {
            error = error + (X[i] - Xviejo[i]) * (X[i] - Xviejo[i]);
        }
        error = sqrt(error);

        double normaX = 0;
        for (int i = 0; i < n; i++) {
            normaX = normaX + X[i] * X[i];
        }
        normaX = sqrt(normaX);
        double error_porcentual = (normaX == 0) ? 0.0 : (error / normaX) * 100;

        // Unica diferencia real con gaussSeidel.cpp: solo se declara
        // divergencia si el error explota de verdad, no ante cualquier
        // repunte puntual (ver comentario del encabezado).
        if (!isfinite(error) || error > 1e10) {
            cout << "El metodo diverge de verdad (error no acotado)." << endl;
            return 1;
        }

        if (iter % 50 == 0 || error <= tolerancia) { // no spamear miles de lineas
            cout << "Iteracion " << iter << ": error = " << error
                 << ", error porcentual = " << error_porcentual << "%" << endl;
        }

    } while (error > tolerancia && iter < max_iter);

    if (error > tolerancia) {
        cout << "Ojo: no se alcanzo la tolerancia pedida en " << max_iter
             << " iteraciones (convergencia demasiado lenta). Resultado parcial:" << endl;
    }

    double normaXFinal = 0;
    for (int i = 0; i < n; i++) {
        normaXFinal = normaXFinal + X[i] * X[i];
    }
    normaXFinal = sqrt(normaXFinal);
    double error_porcentual_final = (normaXFinal == 0) ? 0.0 : (error / normaXFinal) * 100;

    cout << "Solucion:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << X[i] << endl;
    }
    cout << "Error final: " << error << ", error porcentual: " << error_porcentual_final
         << "%, iteraciones: " << iter << endl;

    return 0;
}
