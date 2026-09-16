#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

int main() {

    // ---- Definicion de A y b (desde archivo de texto) ----
    // Mismo formato que eliminacionGausseana.cpp: matriz aumentada [A|b],
    // n filas, n+1 columnas, separadas por espacios.
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

    // ---- 1er aspecto: aii != 0 -> verificar si A es diagonalmente dominante ----
    // Condicion suficiente (no necesaria) para garantizar convergencia.
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

    // ---- 2do aspecto: inicializar el vector X ----
    vector<double> Xviejo(n, 0.0);
    vector<double> Xnuevo(n, 0.0);

    // ---- 3er aspecto: criterio de corte (tolerancia) ----
    double tolerancia;
    cout << "Ingrese la tolerancia: ";
    cin >> tolerancia;

    double error_viejo = 1000;
    double error;
    int iter = 0;

    cout << fixed << setprecision(10); // para que cada iteracion muestre las cifras completas, no 6 sig. figs redondeadas

    do {
        iter++;

        for (int i = 0; i < n; i++) {
            double suma = 0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    suma = suma + A[i][j] * Xviejo[j];
                }
            }
            Xnuevo[i] = (b[i] - suma) / A[i][i];
        }

        error = 0;
        for (int i = 0; i < n; i++) {
            error = error + (Xnuevo[i] - Xviejo[i]) * (Xnuevo[i] - Xviejo[i]);
        }
        error = sqrt(error);

        // Error porcentual: E = e / ||X(k+1)|| * 100
        double normaX = 0;
        for (int i = 0; i < n; i++) {
            normaX = normaX + Xnuevo[i] * Xnuevo[i];
        }
        normaX = sqrt(normaX);
        double error_porcentual = (normaX == 0) ? 0.0 : (error / normaX) * 100;

        if (error > error_viejo) {
            cout << "El metodo no converge." << endl;
            return 1;
        }
        error_viejo = error;

        for (int i = 0; i < n; i++) {
            Xviejo[i] = Xnuevo[i];
        }

        cout << "Iteracion " << iter << ": ";
        for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << Xnuevo[i] << ", ";
        cout << "error = " << error << ", error porcentual = " << error_porcentual << "%" << endl;

    } while (error > tolerancia);

    // ---- Mostrar resultados ----
    double normaXFinal = 0;
    for (int i = 0; i < n; i++) {
        normaXFinal = normaXFinal + Xnuevo[i] * Xnuevo[i];
    }
    normaXFinal = sqrt(normaXFinal);
    double error_porcentual_final = (normaXFinal == 0) ? 0.0 : (error / normaXFinal) * 100;

    cout << "Solucion:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << Xnuevo[i] << endl;
    }
    cout << "Error final: " << error << ", error porcentual: " << error_porcentual_final
         << "%, iteraciones: " << iter << endl;

    return 0;
}
