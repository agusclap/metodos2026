#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

int main() {

    // ---- Definicion de A y b (desde archivo de texto) ----
    // Mismo formato que jacobi.cpp / eliminacionGausseana.cpp: matriz
    // aumentada [A|b], n filas, n+1 columnas, separadas por espacios.
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
    // Un solo array: se actualiza en el lugar (esa es la diferencia con Jacobi).
    vector<double> X(n, 0.0);
    vector<double> Xviejo(n, 0.0);

    // ---- 3er aspecto: criterio de corte (tolerancia) ----
    double tolerancia;
    cout << "Ingrese la tolerancia: ";
    cin >> tolerancia;

    // Coeficiente de relajacion (ver "Gauss-Seidel con Relajacion" del apunte).
    // w=1 -> Gauss-Seidel comun. w<1 -> sub-relajacion (ayuda a converger a
    // un sistema que no converge). w>1 -> sobre-relajacion (acelera uno que
    // ya converge). Si el ejercicio no menciona relajacion, se ingresa w=1.
    double w;
    cout << "Ingrese el coeficiente de relajacion w (1 = Gauss-Seidel sin relajacion): ";
    cin >> w;

    double error_viejo = 1000;
    double error;
    int iter = 0;

    do {
        iter++;
        Xviejo = X; // foto del vector antes de esta pasada, para el error

        for (int i = 0; i < n; i++) {
            double suma = b[i];
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    // j<i ya tiene el valor nuevo de esta pasada; j>i todavia
                    // tiene el valor viejo (no se toco todavia esta vuelta).
                    suma = suma - A[i][j] * X[j];
                }
            }
            double xGS = suma / A[i][i];
            X[i] = w * xGS + (1 - w) * Xviejo[i]; // relajacion
        }

        error = 0;
        for (int i = 0; i < n; i++) {
            error = error + (X[i] - Xviejo[i]) * (X[i] - Xviejo[i]);
        }
        error = sqrt(error);

        // Error porcentual: E = e / ||X(k+1)|| * 100
        double normaX = 0;
        for (int i = 0; i < n; i++) {
            normaX = normaX + X[i] * X[i];
        }
        normaX = sqrt(normaX);
        double error_porcentual = (normaX == 0) ? 0.0 : (error / normaX) * 100;

        if (error > error_viejo) {
            cout << "El metodo no converge." << endl;
            return 1;
        }
        error_viejo = error;

        cout << "Iteracion " << iter << ": error = " << error
             << ", error porcentual = " << error_porcentual << "%" << endl;

    } while (error > tolerancia);

    // ---- Mostrar resultados ----
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
