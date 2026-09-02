#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

void gauss();
void gaussJordan();

// ---- Ingresar A y b (desde archivo de texto) ----
// El archivo debe contener la matriz aumentada [A|b]: n filas, n+1 columnas,
// separadas por espacios, por ejemplo (n=3):
//   1 2 3 4
//   4 5 6 7
//   8 9 1 0
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

    // n se deduce de la cantidad de numeros leidos: total = n*(n+1)
    int total = datos.size();
    n = static_cast<int>((-1 + sqrt(1.0 + 4.0 * total)) / 2);

    if (n <= 0 || n * (n + 1) != total) {
        cerr << "El archivo no tiene el formato esperado (n filas x n+1 columnas)." << endl;
        return false;
    }

    A.assign(n, vector<double>(n));
    b.assign(n, 0.0);
    int idx = 0;
    for (int fila = 0; fila < n; fila++) {
        for (int col = 0; col < n; col++) {
            A[fila][col] = datos[idx++];
        }
        b[fila] = datos[idx++];
    }

    cout << "Matriz A y vector b leidos (n = " << n << "):" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << "| " << b[i] << endl;
    }
    return true;
}

// ---- Triangulacion (comun a Gauss y Gauss-Jordan) ----
// Nota: los indices del pseudocodigo son 1-based (i=1,...,n-1 ; a11 = primer
// elemento). En C++ los arrays son 0-based, asi que todo se corre un lugar:
// i arranca en 0 (no en 1) y det arranca en A[0][0] (no en A[1][1]).
// signoDet: arranca en 1 (afuera) y se multiplica por -1 cada vez que se
// intercambian dos filas de verdad. Cada swap de filas cambia el signo del
// determinante, asi que hay que arrastrarlo para no calcularlo mal.
void triangular(vector<vector<double>> &A, vector<double> &b, int n, int &signoDet) {
    for (int i = 0; i < n - 1; i++) {
        int p = i;
        if (fabs(A[i][i]) < pow(10, -4)) { // Pivoteo: A[i][i] muy chica
            cout << "Pivoteo necesario: |A[" << i+1 << "][" << i+1 << "]| = "
                 << fabs(A[i][i]) << " es muy chico (< 1e-4)." << endl;
            double m = fabs(A[i][i]);
            for (int l = i + 1; l < n; l++) {
                if (fabs(A[l][i]) > m) {
                    m = fabs(A[l][i]);
                    p = l;
                }
            }
            if (p != i) {
                cout << "  -> se intercambia la fila " << i+1 << " con la fila " << p+1 << endl;
                signoDet = -signoDet;
            } else {
                cout << "  -> no se encontro una fila mejor, se sigue con el mismo pivote" << endl;
            }
            for (int l = 0; l < n; l++) {
                double aux = A[p][l];
                A[p][l] = A[i][l];
                A[i][l] = aux;
            }
            double aux = b[p];
            b[p] = b[i];
            b[i] = aux;
        }

        for (int j = i + 1; j < n; j++) {
            double factor = (-A[j][i]) / A[i][i];

            for (int k = i + 1; k < n; k++) {
                A[j][k] = A[i][k] * factor + A[j][k];
            }
            b[j] = b[i] * factor + b[j];
            A[j][i] = 0.0; // queda explicitamente triangular
        }
    }
}

// Norma de Frobenius de A: ||A|| = sqrt( sum_i sum_j A[i][j]^2 ).
// Hay que calcularla ANTES de triangular (triangular pisa A con ceros).
// Sirve para detectar matrices mal condicionadas: si det(A) es chico
// comparado con ||A||, la matriz es "cuasi-singular".
double normaFrobenius(const vector<vector<double>> &A, int n) {
    double suma = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            suma += A[i][j] * A[i][j];
        }
    }
    return sqrt(suma);
}

int main() {
    int opcion;
    do {
        cout << "Menu de opciones:" << endl;
        cout << "1. Gauss (triangulacion + retrosustitucion)" << endl;
        cout << "2. Gauss-Jordan (despeje directo)" << endl;
        cout << "3. Salir" << endl;
        cout << "Seleccione una opcion: ";
        cin >> opcion;

        switch (opcion) {
            case 1:
                gauss();
                break;
            case 2:
                gaussJordan();
                break;
            case 3:
                cout << "Saliendo del programa." << endl;
                break;
            default:
                cout << "Opcion no valida, intente de nuevo." << endl;
        }
    } while (opcion != 3);
    return 0;
}

void gauss() {
    vector<vector<double>> A;
    vector<double> b;
    int n;
    if (!leerSistema(A, b, n)) return;

    double normaA = normaFrobenius(A, n); // antes de triangular, con la matriz original

    int signoDet = 1;
    triangular(A, b, n, signoDet);

    // ---- Calcular det(A) ----
    // signoDet compensa los intercambios de filas hechos durante el pivoteo
    // (cada swap de filas cambia el signo del determinante real).
    double det = signoDet * A[0][0];
    for (int i = 1; i < n; i++) {
        det = det * A[i][i];
    }
    cout << "El valor de la determinante es: " << det << endl;
    cout << "Norma de A (Frobenius): " << normaA << endl;
    if (det == 0) {
        cout << "matriz singular" << endl;
        return;
    }

    // ---- Retrosustitucion ----
    vector<double> x(n);
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];

    for (int i = n - 2; i >= 0; i--) {
        double suma = b[i];
        for (int j = i + 1; j < n; j++) {
            suma = suma - A[i][j] * x[j];
        }
        x[i] = suma / A[i][i];
    }

    cout << "Solucion (Gauss):" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}

void gaussJordan() {
    vector<vector<double>> A;
    vector<double> b;
    int n;
    if (!leerSistema(A, b, n)) return;

    double normaA = normaFrobenius(A, n); // antes de triangular, con la matriz original

    // Misma triangulacion que Gauss (hacia abajo del pivote).
    int signoDet = 1;
    triangular(A, b, n, signoDet);

    // El determinante se calcula ACA, con la diagonal recien triangulada.
    // Ojo: si se calculara despues de normalizar filas (mas abajo), la
    // diagonal ya seria toda 1 y el determinante daria siempre 1 (mal).
    // signoDet compensa los intercambios de filas hechos durante el pivoteo.
    double det = signoDet * A[0][0];
    for (int i = 1; i < n; i++) {
        det = det * A[i][i];
    }
    cout << "El valor de la determinante es: " << det << endl;
    cout << "Norma de A (Frobenius): " << normaA << endl;
    if (det == 0) {
        cout << "matriz singular" << endl;
        return;
    }

    // ---- Eliminacion hacia arriba: lo que suma Gauss-Jordan sobre Gauss ----
    // Se procesa cada columna pivote de atras para adelante (i=n-1,...,1).
    // Para cuando una fila se usa como pivote, ya quedo con un unico valor
    // no nulo en su propia columna (las columnas de mas a la derecha ya se
    // eliminaron en un paso anterior de este mismo loop), asi que alcanza con
    // actualizar b[j] y poner en 0 la entrada de arriba del pivote.
    for (int i = n - 1; i > 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            double factor = (-A[j][i]) / A[i][i];
            b[j] = b[i] * factor + b[j];
            A[j][i] = 0.0;
        }
    }

    // ---- Normalizar: dividir cada fila por su pivote -> queda [I | x] ----
    vector<double> x(n);
    for (int i = 0; i < n; i++) {
        x[i] = b[i] / A[i][i];
    }

    cout << "Solucion (Gauss-Jordan):" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}
