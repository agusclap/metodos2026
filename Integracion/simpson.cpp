#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

void simpsonFuncion();
void simpsonTabla();
vector<double> armarZSpline(const vector<double> &x, const vector<double> &y);
double evaluarSpline(const vector<double> &x, const vector<double> &z, double xhat);

int main() {
    int opcion;
    do {
        cout << "Menu de opciones:" << endl;
        cout << "1. Simpson 1/3 compuesto (funcion definida en el codigo)" << endl;
        cout << "2. Simpson 1/3 compuesto (tabla de datos desde archivo)" << endl;
        cout << "3. Salir" << endl;
        cout << "Seleccione una opcion: ";
        cin >> opcion;

        switch (opcion) {
            case 1:
                simpsonFuncion();
                break;
            case 2:
                simpsonTabla();
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

// ============================================================================
// Opcion 1: pseudocodigo "idem a trapecio, solo se modifica el lazo que
// calcula suma" (asi lo dio la profesora). n=2 -> Simpson 1/3 SIMPLE (el for
// de abajo no ejecuta ninguna vuelta, i=1,...,0 es vacio). n>2 (PAR) ->
// Simpson 1/3 COMPUESTO.
// ============================================================================
void simpsonFuncion() {
    auto f = [](double x) { return sin(2*x) * exp(-x); }; // <-- ACA se elige el problema

    double a, b;
    cout << "Ingrese a y b: ";
    cin >> a >> b;

    int n;
    cout << "Ingrese el numero de subintervalos (n, debe ser PAR; n=2 -> Simpson simple): ";
    cin >> n;
    if (n < 2 || n % 2 != 0) {
        cerr << "n tiene que ser PAR y >= 2." << endl;
        return;
    }

    double h = (b - a) / n;
    double suma = f(a) + f(b);
    for (int i = 1; i <= n/2 - 1; i++) {
        double x = a + 2*i*h;
        suma = suma + 2*f(x) + 4*f(x - h);
    }
    suma = suma + 4*f(b - h);
    double I = (h/3) * suma;

    cout << "La integral aproximada es: I = " << I << endl;
}

// ============================================================================
// Opcion 2: mismo pseudocodigo, pero con tabla de datos. Igual que en
// trapecio.cpp: si los puntos no estan equiespaciados, se genera una tabla
// equiespaciada con el spline cubico natural antes de aplicar Simpson.
// ============================================================================
void simpsonTabla() {
    string nombreArchivo;
    cout << "Ingrese el nombre del archivo con los puntos (x y por linea): ";
    cin >> nombreArchivo;

    ifstream archivo(nombreArchivo);
    if (!archivo.is_open()) {
        cerr << "No se pudo abrir el archivo." << endl;
        return;
    }
    vector<double> x, y;
    double xi, yi;
    while (archivo >> xi >> yi) {
        x.push_back(xi);
        y.push_back(yi);
    }
    archivo.close();

    int npuntos = (int)x.size();
    int n = npuntos - 1; // n subintervalos
    if (n < 2 || n % 2 != 0) {
        cerr << "Simpson 1/3 necesita un numero PAR de subintervalos "
             << "(o sea, una cantidad IMPAR de puntos, >=3)." << endl;
        return;
    }

    double h = (x[n] - x[0]) / n;

    bool equiespaciados = true;
    for (int i = 0; i < n; i++) {
        if (fabs((x[i+1] - x[i]) - h) > 1e-9) {
            equiespaciados = false;
            break;
        }
    }

    vector<double> xp(n + 1), yp(n + 1);
    if (!equiespaciados) {
        cout << "Los puntos no estan equiespaciados -> generando tabla equiespaciada con spline cubico natural." << endl;
        vector<double> z = armarZSpline(x, y);
        for (int i = 0; i <= n; i++) {
            xp[i] = x[0] + i * h;
            yp[i] = evaluarSpline(x, z, xp[i]);
        }
    } else {
        xp = x;
        yp = y;
    }

    // Mismo lazo que en simpsonFuncion(), pero indexando el array yp en vez
    // de evaluar F(x): yp[2i] hace de F(x), yp[2i-1] hace de F(x-h).
    double suma = yp[0] + yp[n];
    for (int i = 1; i <= n/2 - 1; i++) {
        suma = suma + 2*yp[2*i] + 4*yp[2*i - 1];
    }
    suma = suma + 4*yp[n - 1];
    double I = (h/3) * suma;

    cout << "La integral aproximada es: I = " << I << endl;
}

// ============================================================================
// Spline cubico natural (identico a Interpolacion/splineCubico.cpp y a
// Integracion/trapecio.cpp): arma y resuelve el sistema 4n x 4n con Gauss +
// pivoteo, y devuelve z (los 4n coeficientes a_k,b_k,c_k,d_k de cada tramo).
// ============================================================================
vector<double> armarZSpline(const vector<double> &x, const vector<double> &y) {
    int n = (int)x.size() - 1;
    int m = 4 * n;
    vector<vector<double>> A(m, vector<double>(m, 0.0));
    vector<double> b(m, 0.0);

    for (int k = 0; k < n; k++) {
        for (int j = 0; j <= 3; j++) {
            A[2*k][4*k + j]     = pow(x[k],   3 - j);
            A[2*k + 1][4*k + j] = pow(x[k+1], 3 - j);
        }
        b[2*k]     = y[k];
        b[2*k + 1] = y[k+1];
    }
    for (int k = 0; k <= n - 2; k++) {
        int i = 2*n + k;
        for (int j = 0; j <= 2; j++) {
            A[i][4*k + j]     = (3 - j) * pow(x[k+1], 2 - j);
            A[i][4*(k+1) + j] = -(3 - j) * pow(x[k+1], 2 - j);
        }
    }
    for (int k = 0; k <= n - 2; k++) {
        int i = 3*n - 1 + k;
        A[i][4*k]         = 3 * x[k+1];
        A[i][4*k + 1]     = 1.0;
        A[i][4*(k+1)]     = -3 * x[k+1];
        A[i][4*(k+1) + 1] = -1.0;
    }
    A[4*n - 2][0] = 3 * x[0];
    A[4*n - 2][1] = 1.0;
    A[4*n - 1][4*n - 4] = 3 * x[n];
    A[4*n - 1][4*n - 3] = 1.0;

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

double evaluarSpline(const vector<double> &x, const vector<double> &z, double xhat) {
    int n = (int)x.size() - 1;
    for (int i = 0; i < n; i++) {
        // OJO: +1e-9 de tolerancia -- x[i+1] resampleado como x[0]+i*h puede
        // dar un pelito mas grande que el x[n] real de la tabla por redondeo
        // de punto flotante, y una comparacion <= exacta lo dejaria afuera.
        if (xhat >= x[i] - 1e-9 && xhat <= x[i+1] + 1e-9) {
            return z[4*i]*pow(xhat,3) + z[4*i+1]*pow(xhat,2) + z[4*i+2]*xhat + z[4*i+3];
        }
    }
    return NAN;
}
