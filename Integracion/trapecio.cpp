#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

void trapecioFuncion();
void trapecioTabla();
vector<double> armarZSpline(const vector<double> &x, const vector<double> &y);
double evaluarSpline(const vector<double> &x, const vector<double> &z, double xhat);

int main() {
    int opcion;
    do {
        cout << "Menu de opciones:" << endl;
        cout << "1. Trapecio compuesto (funcion definida en el codigo)" << endl;
        cout << "2. Trapecio compuesto (tabla de datos desde archivo)" << endl;
        cout << "3. Salir" << endl;
        cout << "Seleccione una opcion: ";
        cin >> opcion;

        switch (opcion) {
            case 1:
                trapecioFuncion();
                break;
            case 2:
                trapecioTabla();
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
// Opcion 1: pseudocodigo "con funcion". n=1 -> trapecio SIMPLE (el for de abajo
// no ejecuta ninguna vuelta, i=1,...,0 es vacio). n>1 -> trapecio COMPUESTO.
// Es la misma formula para los dos casos, no hace falta codigo separado.
// ============================================================================
void trapecioFuncion() {
    auto f = [](double x) { return sin(2*x) * exp(-x); }; // <-- ACA se elige el problema

    double a, b;
    cout << "Ingrese a y b: ";
    cin >> a >> b;

    int n;
    cout << "Ingrese el numero de subintervalos (n=1 -> trapecio simple, n>1 -> compuesto): ";
    cin >> n;
    if (n < 1) {
        cerr << "n tiene que ser >= 1." << endl;
        return;
    }

    double h = (b - a) / n;
    double suma = f(a) + f(b);
    for (int i = 1; i <= n - 1; i++) {
        double x = a + i * h;
        suma = suma + 2 * f(x);
    }
    double I = (h / 2) * suma;

    cout << "La integral aproximada es: I = " << I << endl;
}

// ============================================================================
// Opcion 2: pseudocodigo "con tabla de datos". Si los puntos del archivo NO
// estan equiespaciados, primero se genera una tabla equiespaciada evaluando
// un spline cubico natural (el de splineCubico.cpp) en n+1 puntos parejos, y
// se aplica el trapecio compuesto sobre esa tabla nueva.
// ============================================================================
void trapecioTabla() {
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
    if (n < 1) {
        cerr << "Hacen falta al menos 2 puntos." << endl;
        return;
    }

    double h = (x[n] - x[0]) / n;

    // Verificar si estan equiespaciados. OJO: comparar con "!=" exacto (como
    // dice el pseudocodigo al pie de la letra) es fragil con numeros de punto
    // flotante leidos de un archivo -- se usa una tolerancia chica en su lugar.
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

    double suma = yp[n] + yp[0];
    for (int i = 1; i <= n - 1; i++) {
        suma = suma + 2 * yp[i];
    }
    double I = (h / 2) * suma;

    cout << "La integral aproximada es: I = " << I << endl;
}

// ============================================================================
// Spline cubico natural (identico a Interpolacion/splineCubico.cpp): arma y
// resuelve el sistema 4n x 4n con Gauss + pivoteo, y devuelve z (los 4n
// coeficientes a_k,b_k,c_k,d_k de cada tramo).
// ============================================================================
vector<double> armarZSpline(const vector<double> &x, const vector<double> &y) {
    int n = (int)x.size() - 1;
    int m = 4 * n;
    vector<vector<double>> A(m, vector<double>(m, 0.0));
    vector<double> b(m, 0.0);

    // Bloque 1: F_k(x_k)=y_k, F_k(x_k+1)=y_k+1
    for (int k = 0; k < n; k++) {
        for (int j = 0; j <= 3; j++) {
            A[2*k][4*k + j]     = pow(x[k],   3 - j);
            A[2*k + 1][4*k + j] = pow(x[k+1], 3 - j);
        }
        b[2*k]     = y[k];
        b[2*k + 1] = y[k+1];
    }
    // Bloque 2: continuidad de F'
    for (int k = 0; k <= n - 2; k++) {
        int i = 2*n + k;
        for (int j = 0; j <= 2; j++) {
            A[i][4*k + j]     = (3 - j) * pow(x[k+1], 2 - j);
            A[i][4*(k+1) + j] = -(3 - j) * pow(x[k+1], 2 - j);
        }
    }
    // Bloque 3: continuidad de F''
    for (int k = 0; k <= n - 2; k++) {
        int i = 3*n - 1 + k;
        A[i][4*k]         = 3 * x[k+1];
        A[i][4*k + 1]     = 1.0;
        A[i][4*(k+1)]     = -3 * x[k+1];
        A[i][4*(k+1) + 1] = -1.0;
    }
    // Bloque 4: spline natural, F''=0 en los extremos
    A[4*n - 2][0] = 3 * x[0];
    A[4*n - 2][1] = 1.0;
    A[4*n - 1][4*n - 4] = 3 * x[n];
    A[4*n - 1][4*n - 3] = 1.0;

    // Eliminacion Gaussiana con pivoteo parcial
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

// Evalua el spline ya resuelto en un punto xhat. A diferencia de
// splineCubico.cpp (que compara con < y > estrictos), aca la comparacion
// incluye los extremos (<=, >=) porque necesitamos evaluar EXACTO en x[0] y
// x[n] (los bordes de la tabla nueva equiespaciada).
double evaluarSpline(const vector<double> &x, const vector<double> &z, double xhat) {
    int n = (int)x.size() - 1;
    for (int i = 0; i < n; i++) {
        if (xhat >= x[i] && xhat <= x[i+1]) {
            return z[4*i]*pow(xhat,3) + z[4*i+1]*pow(xhat,2) + z[4*i+2]*xhat + z[4*i+3];
        }
    }
    return NAN; // xhat fuera del rango de la tabla
}
