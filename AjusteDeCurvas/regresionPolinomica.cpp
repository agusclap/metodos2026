#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

// Regresion polinomica por cuadrados minimos: F(x) = a0 + a1*x + ... + ap*x^p
// Con p=1 esto ES la regresion lineal simple (F(x)=a0+a1*x, "b" y "a" del
// apunte en el orden [a0=ordenada, a1=pendiente]).
//
// Ecuaciones normales: A[l][m] = sum(x_i^(l+m)), b[l] = sum(y_i * x_i^l),
// para l,m = 0,...,p. Se resuelve el sistema (p+1)x(p+1) con Gauss+pivoteo.

int main() {
    // ---- Leer los datos (x,y) desde archivo ----
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
    if (npuntos < 2) {
        cerr << "Hacen falta al menos 2 puntos." << endl;
        return 1;
    }
    cout << "Se leyeron " << npuntos << " puntos." << endl;

    // ---- Ingresar p (grado del polinomio) ----
    int p;
    cout << "Ingrese el grado del polinomio a ajustar (p): ";
    cin >> p;

    // Hacen falta al menos p+1 puntos para determinar los p+1 coeficientes.
    if (npuntos < p + 1) {
        cerr << "No hay suficientes datos (se necesitan al menos " << p + 1
             << " puntos para un polinomio de grado " << p << ")." << endl;
        return 1;
    }

    // ---- Armar A y b (ecuaciones normales) ----
    int m = p + 1; // tamano del sistema
    vector<vector<double>> A(m, vector<double>(m));
    vector<double> b(m);

    for (int l = 0; l < m; l++) {
        double sumaAy = 0;
        for (int i = 0; i < npuntos; i++) {
            sumaAy += pow(x[i], l) * y[i];
        }
        b[l] = sumaAy;

        for (int col = 0; col < m; col++) {
            double sumaAx = 0;
            for (int i = 0; i < npuntos; i++) {
                sumaAx += pow(x[i], l + col);
            }
            A[l][col] = sumaAx;
        }
    }

    // ---- Resolver el sistema: Gauss con pivoteo parcial ----
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

    cout << "\nCoeficientes (F(x) = a0 + a1*x + ... + ap*x^p):" << endl;
    for (int i = 0; i < m; i++) cout << "a" << i << " = " << a[i] << endl;

    // ---- Coeficiente de correlacion ----
    double sumaY = 0;
    for (int i = 0; i < npuntos; i++) sumaY += y[i];
    double yprom = sumaY / npuntos;

    double sumaT = 0, sumaR = 0;
    for (int i = 0; i < npuntos; i++) {
        sumaT += pow(y[i] - yprom, 2);

        double fi = 0;
        for (int j = 0; j < m; j++) fi += a[j] * pow(x[i], j);
        sumaR += pow(y[i] - fi, 2);
    }
    double r = sqrt((sumaT - sumaR) / sumaT);

    // Error cuadratico medio: E = ( (1/N) * sum |F(xi) - yi|^2 )^(1/2).
    // La sumatoria de adentro es exactamente Sr, asi que E = sqrt(Sr/N).
    double ecm = sqrt(sumaR / npuntos);

    cout << "\nSt = " << sumaT << ", Sr = " << sumaR << endl;
    cout << "Coeficiente de correlacion r = " << r << endl;
    cout << "Error cuadratico medio Ecm = " << ecm << endl;

    return 0;
}
