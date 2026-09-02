#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

// Interpolacion de Lagrange: evalua Pn(x) directo en un punto, sin armar
// el polinomio simbolico (esa es la desventaja del metodo).
int main() {
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

    int n = (int)x.size() - 1; // n+1 puntos -> polinomio de grado n
    if (n < 0) {
        cerr << "El archivo no tiene puntos validos." << endl;
        return 1;
    }
    cout << "Se leyeron " << (n + 1) << " puntos (polinomio de grado " << n << ")." << endl;

    double xhat;
    cout << "Ingrese el valor x a interpolar: ";
    cin >> xhat;

    // ---- Calcular Pn(x_hat) ----
    double suma = 0;
    for (int k = 0; k <= n; k++) {
        double prod = 1;
        for (int i = 0; i <= n; i++) {
            if (i != k) {
                prod = prod * (xhat - x[i]) / (x[k] - x[i]);
            }
        }
        suma = suma + y[k] * prod;
    }

    cout << "P(" << xhat << ") = " << suma << endl;

    return 0;
}
