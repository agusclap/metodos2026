#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
using namespace std;

// ============================================================================
// Igual que lagrange.cpp, pero en vez de evaluar Pn(x) en UN punto, lo evalua
// en una grilla fina y lo vuelca a un CSV para graficar.
//
// Se usa cuando los puntos son datos experimentales (vienen de un archivo) y
// NO hay una f(x) exacta contra la cual comparar -- por eso no hay columna de
// error, a diferencia de lagrangeFuncion.cpp.
//
// El rango de la grilla se pide aparte de los datos: poniendo un x_final mayor
// que el ultimo dato se ve la extrapolacion.
// ============================================================================

double lagrange(const vector<double> &x, const vector<double> &y, double xhat) {
    int n = (int)x.size() - 1;
    double suma = 0.0;
    for (int k = 0; k <= n; k++) {
        double prod = 1.0;
        for (int i = 0; i <= n; i++) {
            if (i != k) prod *= (xhat - x[i]) / (x[k] - x[i]);
        }
        suma += y[k] * prod;
    }
    return suma;
}

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

    int n = (int)x.size() - 1;
    if (n < 1) {
        cerr << "Hacen falta al menos 2 puntos." << endl;
        return 1;
    }
    cout << "Se leyeron " << (n + 1) << " puntos (polinomio de grado " << n << ")." << endl;
    cout << "Los datos van de x = " << x[0] << " a x = " << x[n] << endl;

    double xIni, xFin;
    cout << "Ingrese el rango a graficar [x_inicial, x_final]: ";
    cin >> xIni >> xFin;
    if (xFin <= xIni) {
        cerr << "El rango no es valido." << endl;
        return 1;
    }

    string salida;
    cout << "Nombre del archivo CSV a generar: ";
    cin >> salida;

    ofstream out(salida);
    if (!out.is_open()) {
        cerr << "No se pudo crear el archivo." << endl;
        return 1;
    }
    out << setprecision(12);
    out << "x,P_n\n";

    const int npuntos = 400;
    for (int i = 0; i < npuntos; i++) {
        double xh = xIni + i * (xFin - xIni) / (npuntos - 1);
        out << xh << "," << lagrange(x, y, xh) << "\n";
    }
    out.close();

    // Los datos originales aparte, para marcarlos como puntos sobre la curva
    string salidaDatos = "datos_" + salida;
    ofstream outd(salidaDatos);
    outd << setprecision(12);
    outd << "x,y\n";
    for (int i = 0; i <= n; i++) outd << x[i] << "," << y[i] << "\n";
    outd.close();

    cout << "Listo: " << salida << " (curva del polinomio) y "
         << salidaDatos << " (puntos medidos)." << endl;
    return 0;
}
