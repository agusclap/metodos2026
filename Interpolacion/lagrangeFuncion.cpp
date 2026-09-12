#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

// ============================================================================
// Interpolacion de Lagrange aplicada a una funcion conocida.
//
// Para un grado n: genera los n+1 nodos equiespaciados en [x0,xN], evalua la
// funcion exacta en ellos (eso da los puntos a interpolar), arma el polinomio
// de Lagrange y lo evalua en una grilla fina que incluye tanto el rango de
// interpolacion como un margen de extrapolacion afuera de [x0,xN]. Escribe
// todo a un CSV (x, f_exacta, P_n, error) para graficar.
//
// El algoritmo de Lagrange es el mismo de Interpolacion/lagrange.cpp (evaluar
// Pn(x) directo, sin armar el polinomio simbolico).
//
// Lo UNICO que cambia de problema a problema es la funcion de abajo. El resto
// es generico: el dominio, el margen y los n se ingresan por teclado.
// ============================================================================

// ---- ACA se elige el problema: la funcion a interpolar ----
double f(double x) {
    return log(x*x + 1) - sin(x);                  // Problema 2, item 1: f(x) = e^(-x^2)
    // return 4*x*x*x - 3*x*x + 2;       // item 2: f(x) = 4x^3 - 3x^2 + 2
    // return pow(x, x);                 // item 3: f(x) = x^x  (solo definida para x>0)
}

// Evalua el polinomio de Lagrange de grado n (n+1 nodos) en el punto xhat.
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
    string etiqueta;
    cout << "Etiqueta para los archivos de salida (ej: gaussiana): ";
    cin >> etiqueta;

    double x0, xN;
    cout << "Ingrese el rango de interpolacion [x0, xN]: ";
    cin >> x0 >> xN;
    if (xN <= x0) {
        cerr << "El rango no es valido: xN tiene que ser mayor que x0." << endl;
        return 1;
    }

    double margen;
    cout << "Ingrese el margen de extrapolacion (cuanto evaluar antes de x0 y despues de xN): ";
    cin >> margen;
    double xIni = x0 - margen;
    double xFin = xN + margen;

    int cantN;
    cout << "Cuantos valores de n (grado del polinomio) quiere probar: ";
    cin >> cantN;

    const int npuntosGrilla = 400;

    for (int c = 0; c < cantN; c++) {
        int n;
        cout << "Ingrese n (" << (c + 1) << "/" << cantN << "): ";
        cin >> n;
        if (n < 1) {
            cerr << "  n tiene que ser >= 1. Se saltea." << endl;
            continue;
        }

        // ---- Generar n+1 nodos equiespaciados y evaluar f en ellos ----
        vector<double> xNodos(n + 1), yNodos(n + 1);
        bool nodoInvalido = false;
        for (int i = 0; i <= n; i++) {
            xNodos[i] = x0 + i * (xN - x0) / n;
            yNodos[i] = f(xNodos[i]);
            if (!isfinite(yNodos[i])) nodoInvalido = true;
        }
        if (nodoInvalido) {
            cerr << "  La funcion no esta definida en alguno de los nodos de [x0,xN]."
                 << " No se puede interpolar con n = " << n << "." << endl;
            continue;
        }

        // ---- Evaluar P_n(x) y f(x) en una grilla fina (incluye extrapolacion) ----
        string nombreArchivo = etiqueta + "_n" + to_string(n) + ".csv";
        ofstream out(nombreArchivo);
        if (!out.is_open()) {
            cerr << "  No se pudo crear " << nombreArchivo << endl;
            continue;
        }
        out << setprecision(12); // 6 cifras (el default) redondea x y el error al graficar
        out << "x,f_exacta,P_n,error\n";

        double errorMaxInterp = 0.0, errorMaxExtrap = 0.0;
        int sinDefinir = 0;

        for (int i = 0; i < npuntosGrilla; i++) {
            double x = xIni + i * (xFin - xIni) / (npuntosGrilla - 1);
            double fx = f(x);                       // puede no estar definida (log, sqrt, x^x, 1/x...)
            double px = lagrange(xNodos, yNodos, x); // el polinomio SI esta definido en todo x
            double err = fabs(fx - px);

            out << x << "," << (isfinite(fx) ? to_string(fx) : "") << "," << px << ","
                << (isfinite(err) ? to_string(err) : "") << "\n";

            // Los puntos donde la funcion no existe no entran en la estadistica de error
            if (!isfinite(err)) {
                sinDefinir++;
                continue;
            }
            if (x >= x0 && x <= xN) errorMaxInterp = max(errorMaxInterp, err);
            else                    errorMaxExtrap = max(errorMaxExtrap, err);
        }
        out.close();

        cout << "  -> " << nombreArchivo << " generado. Error maximo dentro de [x0,xN] = "
             << errorMaxInterp << ", error maximo en extrapolacion = " << errorMaxExtrap << endl;
        if (sinDefinir > 0) {
            cout << "     (" << sinDefinir << " de " << npuntosGrilla
                 << " puntos de la grilla quedaron fuera del dominio de la funcion"
                 << " y no se contaron en el error)" << endl;
        }
    }

    return 0;
}
