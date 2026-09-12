#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

// Complemento de regresionFuncionesBase.cpp: para cuando los coeficientes
// a0, a1 ya se calcularon aparte (a mano, planteando el sistema, y
// resolviendolo con eliminacionGausseana.cpp u otro Gauss). Este programa NO
// arma ni resuelve el sistema de minimos cuadrados -- solo evalua la funcion
// ya ajustada F(x) = a0*phi0(x) + a1*phi1(x) contra los puntos y calcula
// St, Sr y el coeficiente de correlacion r.

int main() {
    // ---- ACA se eligen las funciones base del problema (igual que en raices.cpp) ----
    auto phi0 = [](double x) { return sin(x); };
    auto phi1 = [](double x) { return cos(x); };

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
    cout << "Se leyeron " << npuntos << " puntos." << endl;

    double a0, a1;
    cout << "Ingrese los coeficientes a0 y a1 (ya calculados a mano/con Gauss): ";
    cin >> a0 >> a1;

    double sumaY = 0;
    for (int i = 0; i < npuntos; i++) sumaY += y[i];
    double yprom = sumaY / npuntos;

    double sumaT = 0, sumaR = 0;
    for (int i = 0; i < npuntos; i++) {
        double fi = a0 * phi0(x[i]) + a1 * phi1(x[i]);
        sumaT += pow(y[i] - yprom, 2);
        sumaR += pow(y[i] - fi, 2);
    }
    double r = sqrt((sumaT - sumaR) / sumaT);

    cout << "\nSt = " << sumaT << ", Sr = " << sumaR << endl;
    cout << "Coeficiente de correlacion r = " << r << endl;

    return 0;
}
