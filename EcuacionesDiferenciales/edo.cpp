#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

void euler();
void heun();
void puntoMedio();

// Solucion exacta del problema (si se conoce -- no siempre es posible, pero
// cuando la EDO es separable, se puede resolver a mano como en el apunte y
// poner el resultado aca). Se usa solo para comparar/graficar; si no se
// conoce para el problema que estes resolviendo, dejarla devolviendo NAN.
double yExacta(double x) { return exp(-x * x); } // <-- ACA se elige el problema (o NAN si no se conoce)

int main() {
    int opcion;
    do {
        cout << "Menu de opciones:" << endl;
        cout << "1. Metodo de Euler" << endl;
        cout << "2. Metodo de Heun" << endl;
        cout << "3. Metodo del punto medio" << endl;
        cout << "4. Salir" << endl;
        cout << "Seleccione una opcion: ";
        cin >> opcion;

        switch (opcion) {
            case 1:
                euler();
                break;
            case 2:
                heun();
                break;
            case 3:
                puntoMedio();
                break;
            case 4:
                cout << "Saliendo del programa." << endl;
                break;
            default:
                cout << "Opcion no valida, intente de nuevo." << endl;
        }
    } while (opcion != 4);
    return 0;
}

// Pide los datos comunes a los 3 metodos (resolver dy/dx=F(x,y) en [x0,xF]
// con y(x0)=y0, dividido en n subintervalos de ancho h).
void leerDatosComunes(double &x0, double &y0, double &xF, int &n, double &h) {
    cout << "Ingrese x0 e y0: ";
    cin >> x0 >> y0;
    cout << "Ingrese xF: ";
    cin >> xF;
    cout << "Ingrese el numero de subintervalos (n): ";
    cin >> n;
    h = (xF - x0) / n;
}

// Imprime la tabla por consola (aproximada, exacta si se conoce -yExacta()-,
// y el error entre las dos) y la exporta a un CSV con las mismas columnas,
// listo para graficar la exacta y la aproximada juntas en un mismo grafico.
void mostrarYExportar(const vector<double> &x, const vector<double> &y, const string &metodo) {
    bool hayExacta = !isnan(yExacta(x[0]));

    cout << fixed << setprecision(6);
    for (size_t i = 0; i < x.size(); i++) {
        cout << "x[" << i << "] = " << x[i] << ", y_aprox[" << i << "] = " << y[i];
        if (hayExacta) {
            double exacta = yExacta(x[i]);
            cout << ", y_exacta = " << exacta << ", error = " << fabs(exacta - y[i]);
        }
        cout << endl;
    }
    cout << "Resultado final: y(" << x.back() << ") = " << y.back();
    if (hayExacta) cout << "  (exacto: " << yExacta(x.back()) << ")";
    cout << endl;
    cout.unsetf(ios::fixed);
    cout << setprecision(6);

    string nombreArchivo;
    cout << "Nombre del archivo CSV para graficar (dejar vacio y tocar enter para no exportar): ";
    cin.ignore();
    getline(cin, nombreArchivo);
    if (nombreArchivo.empty()) return;

    ofstream out(nombreArchivo);
    if (!out.is_open()) {
        cerr << "No se pudo crear " << nombreArchivo << endl;
        return;
    }
    out << setprecision(12);
    out << "x," << metodo << (hayExacta ? ",exacta,error\n" : "\n");
    for (size_t i = 0; i < x.size(); i++) {
        out << x[i] << "," << y[i];
        if (hayExacta) {
            double exacta = yExacta(x[i]);
            out << "," << exacta << "," << fabs(exacta - y[i]);
        }
        out << "\n";
    }
    out.close();
    cout << "Listo: " << nombreArchivo << " generado"
         << (hayExacta ? " (con la exacta y el error, listo para graficar las dos juntas)." : ".") << endl;
}

// ============================================================================
// Metodo de Euler: y[i+1] = y[i] + h*F(x[i],y[i])
// Error de truncamiento local ~ O(h^2) por paso.
// ============================================================================
void euler() {
    auto F = [](double x, double y) { return -2 * x * y; }; // <-- ACA se elige el problema

    double x0, y0, xF, h;
    int n;
    leerDatosComunes(x0, y0, xF, n, h);

    vector<double> x(n + 1), y(n + 1);
    x[0] = x0;
    y[0] = y0;
    for (int i = 0; i < n; i++) {
        x[i + 1] = x[i] + h;
        y[i + 1] = y[i] + h * F(x[i], y[i]);
    }

    mostrarYExportar(x, y, "euler");
}

// ============================================================================
// Metodo de Heun (predictor-corrector): primero se estima y[i+1] con Euler
// (yc), y despues se corrige promediando la pendiente en x[i] con la
// pendiente en x[i+1] (evaluada en esa estimacion yc).
// Error de truncamiento local ~ O(h^3) por paso -- un orden mas que Euler.
// ============================================================================
void heun() {
    auto F = [](double x, double y) { return -2 * x * y; }; // <-- ACA se elige el problema

    double x0, y0, xF, h;
    int n;
    leerDatosComunes(x0, y0, xF, n, h);

    vector<double> x(n + 1), y(n + 1);
    x[0] = x0;
    y[0] = y0;
    for (int i = 0; i < n; i++) {
        x[i + 1] = x[i] + h;
        double yc = y[i] + h * F(x[i], y[i]);                       // prediccion (Euler)
        y[i + 1] = y[i] + (h / 2.0) * (F(x[i], y[i]) + F(x[i + 1], yc)); // correccion
    }

    mostrarYExportar(x, y, "heun");
}

// ============================================================================
// Metodo del punto medio: la pendiente se evalua en el punto medio del
// intervalo [x[i],x[i+1]], no en x[i] como Euler. Como no se conoce y en ese
// punto medio, se estima con medio paso de Euler (h/2).
// Error de truncamiento local ~ O(h^3), igual que Heun.
// ============================================================================
void puntoMedio() {
    auto F = [](double x, double y) { return -2 * x * y; }; // <-- ACA se elige el problema

    double x0, y0, xF, h;
    int n;
    leerDatosComunes(x0, y0, xF, n, h);

    vector<double> x(n + 1), y(n + 1);
    x[0] = x0;
    y[0] = y0;
    for (int i = 0; i < n; i++) {
        x[i + 1] = x[i] + h;
        double xm = (x[i] + x[i + 1]) / 2.0;
        double ym = y[i] + (h / 2.0) * F(x[i], y[i]); // Euler con medio paso
        y[i + 1] = y[i] + h * F(xm, ym);
    }

    mostrarYExportar(x, y, "puntoMedio");
}
