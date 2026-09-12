#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

// Newton-Raphson para el Problema 1 del parcial modelo:
// f(x) = (x+1)/(x+4) - 0.25x , x0 = 1.0
// f'(x) = 3/(x+4)^2 - 0.25   (derivada del cociente + de -0.25x)
// Error segun el mismo criterio (#) de la biseccion: x_viejo en la
// primera iteracion es x0 (no hay "a" y "b" aca, es metodo abierto).
int main() {
    double x_viejo, x_nuevo, error;
    int iter = 0, max_iter = 100;
    double tol = 0.01; // mismo 0.01% que en el enunciado de biseccion

    auto f = [](double x) { return ((x+1)/(x+4)) - 0.25*x; };
    auto df = [](double x) { return 3.0/((x+4)*(x+4)) - 0.25; };

    cout << "Ingrese el valor inicial x0: ";
    cin >> x_viejo;

    do {
        if (fabs(df(x_viejo)) < 1e-4) {
            cerr << "La derivada es muy chica. No se puede continuar." << endl;
            return 1;
        }
        x_nuevo = x_viejo - f(x_viejo) / df(x_viejo);
        error = (fabs(x_nuevo - x_viejo) / (0.5 * fabs(x_nuevo + x_viejo))) * 100;

        iter++;
        cout << "Iteracion " << iter << ": x = " << x_nuevo << ", error = " << error << endl;

        x_viejo = x_nuevo;
    } while (error > tol && iter < max_iter);

    cout << fixed << setprecision(10);
    cout << "La raiz es: " << x_nuevo << ", el error estimado es: " << error
         << "% y el numero de iteraciones es: " << iter << endl;
    cout.unsetf(ios::fixed);
    cout << setprecision(6);
    return 0;
}
