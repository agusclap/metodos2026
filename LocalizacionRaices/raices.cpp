#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

void biseccion();
void falsaposicion();
void puntofijo();
void newtonRaphson();
void secante();
void raizNesima();

int main (void) {
    int opcion;
    do {
        cout << "Menu de opciones:" << endl;
        cout << "1. Biseccion [Metodo Cerrado]" << endl;
        cout << "2. Falsa posicion [Metodo Cerrado]" << endl;
        cout << "3. Punto fijo [Metodo Abierto]" << endl;
        cout << "4. Newton-Raphson [Metodo Abierto]" << endl;
        cout<< "5. Secante [Metodo Abierto]" << endl;
        cout << "6. Raiz n-esima (cuadrada/cubica) [Newton-Raphson especializado]" << endl;
        cout<< "7. Salir" << endl;
        cout << "Seleccione una opcion: ";
        cin >> opcion;

        switch (opcion) {
            case 1:
                biseccion();
                break;
            case 2:
                falsaposicion();
                break;
            case 3:
                puntofijo();
                break;
            case 4:
                newtonRaphson();
                break;
            case 5:
                secante();
                break;
            case 6:
                raizNesima();
                break;
            case 7:
                cout << "Saliendo del programa." << endl;
                break;
            default:
                cout << "Opcion no valida, intente de nuevo." << endl;
        }
    } while (opcion != 7);
    return 0;
}

void biseccion() {
    double a, b, c, c_viejo, fa, fb, fc, error, error_porcentual;
    int iter = 0, max_iter = 100;
    double tol = 1e-8;
    //auto f = [](double x) { return exp(x*x) - 2; };  // punto d: modelo ajustado (a=1,b=-2)
    auto f = [](double x) { return -15.3315 + 58.6036*x - 74.2202*x*x + 30.9504*x*x*x; };  // punto e: polinomio de Lagrange del punto b
    cout << "Ingrese el intervalo [a, b]: ";
    cin >> a >> b;

    fa = f(a);
    fb = f(b);

    if (fa * fb > 0) {
        cerr << "El intervalo no es valido." << endl;
        return;
    }
    c_viejo = a;
    //c_viejo = b; // el enunciado pide x_viejo=3.2=b en la primera iteracion, no a
    do {
        c = (a + b) / 2;
        fc = f(c);

        error = fabs(c - c_viejo);
        // OJO: (1/2) en C++ es division ENTERA (da 0, no 0.5) -> denominador
        // quedaba en 0 -> error=inf. Con 0.5 (o 1.0/2) se fuerza a punto flotante.
        //error = (fabs(c - c_viejo) / (0.5 * fabs(c + c_viejo))) * 100;
        error_porcentual = fabs(error / c) * 100; // solo informativo, no se usa para cortar
        c_viejo = c;
        cout << "Iteracion " << iter + 1 << ": a = " << a << ", b = " << b << ", c = " << c
             << ", error = " << error << ", error porcentual = " << error_porcentual << "%" << endl;

        if (fc == 0.0) {
            error = 0.0;
            error_porcentual = 0.0;
            break;
        }
        if (error < tol) {
            break;
        }

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
        iter++;
    } while (error > tol && iter < max_iter);

    cout << fixed << setprecision(10);
    cout << "La raiz es: " << c << ", el error estimado es: " << error
         << "% y el numero de iteraciones es: " << iter <<" error porcentual: "<<error_porcentual<<endl;
    cout.unsetf(ios::fixed);
    cout << setprecision(6);
}

void falsaposicion() {
    double a, b, c, c_viejo, fa, fb, fc, error = 0.0, error_porcentual;
    int iter = 0, max_iter = 2000;
    double tol = 1e-8;
    auto f = [](double x) { return -0.952854 - 0.206982*x + 2.12295*x*x; };
    cout << "Ingrese el intervalo [a, b]: ";
    cin >> a >> b;

    fa = f(a);
    fb = f(b);

    if (fa * fb > 0) {
        cerr << "El intervalo no es valido." << endl;
        return;
    }
    c_viejo = a;
    c = a; // Inicializar c para evitar uso de variable no inicializada
    do {
        c = (a * fb - b * fa) / (fb - fa);
        fc = f(c);
        if (fc == 0.0) {
            error = 0.0; // Raiz exacta encontrada: no hay diferencia que medir
            break;
        }

        if (fa* fc > 0){
            a = c; // Actualizar a si el signo es positivo
            fa = fc;
        } else {
            b = c; // Actualizar b si el signo es negativo
            fb = fc;
        }

        error = fabs(c - c_viejo);
        c_viejo = c;


        cout << "Iteracion " << iter + 1 << ": a = " << a << ", b = " << b << ", c = " << c
             << ", error = " << error << endl;

        if (error < tol && fc == 0.0) {
            break;
        }

        iter++;
    } while (error > tol && iter < max_iter);
    error_porcentual = ( error / fabs(c) ) * 100;
    cout << "La raiz es: " << c << " y el error es: "<<error<< " y el error porcentual es: "<<error_porcentual<<endl;
}



void puntofijo(){
    double x0, x1, error, gp;
    int iter = 0, max_iter = 10000;
    double tol = 1e-6;
    //auto g = [] (double x) {return sqrt((x+5)/2);};  // problema anterior: puntos fijos de g
    auto g = [] (double x) {return (sin(3*x) - log10(x))/2;};  // 3x+sin(x)-e^x=0 -> x=(e^x-sin(x))/3
    
    cout << "Ingrese el valor inicial x0: ";
    cin >> x0;
    gp = (g(x0 + 0.001) - g(x0))/ 0.001;
    do {
        if(abs(gp) >= 1){
            printf("El metodo no converge");
            return;
        }
        x1 = g(x0);
        error = fabs(x1 - x0);
        x0 = x1;
        cout << "Iteracion " << iter + 1 << ": x0 = " << x0 << ", x1 = " << x1
             << ", error = " << error << endl;

        
        iter++;
    } while (iter < max_iter && error > tol);

    // tol=1e-12 es mas estricta que los 6 digitos que cout muestra por default,
    // asi que hace falta pedir mas precision para que se vea de verdad.
    cout << fixed << setprecision(12);
    cout << "La raiz es: " << x1 << ", el error es: " << error << " y la cantidad total de iteraciones es: " << iter << endl;
    cout.unsetf(ios::fixed);
    cout << setprecision(6);
}


void newtonRaphson(){
    double x0, x1, error;
    int iter = 0, max_iter = 1000;
    double tol = 1e-6;
    auto f = [](double x) { return 2*x + log10(x) - sin(3*x); };
    //auto df = [](double x) { return pow(x,x) * (log(x) + 1); }; // derivada analitica (item a)
    auto df = [f](double x) { return ((3*f(x) - 4*f(x-0.01) + f(x-2*0.01))/(2*0.01)); }; // derivada aproximada, diferencia finita (item c: h=0.5, cambiar por 0.01 para el segundo caso)
    cout << "Ingrese el valor inicial x0: ";
    cin >> x0;

    do {
        if(fabs(df(x0)) < 1e-4) {
            cerr << "La derivada es pequeña. No se puede continuar." << endl;
            return;
        }
        x1 = x0 - f(x0) / df(x0);
        error = fabs(x1 - x0);
        //error = (fabs(x1 - x0) / (0.5 * fabs(x1 + x0))) * 100;
        cout << "Iteracion " << iter + 1 << ": x0 = " << x0 << ", x1 = " << x1
             << ", error = " << error << endl;

        x0 = x1;
        iter++;
    } while (iter < max_iter && error > tol);
    cout << "La raiz es: " << x1 << ", el error es: " << error << " y la cantidad total de iteraciones es: " << iter << endl;
    cout.unsetf(ios::fixed);
    cout << setprecision(6);
}

void secante() {
    double x0, x1, x2, error;
    int iter = 0, max_iter = 100;
    double tol = 1e-6;
    auto f = [](double x) { return x * x - 9 * x + 2; };
    cout << "Ingrese los valores iniciales x0 y x1: ";
    cin >> x0 >> x1;

    do {
        if(fabs(f(x1) - f(x0)) < 1e-6) {
            cerr << "Denominador pequenio. No se puede continuar." << endl;
            return;
        }
        x2 = x1 - ((f(x1) * (x1 - x0)) / (f(x1) - f(x0)));
        error = fabs(x2 - x1);
        cout << "Iteracion " << iter + 1 << ": x0 = " << x0 << ", x1 = " << x1 << ", x2 = " << x2
             << ", error = " << error << endl;

        x0 = x1;
        x1 = x2;
        iter++;
    } while (iter < max_iter && error > tol);

    cout << "La raiz es: " << x1 << ", el error es: " << error << " y la cantidad total de iteraciones es: " << iter << endl;

}

// Algoritmo de la raiz cuadrada y cubica (Newton-Raphson aplicado a f(x)=x^n - A).
// Formula (1), n=2:  p_k = (1/2)*(p_{k-1} + A/p_{k-1})
// Formula (2), n=3:  p_k = (2*p_{k-1} + A/p_{k-1}^2) / 3
void raizNesima() {
    double A, p, p_viejo, error;
    int n, iter = 0, max_iter = 100;
    double tol = 1e-12;

    cout << "Ingrese el numero A (se calcula su raiz n-esima): ";
    cin >> A;
    cout << "Ingrese n (2 = raiz cuadrada, 3 = raiz cubica): ";
    cin >> n;

    if (n != 2 && n != 3) {
        cerr << "n debe ser 2 o 3." << endl;
        return;
    }
    if (n == 2 && A < 0) {
        cerr << "A no puede ser negativo para raiz cuadrada (raiz no real)." << endl;
        return;
    }

    cout << "Ingrese el valor inicial p0: ";
    cin >> p;

    if (p == 0) {
        cerr << "p0 no puede ser 0 (se divide por p0)." << endl;
        return;
    }

    do {
        p_viejo = p;
        if (n == 2) {
            p = 0.5 * (p_viejo + A / p_viejo);
        } else { // n == 3
            p = (2 * p_viejo + A / (p_viejo * p_viejo)) / 3;
        }
        error = fabs(p - p_viejo);
        iter++;
        cout << "Iteracion " << iter << ": p = " << p << ", error = " << error << endl;
    } while (error > tol && iter < max_iter);

    cout << fixed << setprecision(12);
    cout << "La raiz " << n << "-esima de " << A << " es: " << p
         << ", el error es: " << error << " y la cantidad total de iteraciones es: " << iter << endl;
    cout.unsetf(ios::fixed);
    cout << setprecision(6);
}