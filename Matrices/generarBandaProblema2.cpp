#include <iostream>
#include <fstream>
using namespace std;

// Genera el sistema tridiagonal (banda p=1) del Problema 2:
// A[i][i]=2, A[i][i-1]=A[i][i+1]=1 (si existen), b[i]=6 salvo en los
// extremos (i=1 o i=n) donde b[i]=4.5. Lo escribe como matriz aumentada
// [A|b] en un archivo de texto, listo para leer con leerSistema().
int main() {
    int n;
    cout << "Ingrese n (tamano del sistema, el del enunciado es 100): ";
    cin >> n;

    string nombreArchivo;
    cout << "Nombre del archivo a generar: ";
    cin >> nombreArchivo;

    ofstream archivo(nombreArchivo);
    if (!archivo.is_open()) {
        cerr << "No se pudo crear el archivo." << endl;
        return 1;
    }

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            double valor = 0;
            if (j == i)                  valor = 2;
            else if (j == i - 1 || j == i + 1) valor = 1;
            archivo << valor << " ";
        }
        double b = (i == 1 || i == n) ? 4.5 : 6.0;
        archivo << b << endl;
    }

    archivo.close();
    cout << "Listo: " << nombreArchivo << " generado (n=" << n << ")." << endl;
    return 0;
}
